"""
Shared ML pipeline functions for the PD/AD/HC multiclass classifier.

Used by ml_01_all_taxa.ipynb, ml_02_significant_taxa.ipynb,
ml_03_significant_taxa_clinical.ipynb, ml_04_significant_taxa_clinical_pathway.ipynb,
and ml_02_significant_taxa_two_models.ipynb (the platform-matched binary variant).

Keep this file in the same folder as those notebooks (or add its folder to
sys.path) - it is imported, not copy-pasted, so that a fix here applies to
all configs at once instead of drifting between separately-edited copies.
"""

import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import StratifiedGroupKFold, GridSearchCV, LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import (
    roc_auc_score,
    f1_score,
    accuracy_score,
    confusion_matrix,
    classification_report,
)

import matplotlib.pyplot as plt
import joblib

RANDOM_STATE = 42

DEFAULT_PARAM_GRID = {
    "n_estimators": [100, 300, 500],
    "max_depth": [None, 10, 20],
    "min_samples_leaf": [1, 5, 10],
}


def match_taxon_columns(taxa, all_taxa_columns, label="taxa", prefix="taxon__", strip_prefix="g__"):
    """
    Map a list of bare feature names (e.g. 'g__Blautia', 'g__GGB9350' for
    genera, or 'P164-PWY' for MetaCyc pathway IDs) onto the matching
    prefixed columns in the master feature table (e.g. 'taxon__g__Blautia',
    'pathway__P164-PWY'). Handles the prefix, GGB codes, and bracketed names
    (e.g. 'g__[Eubacterium]_xylanophilum_group') exactly as they appear in
    the master table, so a feature list can be used directly without
    renaming.

    `prefix`/`strip_prefix` default to the taxon case ('taxon__'/'g__') for
    backward compatibility - pass prefix="pathway__", strip_prefix="" for
    pathway IDs, which have no 'g__'-style sub-prefix to strip.

    Warns (does not fail) on any feature with no matching column - this can
    legitimately happen if a genus/pathway was significant in ANCOM-BC2 or
    MaAsLin2 but didn't survive ml_00's prevalence filter (e.g.
    g__Granulicatella in the PD-specific taxa list, missing from the master
    table entirely).
    """
    matched_cols = []
    missing = []
    for feat in taxa:
        clean_name = feat.replace(strip_prefix, "") if strip_prefix else feat
        candidates = [
            col for col in all_taxa_columns
            if col == feat
            or col == f"{prefix}{feat}"
            or col == f"{prefix}{clean_name}"
            or col == clean_name
        ]
        if candidates:
            matched_cols.append(candidates[0])
        else:
            missing.append(feat)

    if missing:
        print(f"WARNING: {len(missing)}/{len(taxa)} {label} not found in the "
              f"master feature table (likely filtered out upstream in ml_00 "
              f"for low prevalence): {missing}")

    return matched_cols


def check_platform_leakage(df, feature_cols, disease_col="disease_group",
                            group_col="study_id", classes=("PD", "AD"), n_splits=5):
    """
    Study-grouped AUC for distinguishing `classes` using only `feature_cols`.

    This project has a documented platform confound (PD = 100% shotgun,
    AD = 100% 16S). A near-perfect AUC here (~0.95+) is a red flag that the
    feature set is encoding platform, not disease biology - treat the main
    multiclass result with real suspicion if this comes back high, rather
    than treating it as a bonus finding.

    NOTE: every study in this project is single-arm (a PD study or an AD
    study - never both), so a true leave-one-study-out split is structurally
    impossible here: any single held-out study is 100% one class, AUC is
    undefined for a fold with only one class present, and the check would
    return nan every time regardless of the features used. This uses
    StratifiedGroupKFold instead - studies are still never split across
    train/test (so this is still an unseen-study generalisation check, not
    a leaky within-study one), but several studies are pooled per fold so
    both classes are actually present in every test fold.
    """
    subset = df[df[disease_col].isin(classes)].reset_index(drop=True)
    X = subset[feature_cols].values
    y = (subset[disease_col] == classes[0]).astype(int).values
    groups = subset[group_col].values

    n_studies_per_class = subset.groupby(disease_col)[group_col].nunique()
    splits = min(n_splits, int(n_studies_per_class.min()))
    if splits < 2:
        print(f"Cannot run platform-leakage check: the smaller class has only "
              f"{int(n_studies_per_class.min())} study/studies (need >= 2 for a "
              f"study-grouped split). Skipping.")
        return np.nan, np.nan

    cv = StratifiedGroupKFold(n_splits=splits, shuffle=True, random_state=RANDOM_STATE)
    aucs = []
    for train_idx, test_idx in cv.split(X, y, groups=groups):
        if len(set(y[test_idx])) < 2:
            continue  # shouldn't happen with a stratified split, but guard anyway
        clf = RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                      random_state=RANDOM_STATE, n_jobs=-1)
        clf.fit(X[train_idx], y[train_idx])
        proba = clf.predict_proba(X[test_idx])[:, 1]
        aucs.append(roc_auc_score(y[test_idx], proba))

    if not aucs:
        print("WARNING: no fold had both classes present after stratified grouping - "
              "platform-leakage check inconclusive. Try lowering n_splits.")
        return np.nan, np.nan

    mean_auc, std_auc = np.mean(aucs), np.std(aucs)
    print(f"{classes[0]}-vs-{classes[1]} study-grouped {splits}-fold AUC using "
          f"{len(feature_cols)} features: {mean_auc:.3f} +/- {std_auc:.3f}")
    if mean_auc >= 0.95:
        print("WARNING: this is very high - likely platform, not disease biology, "
              "is driving separation. Flag this explicitly if reporting the main model.")
    return mean_auc, std_auc


def run_hyperparameter_search(X_scaled, y, groups, param_grid=None, n_splits=5):
    """StratifiedGroupKFold-tuned RandomForest via GridSearchCV, scored on f1_macro."""
    param_grid = param_grid or DEFAULT_PARAM_GRID
    cv_inner = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_STATE)
    rf = RandomForestClassifier(class_weight="balanced", random_state=RANDOM_STATE, n_jobs=-1)
    grid = GridSearchCV(rf, param_grid, cv=cv_inner, scoring="f1_macro", n_jobs=-1, verbose=1)
    grid.fit(X_scaled, y, groups=groups)
    print("Best params:", grid.best_params_)
    print("Best CV f1_macro:", grid.best_score_)
    return grid

def select_features_lasso(
    X_train,
    y_train,
    feature_cols,
    random_state=RANDOM_STATE,
    min_features=1,
):
    """
    L1-regularised Logistic Regression feature selection.

    IMPORTANT:
    This function receives TRAINING DATA ONLY when used inside LOSO.

    Returns:
        selected_indices
        selected_feature_names
        selector_scaler
        lasso_model
    """

    # L1 regularisation is scale-sensitive, so scale only for LASSO
    selector_scaler = StandardScaler()

    X_train_scaled = selector_scaler.fit_transform(X_train)

    n_classes = len(np.unique(y_train))

    if n_classes < 2:
        raise ValueError(
            "LASSO feature selection requires at least two classes "
            "in the training data."
        )

    # Multiclass/binary L1 logistic regression
    lasso_model = LogisticRegressionCV(
        Cs=10,
        cv=5,
        penalty="l1",
        solver="saga",
        scoring="f1_macro",
        class_weight="balanced",
        max_iter=10000,
        n_jobs=-1,
        random_state=random_state,
    )

    lasso_model.fit(X_train_scaled, y_train)

    # Binary case: coef_.shape = (1, n_features)
    # Multiclass case: coef_.shape = (n_classes, n_features)
    #
    # Keep a feature if it has a non-zero coefficient
    # for at least one class.
    importance = np.max(
        np.abs(lasso_model.coef_),
        axis=0
    )

    selected_indices = np.where(
        importance > 1e-8
    )[0]

    # Safety fallback
    if len(selected_indices) < min_features:

        n_keep = min(
            min_features,
            len(feature_cols)
        )

        selected_indices = np.argsort(
            importance
        )[-n_keep:]

        print(
            f"WARNING: LASSO selected fewer than "
            f"{min_features} features. "
            f"Keeping top {n_keep} coefficient-ranked features."
        )

    selected_indices = np.sort(
        selected_indices
    )

    selected_feature_names = [
        feature_cols[i]
        for i in selected_indices
    ]

    return (
        selected_indices,
        selected_feature_names,
        selector_scaler,
        lasso_model,
    )

def run_loso_evaluation(X_scaled, y, groups, best_params, label_encoder, output_dir):
    """
    Leave-One-Study-Out evaluation. This is the number that should go in the
    thesis - not the GridSearchCV score above, which is in-sample-ish relative
    to this held-out-by-study evaluation.

    IMPORTANT (multiclass case, n_classes > 2): every study in this dataset is
    single-arm (PD+HC or AD+HC only - no study contains all three classes), so
    no individual fold can ever contain all 3 classes in y_true. A per-fold
    multiclass OVR AUC is therefore undefined for every fold, and averaging
    per-fold AUCs always returns nan - this is structural, not a bug you can
    tune away by fold count or class weighting. The fix is to pool the
    out-of-fold predicted probabilities across ALL folds first (where every
    class is guaranteed to appear at least once across the full study sweep)
    and compute ONE global AUC from the pooled predictions. Per-fold F1_macro
    is still reported per fold as a diagnostic, but treat it as approximate
    for any fold where n_classes_present < n_classes - it isn't averaged into
    the headline numbers below.

    BINARY case (n_classes == 2, e.g. a platform-matched PD-vs-HC or
    AD-vs-HC model): sklearn's roc_auc_score rejects the multi_class="ovr"
    argument for a 2-class target (it expects a strictly multiclass target
    when that argument is supplied), so the pooled AUC is computed the plain
    binary way instead, from the positive-class column of the pooled
    out-of-fold probabilities. Everything else below (pooling logic, per-fold
    F1_macro, saved files) is identical between the two cases.
    """
    logo = LeaveOneGroupOut()
    fold_results = []
    all_y_true, all_y_pred, all_y_proba = [], [], []
    n_classes = len(label_encoder.classes_)
    class_labels = np.arange(n_classes)

    for fold_i, (train_idx, test_idx) in enumerate(logo.split(X_scaled, y, groups=groups)):
        held_out_study = groups[test_idx][0]
        model = RandomForestClassifier(**best_params, class_weight="balanced",
                                        random_state=RANDOM_STATE, n_jobs=-1)
        model.fit(X_scaled[train_idx], y[train_idx])
        y_pred = model.predict(X_scaled[test_idx])
        y_proba = model.predict_proba(X_scaled[test_idx])

        all_y_true.extend(y[test_idx])
        all_y_pred.extend(y_pred)
        all_y_proba.append(y_proba)

        present_classes = np.unique(y[test_idx])
        f1_macro = f1_score(y[test_idx], y_pred, average="macro",
                             labels=class_labels, zero_division=0)
        fold_results.append({
            "held_out_study": held_out_study,
            "n_test_samples": len(test_idx),
            "n_classes_present": len(present_classes),
            "f1_macro": f1_macro,
        })
        print(f"Fold {fold_i} (held out {held_out_study}): "
              f"n={len(test_idx)}, classes_present={len(present_classes)}/{n_classes}, "
              f"F1_macro={f1_macro:.3f}")

    loso_df = pd.DataFrame(fold_results)
    loso_df.to_csv(Path(output_dir) / "loso_fold_results.csv", index=False)

    # --- Pooled headline metrics (the numbers for the thesis) ---
    all_y_true_arr = np.array(all_y_true)
    all_y_proba_arr = np.vstack(all_y_proba)

    if n_classes == 2:
        # Binary target: sklearn's roc_auc_score does not accept multi_class
        # for a 2-class problem, so score directly off the positive-class
        # (second column) probability instead of the OVR path used below.
        pooled_auc = roc_auc_score(all_y_true_arr, all_y_proba_arr[:, 1])
    else:
        pooled_auc = roc_auc_score(all_y_true_arr, all_y_proba_arr, multi_class="ovr",
                                    average="macro", labels=class_labels)
    pooled_f1 = f1_score(all_y_true_arr, all_y_pred, average="macro",
                          labels=class_labels, zero_division=0)

    print("\n--- LOSO summary (pooled out-of-fold predictions across all folds) ---")
    print(f"Pooled AUC ({'binary' if n_classes == 2 else 'ovr, macro'}): {pooled_auc:.3f}")
    print(f"Pooled F1_macro: {pooled_f1:.3f}")
    print(f"(Per-fold F1_macro, reference only, many folds missing a class by design: "
          f"{loso_df['f1_macro'].mean():.3f} +/- {loso_df['f1_macro'].std():.3f})")

    with open(Path(output_dir) / "loso_pooled_summary.txt", "w") as f:
        f.write(f"Pooled AUC ({'binary' if n_classes == 2 else 'ovr, macro'}): {pooled_auc:.4f}\n")
        f.write(f"Pooled F1_macro: {pooled_f1:.4f}\n")
        f.write(f"N folds: {len(loso_df)}\n")
        f.write(f"N samples: {len(all_y_true_arr)}\n")

    return loso_df, all_y_true, all_y_pred, all_y_proba

def run_loso_evaluation_with_lasso(X,y,groups,feature_cols,best_params,label_encoder,output_dir,min_features=1,):
    """
    Leave-One-Study-Out evaluation with fold-wise LASSO feature selection.
    Workflow for each LOSO fold:
        training studies only
            -> LASSO feature selection
            -> Random Forest using selected features
            -> predict held-out study
    Returns:
        loso_df
        all_y_true
        all_y_pred
        all_y_proba
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True,exist_ok=True)

    logo = LeaveOneGroupOut()

    all_y_true = []
    all_y_pred = []
    all_y_proba = []

    fold_results = []

    n_classes = len(label_encoder.classes_)
    class_labels = np.arange(n_classes)

    for fold_i, (train_idx,test_idx) in enumerate(logo.split(X,y,groups=groups)):
        held_out_study = groups[test_idx][0]
        X_train = X[train_idx]
        X_test = X[test_idx]
        y_train = y[train_idx]
        y_test = y[test_idx]

        # ==================================================
        # 1. LASSO FEATURE SELECTION
        # TRAINING DATA ONLY
        # ==================================================

        (selected_indices,selected_feature_names,selector_scaler,lasso_model,) = select_features_lasso(X_train=X_train, y_train=y_train,feature_cols=feature_cols,min_features=min_features,)

        print(f"\nFold {fold_i} (held out {held_out_study})")
        print(f"LASSO selected {len(selected_feature_names)}/{len(feature_cols)} features")

        # ==================================================
        # 2. RESTRICT TO SELECTED FEATURES
        # RF does not need scaling, so use original feature values after LASSO determines which columns to keep.
        # ==================================================
        X_train_selected = X_train[:,selected_indices]
        X_test_selected = X_test[:,selected_indices]

        # ==================================================
        # 3. RANDOM FOREST
        # ==================================================
        model = RandomForestClassifier(**best_params,class_weight="balanced",random_state=RANDOM_STATE,n_jobs=-1,)
        model.fit(X_train_selected,y_train)
        y_pred = model.predict(X_test_selected)
        y_proba = model.predict_proba(X_test_selected)

        # ==================================================
        # 4. STORE POOLED PREDICTIONS
        # ==================================================
        all_y_true.extend(y_test)
        all_y_pred.extend(y_pred)
        all_y_proba.append(y_proba)

        # ==================================================
        # 5. FOLD METRICS
        # ==================================================
        accuracy = accuracy_score(y_test,y_pred)
        f1_macro = f1_score(y_test,y_pred,average="macro",labels=class_labels,zero_division=0,)
        present_classes = np.unique(y_test)

        # Per-fold AUC
        # Multiclass AUC is only valid if all classes are represented in that held-out fold.
        # Binary AUC is valid only if both classes occur.
        fold_auc = np.nan
        try:
            if n_classes == 2:
                if len(present_classes) == 2:
                    fold_auc = roc_auc_score(y_test,y_proba[:, 1])
            else:
                if len(present_classes) == n_classes:
                    fold_auc = roc_auc_score(y_test,y_proba,multi_class="ovr",average="macro",labels=class_labels,)
        except ValueError:
            fold_auc = np.nan

        # ==================================================
        # 6. ACTUAL / PREDICTED CLASS COUNTS
        # ==================================================
        actual_counts = {}
        predicted_counts = {}
        for class_idx, class_name in enumerate(label_encoder.classes_):
            actual_counts[f"actual_{class_name}"] = int(np.sum(y_test == class_idx))
            predicted_counts[f"predicted_{class_name}"] = int(np.sum(y_pred == class_idx))

        # ==================================================
        # 7. SAVE FOLD RESULT
        # ==================================================
        fold_result = {"held_out_study":held_out_study,
                       "n_test_samples":len(test_idx),
                       "n_classes_present":len(present_classes),
                       "n_selected_features":len(selected_feature_names),
                       "selected_features":";".join(selected_feature_names),
                       "accuracy":accuracy,
                       "f1_macro":f1_macro,
                       "auc":fold_auc,
                       }
        fold_result.update(actual_counts)
        fold_result.update(predicted_counts)
        fold_results.append(fold_result)
        print(f"n={len(test_idx)}, classes_present={len(present_classes)}/{n_classes}, "
              f"Accuracy={accuracy:.3f}, F1_macro={f1_macro:.3f}, AUC={fold_auc:.3f}"
            if not np.isnan(fold_auc)
            else
            f"n={len(test_idx)}, "
            f"classes_present="
            f"{len(present_classes)}/{n_classes}, "
            f"Accuracy={accuracy:.3f}, "
            f"F1_macro={f1_macro:.3f}, "
            f"AUC=NA"
        )

    # ======================================================
    # STUDY-LEVEL RESULTS
    # ======================================================
    loso_df = pd.DataFrame(fold_results)
    loso_df.to_csv(output_dir / "loso_study_level_metrics_with_lasso.csv",index=False,)

    # ======================================================
    # POOLED RESULTS
    # ======================================================
    all_y_true_arr = np.array(all_y_true)
    all_y_pred_arr = np.array(all_y_pred)
    all_y_proba_arr = np.vstack(all_y_proba)
    pooled_accuracy = accuracy_score(all_y_true_arr,all_y_pred_arr)
    pooled_f1 = f1_score(all_y_true_arr,all_y_pred_arr,average="macro",labels=class_labels,zero_division=0,)

    if n_classes == 2:
        pooled_auc = roc_auc_score(all_y_true_arr,all_y_proba_arr[:, 1])
    else:
        pooled_auc = roc_auc_score(all_y_true_arr,all_y_proba_arr,multi_class="ovr",average="macro",labels=class_labels,)
    print("\n--- LOSO summary (pooled out-of-fold predictions across all folds) ---")

    print(f"Pooled AUC ({'binary' if n_classes == 2 else 'ovr, macro'}): {pooled_auc:.3f}")
    print(f"Pooled Accuracy: {pooled_accuracy:.3f}")
    print(f"Pooled F1_macro: {pooled_f1:.3f}")
    print(f"Per-fold F1_macro: {loso_df['f1_macro'].mean():.3f} +/- {loso_df['f1_macro'].std():.3f}")
    with open(output_dir / "loso_pooled_summary_with_lasso.txt","w") as f:
        f.write(f"Pooled AUC: {pooled_auc:.4f}\n")
        f.write(f"Pooled Accuracy: {pooled_accuracy:.4f}\n")
        f.write(f"Pooled F1_macro: {pooled_f1:.4f}\n")
        f.write(f"N folds: {len(loso_df)}\n")
        f.write(f"N samples: {len(all_y_true_arr)}\n")
    return (loso_df,all_y_true,all_y_pred,all_y_proba)

def plot_confusion_matrix(all_y_true, all_y_pred, label_encoder, output_dir):
    """Confusion matrix from pooled LOSO out-of-fold predictions (not final-model self-predictions)."""
    cm = confusion_matrix(all_y_true, all_y_pred, labels=np.arange(len(label_encoder.classes_)))
    cm_df = pd.DataFrame(cm, index=label_encoder.classes_, columns=label_encoder.classes_)
    print(cm_df)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.imshow(cm, cmap="Blues")
    ax.set_xticks(range(len(label_encoder.classes_)))
    ax.set_yticks(range(len(label_encoder.classes_)))
    ax.set_xticklabels(label_encoder.classes_)
    ax.set_yticklabels(label_encoder.classes_)
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_title("Confusion matrix (pooled LOSO out-of-fold predictions)")
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, cm[i, j], ha="center", va="center")
    plt.tight_layout()
    plt.savefig(Path(output_dir) / "confusion_matrix_loso.png", dpi=300)
    plt.show()

    print(classification_report(all_y_true, all_y_pred, target_names=label_encoder.classes_))
    return cm_df


def run_shap_analysis(final_model, X_scaled, feature_cols, label_encoder, output_dir, max_display=None):
    """
    SHAP TreeExplainer, per-class beeswarm. `final_model` should be fit on ALL
    data (not a LOSO fold) - this is for interpretation, not a performance claim.
    """
    import shap
    explainer = shap.TreeExplainer(final_model)
    shap_values = explainer.shap_values(X_scaled)

    # Newer shap returns a single (n_samples, n_features, n_classes) array;
    # older shap returns a list of (n_samples, n_features) arrays, one per class.
    if isinstance(shap_values, list):
        get_class_values = lambda idx: shap_values[idx]
    else:
        get_class_values = lambda idx: shap_values[:, :, idx]

    for class_idx, class_name in enumerate(label_encoder.classes_):
        print(f"Generating SHAP beeswarm for class: {class_name}")
        kwargs = dict(feature_names=feature_cols, show=False)
        if max_display:
            kwargs["max_display"] = max_display
        shap.summary_plot(get_class_values(class_idx), X_scaled, **kwargs)
        plt.title(f"SHAP summary — class: {class_name}")
        plt.tight_layout()
        plt.savefig(Path(output_dir) / f"shap_beeswarm_{class_name}.png", dpi=300, bbox_inches="tight")
        plt.show()
        plt.clf()

    return shap_values


def save_artifacts(model, scaler, label_encoder, feature_cols, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, output_dir / "rf_final_model.joblib")
    joblib.dump(scaler, output_dir / "scaler.joblib")
    joblib.dump(label_encoder, output_dir / "label_encoder.joblib")
    with open(output_dir / "feature_list.txt", "w") as f:
        f.write("\n".join(feature_cols))
    print("Saved model, scaler, label encoder, and feature list to:", output_dir)

def scale_within_platform(df, feature_cols, platform_col="platform"):
    """
    Z-score each feature separately within each platform group, instead of
    globally. Returns (X_scaled, scalers) where `scalers` is a dict of
    {platform_name: fitted StandardScaler} - keep this dict if you need to
    transform new samples later, since there's no single global scaler.
    """
    X_scaled = np.zeros((len(df), len(feature_cols)), dtype=float)
    scalers = {}
    for platform, group_df in df.groupby(platform_col):
        idx = df.index.get_indexer(group_df.index)
        s = StandardScaler()
        X_scaled[idx] = s.fit_transform(group_df[feature_cols].values)
        scalers[platform] = s
    return X_scaled, scalers

def fit_final_lasso_selector(X,y,feature_cols,min_features=1,):
    """
    Fit LASSO feature selection on ALL available data.
    This is used AFTER LOSO validation to define the final
    feature set for the final exported model and SHAP analysis.
    """
    (selected_indices,selected_feature_names,selector_scaler,lasso_model,) = select_features_lasso(X_train=X,y_train=y, feature_cols=feature_cols,min_features=min_features,)
    print("\nFinal full-data LASSO selection:")
    print(f"{len(selected_feature_names)} / {len(feature_cols)} features retained")
    return (selected_indices,selected_feature_names,selector_scaler,lasso_model,)
