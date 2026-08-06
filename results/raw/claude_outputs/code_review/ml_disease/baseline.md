# SeedBench-Bio Model Response

## Metadata

- **Prompt variant:** baseline
- **Prompt file:** prompts/code/baseline.md
- **Model:**
provider: Anthropic
name: Claude Sonnet
thinking: medium
- **Task file:** benchmark_tasks/code_review/ml_disease_prediction.py
- **Ground truth file:** ground_truth/code_review/ml_disease_prediction.json
- **Timestamp:** 2026-08-05 19:23:00
- **Temperature:** 0

---
## Prompt sent to model

Review the following bioinformatics code.

Identify any objective scientific, statistical, or computational errors.

For each error, explain why it is incorrect and suggest a correction.

Ignore coding style, formatting, and personal preferences.
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score

# ------------------------------------------------------------------------
# 1. Load data
#    expression_features.csv: samples x genes (rows = samples, cols = genes)
#    labels.csv: sample_id, disease_status (0 = control, 1 = case)
# ------------------------------------------------------------------------
X = pd.read_csv("expression_features.csv", index_col=0)
y = pd.read_csv("labels.csv", index_col=0)["disease_status"]

# ------------------------------------------------------------------------
# 2. Scale features
#    Standardize the expression matrix before modeling so that all genes
#    are on a comparable scale.
# ------------------------------------------------------------------------
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)   # fit on the FULL dataset

# ------------------------------------------------------------------------
# 3. Feature selection
#    Reduce dimensionality by keeping the top 50 most informative genes
#    before splitting, so training and testing use the same feature set.
# ------------------------------------------------------------------------
selector = SelectKBest(score_func=f_classif, k=50)
X_selected = selector.fit_transform(X_scaled, y)   # fit on the FULL dataset

# ------------------------------------------------------------------------
# 4. Train/test split
# ------------------------------------------------------------------------
X_train, X_test, y_train, y_test = train_test_split(
    X_selected, y, test_size=0.25
)

# ------------------------------------------------------------------------
# 5. Train model
# ------------------------------------------------------------------------
clf = LogisticRegression(max_iter=1000)
clf.fit(X_train, y_train)

# ------------------------------------------------------------------------
# 6. Evaluate model performance
# ------------------------------------------------------------------------
y_pred = clf.predict(X_train)
acc = accuracy_score(y_train, y_pred)

print(f"Model accuracy: {acc:.3f}")
print("Model training complete.")

---

## Model response

Here's a review of the pipeline. There are several serious methodological errors — most of them forms of **data leakage** and **evaluation errors** that will make the reported accuracy meaningless (almost certainly inflated).

## 1. Scaling fit on the full dataset before the train/test split (data leakage)
```python
X_scaled = scaler.fit_transform(X)   # fit on the FULL dataset
```
`StandardScaler` computes per-gene mean/SD. Fitting it on the full dataset lets information from the test samples (their contribution to the mean/SD) leak into the "training" transformation. The scaler must be fit **only on the training set**, then used to transform both train and test:
```python
X_train_raw, X_test_raw, y_train, y_test = train_test_split(X, y, test_size=0.25, stratify=y, random_state=0)
scaler = StandardScaler().fit(X_train_raw)
X_train = scaler.transform(X_train_raw)
X_test = scaler.transform(X_test_raw)
```

## 2. Feature selection fit on the full dataset, using labels (severe data leakage)
```python
selector = SelectKBest(score_func=f_classif, k=50)
X_selected = selector.fit_transform(X_scaled, y)   # fit on the FULL dataset
```
This is the most serious error. `SelectKBest` uses `y` to score every gene by its association with disease status across **all** samples, including the ones that will later become the "test" set. The 50 genes chosen are therefore selected partly *because* they discriminate well on the test samples too. This is a classic and well-documented source of massively overoptimistic accuracy in omics pipelines (this exact mistake — selecting features on the full dataset before CV/split — is one of the most common reported errors in microarray/RNA-seq classification papers). Feature selection must be done **inside** the training fold only:
```python
selector = SelectKBest(score_func=f_classif, k=50).fit(X_train, y_train)
X_train_sel = selector.transform(X_train)
X_test_sel = selector.transform(X_test)
```
If you eventually want to use cross-validation (recommended, see point 6), scaling and feature selection should be wrapped in a `Pipeline` so they're refit within each fold rather than done once globally.

## 3. Model evaluated on the training set, not the test set
```python
y_pred = clf.predict(X_train)
acc = accuracy_score(y_train, y_pred)
```
This defeats the entire purpose of doing a train/test split — you're measuring how well the model fits data it already saw, which is not a valid estimate of generalization performance. It should be:
```python
y_pred = clf.predict(X_test)
acc = accuracy_score(y_test, y_pred)
```

## 4. No stratification of the train/test split
```python
X_train, X_test, y_train, y_test = train_test_split(X_selected, y, test_size=0.25)
```
Disease case/control cohorts are frequently imbalanced. Without `stratify=y`, a random split can produce a training or test set with a very different case/control ratio than the full dataset (or, with small n, a test set with very few/no cases of one class), which biases and destabilizes the accuracy estimate. Use `train_test_split(..., stratify=y)`.

## 5. Accuracy alone is an inappropriate metric here
For case/control classification, especially with any class imbalance (common in disease cohorts), accuracy can look high simply by predicting the majority class. This isn't a coding bug but a statistical validity issue: report balanced accuracy, ROC-AUC, precision/recall, or F1 alongside (or instead of) raw accuracy.

## 6. Single train/test split is a weak/high-variance estimator for typical expression-data sample sizes
Gene expression studies often have small sample sizes (tens to low hundreds of samples) and high dimensionality. A single 75/25 split gives a highly variable, easily "lucky"/"unlucky" performance estimate. The standard approach is k-fold (or stratified k-fold / repeated k-fold) cross-validation, with scaling and feature selection refit inside each fold to avoid the leakage in points 1–2.

## 7. No verification that `X` and `y` are sample-aligned
```python
X = pd.read_csv("expression_features.csv", index_col=0)
y = pd.read_csv("labels.csv", index_col=0)["disease_status"]
...
X_scaled = scaler.fit_transform(X)   # returns a numpy array — index is dropped
```
Once `X` is converted to a numpy array, all positional alignment with `y` depends on the two files having identical row order. If `labels.csv` and `expression_features.csv` are not already sorted identically by `sample_id`, samples will be silently mismatched with the wrong labels, corrupting everything downstream. This should be checked/enforced explicitly, e.g.:
```python
X, y = X.align(y, join="inner", axis=0)
assert (X.index == y.index).all()
```

---

### Summary of required fixes (in order)
1. Split data into train/test (stratified) **first**.
2. Fit `StandardScaler` on train only; transform both.
3. Fit `SelectKBest` on train only (post-scaling); transform both.
4. Evaluate on `X_test`/`y_test`, not `X_train`/`y_train`.
5. Report more than accuracy (AUC/F1/balanced accuracy).
6. Prefer cross-validation (ideally via a `Pipeline` wrapping scaling + selection + classifier) over a single split, given typical expression-data sample sizes.
7. Explicitly verify sample-ID alignment between `X` and `y` before modeling.

Errors 1, 2, and 3 (leakage in scaling, leakage in feature selection, and testing on training data) are the most severe — each independently would produce an accuracy number that overstates real-world predictive performance, and together they compound.

---

## Human notes 

- Errors identified:E1-E5
- Missed errors: E6
- Incorrect flags: 0