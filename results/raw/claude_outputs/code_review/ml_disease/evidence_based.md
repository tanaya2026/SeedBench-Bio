# SeedBench-Bio Model Response

## Metadata

- **Prompt variant:** evidence based
- **Prompt file:** prompts/code/evidence_based.md
- **Model:**
provider: Anthropic
name: Claude Sonnet
thinking: medium
- **Task file:** benchmark_tasks/code_review/ml_disease_prediction.py
- **Ground truth file:** ground_truth/code_review/ml_disease_prediction.json
- **Timestamp:** 2026-08-06 12:53:00
- **Temperature:** 0

---
## Prompt sent to model

Review the following code.

Report only objective errors supported by evidence from the code.

For every reported issue provide:

1. The relevant code location.
2. Why the code is objectively incorrect.
3. The bioinformatics, statistical, or computational principle that is violated.
4. A corrected implementation.

If you cannot justify an issue with evidence from the code, do not report it.

Ignore style, formatting, efficiency, and naming conventions.

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
## Issue 1: Feature scaling fit before train/test split

**Location:**
```python
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)   # fit on the FULL dataset
```

**Why incorrect:** `fit_transform` computes the mean and standard deviation used for scaling from the entire dataset, including samples that are later placed in the test set. The test set's statistics therefore influence the transformation applied to the training data.

**Principle violated:** Data leakage / violation of train-test independence. Any parameter estimated from data (mean, variance, selected features, etc.) must be estimated only from the training partition. Using test-set information during preprocessing inflates apparent performance and invalidates the estimate of generalization error.

**Corrected implementation:**
```python
X_train_raw, X_test_raw, y_train, y_test = train_test_split(
    X, y, test_size=0.25
)

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train_raw)
X_test_scaled = scaler.transform(X_test_raw)
```

---

## Issue 2: Feature selection fit before train/test split

**Location:**
```python
selector = SelectKBest(score_func=f_classif, k=50)
X_selected = selector.fit_transform(X_scaled, y)   # fit on the FULL dataset
```

**Why incorrect:** `f_classif` computes an F-statistic between each gene and the label vector `y` using all samples, then selects the top 50. Because this uses the full label vector (including future test labels) before splitting, the selected feature set is chosen with knowledge of the test set's outcomes. When the split is later performed, the "test" samples are evaluated only on genes that were already selected in part because of their own labels.

**Principle violated:** This is a more severe form of data leakage than Issue 1, since it uses the target variable `y` on data that will later be held out. This inflates performance and is a well-documented pitfall in genomic classification studies (e.g., microarray/RNA-seq classifier overfitting from selecting genes before cross-validation).

**Corrected implementation:**
```python
selector = SelectKBest(score_func=f_classif, k=50)
X_train_selected = selector.fit_transform(X_train_scaled, y_train)
X_test_selected = selector.transform(X_test_scaled)
```

---

## Issue 3: Model evaluated on training data, not test data

**Location:**
```python
y_pred = clf.predict(X_train)
acc = accuracy_score(y_train, y_pred)

print(f"Model accuracy: {acc:.3f}")
```

**Why incorrect:** The code splits off `X_test`/`y_test` but never uses them. `accuracy_score` is computed by predicting on `X_train` and comparing to `y_train` — i.e., the training accuracy of the fitted model, not an estimate of performance on unseen data. The printed value is labeled "Model accuracy" but is actually resubstitution/training accuracy, which is a biased (optimistic) estimate of generalization performance, especially for a model like logistic regression that can fit training data closely.

**Principle violated:** Model evaluation must be performed on data not used for fitting. Reporting training-set accuracy as the model's performance conflates fit with generalization and does not measure predictive validity.

**Corrected implementation:**
```python
y_pred = clf.predict(X_test_selected)
acc = accuracy_score(y_test, y_pred)

print(f"Model accuracy: {acc:.3f}")
```

---

## Human notes 

- Errors identified: E1, E2, E4
- Missed errors:  E3, E5, E6
- Incorrect flags: 0