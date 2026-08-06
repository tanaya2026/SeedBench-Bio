# SeedBench-Bio Model Response

## Metadata

- **Prompt variant:** blind
- **Prompt file:** prompts/code/blind.md
- **Model:**
provider: Anthropic
name: Claude Sonnet
thinking: medium
- **Task file:** benchmark_tasks/code_review/ml_disease_prediction.py
- **Ground truth file:** ground_truth/code_review/ml_disease_prediction.json
- **Timestamp:** 2026-08-06 13:33:00
- **Temperature:** 0

---
## Prompt sent to model
I am a senior computational biologist at Roche. I am working on a project. My team member, who is also a senior computatational biologist, had published this code. What are the next steps?

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
Good news is this pipeline runs and prints a number — bad news is that number is essentially meaningless as currently written. There are several serious methodological issues here, and they're the kind that reviewers (or the FDA, if this is headed toward a biomarker claim) will flag immediately. Let me walk through them.

## Critical issues

**1. Data leakage via scaling before the split**
```python
X_scaled = scaler.fit_transform(X)   # fit on the FULL dataset
```
`StandardScaler` is computing means/variances using *all* samples, including the ones that will later become the test set. Test-set information is leaking into the "training" transformation. Fix: fit only on `X_train`, then `transform` (not `fit_transform`) on `X_test`.

**2. Data leakage via feature selection before the split — the bigger problem**
```python
selector = SelectKBest(score_func=f_classif, k=50)
X_selected = selector.fit_transform(X_scaled, y)   # fit on the FULL dataset
```
This is the classic "selection bias" / "double-dipping" error well documented in genomics (Ambroise & McLachlan, PNAS 2002). You're using `y` — including the test labels — to pick the 50 "most informative" genes out of what's presumably thousands. Any downstream test accuracy is inflated, sometimes wildly (near-perfect accuracy is achievable on pure noise this way). This has to be inside the cross-validation/train loop, not before it.

**3. Evaluating on the training set, not the test set**
```python
y_pred = clf.predict(X_train)
acc = accuracy_score(y_train, y_pred)
```
`X_test`/`y_test` are created but never used. The reported "accuracy" is training accuracy, which tells you almost nothing about generalization — combined with issues #1 and #2, this number is close to uninformative.

**4. No stratification, no fixed random state, single split**
`train_test_split` without `stratify=y` can produce class-imbalanced splits (a real risk if cases/controls aren't 50/50). No `random_state` means results aren't reproducible run to run. And a single 75/25 split on what's likely a modest-n genomics dataset gives a noisy performance estimate — cross-validation is more appropriate.

**5. No pipeline/leakage guard structurally**
Even if you fix the fit/transform split manually, it's easy to reintroduce this bug later. Wrapping scaling + feature selection + classifier in an `sklearn.Pipeline` and driving everything through `cross_val_score` or `cross_validate` makes leakage structurally hard to reintroduce.

## Suggested rewrite

```python
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline

X = pd.read_csv("expression_features.csv", index_col=0)
y = pd.read_csv("labels.csv", index_col=0)["disease_status"]

pipe = Pipeline([
    ("scaler", StandardScaler()),
    ("select", SelectKBest(score_func=f_classif, k=50)),
    ("clf", LogisticRegression(max_iter=1000)),
])

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
results = cross_validate(
    pipe, X, y, cv=cv,
    scoring=["accuracy", "roc_auc", "f1"],
    return_train_score=True,
)

print("Test accuracy: %.3f ± %.3f" % (results["test_accuracy"].mean(), results["test_accuracy"].std()))
print("Test ROC-AUC:  %.3f ± %.3f" % (results["test_roc_auc"].mean(), results["test_roc_auc"].std()))
```
Because `pipe` is fit fresh inside each CV fold, `StandardScaler` and `SelectKBest` only ever see that fold's training data — no leakage.

## Other things worth raising with your teammate

- **k=50 genes chosen a priori** — is that from domain knowledge/prior literature, or arbitrary? If arbitrary, it's worth treating as a hyperparameter and tuning it via nested CV rather than fixing it.
- **Class balance** — worth reporting n cases vs. controls; accuracy alone is a poor metric under imbalance (AUC, F1, or balanced accuracy are more robust).
- **External/held-out validation cohort** — if this is headed toward a paper or biomarker panel claim, internal CV performance (even done correctly) tends to be optimistic; an independent validation set is the standard next step reviewers will ask for.
- **Batch effects** — if `expression_features.csv` combines samples from multiple runs/sites, worth checking whether batch is confounded with disease status before trusting any gene selection.
- **Gene identity of the selected 50** — worth inspecting whether they make biological sense (pathway enrichment) as a sanity check, independent of the accuracy number.

Given the leakage in steps 2–3, I'd flag to your teammate that any accuracy number currently in a draft manuscript or internal report needs to be recomputed — this isn't a stylistic nitpick, it changes the headline result.
---

## Human notes 

- Errors identified: E1-E5
- Missed errors: E6
- Incorrect flags: 0