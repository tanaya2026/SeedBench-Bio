#!/usr/bin/env python
# ==============================================================================
# ml_model_disease.py
#
# Pipeline: predict disease status (case/control) from a gene expression
# feature matrix (samples x genes).
#
# NOTE FOR BENCHMARK USERS: this script contains INTENTIONALLY SEEDED ERRORS.
# Do not "fix" this file — it is a fixed stimulus for SeedBench-Bio.
# See ground_truth/ml_disease_prediction.json for the answer key.
# ==============================================================================

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