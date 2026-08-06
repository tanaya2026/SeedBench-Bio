| task    | prompt   | TP | FP | FN | Precision | Recall |    F1 | Critical Error Recall | Reasoning Quality |
| ------- | -------- | -: | -: | -: | --------: | -----: | ----: | --------------------: ||
| bulkRNA | baseline |  5 |  0 |  3 |      1.00 |  0.625 | 0.769 | 1                     ||
| seurat_qc | baseline |  5 |  0 |  0 |      1.00 |  1.00 | 1.00 | 1                     ||
| ml_disease | baseline |  5 |  0 |  1 |      1.00 |  0.833 | 0.907 | 1                  ||
| claims | baseline | 7 | 0 | 0| 1.00 | 1.00 | 1.00 ||1.571|
| claims | group_review | 7 | 0 | 0| 1.00 | 1.00 | 1.00 ||1.428|
| claims | senior_review | 7 | 0 | 0| 1.00 | 1.00 | 1.00 ||1.428|
|bulkRNA | evidence_based | 3 |0 | 5 | 1.00 | 0.375 | 0.5454 | 0.6||
|ml_disease| evidence_based | 3| 0| 3 | 1.00 | 0.5 | 0.666|0.75||
| seurat_qc | evidence_based | 4| 0 |1 |1.00 | 0.8 |0.88 |1.00 ||
|bulkRNA| blind | 7| 0 | 1 | 1.00 | 0.875 | 0.933| 0.8 ||
|ml_disease| blind | 5 | 0| 1 | 1.00 | 0.833 | 0.907 | 1| ||
| seurat_qc | blind | 4| 0 |1 |1.00 | 0.8 |0.88 |1.00 ||
| claims | blind | 7 | 0 | 0| 1.00 | 1.00 | 1.00 ||1.571|