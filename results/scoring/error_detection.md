| run_id | task_id | prompt    | error_id | category                        | severity | found | Reasoning Quality |
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- |
| R001   | bulkRNA | baseline  | E1       | reference_mismatch              | critical | TRUE  ||
| R001   | bulkRNA | baseline  | E2       | missing_qc                      | moderate | FALSE ||
| R001   | bulkRNA | baseline  | E3       | missing_preprocessing           | moderate | FALSE ||
| R001   | bulkRNA | baseline  | E4       | inappropriate_deduplication     | critical | TRUE  ||
| R001   | bulkRNA | baseline  | E5       | strandedness_ignored            | critical | TRUE  ||
| R001   | bulkRNA | baseline  | E6       | missing_pairing_validation      | minor    | FALSE ||
| R001   | bulkRNA | baseline  | E7       | sample_identity_mismatch_risk   | critical | TRUE  ||
| R001   | bulkRNA | baseline  | E8       | missing_normalization           | critical | TRUE  ||
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- ||
| R002 | seurat_qc | baseline | E1 | preprocessing            | critical | TRUE ||
| R002 | seurat_qc | baseline | E2 | species_annotation            | critical | TRUE ||
| R002 | seurat_qc | baseline | E3 | statistics            | critical | TRUE ||
| R002 | seurat_qc | baseline | E4 | biological_annotation            | critical | TRUE ||
| R002 | seurat_qc | baseline | E5 | workflow_logic            | critical | TRUE ||
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- ||
| R003 | ml_disease | baseline | E1 | data_leakage_scaling            | critical | TRUE ||
| R003 | ml_disease | baseline | E2 | data_leakage_feature_selection            | critical | TRUE ||
| R003 | ml_disease | baseline | E3 | split_reproducibility_and_balance            | moderate | TRUE ||
| R003 | ml_disease | baseline | E4 | testing_on_training_data            | critical | TRUE ||
| R003 | ml_disease | baseline | E5 | wrong_metric            | critical | TRUE ||
| R003 | ml_disease | baseline | E6 | class_imbalance_unhandled            | moderate | TRUE ||
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- ||
| R004 | claims | baseline | E1 | small_n_enrichment_overreach | | TRUE| 2 |
| R004 | claims | baseline | E2 | unsupervised_result_given_unearned_label | | TRUE| 2 |
| R004 | claims | baseline | E3 | gene_set_overrepresentation_misread_as_universal | | TRUE| 2|
| R004 | claims | baseline | E4 | statistical_significance_conflated_with_expression_magnitude | | TRUE| 1 |
| R004 | claims | baseline | E5 | batch_effect_misread_as_biology | | TRUE| 2 |
| R004 | claims | baseline | E6 | correlation_causation | | TRUE|1 |
| R004 | claims | baseline | E7 | multiple_testing_ignored | | TRUE|1 |
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- |
| R005 | claims | group_review | E1 | small_n_enrichment_overreach | | TRUE| 1 |
| R005 | claims | group_review | E2 | unsupervised_result_given_unearned_label | | TRUE| 2 |
| R005 | claims | group_review | E3 | gene_set_overrepresentation_misread_as_universal | | TRUE| 1|
| R005 | claims | group_review | E4 | statistical_significance_conflated_with_expression_magnitude | | TRUE| 1 |
| R005 | claims | group_review | E5 | batch_effect_misread_as_biology | | TRUE| 2 |
| R005 | claims | group_review | E6 | correlation_causation | | TRUE|1 |
| R005 | claims | group_review | E7 | multiple_testing_ignored | | TRUE|2 |
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- |
| R006 | claims | senior_review | E1 | small_n_enrichment_overreach | | TRUE| 1 |
| R006 | claims | senior_review | E2 | unsupervised_result_given_unearned_label | | TRUE| 2 |
| R006 | claims | senior_review | E3 | gene_set_overrepresentation_misread_as_universal | | TRUE| 1|
| R006 | claims | senior_review | E4 | statistical_significance_conflated_with_expression_magnitude | | TRUE| 1 |
| R006 | claims | senior_review | E5 | batch_effect_misread_as_biology | | TRUE| 2 |
| R006 | claims | senior_review | E6 | correlation_causation | | TRUE|1 |
| R006 | claims | senior_review | E7 | multiple_testing_ignored | | TRUE|2 |
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- |
| R007   | bulkRNA | evidence_based  | E1       | reference_mismatch              | critical | FALSE ||
| R007   | bulkRNA | evidence_based  | E2       | missing_qc                      | moderate | FALSE ||
| R007   | bulkRNA | evidence_based  | E3       | missing_preprocessing           | moderate | FALSE ||
| R007   | bulkRNA | evidence_based  | E4       | inappropriate_deduplication     | critical | TRUE  ||
| R007   | bulkRNA | evidence_based  | E5       | strandedness_ignored            | critical | TRUE  ||
| R007   | bulkRNA | evidence_based  | E6       | missing_pairing_validation      | minor    | FALSE ||
| R007   | bulkRNA | evidence_based  | E7       | sample_identity_mismatch_risk   | critical | TRUE  ||
| R007   | bulkRNA | evidence_based  | E8       | missing_normalization           | critical | FALSE ||
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- ||
| R008 | ml_disease | evidence_based | E1 | data_leakage_scaling            | critical | TRUE ||
| R008 | ml_disease | evidence_based | E2 | data_leakage_feature_selection            | critical | TRUE ||
| R008 | ml_disease | evidence_based | E3 | split_reproducibility_and_balance            | moderate | FALSE ||
| R008 | ml_disease | evidence_based | E4 | testing_on_training_data            | critical | TRUE ||
| R008 | ml_disease | evidence_based | E5 | wrong_metric            | critical | FALSE ||
| R008 | ml_disease | evidence_based | E6 | class_imbalance_unhandled            | moderate | FALSE ||
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- ||
| R009 | seurat_qc | evidence_based | E1 | preprocessing            | critical | TRUE ||
| R009 | seurat_qc | evidence_based | E2 | species_annotation            | critical | TRUE ||
| R009 | seurat_qc | evidence_based | E3 | statistics            | critical | TRUE ||
| R009 | seurat_qc | evidence_based | E4 | biological_annotation            | critical | FALSE ||
| R009 | seurat_qc | evidence_based | E5 | workflow_logic            | critical | TRUE ||
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- ||
| R0010   | bulkRNA | blind  | E1       | reference_mismatch              | critical | FALSE ||
| R0010  | bulkRNA | blind   | E2       | missing_qc                      | moderate | TRUE ||
| R0010   | bulkRNA | blind   | E3       | missing_preprocessing           | moderate | TRUE ||
| R0010   | bulkRNA | blind   | E4       | inappropriate_deduplication     | critical | TRUE  ||
| R0010   | bulkRNA | blind   | E5       | strandedness_ignored            | critical | TRUE  ||
| R0010   | bulkRNA | blind   | E6       | missing_pairing_validation      | minor    | TRUE ||
| R0010   | bulkRNA | blind   | E7       | sample_identity_mismatch_risk   | critical | TRUE  ||
| R0010   | bulkRNA | blind   | E8       | missing_normalization           | critical | TRUE ||
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- ||
| R0011 | ml_disease | blind | E1 | data_leakage_scaling            | critical | TRUE ||
| R0011| ml_disease | blind | E2 | data_leakage_feature_selection            | critical | TRUE ||
| R0011 | ml_disease | blind | E3 | split_reproducibility_and_balance            | moderate | TRUE ||
| R0011 | ml_disease | blind | E4 | testing_on_training_data            | critical | TRUE ||
| R0011 | ml_disease | blind | E5 | wrong_metric            | critical | TRUE ||
| R0011 | ml_disease | blind | E6 | class_imbalance_unhandled            | moderate | FALSE ||
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- ||
| R0012 | seurat_qc | blind | E1 | preprocessing            | critical | TRUE ||
| R0012 | seurat_qc | blind | E2 | species_annotation            | critical | TRUE ||
| R0012 | seurat_qc | blind | E3 | statistics            | critical | TRUE ||
| R0012 | seurat_qc | blind | E4 | biological_annotation            | critical | FALSE ||
| R0012 | seurat_qc | blind | E5 | workflow_logic            | critical | TRUE ||
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- ||
| R0013 | claims | blind | E1 | small_n_enrichment_overreach | | TRUE| 2 |
| R0013 | claims | blind | E2 | unsupervised_result_given_unearned_label | | TRUE| 2 |
| R0013 | claims | blind | E3 | gene_set_overrepresentation_misread_as_universal | | TRUE| 1|
| R0013 | claims | blind | E4 | statistical_significance_conflated_with_expression_magnitude | | TRUE| 1 |
| R0013 | claims | blind | E5 | batch_effect_misread_as_biology | | TRUE| 2 |
| R0013 | claims | blind | E6 | correlation_causation | | TRUE|2 |
| R0013 | claims | blind | E7 | multiple_testing_ignored | | TRUE|1 |
| ------ | ------- | --------- | -------- | --------------------            | -------- | ----- |


