# QUEST-HPV v3.1.3

[English](README.md) | [简体中文](README.zh-CN.md)

**QUEST-HPV** is a patient-disjoint, clone-first deep-learning workflow for identifying HPV-specific tumor-infiltrating T-cell clonotypes from single-cell RNA expression, paired CDR3α/CDR3β sequences, and V(D)J annotations.

> Research-use software. It is not a clinical diagnostic device and should not be used as the sole basis for patient-care decisions.

## What is included

The repository consolidates the validated v3.1 training workflow and the complete unknown-cohort deployment workflow into one installable package. It retains:

- deterministic HPV biological eligibility gating;
- patient-specific paired CDR3α+CDR3β biological clonotypes;
- learned RNA cell attention with optional curated T-cell program scores;
- a cross-chain paired-TCR Transformer and V(D)J embeddings;
- rare-positive ranking objectives, online hard-negative mining, small-clone-focused training, EMA checkpoints, and multi-seed ensembles;
- Validation-only checkpoint selection, calibration, and thresholding on all pre-gate conditional scores;
- independent-test evaluation, patient-cluster bootstrap, interpretation, split-search, mixed-patient CV planning, and retrained architecture ablations;
- inference for unlabeled/unknown cohorts without re-fitting preprocessing, vocabularies, calibration, or thresholds.

## v3.1.3 HPV-negative control learning and cell-level deployment output

v3.1.3 fixes the supervised-domain gap in biological-gate mode:

- `specificity_train_hpv_positive_only=false` includes HPV-negative tumour TIL clonotypes as explicit `nonNeo` negatives in the main classifier, modality auxiliary losses, ranking objectives, and sampler;
- online hard-negative mining searches all supervised negatives using the **pre-gate conditional specificity score**, so high-scoring HPV-negative controls are actively revisited;
- `selection_domain="all_conditional"` uses all Validation clonotypes and their learned pre-gate probabilities for checkpoint selection, temperature scaling, and threshold selection;
- the biological gate remains a deployment rule: known HPV-negative cases have final `probability=0`, but their conditional logits genuinely affect fitting and model selection;
- `patient_adversary_hpv_positive_only=false` extends patient-adversarial regularization to every sampled patient to reduce cohort/patient shortcut learning;
- training histories now report HPV-negative hard-negative counts, HPV-negative conditional-score summaries, and both all-conditional and all-gated AP.

Unknown-cohort delivery is now cell-level. The model still aggregates cells into patient-specific paired-TCR bags internally, but it writes only per-cell `cell_id`/`Label` and detailed cell probabilities; clone-level prediction CSVs are no longer produced.

## v3.1.2 ablation fix

In v3.1.1, post-hoc perturbations and retrained architecture ablations were separate code paths that were not fully connected to the default workflow: `ablation.enabled` was not consumed by the primary Python API, the recommended config and notebook disabled the retrained suite, and post-hoc files were nested under `interpretability/` and disappeared when `--skip-interpretability` was used. A completed primary run therefore did not imply that architecture ablations had run.

v3.1.2 connects these paths and adds automatic primary-run reuse, a saved-bundle `run-posthoc-ablations` command, top-level/subworkflow manifests, reference-minus-variant tables, patient-cluster paired bootstrap intervals, and robust calibration histograms for tiny or constant-score datasets.

## Repository layout

```text
QUEST-HPV-v3.1.3/
├── src/quest_hpv/
│   ├── config.py            # dataclass configuration and biological/model constants
│   ├── data.py              # validation, clone construction, preprocessing and sampling
│   ├── models.py            # RNA, paired-TCR, V(D)J and fusion networks
│   ├── losses.py            # focal, ranking, contrastive and smooth-AP objectives
│   ├── metrics.py           # calibration, clone/patient metrics and bootstrap helpers
│   ├── training.py          # fitting, EMA, multi-view prediction and ensembles
│   ├── interpretation.py    # attribution, ablation and publication/QC figures
│   ├── workflows.py         # split search, training/evaluation and deployment runners
│   ├── pipeline.py          # backward-compatible public facade
│   ├── cli.py               # unified command-line interface
│   ├── __init__.py
│   └── __main__.py
├── configs/recommended.json
├── notebooks/QUEST_HPV_v3_1_workflow.ipynb
├── tests/test_smoke.py
├── quest_hpv_tcell_pipeline_v3_1.py  # backward-compatible launcher
├── pyproject.toml
└── requirements.txt
```

## Installation

Python 3.9–3.12 is supported; Python 3.10–3.12 is recommended.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e .
```

UMAP is optional; when it is unavailable, latent visualization automatically falls back to PCA.

```bash
python -m pip install -e ".[umap]"
```

Verify the installation:

```bash
quest-hpv --version
quest-hpv --help
```

In Jupyter or when the console script is not on `PATH`, use the equivalent and more robust form:

```bash
python -m quest_hpv --help
```

The legacy invocation remains available:

```bash
python quest_hpv_tcell_pipeline_v3_1.py --help
```

## Input data

### Required columns for labeled Train/Validation/Test matrices

| Column | Meaning |
|---|---|
| `ID` | Unique cell identifier within each partition |
| `Label` | Binary specificity label; defaults are `Neo` and `nonNeo` |
| `patient_id` | Patient identifier |
| `CTgene` | Combined V(D)J annotation containing TRAV/TRAJ/TRBV/TRBD/TRBJ tokens when available |
| `cdr3_aa1` | CDR3α amino-acid sequence |
| `cdr3_aa2` | CDR3β amino-acid sequence |
| `HPV` | HPV status, such as `HPV+` or `HPV-` |
| expression columns | Pre-normalized numeric RNA features |

Optional metadata columns include `cell_names` and `TCR_Frequency`. Clone size and TCR frequency are not used as inference features.

The default automatic expression-column discovery follows the existing matrix convention: expression features lie after `pMT` and before `patient_id`. For a different column order, provide a locked gene list using `--gene-list`.

### Required columns for an unlabeled/unknown cohort

The unknown cohort does not need `Label`. It must contain `ID`, `patient_id`, `CTgene`, `cdr3_aa1`, `cdr3_aa2`, `HPV`, and every expression feature stored in the fitted `preprocessor.joblib`.

If required metadata are stored separately, provide `--metadata` for standalone prediction or `--unknown-metadata` when predicting immediately after training.

## Recommended workflow

### 1. Train, calibrate, evaluate, and run architecture ablations

```bash
quest-hpv train-evaluate \
  --train-data Model_construction_train_data.csv \
  --validation-data Model_construction_validation_data.csv \
  --test-data Model_construction_test_data.csv \
  --metadata Model_construction_patient_metadata.csv \
  --config configs/recommended.json \
  --output-dir QUEST_v3_1_results \
  --run-ablations
```

The workflow audits patient and exact paired-CDR3 overlap, canonicalizes clone labels, fits preprocessing on Train only, trains with HPV-negative TILs as supervised negatives, selects checkpoints/calibration/thresholds from all Validation conditional scores, evaluates the locked Independent Test, runs post-hoc perturbations, and retrains the configured architecture variants. The completed primary run is reused as the reference.

`configs/recommended.json` sets `ablation.enabled=true`; the explicit CLI flag makes the experimental intent visible in logs and scripts. Use `--skip-ablations` to override the config.

For a quick debugging run:

```bash
quest-hpv train-evaluate ... \
  --single-seed \
  --skip-interpretability \
  --skip-ablations
```

### 2. Primary model already exists: run only post-hoc ablations

```bash
quest-hpv run-posthoc-ablations \
  --test-data Model_construction_test_data.csv \
  --metadata Model_construction_patient_metadata.csv \
  --primary-results-dir QUEST_v3_1_results
```

This reads the saved ensemble, locked preprocessor, calibration, threshold, and saved run config; it does not retrain the primary model. It writes modality and individual gene-program perturbation tables/figures to `QUEST_v3_1_results/interpretability_posthoc_rerun/` by default.

### 3. Predict an unknown cohort with the locked ensemble

```bash
quest-hpv predict \
  --input-data unknown_cells.csv \
  --model-bundle QUEST_v3_1_results/quest_hpv_v3_1_ensemble.pt \
  --preprocessor QUEST_v3_1_results/preprocessor.joblib \
  --output-dir QUEST_v3_1_unknown_predictions \
  --eval-bag-views 4 \
  --max-cells-per-clone 128
```

`predict-unknown` is an alias of `predict`.

Prediction semantics (one row per cell):

- `cell_id`: copied from the configured cell-ID column (default `ID`);
- `Label`: string call, defaulting to `Neo` / `nonNeo`;
- `predicted_label_numeric`: the corresponding 1 / 0 call;
- `conditional_probability`: learned specificity probability before the deterministic HPV gate;
- `probability`: final HPV-specific probability after the gate;
- `seed_logit_sd`: disagreement across ensemble seeds;
- `view_logit_sd`: sensitivity to deterministic multi-view sampling of expanded clone bags;
- `prediction_uncertainty`: `seed_logit_sd + view_logit_sd`, provided as a QC indicator.

Cells sharing a patient-specific paired CDR3α/CDR3β receptor inherit the same internal clone-bag prediction; there is no separate per-cell reclassification stage.

### 4. Train and predict an unknown cohort in one command

```bash
quest-hpv train-evaluate \
  --train-data Model_construction_train_data.csv \
  --validation-data Model_construction_validation_data.csv \
  --test-data Model_construction_test_data.csv \
  --metadata Model_construction_patient_metadata.csv \
  --config configs/recommended.json \
  --output-dir QUEST_v3_1_results \
  --unknown-data unknown_cells.csv \
  --run-ablations
```

Unknown-cohort outputs are written to `QUEST_v3_1_results/unknown_prediction/`.

## Additional commands

### Metadata audit

```bash
quest-hpv audit-metadata \
  --metadata Model_construction_patient_metadata.csv \
  --config configs/recommended.json \
  --output-dir QUEST_v3_1_metadata_audit
```

### Development-only patient split search

```bash
quest-hpv optimize-development-split \
  --train-data Model_construction_train_data.csv \
  --validation-data Model_construction_validation_data.csv \
  --test-data Model_construction_test_data.csv \
  --metadata Model_construction_patient_metadata.csv \
  --config configs/recommended.json \
  --output-dir QUEST_v3_1_split_search
```

This is a label-aware development procedure. It may reassign only Train and Validation patients; the Independent Test remains locked and must not be re-selected by performance.

### Mixed HPV-positive patient-held-out CV plan

```bash
quest-hpv make-mixed-patient-cv \
  --metadata Model_construction_patient_metadata.csv \
  --config configs/recommended.json \
  --output mixed_patient_cv_plan_v3_1.csv
```

### Retrained architecture ablations with primary-run reuse

```bash
quest-hpv run-ablations \
  --train-data Model_construction_train_data.csv \
  --validation-data Model_construction_validation_data.csv \
  --test-data Model_construction_test_data.csv \
  --metadata Model_construction_patient_metadata.csv \
  --config configs/recommended.json \
  --primary-results-dir QUEST_v3_1_results \
  --output-dir QUEST_v3_1_ablation
```

`--primary-results-dir` reuses the completed primary `run_summary.json` and clone predictions. All non-reference RNA/TCR/VDJ variants are still trained from scratch on identical patient partitions. Without this option, the reference is retrained as well. This is distinct from post-hoc feature perturbation.

### Write the built-in default configuration

```bash
quest-hpv write-default-config --output quest_v3_1_default_config.json
```

## Principal outputs

Training/evaluation outputs include:

- `quest_hpv_v3_1_ensemble.pt` and `preprocessor.joblib`;
- clone predictions for Train, Validation, and Independent Test;
- final-gated all-patient, ungated all-conditional, and HPV-positive-domain metrics;
- patient-level metrics and patient-cluster bootstrap confidence intervals;
- overlap, label, gene-program, and V(D)J out-of-vocabulary audits;
- PR/ROC, calibration, patient, clone-size, uncertainty, and interpretation figures;
- RNA attention, integrated gradients, TCR residue occlusion, latent embedding, and hard-error tables when enabled.

Ablation outputs are indexed explicitly:

```text
QUEST_v3_1_results/
├── ablation_manifest.json
├── interpretability/
│   ├── posthoc_ablation_manifest.json
│   ├── modality_ablation.csv
│   ├── gene_program_ablation.csv
│   ├── Figure_modality_ablation.{png,pdf,svg}
│   └── Figure_gene_program_ablation.{png,pdf,svg}
└── architecture_ablation/
    ├── architecture_ablation_manifest.json
    ├── architecture_ablation_metrics.csv
    ├── architecture_ablation_metric_deltas.csv
    ├── architecture_ablation_paper_table.csv
    ├── architecture_ablation_paired_bootstrap_summary.csv
    └── Figure_architecture_ablation_*.{png,pdf,svg}
```

Post-hoc ablation asks what happens when an input is removed from the locked model at inference. Architecture ablation asks what happens after removing a component and retraining. They should be reported separately.

Unknown-cohort outputs include:

- `unknown_cell_id_Label.csv`, containing only `cell_id` and `Label`;
- `unknown_cell_predictions.csv`, with one row per cell and its label, probabilities, gate, and uncertainty;
- `unknown_RNA_cell_attention.csv` and top-attention cells;
- `unknown_input_audit.csv` and `unknown_VDJ_OOV.csv`;
- `unknown_prediction_summary.json`;
- `Figure_unknown_prediction_overview.{png,pdf,svg}` at cell level.

No `unknown_clone_predictions.csv` or other clone-level prediction result is written.

## Configuration notes

The supplied `configs/recommended.json` uses three ensemble seeds, four evaluation bag views, balanced-accuracy threshold selection, all-labelled-clone specificity training, all-conditional Validation selection, all-patient adversarial regularization, small-clone losses, interpretation, bootstrap iterations, and publication figure formats.

Expression values are assumed to be normalized before input. Optional train-fitted standardization is feature scaling only and is not batch correction.

## Tests

```bash
python -m pip install -e ".[test]"
pytest
```

The test suite checks standalone and saved-run configuration loading, clone construction, HPV gating, preprocessor/model forward execution, CLI ablation switches, constant-score plotting, saved-bundle inference on a synthetic unlabeled cohort, and a tiny end-to-end Train → Validation calibration → Independent Test → integrated architecture ablation → saved-bundle post-hoc ablation run.

## Reproducibility and release checklist

Before a public release:

1. confirm the repository does not contain patient-level raw data or trained weights that cannot be shared;
2. choose and add the project license—the original files did not specify one, so this package intentionally does not guess it;
3. add author/institution information and the preferred manuscript citation;
4. record the exact environment and data version used for the reported results.
