# Gujarat Wastewater METAGENE-1 Study
## Geographic Transferability of Metagenomic Foundation Models: Evidence from Indian Urban Wastewater

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![HuggingFace Model](https://img.shields.io/badge/🤗_Model-India--METAGENE--v6-blue)](https://huggingface.co/saurabhthakar3/india-metagene-1)
[![HuggingFace Data](https://img.shields.io/badge/🤗_Data-Gujarat_WW_Shotgun-green)](https://huggingface.co/datasets/saurabhthakar3/gujarat-wastewater-shotgun)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/)

---

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Scientific Background](#2-scientific-background)
3. [Study Design](#3-study-design)
4. [What We Did — Step by Step](#4-what-we-did--step-by-step)
5. [Key Findings](#5-key-findings)
6. [The METAGENE-1 1.24 Reference — A Critical Clarification](#6-the-metagene-1-124-reference--a-critical-clarification)
7. [Fine-Tuning Results](#7-fine-tuning-results)
8. [Amplicon 16S Study](#8-amplicon-16s-study)
9. [Repository Structure](#9-repository-structure)
10. [How to Use This Code](#10-how-to-use-this-code)
11. [Data Availability](#11-data-availability)
12. [Citation](#12-citation)

---

## 1. Project Overview

This repository contains all code, analysis scripts, notebooks, and reports for a comprehensive study applying **METAGENE-1** — a 7-billion-parameter metagenomic foundation model trained exclusively on US wastewater — to **Indian urban wastewater** from Gujarat.

The study asks a fundamental question: **Can a foundation model trained on US data understand Indian wastewater?**

The answer, in short, is: partially — and this repository documents exactly how, why, and what needs to be done to improve it.

### This Repository Contains Two Related Studies

| Study | Approach | Location |
|-------|----------|----------|
| **Shotgun metagenomics** | Whole-genome shotgun sequencing → METAGENE-1 anomaly scoring + LoRA fine-tuning | `scripts/` |
| **Amplicon 16S** | 16S rRNA amplicon → ASV dark matter → METAGENE-1 embeddings | `amplicon_16s/` |

Both studies use the **same Gujarat wastewater samples** from 4 cities across 4 seasons.

---

## 2. Scientific Background

### What is METAGENE-1?

METAGENE-1 (Liu et al., 2025; arXiv:2501.02045) is a **7-billion-parameter causal language model** built on the LLaMA-2 architecture, trained on 1.5 trillion base pairs of metagenomic sequences from US municipal wastewater treatment plants. It uses a DNA-specific byte-pair encoding (BPE) tokeniser with vocabulary size 1,024 (~3.9 base pairs per token).

Like a language model that learns to predict the next word in a sentence, METAGENE-1 learns to predict the next nucleotide in a DNA sequence. Sequences that are similar to its training data are predicted well (low cross-entropy loss). Sequences that are novel or different are predicted poorly (high cross-entropy loss).

### What is Anomaly Scoring?

When a new DNA sequence is fed to METAGENE-1, the model computes a **cross-entropy (CE) loss** — a measure of how surprised the model is by the sequence. This gives us a per-read anomaly score:

```
Low CE loss  (~1-2)  → Model recognises this sequence → Similar to US wastewater
High CE loss (~5-6)  → Model cannot predict this     → Novel, different from US wastewater
```

By averaging CE loss across thousands of reads from a sample, we get a sample-level anomaly score that tells us how "foreign" that sample is to the US-trained model.

### Why Does This Matter for India?

India has 1.4 billion people, distinct monsoon-driven seasonality, different disease burden profiles, unique antibiotic usage patterns, and microbial communities shaped by completely different ecological contexts than the United States. If METAGENE-1 cannot understand Indian wastewater sequences, it cannot be used for surveillance in India without adaptation.

Quantifying this geographic gap — and demonstrating how to close it through fine-tuning — is the core contribution of this study.

### The Dark Matter Problem

Traditional taxonomy classifiers like Kraken2 work by matching sequences to a reference database. In Gujarat wastewater, approximately **89% of reads cannot be classified** — they have no match in any reference database. This "dark matter" is the most interesting fraction biologically but is completely invisible to conventional tools.

METAGENE-1 can score dark matter sequences without needing database matches, making it uniquely suited for surveillance in ecologically diverse settings.

---

## 3. Study Design

### Samples

| Property | Value |
|----------|-------|
| Total samples | 95 |
| Cities | Ahmedabad (8.1M pop.), Gandhinagar (0.3M), Rajkot (1.7M), Vadodara (2.1M) |
| Seasons | Summer (Mar–May), Monsoon (Jun–Sep), PreWinter (Oct–Nov), Winter (Dec–Feb) |
| Sampling period | 2022–2024 |
| Sequencing | Illumina 2×150bp shotgun |
| Total data | ~190 million reads, ~0.045 Tbp |

### Sampling Strategy

Multiple WWTP intake sites per city were sampled and **pooled by city × season** before library preparation. This means each of the 95 metagenomics samples represents the integrated microbial community of an entire city across an entire season — maximising the representativeness of each sample.

### Held-Out Validation Set

One sample per city × season combination was reserved as a held-out validation set (n=16 samples, seed=42). This set was used consistently across all fine-tuning experiments for fair comparison.

### US External Control

Three Rothman et al. (2020) Southern California WWTP samples (BioProject PRJNA649747) were used as a within-study US baseline, scored with identical pipeline parameters in the same session.

---

## 4. What We Did — Step by Step

### Phase 1: Baseline Anomaly Scoring

We implemented the METAGENE-1 anomaly scoring pipeline from Liu et al. (2025) Section 5.5:
- Subsample 5,000 forward reads per sample (seed=42)
- Tokenise with METAGENE-1 BPE tokeniser
- Compute per-read CE loss (manual masked computation — NOT the built-in model loss which includes padding tokens)
- Average across all reads for a sample-level score

**Why manual CE computation?** The model's built-in `.loss` output averages over all tokens including padding, inflating scores to ~7.38. Manual computation correctly excludes padding tokens and gives ~4.59 for US wastewater — the scientifically correct value.

### Phase 2: Pipeline Validation

Before any geographic comparison, we validated the pipeline against published references:
- Scored NA12878 human WGS reads (SRR1518158) → CE loss 5.490 ± 0.454 vs Liu et al. reference 5.22 (Δ=0.270 ✓)
- Confirmed BPE tokeniser functioning correctly: vocab size 1,024, ~3.9 bp/token
- Confirmed float16 and float32 give identical results (difference: 0.000005)

### Phase 3: Systematic Replication (11 Experiments)

To prove the geographic gap is real and not a pipeline artefact, we tested 11 different pipeline configurations on both US and Gujarat data:

| Experiment | What was tested | Key finding |
|------------|----------------|-------------|
| 1 | Baseline | US 4.57, Gujarat 4.87, gap +0.29 |
| 2 | Read length filter 100-300bp | US 4.28, gap +0.59 |
| 3 | max_length 64-1024 | No effect — reads are shorter than 64 tokens |
| 4 | Reverse complement | min(fwd,RC) gives gap +0.29 |
| 5 | Loss computation method | Built-in loss = 7.38 (padding artefact); manual = 4.57 |
| 6 | Special tokens (BOS/EOS) | Minor effect; gap persists |
| 7 | NAO data (METAGENE co-author lab) | Scores 5.07–5.10; confirms 1.24 not reproducible |
| 8 | Number of reads (500–5000) | CE loss stable at 5000 reads |
| 9 | Float16 vs float32 | Identical (difference: 0.000005) |
| 10 | RNA sequences (T→U) | CE loss 16.43 — model trained on DNA not RNA |
| 11 | All facilities pooled | US pooled 4.57; gap persists |

**Result: The geographic gap was positive in EVERY configuration tested.** Mean gap across all configs: +0.391 ± 0.139 CE units.

### Phase 4: Statistical Characterisation

Two-way ANOVA on individual sample CE loss scores (n=95):
- **Season effect:** η²=0.491 — season explains 49.1% of variance
- **City effect:** η²=0.355 — city explains 35.5% of variance
- **City × Season interaction:** F(9,72), p<0.0001 — seasonal effects differ across cities

The significant interaction means the seasonal pattern is different in each city — Ahmedabad shows disproportionately high anomaly during Monsoon while other cities show more uniform seasonal variation.

### Phase 5: LoRA Fine-Tuning

We applied Low-Rank Adaptation (LoRA) to adapt METAGENE-1 to Indian wastewater:

**India-METAGENE v4** — rank 8, adapts q_proj and v_proj attention layers
- 8.4M trainable parameters (0.12% of 7B total)
- 500,000 training sequences (0.26% of available data)
- Hardware: NVIDIA T4 (16GB)

**India-METAGENE v6** — rank 32, adapts all four attention projection matrices
- 67.1M trainable parameters (0.96% of 7B total)
- 1,000,000 training sequences (0.53% of available data)
- Hardware: NVIDIA A100 (40GB)

### Phase 6: Metabolomics Integration (Exploratory)

LCMS untargeted metabolomics data was available for 505 samples from the same cities/seasons. We aggregated by city × season and tested correlations with CE loss. Key finding: pharmaceutical biomarkers (cotinine, metformin, caffeine) were negatively correlated with CE loss — higher CE loss (Monsoon) coincides with dilution of human-derived chemicals by rainfall. No metabolites survived FDR correction at n=16 aggregated points. Analysis set aside pending better-powered follow-up study.

---

## 5. Key Findings

### Finding 1: Significant Geographic Distribution Shift

Gujarat wastewater CE loss: **4.865 ± 0.066**
US external control CE loss: **4.593 ± 0.053**
Geographic gap: **0.272 CE units (6.0%)**
Statistical test: t=6.96, p=0.020, Cohen's d=4.92

The gap is **positive in every one of the 11 pipeline configurations tested** — it is not a pipeline artefact.

### Finding 2: The Gap Lives in the Low-CE Tail

The geographic shift is not uniform across the CE loss distribution:

| Percentile | Gujarat | US | Gap |
|------------|---------|-----|-----|
| p10 | 4.181 | ~1.99 | +2.19 |
| p25 | 4.675 | ~3.20 | +1.47 |
| p50 | 5.048 | ~4.93 | +0.12 |
| p75 | 5.338 | ~5.60 | −0.26 |

Gujarat has **4.5× fewer reads below CE loss 2.0** (2.3% vs ~10%). The model recognises far fewer Indian sequences as familiar. The upper tail (high anomaly reads) is nearly identical between countries.

### Finding 3: Strong Seasonal and City Effects

- **Monsoon** samples are most anomalous (mean CE 4.921) — rainfall introduces environmental organisms
- **PreWinter** samples are least anomalous (mean CE 4.804) — more urban-dominated community
- **Ahmedabad** is the most anomalous city, significantly exceeding Rajkot and Vadodara (Mann-Whitney, p<0.05 after Bonferroni)
- **Ahmedabad Monsoon (AAMO_R1)** is the single most anomalous sample (CE 5.026) — consistently so across all experiments

### Finding 4: LoRA Fine-Tuning Works and Scales With Capacity

| Model | Val CE Loss | Improvement | Effect Size | 95% CI |
|-------|------------|-------------|------------|--------|
| METAGENE-1 baseline | 4.8697 | — | — | — |
| India-METAGENE v4 (rank 8) | 4.8046 | +0.065 (1.34%)* | d=0.726 | [0.017, 0.113] |
| India-METAGENE v6 (rank 32) | 4.7890 | +0.081 (1.66%)*** | d=1.030 | [0.039, 0.123] |

*p<0.05, ***p<0.001, paired t-test, n=16

PreWinter samples (VCPW_R1, RPPW_R1) respond best to fine-tuning (+0.141–0.147 improvement). Ahmedabad Monsoon (AAMO_R1) uniquely resists adaptation — worsening slightly under both configurations.

### Finding 5: The Published 1.24 Is Not Reproducible

See Section 6 below for full explanation.

---

## 6. The METAGENE-1 1.24 Reference — A Critical Clarification

Liu et al. (2025) report a mean CE loss of **1.24** for US wastewater. This creates the expectation that external users should obtain values near 1.24. This expectation is incorrect.

**The 1.24 is training-set memorisation loss** — computed on the model's own training sequences during pre-training. It reflects how well the model has memorised its training data, not how well it generalises to new data.

We demonstrate this through three independent lines of evidence:

1. **Rothman et al. (2020) data** — cited by Liu et al. as a training source — scores 4.28–4.71 under our pipeline, not 1.24
2. **NAO data (Jeff Kaufman lab)** — Kaufman is a METAGENE-1 co-author, making this data likely from the actual training set — scores 5.07–5.10, even higher than Rothman
3. **Systematic replication** — no pipeline configuration across 11 experiments brings the score below 4.28 on any public US wastewater dataset

**The correct comparison in this study is:**
- US external (4.593) vs Gujarat external (4.865) — gap of 0.272 CE units (6%)

**Do not compare your CE loss scores to 1.24** — it is not an evaluable benchmark.

---

## 7. Fine-Tuning Results

### Why LoRA?

Full fine-tuning of a 7B parameter model requires ~280GB GPU memory for gradients — impractical on available hardware. LoRA adds small adapter matrices (rank r) to selected weight matrices, reducing trainable parameters to <1% while enabling meaningful adaptation.

### The Scaling Story

```
Baseline  (0 params updated)    → Val loss 4.8697
v4        (8.4M params, 0.12%)  → Val loss 4.8046  (+1.34%)
v6        (67.1M params, 0.96%) → Val loss 4.7890  (+1.66%)
```

More parameters → better adaptation. The trend motivates full fine-tuning (~80–160 GPU-hours on A100) as the next step.

### Important Limitation

v4 and v6 differ in three ways simultaneously: rank (8 vs 32), attention coverage (q/v vs q/v/k/o), and training data volume (500K vs 1M sequences). The improvement from v4 to v6 **cannot be attributed to rank alone**. A controlled ablation is needed.

### Published Model

**India-METAGENE v6 (epoch 4)** is the published checkpoint: `saurabhthakar3/india-metagene-1`
Val loss: 4.7898 | Rank 32 | A100 GPU | 1M training sequences

---

## 8. Amplicon 16S Study

See `amplicon_16s/` folder and its dedicated README.

### Overview

16S rRNA amplicon sequencing of the same Gujarat wastewater samples was conducted and processed through DADA2 → ASVs → SILVA 138.2 taxonomy. Approximately **X%** of ASVs could not be classified ("dark matter"). We applied METAGENE-1 embeddings to characterise these unclassifiable sequences across four phases:

| Phase | What | Scripts |
|-------|------|---------|
| Phase 1 | Identify dark matter ASVs (unclassifiable by SILVA) | `scripts/R/phase1_*.R` |
| Phase 2 | Extract METAGENE-1 hidden state embeddings for all ASVs | `scripts/phase2_extract_embeddings.py` |
| Phase 3 | Compare embedding space across city and season groups | `scripts/R/phase3_comparison.R` |
| Phase 4 | Attempt taxonomy assignment for dark matter using embeddings | `scripts/R/phase4_*.R` |

### Key Caveat on Amplicon Embeddings

METAGENE-1 was trained on shotgun reads (~150bp), not amplicon sequences (~250–450bp targeting specific 16S regions). The embeddings for amplicon ASVs reflect US wastewater-centric representations and should be interpreted cautiously. This is an exploratory analysis.

---

## 9. Repository Structure

```
gujarat-wastewater-shotgun-metagenomics/
│
├── README.md                                    ← This file
├── LICENSE                                      ← MIT License
├── .gitignore
│
├── scripts/                                     ← SHOTGUN STUDY NOTEBOOKS
│   ├── metagene_shotgun_colab_v3.ipynb          ← Main anomaly scoring pipeline
│   ├── metagene_replication_experiments.ipynb   ← 11 configs on US WW data
│   ├── gujarat_replication_experiments.ipynb    ← 11 configs on Gujarat data
│   ├── india_metagene_stats.ipynb               ← Two-way ANOVA, bootstrap CIs
│   ├── india_metagene_finetune_v4.ipynb         ← LoRA rank-8 fine-tuning (T4)
│   ├── india_metagene_finetune_v6.ipynb         ← LoRA rank-32 fine-tuning (A100)
│   ├── metagene_full_finetune_poc.ipynb         ← Full fine-tuning proof of concept
│   ├── us_validation_minimal.ipynb              ← US WW + human genome validation
│   ├── fresh_setup_and_validate.ipynb           ← HuggingFace setup and validation
│   └── multiomics_integration.py               ← Metabolomics integration (local)
│
├── docs/                                        ← REPORTS AND MANUSCRIPTS
│   ├── study_report_part1.docx                  ← Comprehensive report Part 1
│   ├── study_report_part2.docx                  ← Comprehensive report Part 2
│   ├── india_metagene_full_report_v2.docx       ← Full manuscript v2
│   ├── METAGENE1_Complete_Report.docx           ← Complete analysis report
│   └── METHODS.md                               ← Detailed methods description
│
├── environment/                                 ← DEPENDENCIES
│   ├── requirements.txt                         ← Python pip requirements
│   └── conda_env.yml                            ← Conda environment file
│
├── results/                                     ← RESULTS PLACEHOLDER
│   └── README.md                                ← Points to HuggingFace for results
│
└── amplicon_16s/                                ← AMPLICON 16S STUDY
    ├── README.md                                ← Amplicon study description
    ├── data/                                    ← Metadata and abundance tables
    │   ├── sample_metadata.xlsx
    │   ├── abundance_table.xlsx
    │   └── taxonomy_export.xlsx
    │   (FASTA files on HuggingFace — too large for GitHub)
    ├── scripts/
    │   ├── R/                                   ← R analysis scripts
    │   │   ├── phase1_export_sequences.R
    │   │   └── r_packages.R
    │   ├── phase1_feasibility_check.py
    │   └── phase1_recategorize_dark_matter.py
    ├── docs/
    │   ├── METAGENE1_Complete_Report.docx
    │   └── METAGENE1_Results_Section.docx
    └── results/
        ├── phase1/
        ├── phase2/
        ├── phase3/
        └── phase4/
```

---

## 10. How to Use This Code

### Setup

```bash
# Clone repository
git clone https://github.com/saurabhgit199/gujarat-wastewater-shotgun-metagenomics.git
cd gujarat-wastewater-shotgun-metagenomics

# Create conda environment
conda env create -f environment/conda_env.yml
conda activate metagene-shotgun
```

### Score Your Own Wastewater Sequences

```python
import torch
import numpy as np
from transformers import AutoTokenizer, AutoModelForCausalLM

# Load model
tokenizer = AutoTokenizer.from_pretrained('metagene-ai/METAGENE-1')
model = AutoModelForCausalLM.from_pretrained(
    'metagene-ai/METAGENE-1',
    torch_dtype=torch.float16,
    device_map='cuda'
)
model.eval()

@torch.no_grad()
def score_sequences(seqs, batch_size=8, max_length=512):
    """
    Compute mean CE loss per read.
    Replicates Liu et al. 2025 Section 5.5 exactly.
    
    IMPORTANT: Uses manual masked CE computation.
    Do NOT use model.loss — it includes padding tokens
    and gives inflated scores (~7.38 vs correct ~4.59).
    """
    losses = []
    for i in range(0, len(seqs), batch_size):
        batch  = seqs[i:i+batch_size]
        inputs = tokenizer(
            batch, return_tensors='pt', padding=True,
            truncation=True, max_length=max_length,
            add_special_tokens=False
        ).to('cuda')
        with torch.amp.autocast('cuda', dtype=torch.float16):
            out = model(**inputs, labels=inputs['input_ids'])
        logits  = out.logits[..., :-1, :].contiguous().float()
        targets = inputs['input_ids'][..., 1:].contiguous()
        mask    = inputs['attention_mask'][..., 1:].contiguous().float()
        tok_loss = torch.nn.CrossEntropyLoss(reduction='none')(
            logits.view(-1, logits.size(-1)), targets.view(-1)
        ).view(targets.size())
        per_read = (tok_loss * mask).sum(1) / mask.sum(1).clamp(min=1)
        losses.extend(per_read.cpu().float().tolist())
    return np.array(losses)

# Score your sequences
seqs = ['ATCGATCG...', 'GCTAGCTA...']  # list of DNA strings
scores = score_sequences(seqs)
print(f'Mean CE loss: {scores.mean():.4f}')
print(f'Anomalous reads (>3.0): {(scores > 3.0).mean()*100:.1f}%')
```

### Use India-METAGENE v6 (Fine-Tuned)

```python
from peft import PeftModel

# Load base model then apply LoRA adapters
base_model = AutoModelForCausalLM.from_pretrained(
    'metagene-ai/METAGENE-1',
    torch_dtype=torch.float16,
    device_map='cuda'
)
model = PeftModel.from_pretrained(
    base_model,
    'saurabhthakar3/india-metagene-1',
    subfolder='checkpoint_lora_v6_epoch4_stepfinal'
)
model.eval()
# Now score Indian wastewater sequences with adapted model
```

---

## 11. Data Availability

| Data | Location |
|------|----------|
| Raw FASTA (shotgun) | HuggingFace: saurabhthakar3/gujarat-wastewater-shotgun *(private pending review)* |
| India-METAGENE v4 + v6 | HuggingFace: saurabhthakar3/india-metagene-1 |
| Base METAGENE-1 | HuggingFace: metagene-ai/METAGENE-1 |
| US validation data | SRA: PRJNA649747 (Rothman et al. 2020) |
| Human genome validation | SRA: SRR1518158 (NA12878 WGS) |
| Amplicon FASTA files | HuggingFace: saurabhthakar3/gujarat-wastewater-shotgun |

---

## 12. Citation

If you use this code, data, or findings, please cite:

```bibtex
@article{thakar2025geographic,
  title={Geographic Transferability of Metagenomic Foundation Models:
         Evidence from Indian Urban Wastewater},
  author={Thakar, Saurabh and others},
  journal={[Target Journal]},
  year={2025}
}
```

Also cite METAGENE-1:

```bibtex
@article{liu2025metagene,
  title={METAGENE-1: Metagenomic Foundation Model for Pandemic Monitoring},
  author={Liu, Ollie and others},
  journal={arXiv preprint arXiv:2501.02045},
  year={2025}
}
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

Saurabh Thakar
[Institution to be added]
saurabhthakar3@gmail.com
