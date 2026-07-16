# Evaluating METAGENE-1 for 16S rRNA V4 Amplicon Analysis of Indian Wastewater

> **A rigorous, phased evaluation of a 7B-parameter metagenomic foundation model on targeted amplicon data — with honest negative results.**

[![Paper](https://img.shields.io/badge/METAGENE--1-arXiv:2501.02045-red)](https://arxiv.org/abs/2501.02045)
[![Model](https://img.shields.io/badge/Model-HuggingFace-yellow)](https://huggingface.co/metagene-ai/metagene-1)

## Summary

This project systematically evaluates whether [METAGENE-1](https://arxiv.org/abs/2501.02045), a 7-billion parameter autoregressive transformer pretrained on >1.5 trillion base pairs of US wastewater shotgun metagenomics, can add value to standard 16S rRNA V4 amplicon analysis of Indian wastewater.

**Key finding:** METAGENE-1's primary utility for 16S amplicon data is **sequence quality control** (identifying 262K non-target sequences) and **phylum-level dark matter classification** (2.2% of reads newly classified), not community-level ecological analysis where standard methods (Bray-Curtis, SILVA) outperform embeddings.

## Study Design

| Parameter | Value |
|-----------|-------|
| Cities | Ahmedabad, Gandhinagar, Rajkot, Vadodara (Gujarat, India) |
| Sites | 24 (6 per city) |
| Seasons | Pre-Winter, Winter, Summer, Monsoon |
| Samples | 505 |
| Target | 16S rRNA V4 (515F/806R, ~253 bp) |
| Total ASVs | 371,775 |
| Total reads | 107,962,810 |
| Taxonomy | SILVA 138.2 |

## Evaluation Pipeline

The evaluation followed a phased approach, with each phase serving as a go/no-go checkpoint:

```
Phase 1: Feasibility ──→ Phase 2: Embeddings ──→ Phase 3: Ecological ──→ Phase 4: Taxonomy
  (Is data in-dist?)     (Extract repr.)        (Compare methods)      (Classify dark matter)
```

### Phase 1: Feasibility Assessment

**Question:** Can METAGENE-1 process 16S V4 amplicon sequences from Indian wastewater?

| Test | Result | Conclusion |
|------|--------|------------|
| Classified ASVs (n=200) | Mean loss 1.465, 94% in-distribution | ✅ Recognized by model |
| Dark matter ASVs (n=200) | Mean loss 4.773, 11.5% in-distribution | ❌ Initially appeared OOD |
| Primer-controlled validation | Trimmed DM: loss 1.729, 94% in-dist | ✅ OOD driven by non-16S sequences |
| Dark matter profiling (n=50K) | 89.4% lack primer binding sites | 262K sequences are non-target amplification |
| 6-category re-categorization | 22K strict-match DM embeddable | Only 7.9% of dark matter is real 16S |

**Key insight:** The primer-controlled 4-group experiment overturned the initial dark matter OOD result, demonstrating that the high loss was caused by non-16S sequences lacking primer binding sites, not by biological novelty.

### Phase 2: Embedding Extraction

- Extracted 4,096-dim mean-pooled last hidden state embeddings for 8,243 ASVs (95% of reads)
- Processing: 7.3 hours on CPU (float32), batch size 14
- Sanity check: no zero embeddings, norm 14.69 ± 1.73

### Phase 3: Ecological Analysis Comparison

**Question:** Do embedding-based distances improve ordination and diversity analysis?

| Metric | Bray-Curtis | Embedding | Winner |
|--------|------------|-----------|--------|
| PERMANOVA R² (City) | 0.0215 | 0.0114 | Bray-Curtis |
| PERMANOVA R² (Season) | 0.0896 | 0.0297 | Bray-Curtis |
| PERMANOVA R² (Site) | 0.0780 | 0.0408 | Bray-Curtis |
| RF City accuracy | 58.4% | 50.3% | Bray-Curtis |
| RF Season accuracy | 81.0% | 62.8% | Bray-Curtis |
| Mantel (emb vs taxonomy) | — | r=0.309, p=0.001 | Embeddings recover taxonomy |

**Conclusion:** Bray-Curtis outperforms embedding-based distances on all community-level metrics. Embeddings capture individual sequence-level taxonomic structure (Mantel r=0.31) but lose compositional information when averaged to sample level.

### Phase 4: Taxonomy Assignment

**Validation (genus-holdout cross-validation):**

| Level | Accuracy | Reliability |
|-------|----------|-------------|
| Phylum | 89.0% | Reliable |
| Class | 86.5% | Reliable |
| Order | 70.9% | Putative |
| Family | 47.2% | Putative |
| Genus | Not validated | Speculative |

**Dark matter assignment:** 22,173 ASVs assigned phylum-level taxonomy via K=5 nearest neighbor in embedding space. Top predicted phyla: Firmicutes (22.6%), Proteobacteria (18.8%), Bacteroidota (12.4%).

**SILVA vs METAGENE-1 agreement on classified ASVs:**

| Level | Agreement |
|-------|-----------|
| Kingdom | 99.2% |
| Phylum | 74.0% |
| Class | 71.7% |
| Family | 96.9%* |
| Genus | 94.3%* |

*\*Inflated by near-duplicate sequences (83% of ASV pairs have cosine distance < 0.05)*

### Final Read Accounting

| Category | Reads | Percentage |
|----------|-------|------------|
| SILVA classified to Genus | 90,638,544 | 84.0% |
| SILVA partial (Phylum, no Genus) | 14,774,005 | 13.7% |
| METAGENE-1 dark matter (Phylum) | 2,385,956 | 2.2% |
| Unclassified by both methods | 164,305 | 0.2% |

## What Worked vs What Didn't

### ✅ Worked Well
- **Sequence QC:** Loss-based OOD detection independently identified 262K non-target sequences, confirmed by primer analysis
- **Cross-domain transfer:** Model trained on US shotgun metagenomics recognizes Indian 16S V4 sequences
- **Phylum-level classification:** 89% validated accuracy for dark matter ASVs
- **Community composition consistency:** SILVA and METAGENE-1 agree on overall composition patterns across cities and seasons

### ❌ Didn't Work
- **Ecological ordination:** Bray-Curtis beats embeddings on every metric
- **Fine-grained taxonomy:** Family-level accuracy only 47% in holdout validation
- **Genus prediction:** Not reliable due to high sequence conservation in V4 region

### ⚠️ Important Caveats
- Standard cross-validation (98% genus accuracy) was misleading — genus-holdout revealed true accuracy of 47% at family level
- 83% of ASV pairs have cosine distance < 0.05, limiting discriminative power
- The 16S V4 region (~253 bp of a single conserved gene) provides insufficient sequence diversity for foundation model embeddings to outperform purpose-built tools

## Repository Structure

```
├── README.md
├── scripts/
│   ├── phase1_feasibility_check.py       # OOD loss computation
│   ├── phase1_dark_matter_validation.py  # 4-group primer-controlled test
│   ├── profile_no_primer_dark_matter.py  # Sequence profiling
│   ├── phase1_recategorize_dark_matter.py # 6-category fuzzy matching
│   ├── phase2_extract_embeddings.py      # Embedding extraction with checkpointing
│   ├── phase3_comparison.py              # PERMANOVA, classification, ordination
│   ├── phase4_validation.py              # Standard CV (baseline)
│   ├── phase4_validation_v2.py           # Genus/family holdout (rigorous)
│   ├── phase4_taxonomy_assignment.py     # Dark matter classification
│   ├── compare_three_phyloseq.R          # SILVA vs METAGENE-1 comparison
│   └── compare_by_groups.R              # City/season composition plots
├── results/
│   ├── phase1/                           # Loss distributions, validation plots
│   ├── phase2/                           # Embedding files
│   ├── phase3/                           # Ordination plots, PERMANOVA results
│   └── phase4/                           # Taxonomy assignments, comparison plots
└── docs/
    ├── METAGENE1_Results_Section.docx    # Publication-ready results
    └── METAGENE1_Complete_Report.docx    # Detailed technical report
```

## Reproducibility

### Requirements
```bash
# Python
pip install torch transformers numpy scipy pandas matplotlib seaborn scikit-learn

# R
install.packages(c("phyloseq", "ggplot2", "gridExtra"))
```

### Hardware Used
- CPU: 16 cores
- RAM: 251 GB
- GPU: NVIDIA T1000 8GB (insufficient for 7B model; all inference on CPU)
- Total compute time: ~40 hours (Phase 2 embeddings: 7h, Phase 4 dark matter: 20h, misc: 13h)

### Data Requirements
- Phyloseq object with SILVA 138.2 taxonomy
- FASTA files for classified and dark matter ASVs
- METAGENE-1 model weights from [HuggingFace](https://huggingface.co/metagene-ai/metagene-1)

## Citation

If you use this evaluation framework:
```
METAGENE-1 paper:
Liu et al. (2025). METAGENE-1: Metagenomic Foundation Model for Pandemic Monitoring.
arXiv:2501.02045
```

## Lessons Learned

1. **Always validate at the right granularity.** Standard cross-validation gave 98% genus accuracy; genus-holdout revealed the true accuracy was 47% at family level. Near-duplicate sequences inflate naive metrics.

2. **Negative results are results.** Showing that embedding-based ordination doesn't outperform Bray-Curtis for 16S data is as informative as showing it does — it defines the boundary of the tool's applicability.

3. **Foundation models aren't universal.** METAGENE-1 excels at what it was trained for (shotgun metagenomics QC, anomaly detection). Applying it to a different data type (targeted amplicons) reveals its limitations.

4. **Primer analysis is undervalued.** The most impactful finding wasn't from the neural network — it was discovering that 89% of dark matter ASVs lack primer binding sites, indicating they're non-target amplification products.

5. **Phased evaluation saves time.** Each phase served as a checkpoint. If Phase 1 had shown all sequences were OOD, we would have stopped — saving 30+ hours of compute.

## License

Code in this repository is available under MIT License. The METAGENE-1 model has its own license — see the [HuggingFace page](https://huggingface.co/metagene-ai/metagene-1).
