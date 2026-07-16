# Detailed Methodology

## Model Details

**METAGENE-1** (Liu et al., 2025)
- Architecture: Llama-2-7B style decoder-only transformer
- Parameters: 7 billion
- Training data: >1.5 trillion base pairs from US wastewater shotgun metagenomics
- Tokenizer: BPE with 1,024 vocabulary
- Context length: 512 tokens
- Embedding extraction: Mean pooling of last hidden state (4,096 dimensions)

## Phase 1: In-Distribution Assessment

### Methodology
Following METAGENE-1 paper Section 5.5, we computed length-normalized cross-entropy loss for each ASV sequence. The model processes each sequence autoregressively, and sequences that receive low loss are considered "in-distribution" — the model has learned meaningful representations for similar sequences during pretraining.

### Key Parameters
- OOD threshold: 3.0 (from paper's anomaly detection study)
- Paper reference values: metagenomics = 1.24, human genome = 5.22, random = 5.83
- Sequence preparation: prepend '_' character per paper Appendix B
- Primer trimming: IUPAC-aware regex for 515F/806R

### Primer Trimming Details
- Forward: 515F = `GTGYCAGCMGCCGCGGTAA` (19 bp)
- Reverse complement of 806R = `ATTAGAWACCCBNGTAGTCC` (20 bp)
- IUPAC degenerate bases resolved via regex: Y=[CT], M=[AC], W=[AT], V=[ACG], B=[CGT], N=[ACGT]

### Controlled Validation Design
The initial dark matter OOD result was challenged by a 4-group experiment:
1. Classified ASVs with both primers found and trimmed
2. Dark matter ASVs with both primers found and trimmed
3. Dark matter ASVs where neither primer was found
4. Random DNA sequences (length-matched)

This design isolates the effect of primer status from biological sequence content.

### Fuzzy Primer Re-categorization
For sequences where strict primer matching failed, we applied:
- Hamming distance matching allowing ≤3 mismatches in primer binding region
- Reverse orientation detection (sequences starting with 806R forward)
- Reverse complementation followed by 515F matching

## Phase 2: Embedding Extraction

### Method
Mean-pooled last hidden state embeddings, following Gene-MTEB methodology (paper Section 5.3).

```python
outputs = model(input_ids, attention_mask, output_hidden_states=True)
last_hidden = outputs.hidden_states[-1]  # (batch, seq_len, 4096)
mask_expanded = attention_mask.unsqueeze(-1).float()
mean_pooled = (last_hidden * mask_expanded).sum(dim=1) / mask_expanded.sum(dim=1)
```

### ASV Selection
Top 8,243 ASVs by total read abundance (covering 95% of all reads), determined by cumulative abundance ranking across all 505 samples.

## Phase 3: Ecological Comparison

### Distance Matrices
1. **Bray-Curtis:** Standard dissimilarity on relative abundance profiles
2. **Embedding-based:** Abundance-weighted average of ASV embeddings per sample, then pairwise cosine distance

### PERMANOVA
Simplified implementation with 999 permutations. R² computed as between-group sum of squares divided by total sum of squares on the distance matrix.

### Classification
Random Forest (100 trees, 5-fold stratified CV) trained on three feature sets:
- Raw relative abundance (Bray-Curtis features)
- Sample embedding vectors (4,096 dimensions)
- Genus-level aggregated abundance

## Phase 4: Taxonomy Assignment

### Validation Design
**Standard CV (baseline):** 5-fold stratified on phylum, random ASV holdout. This was found to be inflated due to near-duplicate sequences.

**Genus-holdout (rigorous):** Each of 377 genera (≥3 ASVs) sequentially held out. K=5 nearest neighbors from remaining ASVs predict higher-level taxonomy. This tests generalization across genera, not memorization of near-duplicates.

**Family-holdout:** Same approach at family level (169 families).

### KNN Assignment
- Distance metric: Cosine distance in 4,096-dim embedding space
- K = 5 nearest neighbors
- Prediction: Majority vote
- Confidence calibration based on validated accuracy at each taxonomic level:
  - Phylum/Class: "reliable" (85-89%)
  - Order/Family: "putative" (45-71%)
  - Genus: "speculative" (not validated via holdout)

### Distance-Based Confidence
Correct predictions had significantly shorter nearest-neighbor distances than incorrect ones (Family level: mean 0.0041 vs 0.0082, Mann-Whitney U p < 10⁻¹⁰), supporting distance as a confidence metric.
