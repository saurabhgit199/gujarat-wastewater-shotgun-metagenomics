#!/usr/bin/env python3
"""
Metagenomics × Metabolomics Integration
Gujarat Urban Wastewater — Local Version

Run from the directory containing your MASTERSHEET.xlsx:
    python multiomics_integration.py

Or specify the file:
    python multiomics_integration.py --lcms path/to/MASTERSHEET.xlsx
"""

import argparse
import json
import re
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # non-interactive backend — saves figures to disk
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings('ignore')

# ── Output directory ─────────────────────────────────────────────
OUT = Path('multiomics_results')
OUT.mkdir(exist_ok=True)

# ════════════════════════════════════════════════════════════════
# 1. CONFIGURATION
# ════════════════════════════════════════════════════════════════

# Time point → Season mapping (confirmed from phyloseq metadata)
TP_TO_SEASON = {
    '03': 'PreWinter', '04': 'PreWinter',
    '05': 'PreWinter', '06': 'PreWinter',
    '07': 'Winter',    '08': 'Winter',
    '09': 'Winter',    '10': 'Winter',
    '11': 'Winter',    '12': 'Winter',
    '13': 'Summer',    '14': 'Summer',
    '15': 'Summer',    '16': 'Summer',
    '17': 'Summer',    '18': 'Summer',
    '19': 'Summer',    # Jun early = Summer
    '20': 'Monsoon',   # Jun late  = Monsoon
    '21': 'Monsoon',   '22': 'Monsoon',
    '23': 'Monsoon',   '24': 'Monsoon',
    '25': 'Monsoon',   '26': 'Monsoon',
}

CITY_MAP = {
    'A': 'Ahmedabad',
    'G': 'Gandhinagar',
    'R': 'Rajkot',
    'V': 'Vadodara',
}

# Metagenomics CE loss — from validated METAGENE-1 systematic experiments
# (Baseline_allreads configuration, 16 city × season validation samples)
META_CE_LOSS = {
    ('Ahmedabad',   'Monsoon')   : 5.0262,
    ('Ahmedabad',   'Summer')    : 4.8895,
    ('Ahmedabad',   'Winter')    : 4.8577,
    ('Ahmedabad',   'PreWinter') : 4.8956,
    ('Gandhinagar', 'PreWinter') : 4.8595,
    ('Gandhinagar', 'Monsoon')   : 4.8725,
    ('Gandhinagar', 'Summer')    : 4.8672,
    ('Gandhinagar', 'Winter')    : 4.8754,
    ('Rajkot',      'Summer')    : 4.9529,
    ('Rajkot',      'Winter')    : 4.8309,
    ('Rajkot',      'PreWinter') : 4.7817,
    ('Rajkot',      'Monsoon')   : 4.8600,
    ('Vadodara',    'Monsoon')   : 4.9275,
    ('Vadodara',    'Summer')    : 4.8758,
    ('Vadodara',    'PreWinter') : 4.7227,
    ('Vadodara',    'Winter')    : 4.8651,
}

CITY_COLORS = {
    'Ahmedabad':   '#1f77b4',
    'Gandhinagar': '#ff7f0e',
    'Rajkot':      '#2ca02c',
    'Vadodara':    '#d62728',
}
SEASON_MARKERS = {
    'Summer': 'o', 'Monsoon': 's',
    'PreWinter': '^', 'Winter': 'D',
}


# ════════════════════════════════════════════════════════════════
# 2. LOAD LCMS DATA
# ════════════════════════════════════════════════════════════════

def find_lcms_file(path_hint=None):
    """Find LCMS master sheet automatically."""
    candidates = [
        path_hint,
        'MASTERSHEET.xlsx',
        'Sorted_LCMS_Data.xlsx',
        'sorted_LCMS_Data.xlsx',
        'SITE_WISE.xlsx',
        'output.xlsx',
    ]
    for c in candidates:
        if c and Path(c).exists():
            return Path(c)
    # Search common locations
    home = Path.home()
    for pattern in [
        '**/MASTERSHEET.xlsx',
        '**/Sorted_LCMS_Data.xlsx',
        '**/LCMS_sorted_data_AUC/MASTERSHEET.xlsx',
    ]:
        matches = list(home.glob(pattern))
        if matches:
            return matches[0]
    return None


def load_lcms(filepath):
    """Load and orient the LCMS matrix as samples × metabolites."""
    print(f'\nLoading: {filepath}')
    xl  = pd.ExcelFile(filepath)
    print(f'  Sheets: {xl.sheet_names}')
    raw = pd.read_excel(filepath, sheet_name=0, index_col=0)
    print(f'  Raw shape: {raw.shape}')

    # Auto-detect orientation
    pat = re.compile(r'^[AGRVagrv][A-Za-z]\d{2}', re.IGNORECASE)
    n_col = sum(bool(pat.match(str(c))) for c in raw.columns)
    n_row = sum(bool(pat.match(str(r))) for r in raw.index)

    if n_col > n_row:
        df = raw.T
        print(f'  Transposed → samples × metabolites')
    else:
        df = raw.copy()
        print(f'  Already → samples × metabolites')

    df = df.apply(pd.to_numeric, errors='coerce')
    print(f'  Shape after orient: {df.shape}')
    print(f'  Sample IDs (first 8): {df.index[:8].tolist()}')
    print(f'  Missing: {df.isna().mean().mean()*100:.1f}%')
    return df


# ════════════════════════════════════════════════════════════════
# 3. DECODE METADATA + AGGREGATE
# ════════════════════════════════════════════════════════════════

def decode_sample_id(sid):
    """Extract city, season, city_season from LCMS sample ID."""
    sid = str(sid).strip()
    m   = re.match(r'^([AGRVagrv])([A-Za-z])(\d{2})', sid)
    if not m:
        return None
    city   = CITY_MAP.get(m.group(1).upper(), None)
    tp     = m.group(3)
    season = TP_TO_SEASON.get(tp, None)
    if not city or not season:
        return None
    return {
        'lcms_id':     sid,
        'city':        city,
        'season':      season,
        'timepoint':   tp,
        'city_season': f'{city}_{season}',
    }


def aggregate_by_city_season(df_lcms):
    """
    Aggregate LCMS matrix by city × season.
    Mirrors the shotgun pooling strategy exactly.
    Returns aggregated dataframe + sample counts per group.
    """
    records = []
    for sid in df_lcms.index:
        d = decode_sample_id(sid)
        if d:
            records.append(d)

    df_meta_lcms = pd.DataFrame(records).set_index('lcms_id')
    print(f'\nDecoded {len(df_meta_lcms)} / {len(df_lcms)} sample IDs')

    # Join city_season
    df_work             = df_lcms.copy()
    df_work['__cs__']   = df_meta_lcms.reindex(df_work.index)['city_season']
    df_work             = df_work[df_work['__cs__'].notna()]

    met_cols = [c for c in df_work.columns if c != '__cs__']
    df_agg   = df_work.groupby('__cs__')[met_cols].mean()
    n_per    = df_work.groupby('__cs__').size()

    print(f'\nAggregated matrix: {df_agg.shape[0]} city-seasons × {df_agg.shape[1]} metabolites')
    print(f'Samples per group (mean {n_per.mean():.1f}, min {n_per.min()}):')
    for cs, n in sorted(n_per.items()):
        print(f'  {cs:<32} n={n}')

    return df_agg, n_per, df_meta_lcms


# ════════════════════════════════════════════════════════════════
# 4. MATCH + PREPROCESS
# ════════════════════════════════════════════════════════════════

def build_meta_df():
    """Build metagenomics dataframe."""
    rows = []
    for (city, season), ce in META_CE_LOSS.items():
        rows.append({
            'city_season': f'{city}_{season}',
            'city':        city,
            'season':      season,
            'ce_loss':     ce,
        })
    return pd.DataFrame(rows).set_index('city_season')


def preprocess(df_lcms_matched, presence_threshold=0.5):
    """Log10 + autoscale LCMS data."""
    # Filter: keep metabolites present in ≥50% of samples
    keep   = df_lcms_matched.notna().mean() >= presence_threshold
    df_f   = df_lcms_matched.loc[:, keep].copy()
    print(f'After ≥{presence_threshold*100:.0f}% presence filter: {df_f.shape[1]} metabolites')

    # Impute: half-minimum
    for col in df_f.columns:
        half_min = df_f[col].min(skipna=True) / 2
        df_f[col].fillna(half_min, inplace=True)

    # Log10
    df_log = np.log10(df_f.clip(lower=1e-10) + 1)

    # Autoscale
    sc      = StandardScaler()
    df_sc   = pd.DataFrame(
        sc.fit_transform(df_log),
        index=df_log.index, columns=df_log.columns
    )
    return df_sc, df_f


# ════════════════════════════════════════════════════════════════
# 5. ANALYSES
# ════════════════════════════════════════════════════════════════

def analysis_diversity(df_met_raw, ce_loss, city, season):
    """CE loss vs metabolite richness and diversity."""
    print('\n' + '='*60)
    print('ANALYSIS 1: CE loss vs metabolite richness / diversity')
    print('='*60)

    ce_arr = ce_loss.values

    richness = (df_met_raw > df_met_raw.min() * 0.6).sum(axis=1)
    r1, p1   = stats.spearmanr(ce_arr, richness.reindex(ce_loss.index).values)

    def shannon(row):
        v = row.dropna().clip(lower=0).values
        v = v[v > 0]
        if not len(v): return 0
        p = v / v.sum()
        return float(-np.sum(p * np.log(p + 1e-12)))

    diversity = df_met_raw.apply(shannon, axis=1)
    r2, p2    = stats.spearmanr(ce_arr, diversity.reindex(ce_loss.index).values)

    total_auc = df_met_raw.sum(axis=1)
    r3, p3    = stats.spearmanr(ce_arr, total_auc.reindex(ce_loss.index).values)

    print(f'CE loss vs richness   : r={r1:.3f}, p={p1:.4f}')
    print(f'CE loss vs diversity  : r={r2:.3f}, p={p2:.4f}')
    print(f'CE loss vs total AUC  : r={r3:.3f}, p={p3:.4f}')

    print(f'\n{"City-Season":<32} {"CE Loss":>9} {"Richness":>10} {"Shannon":>10}')
    print('-' * 65)
    for cs in ce_loss.sort_values(ascending=False).index:
        print(f'{cs:<32} {ce_loss[cs]:>9.4f} '
              f'{richness.get(cs, 0):>10.0f} '
              f'{diversity.get(cs, 0):>10.3f}')
    return {'r_richness': r1, 'p_richness': p1,
            'r_shannon': r2, 'p_shannon': p2}


def analysis_metabolite_correlations(df_met_scaled, ce_loss):
    """Spearman correlation of each metabolite with CE loss."""
    print('\n' + '='*60)
    print('ANALYSIS 2: Metabolite–CE loss correlations')
    print('='*60)

    ce_arr = ce_loss.values
    corrs  = []
    for col in df_met_scaled.columns:
        arr = df_met_scaled[col].reindex(ce_loss.index).values
        if np.std(arr) < 1e-6: continue
        r, p = stats.spearmanr(ce_arr, arr)
        corrs.append({'metabolite': str(col), 'r': r, 'p': p})

    df_c = pd.DataFrame(corrs).sort_values('r', key=abs, ascending=False)

    # FDR correction
    try:
        from statsmodels.stats.multitest import multipletests
        p_vals = df_c['p'].fillna(1.0).values
        _, p_fdr, _, _ = multipletests(p_vals, method='fdr_bh')
        df_c['p_fdr'] = p_fdr
    except ImportError:
        df_c['p_fdr'] = (df_c['p'].fillna(1.0) * len(df_c)).clip(upper=1.0)

    print(f'\nTop 20 metabolites by |r|:')
    print(f'{"Metabolite":<40} {"r":>7} {"p":>9} {"p_FDR":>9}')
    print('-' * 70)
    for _, row in df_c.head(20).iterrows():
        sig = ('***' if row['p_fdr']<0.001 else
               ('**'  if row['p_fdr']<0.01  else
                ('*'   if row['p_fdr']<0.05  else '')))
        direction = '↑' if row['r'] > 0 else '↓'
        print(f'{str(row["metabolite"])[:39]:<40} '
              f'{row["r"]:>7.3f} {row["p"]:>9.4f} '
              f'{row["p_fdr"]:>9.4f} {sig}{direction}')

    sig = df_c[df_c['p_fdr'] < 0.05]
    print(f'\nFDR-significant: {len(sig)} metabolites '
          f'(pos: {(sig["r"]>0).sum()}, neg: {(sig["r"]<0).sum()})')
    return df_c


def analysis_seasonal_validation(df_met_scaled, ce_loss, city, season):
    """Does metabolomics PC1 agree with CE loss seasonal/city pattern?"""
    print('\n' + '='*60)
    print('ANALYSIS 3: Seasonal/city validation')
    print('='*60)

    # Fill remaining NaN before PCA
    X = df_met_scaled.reindex(ce_loss.index).copy()
    X = X.fillna(X.mean())
    pca = PCA(n_components=1)
    pc1 = pca.fit_transform(X).flatten()
    pc1_ser = pd.Series(pc1, index=ce_loss.index)

    r, p = stats.spearmanr(ce_loss.values, pc1_ser.values)
    print(f'Overall Spearman r (CE vs Met PC1): r={r:.3f}, p={p:.4f}')

    print(f'\n{"City-Season":<32} {"CE Loss":>9} {"Met PC1":>10}')
    print('-' * 54)
    for cs in ce_loss.sort_values(ascending=False).index:
        print(f'{cs:<32} {ce_loss[cs]:>9.4f} {pc1_ser[cs]:>10.3f}')

    seasons = ['Summer', 'Monsoon', 'PreWinter', 'Winter']
    print(f'\nSeasonal means:')
    print(f'{"Season":<12} {"CE mean":>10} {"Met PC1":>10}')
    for s in seasons:
        keys_s = [k for k in ce_loss.index if k.endswith(s)]
        if keys_s:
            print(f'{s:<12} {ce_loss[keys_s].mean():>10.4f} '
                  f'{pc1_ser[keys_s].mean():>10.3f}')

    if p < 0.05:
        print('\n✓ Significant agreement — both methods capture same variation')
    elif abs(r) > 0.4:
        print('\n~ Moderate trend (low n=16 limits power)')
    else:
        print('\n→ Complementary — methods capture different biology')
    return r, p, pc1_ser


# ════════════════════════════════════════════════════════════════
# 6. FIGURES
# ════════════════════════════════════════════════════════════════

def make_figures(df_met_scaled, ce_loss, city, season, df_corr, pc1_ser):
    """Generate and save all figures."""
    print('\nGenerating figures...')

    X_fig = df_met_scaled.reindex(ce_loss.index).copy().fillna(0)
    pca = PCA(n_components=min(3, len(ce_loss) - 1))
    pcs = pca.fit_transform(X_fig)
    evr = pca.explained_variance_ratio_

    # ── Figure 1: Three-panel ─────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel 1: PCA coloured by CE loss
    ax   = axes[0]
    norm = mcolors.Normalize(vmin=ce_loss.min(), vmax=ce_loss.max())
    cmap = plt.cm.RdYlGn_r
    for i, cs in enumerate(ce_loss.index):
        mk  = SEASON_MARKERS.get(season[cs], 'o')
        col = cmap(norm(ce_loss[cs]))
        ax.scatter(pcs[i, 0], pcs[i, 1], c=[col], marker=mk,
                   s=200, edgecolor='black', linewidth=0.8, zorder=5)
        ax.annotate(cs.split('_')[0][:3], (pcs[i, 0], pcs[i, 1]),
                    fontsize=7, ha='center', va='bottom',
                    xytext=(0, 5), textcoords='offset points')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    plt.colorbar(sm, ax=ax, label='CE Loss', shrink=0.8)
    ax.set_xlabel(f'PC1 ({evr[0]*100:.1f}%)', fontsize=11)
    ax.set_ylabel(f'PC2 ({evr[1]*100:.1f}%)', fontsize=11)
    ax.set_title('Metabolomics PCA\nColoured by CE Loss',
                 fontweight='bold', fontsize=11)

    # Panel 2: PCA coloured by city
    ax = axes[1]
    for i, cs in enumerate(ce_loss.index):
        col = CITY_COLORS.get(city[cs], 'grey')
        mk  = SEASON_MARKERS.get(season[cs], 'o')
        ax.scatter(pcs[i, 0], pcs[i, 1], c=col, marker=mk,
                   s=200, edgecolor='black', linewidth=0.8, zorder=5)
    ch = [Line2D([0],[0], marker='o', color='w',
                 markerfacecolor=c, markersize=10, label=n)
          for n, c in CITY_COLORS.items()]
    sh = [Line2D([0],[0], marker=m, color='grey',
                 markersize=10, linewidth=0, label=s)
          for s, m in SEASON_MARKERS.items()]
    ax.legend(handles=ch+sh, fontsize=7, ncol=2)
    ax.set_xlabel(f'PC1 ({evr[0]*100:.1f}%)', fontsize=11)
    ax.set_ylabel(f'PC2 ({evr[1]*100:.1f}%)', fontsize=11)
    ax.set_title('Metabolomics PCA\nCity (colour) × Season (marker)',
                 fontweight='bold', fontsize=11)

    # Panel 3: CE loss vs top metabolite
    ax      = axes[2]
    top_met = df_corr.iloc[0]['metabolite']
    top_r   = df_corr.iloc[0]['r']
    top_p   = df_corr.iloc[0]['p']
    met_v   = df_met_scaled[top_met].reindex(ce_loss.index)
    for cs in ce_loss.index:
        col = CITY_COLORS.get(city[cs], 'grey')
        mk  = SEASON_MARKERS.get(season[cs], 'o')
        ax.scatter(ce_loss[cs], met_v[cs], c=col, marker=mk,
                   s=200, edgecolor='black', linewidth=0.8, zorder=5)
        ax.annotate(cs.split('_')[0][:3], (ce_loss[cs], met_v[cs]),
                    fontsize=7, ha='left',
                    xytext=(4, 0), textcoords='offset points')
    m_fit, b_fit = np.polyfit(ce_loss.values, met_v.values, 1)
    xl = np.linspace(ce_loss.min(), ce_loss.max(), 100)
    ax.plot(xl, m_fit*xl+b_fit, 'k--', lw=1.5, alpha=0.7)
    ax.set_xlabel('METAGENE-1 CE Loss', fontsize=11)
    ax.set_ylabel(f'{str(top_met)[:25]}', fontsize=9)
    ax.set_title(f'Top correlated metabolite\nr={top_r:.3f}, p={top_p:.3f}',
                 fontweight='bold', fontsize=11)

    plt.suptitle(
        'Gujarat Wastewater: Metagenomics × Metabolomics Integration\n'
        f'({len(ce_loss)} city × season combinations)',
        fontsize=11, fontweight='bold'
    )
    plt.tight_layout()
    fig.savefig(OUT / 'multiomics_pca_figure.pdf', dpi=150, bbox_inches='tight')
    fig.savefig(OUT / 'multiomics_pca_figure.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  ✓ {OUT}/multiomics_pca_figure.pdf')

    # ── Figure 2: Heatmap ─────────────────────────────────────
    top30   = df_corr.head(min(30, len(df_corr)))['metabolite'].tolist()
    df_heat = df_met_scaled.reindex(ce_loss.index)[top30]

    row_colors = pd.DataFrame({
        'City': pd.Series(
            [CITY_COLORS.get(city[cs], 'grey') for cs in df_heat.index],
            index=df_heat.index
        )
    })
    try:
        g = sns.clustermap(
            df_heat, row_colors=row_colors,
            cmap='RdBu_r', figsize=(14, 8),
            xticklabels=True, yticklabels=True,
            vmin=-3, vmax=3,
            dendrogram_ratio=0.12,
        )
        g.ax_heatmap.set_xticklabels(
            g.ax_heatmap.get_xmajorticklabels(),
            fontsize=7, rotation=45, ha='right'
        )
        g.ax_heatmap.set_yticklabels(
            g.ax_heatmap.get_ymajorticklabels(), fontsize=8
        )
        g.fig.suptitle(
            'Top Metabolites Correlated with METAGENE-1 CE Loss\n'
            'Gujarat urban wastewater — city × season aggregated',
            fontsize=10, fontweight='bold', y=1.02
        )
        g.fig.savefig(OUT / 'metabolite_heatmap.pdf', dpi=150,
                      bbox_inches='tight')
        g.fig.savefig(OUT / 'metabolite_heatmap.png', dpi=150,
                      bbox_inches='tight')
        plt.close('all')
        print(f'  ✓ {OUT}/metabolite_heatmap.pdf')
    except Exception as e:
        print(f'  Heatmap failed: {e}')

    # ── Figure 3: CE loss vs Met PC1 scatter ─────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    for cs in ce_loss.index:
        col = CITY_COLORS.get(city[cs], 'grey')
        mk  = SEASON_MARKERS.get(season[cs], 'o')
        ax.scatter(ce_loss[cs], pc1_ser[cs], c=col, marker=mk,
                   s=180, edgecolor='black', linewidth=0.8, zorder=5)
        ax.annotate(cs.replace('_', '\n'), (ce_loss[cs], pc1_ser[cs]),
                    fontsize=6, ha='center', va='bottom',
                    xytext=(0, 6), textcoords='offset points')
    r_v, p_v = stats.spearmanr(ce_loss.values, pc1_ser.values)
    m_f, b_f = np.polyfit(ce_loss.values, pc1_ser.values, 1)
    xl        = np.linspace(ce_loss.min(), ce_loss.max(), 100)
    ax.plot(xl, m_f*xl+b_f, 'k--', lw=1.5, alpha=0.7)
    ax.set_xlabel('METAGENE-1 CE Loss', fontsize=12)
    ax.set_ylabel('Metabolomics PC1', fontsize=12)
    ax.set_title(
        f'Metagenomics vs Metabolomics\n'
        f'Spearman r = {r_v:.3f}, p = {p_v:.4f}',
        fontweight='bold', fontsize=11
    )
    ch2 = [Line2D([0],[0], marker='o', color='w',
                  markerfacecolor=c, markersize=10, label=n)
           for n, c in CITY_COLORS.items()]
    sh2 = [Line2D([0],[0], marker=m, color='grey',
                  markersize=10, linewidth=0, label=s)
           for s, m in SEASON_MARKERS.items()]
    ax.legend(handles=ch2+sh2, fontsize=8, ncol=2)
    plt.tight_layout()
    fig.savefig(OUT / 'ce_vs_metabolomics_pc1.pdf', dpi=150,
                bbox_inches='tight')
    fig.savefig(OUT / 'ce_vs_metabolomics_pc1.png', dpi=150,
                bbox_inches='tight')
    plt.close(fig)
    print(f'  ✓ {OUT}/ce_vs_metabolomics_pc1.pdf')


# ════════════════════════════════════════════════════════════════
# 7. MAIN
# ════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Metagenomics × Metabolomics integration')
    parser.add_argument('--lcms', type=str, default=None,
                        help='Path to LCMS master sheet (xlsx)')
    args = parser.parse_args()

    # ── Find and load LCMS data ───────────────────────────────
    lcms_path = find_lcms_file(args.lcms)
    if lcms_path is None:
        print('ERROR: Could not find LCMS file.')
        print('Run as:  python multiomics_integration.py --lcms path/to/MASTERSHEET.xlsx')
        sys.exit(1)

    df_lcms = load_lcms(lcms_path)

    # ── Aggregate by city × season ────────────────────────────
    df_agg, n_per, df_meta_lcms = aggregate_by_city_season(df_lcms)

    # ── Build metagenomics dataframe ──────────────────────────
    df_meta = build_meta_df()

    # ── Match ─────────────────────────────────────────────────
    common = sorted(set(df_meta.index) & set(df_agg.index))
    print(f'\nMatched city × season pairs: {len(common)}')
    if not common:
        print('ERROR: No matching city-season keys.')
        print('Check that TP_TO_SEASON mapping produces correct season names.')
        sys.exit(1)

    df_lcms_matched = df_agg.loc[common]
    df_meta_matched = df_meta.loc[common]
    ce_loss = df_meta_matched['ce_loss']
    city    = df_meta_matched['city']
    season  = df_meta_matched['season']

    print(f'\nCE loss summary:')
    print(f'  Most anomalous  : {ce_loss.idxmax()} ({ce_loss.max():.4f})')
    print(f'  Least anomalous : {ce_loss.idxmin()} ({ce_loss.min():.4f})')

    # ── Preprocess metabolomics ───────────────────────────────
    df_met_scaled, df_met_raw = preprocess(df_lcms_matched)

    # ── Analyses ──────────────────────────────────────────────
    div_results = analysis_diversity(df_met_raw, ce_loss, city, season)
    df_corr     = analysis_metabolite_correlations(df_met_scaled, ce_loss)
    r_v, p_v, pc1_ser = analysis_seasonal_validation(
        df_met_scaled, ce_loss, city, season)

    # ── Figures ───────────────────────────────────────────────
    make_figures(df_met_scaled, ce_loss, city, season, df_corr, pc1_ser)

    # ── Save results ──────────────────────────────────────────
    print(f'\nSaving results to {OUT}/')
    df_met_scaled.to_csv(OUT / 'metabolomics_city_season_scaled.csv')
    df_agg.loc[common].to_csv(OUT / 'metabolomics_city_season_aggregated.csv')
    df_corr.to_csv(OUT / 'metabolite_ce_correlations.csv', index=False)
    df_meta_matched.to_csv(OUT / 'metagenomics_ce_loss_matched.csv')

    sig_fdr = df_corr[df_corr['p_fdr'] < 0.05]
    summary = {
        'n_city_season_pairs':         len(common),
        'n_metabolites_raw':           df_lcms_matched.shape[1],
        'n_metabolites_after_filter':  df_met_scaled.shape[1],
        'ce_loss_most_anomalous':      str(ce_loss.idxmax()),
        'ce_loss_least_anomalous':     str(ce_loss.idxmin()),
        'spearman_r_met_pc1_vs_ce':    float(r_v),
        'spearman_p_met_pc1_vs_ce':    float(p_v),
        'n_fdr_sig_metabolites':       int(len(sig_fdr)),
        'diversity_r_richness':        float(div_results['r_richness']),
        'diversity_p_richness':        float(div_results['p_richness']),
        'top5_positive_metabolites':   df_corr[df_corr['r']>0].head(5)['metabolite'].tolist(),
        'top5_negative_metabolites':   df_corr[df_corr['r']<0].head(5)['metabolite'].tolist(),
    }
    with open(OUT / 'summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    print(f'\n{"="*60}')
    print(f'DONE — results in {OUT}/')
    print(f'{"="*60}')
    for f in sorted(OUT.iterdir()):
        print(f'  {f.name}')

    print(f'\nKEY RESULTS:')
    print(f'  Matched pairs              : {len(common)}')
    print(f'  Metabolites retained       : {df_met_scaled.shape[1]}')
    print(f'  FDR-sig metabolites        : {len(sig_fdr)}')
    print(f'  Met PC1 vs CE loss r       : {r_v:.3f} (p={p_v:.4f})')
    print(f'  Most anomalous city-season : {ce_loss.idxmax()} ({ce_loss.max():.4f})')


if __name__ == '__main__':
    main()
