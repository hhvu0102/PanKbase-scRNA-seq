#!/usr/bin/env python3
# coding: utf-8

import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import tables
import anndata
from typing import Dict, Optional
import numpy as np
import scipy.sparse as sp
from scipy import io
import glob
import os
import upsetplot
from scipy.io import mmread
import csv

import argparse

parser = argparse.ArgumentParser("Plot QC metrics per sample")
parser.add_argument("-sample", "--sample", help="Donor ID.", type=str)
parser.add_argument("-v2bc", "--v2bc", help="V2 chemistry whitelist barcodes.", type=str)
parser.add_argument("-v3bc", "--v2bc", help="V3 chemistry whitelist barcodes.", type=str)
parser.add_argument("-knee", "--knee", help="File with knee, infection and end cliff points.", type=str)
parser.add_argument("-passBC", "--passBC", help="File with high qual cells by emptyDrops.", type=str)
parser.add_argument("-cellbender", "--cellbender", help="Path to cellbender h5 file.", type=str)
parser.add_argument("-rnametrics", "--rnametrics", help="Path to metrics file, typically in qc directory.", type=str)
parser.add_argument("-starsolo", "--starsolo", help="Path to Starsolo/*Solo.out results.", type=str)
parser.add_argument("-qcPlot", "--qcPlot", help="Path to save qcPlot plots.", type=str)
parser.add_argument("-upsetPlot", "--upsetPlot", help="Path to save upset plots.", type=str)
parser.add_argument("-outmetrics", "--outmetrics", help="Path to save all metrics results.", type=str)

args = parser.parse_args()

sample = args.sample
passQC = args.passBC
knee = args.knee

with open('/nfs/turbo/umms-scjp-pank/4_integration/data/202503_freeze/20250314_meta_runs.txt', 'r') as file: #hardcoded
    reader = csv.reader(file, delimiter='\t')
    for row in reader:
        if row[1] == sample:
            metadata = row

if metadata[4] == 'V2':
    RNA_BARCODE_WHITELIST = args.v2bc
else:
    RNA_BARCODE_WHITELIST = args.v3bc


CELLBENDER = args.cellbender
RNA_METRICS = args.rnametrics
GENE_FULL_EXON_OVER_INTRON_COUNTS = args.starsolo + '/GeneFull_ExonOverIntron/raw'
GENE_COUNTS = args.starsolo + '/Gene/raw'

THRESHOLD_CELLBENDER_MIN_CELL_PROBABILITY = 0.99

#### FUNCTIONS FROM CELLBENDER
def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d


def anndata_from_h5(file: str,
                    analyzed_barcodes_only: bool = True) -> 'anndata.AnnData':
    """Load an output h5 file into an AnnData object for downstream work.
    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count matrix.
            True to load a limited set of barcodes: only those analyzed by the
            algorithm. This allows relevant latent variables to be loaded
            properly into adata.obs and adata.obsm, rather than adata.uns.
    Returns:
        adata: The anndata object, populated with inferred latent variables
            and metadata.
    """

    d = dict_from_h5(file)
    X = sp.csc_matrix((d.pop('data'), d.pop('indices'), d.pop('indptr')),
                      shape=d.pop('shape')).transpose().tocsr()

    # check and see if we have barcode index annotations, and if the file is filtered
    barcode_key = [k for k in d.keys() if (('barcode' in k) and ('ind' in k))]
    if len(barcode_key) > 0:
        max_barcode_ind = d[barcode_key[0]].max()
        filtered_file = (max_barcode_ind >= X.shape[0])
    else:
        filtered_file = True

    if analyzed_barcodes_only:
        if filtered_file:
            # filtered file being read, so we don't need to subset
            print('Assuming we are loading a "filtered" file that contains only cells.')
            pass
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        elif 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        else:
            print('Warning: analyzed_barcodes_only=True, but the key '
                  '"barcodes_analyzed_inds" or "barcode_indices_for_latents" '
                  'is missing from the h5 file. '
                  'Will output all barcodes, and proceed as if '
                  'analyzed_barcodes_only=False')

    # Construct the anndata object.
    adata = anndata.AnnData(X=X,
                            obs={'barcode': d.pop('barcodes').astype(str)},
                            var={'gene_name': (d.pop('gene_names') if 'gene_names' in d.keys()
                                               else d.pop('name')).astype(str)},
                            dtype=X.dtype)
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_name', inplace=True)

    # For CellRanger v2 legacy format, "gene_ids" was called "genes"... rename this
    if 'genes' in d.keys():
        d['id'] = d.pop('genes')

    # For purely aesthetic purposes, rename "id" to "gene_id"
    if 'id' in d.keys():
        d['gene_id'] = d.pop('id')

    # If genomes are empty, try to guess them based on gene_id
    if 'genome' in d.keys():
        if np.array([s.decode() == '' for s in d['genome']]).all():
            if '_' in d['gene_id'][0].decode():
                print('Genome field blank, so attempting to guess genomes based on gene_id prefixes')
                d['genome'] = np.array([s.decode().split('_')[0] for s in d['gene_id']], dtype=str)

    # Add other information to the anndata object in the appropriate slot.
    _fill_adata_slots_automatically(adata, d)

    # Add a special additional field to .var if it exists.
    if 'features_analyzed_inds' in adata.uns.keys():
        adata.var['cellbender_analyzed'] = [True if (i in adata.uns['features_analyzed_inds'])
                                            else False for i in range(adata.shape[1])]

    if analyzed_barcodes_only:
        for col in adata.obs.columns[adata.obs.columns.str.startswith('barcodes_analyzed')
                                     | adata.obs.columns.str.startswith('barcode_indices')]:
            try:
                del adata.obs[col]
            except Exception:
                pass
    else:
        # Add a special additional field to .obs if all barcodes are included.
        if 'barcodes_analyzed_inds' in adata.uns.keys():
            adata.obs['cellbender_analyzed'] = [True if (i in adata.uns['barcodes_analyzed_inds'])
                                                else False for i in range(adata.shape[0])]

    return adata


def _fill_adata_slots_automatically(adata, d):
    """Add other information to the adata object in the appropriate slot."""

    for key, value in d.items():
        try:
            if value is None:
                continue
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == adata.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == adata.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print('Unable to load data into AnnData: ', key, value, type(value))


#### END FUNCTIONS FROM CELLBENDER

def cellbender_anndata_to_cell_probability(a):
    return a.obs.cell_probability


def cellbender_anndata_to_sparse_matrix(adata, min_cell_probability=0):
    barcodes = adata.obs[adata.obs.cell_probability>=min_cell_probability].index.to_list()
    features = adata.var.gene_id.to_list()
    matrix = adata[adata.obs.cell_probability>=min_cell_probability].X.transpose()
    return {'features': features, 'barcodes': barcodes, 'matrix': matrix}


def umi_count_after_decontamination(adata):
    x = cellbender_anndata_to_sparse_matrix(adata)
    return dict(zip(x['barcodes'], x['matrix'].sum(axis=0).tolist()[0]))


def barcode_rank_plot(metrics, ax):
    df = metrics.sort_values('rna_umis', ascending=False)
    df['barcode_rank'] = range(1, len(df) + 1)
    sns.scatterplot(x='barcode_rank', y='rna_umis', data=df, ax=ax, hue='pass_all_filters', palette={True: 'red', False: 'black'}, edgecolor=None, alpha=0.2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Barcode rank')
    ax.set_ylabel('UMIs')
    return ax


def rna_umis_vs_rna_mito_plot(metrics, ax):
    sns.scatterplot(x='rna_umis', y='rna_fraction_mitochondrial', data=metrics, ax=ax, hue='pass_all_filters', palette={True: 'red', False: 'black'}, edgecolor=None, alpha=0.02, s=3)
    ax.set_xscale('log')
    ax.set_xlabel('UMIs')
    ax.set_ylabel('Fraction mito. (RNA)')
    return ax


def rna_umis_vs_exon_to_full_gene_body_ratio(metrics, ax):
    sns.scatterplot(x='rna_umis', y='rna_exon_to_full_gene_body_ratio', data=metrics, ax=ax, hue='pass_all_filters', palette={True: 'red', False: 'black'}, edgecolor=None, alpha=0.02, s=3)
    ax.set_xscale('log')
    ax.set_xlabel('UMIs')
    ax.set_ylabel('Exon/full-gene-body ratio (RNA)')
    return ax


def cellbender_fraction_removed(metrics, ax):
    sns.scatterplot(x='rna_umis', y='fraction_cellbender_removed', data=metrics, ax=ax, hue='pass_all_filters', palette={True: 'red', False: 'black'}, edgecolor=None, alpha=0.05)
    ax.set_xscale('log')
    ax.set_xlabel('UMIs')
    ax.set_ylabel('Fraction ambient')
    return ax


def cellbender_cell_probabilities(metrics, ax):
    sns.histplot(x='cell_probability', data=metrics[(metrics.filter_rna_emptyDrops) & (metrics.filter_rna_max_mito)], ax=ax, bins=20)
    ax.set_xlabel('Cellbender cell prob.\nfor cells by EmptyDrops and mito. thresholds')
    return ax


def rna_umis_vs_atac_hqaa_plot(metrics, ax):
    sns.scatterplot(x='rna_umis', y='atac_hqaa', data=metrics, ax=ax, hue='pass_all_filters', palette={True: 'red', False: 'black'}, edgecolor=None, alpha=0.02, s=3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('UMIs (RNA)')
    ax.set_ylabel('Pass filter reads (ATAC)')
    return ax


def atac_hqaa_vs_atac_tss_enrichment_plot(metrics, ax):
    sns.scatterplot(x='atac_hqaa', y='atac_tss_enrichment', data=metrics, ax=ax, hue='pass_all_filters', palette={True: 'red', False: 'black'}, edgecolor=None, alpha=0.02, s=3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Pass filter reads (ATAC)')
    ax.set_ylabel('TSS enrichment')
    return ax

# Multi-Otsu function
from skimage.filters import threshold_multiotsu
def estimate_threshold(x, classes=3):
    # do on logscale
    values = np.log10(x).values
    values = values.reshape((len(values),1))
    thresholds = threshold_multiotsu(image=values, classes=classes, nbins=256)
    # convert back to linear scale
    thresholds = [pow(10, i) for i in thresholds]
    UMI_THRESHOLD = round(thresholds[classes - 2])
    return UMI_THRESHOLD

rna_barcodes = pd.read_csv(RNA_BARCODE_WHITELIST, header=None)[0].to_list()

#load metrics df
adata = anndata_from_h5(CELLBENDER, analyzed_barcodes_only=True)
rna_metrics = pd.read_csv(RNA_METRICS, sep='\t')
rna_metrics = rna_metrics[rna_metrics.barcode!='-']
metrics = rna_metrics.set_index('barcode').rename(columns=lambda x: 'rna_' + x)

metrics = metrics.reset_index()
cell_probability = cellbender_anndata_to_cell_probability(adata)
post_cellbender_umis = umi_count_after_decontamination(adata)

metrics['cell_probability'] = metrics.barcode.map(lambda x: cell_probability[x] if x in cell_probability else np.nan)
metrics['post_cellbender_umis'] = metrics.barcode.map(lambda x: post_cellbender_umis[x] if x in post_cellbender_umis else np.nan)
metrics['fraction_cellbender_removed'] = (metrics.rna_umis - metrics.post_cellbender_umis) / metrics.rna_umis
metrics['mt_pct'] = metrics.rna_fraction_mitochondrial * 100
metrics['pct_cellbender_removed'] = metrics.fraction_cellbender_removed * 100


bc = pd.read_csv(passQC, header=0, delim_whitespace="\t") 
metrics['filter_rna_emptyDrops'] = metrics['barcode'].isin(bc.barcode)


KNEE_FILE = knee
with open(KNEE_FILE, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    next(reader, None)
    for row in reader:
        knee = round(float(row[0]))
        inflection = round(float(row[1]))
        inflection_rank = round(float(row[2]))
        knee_rank = round(float(row[3]))
        endCliff = round(float(row[4]))
        end_cliff_rank = round(float(row[5]))
        plateau = round(float(row[6]))

### 
### get THRESHOLD_RNA_MAX_MITO
THRESHOLD_RNA_MAX_MITO = estimate_threshold(metrics[((metrics.rna_umis < endCliff) |
                                                    (metrics.filter_rna_emptyDrops == True)) &
                                                    (metrics.mt_pct>0) & (metrics.mt_pct<30)].mt_pct.astype(float),
                                           classes = 3)

THRESHOLD_FRACTION_CB_REMOVED = estimate_threshold(metrics[(metrics.filter_rna_emptyDrops == True) &
                                                         (metrics.pct_cellbender_removed>0) &
                                                         (np.isnan(metrics.pct_cellbender_removed) == False)].pct_cellbender_removed.astype(float),
                                                  classes = 2)


### get cells that passed all thresholds
metrics['filter_cellbender_cell_probability'] = metrics.cell_probability >= THRESHOLD_CELLBENDER_MIN_CELL_PROBABILITY
metrics['filter_pct_cellbender_removed'] = metrics.pct_cellbender_removed <= THRESHOLD_FRACTION_CB_REMOVED
metrics['filter_rna_max_mito'] = metrics.rna_fraction_mitochondrial <= THRESHOLD_RNA_MAX_MITO/100

metrics['pass_all_filters'] = metrics.filter(like='filter_').all(axis=1)


# List of pass-QC barcodes
pass_qc_nuclei = list(sorted(metrics[metrics.pass_all_filters].barcode.to_list()))

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(5*3, 8))

ax = axs[0, 0]
barcode_rank_plot(metrics, ax)
ax.axhline(knee, color='red', ls='--', label='Knee')
ax.axhline(inflection, color='green', ls='--', label='Inflection')
ax.axhline(endCliff, color='blue', ls='--', label='End cliff')
ax.axhline(plateau, color='orange', ls='--', label='Plateau')
ax.legend()
ax.set_title('(A)')

ax = axs[0, 1]
rna_umis_vs_rna_mito_plot(metrics, ax)
ax.axhline(THRESHOLD_RNA_MAX_MITO/100, color='blue', ls='--')
ax.axvline(inflection, color='green', ls='--')
ax.set_title('(B)')

ax = axs[0, 2] #this plot is subjected to changes when optimizing for %MT thresholds
sns.histplot(x='rna_fraction_mitochondrial', data=metrics[((metrics.rna_umis < endCliff) | (metrics.filter_rna_emptyDrops == True)) & (metrics.mt_pct>0) & (metrics.mt_pct<30)], ax=ax, log_scale=True)
ax.axvline(THRESHOLD_RNA_MAX_MITO/100, color='blue', ls='--', label='Fraction mito. threshold Multi-otsu= {:,}'.format(THRESHOLD_RNA_MAX_MITO/100))
ax.legend()
ax.set_xlabel('Fraction mito. across cells detected by\nEmptyDrops or UMIs < end cliff\nand 0 < fraction mito. < 0.3')
ax.set_title('(C)')

ax = axs[1, 0]
cellbender_fraction_removed(metrics, ax)
ax.axhline(THRESHOLD_FRACTION_CB_REMOVED/100, color='red', ls='--')
ax.set_title('(D)')

ax = axs[1, 1]
cellbender_cell_probabilities(metrics, ax)
ax.set_title('(E)')

ax = axs[1, 2]
sns.histplot(x='fraction_cellbender_removed', data=metrics[(metrics.barcode!='-')], ax=ax, log_scale=True)
ax.axvline(THRESHOLD_FRACTION_CB_REMOVED/100, color='red', ls='--', label='Fraction ambient removed thres. Multi-otsu = {:,}'.format(THRESHOLD_FRACTION_CB_REMOVED))
ax.legend()
ax.set_xlabel('Fraction ambient removed')
ax.set_title('(F)')

fig.suptitle('{:,} pass QC cells'.format(len(pass_qc_cells)) + " " + donor)

fig.tight_layout()
fig.savefig(args.qcPlot, bbox_inches='tight', dpi=300)

# Plot the number of cells passing each filter
fig, ax = plt.subplots(figsize=(7, 6))
ax.remove()

for_upset = metrics.filter(like='filter_').rename(columns=lambda x: 'pass_' + x)
for_upset = for_upset.groupby(for_upset.columns.to_list()).size()
upsetplot.plot(for_upset, fig=fig, sort_by='cardinality', show_counts=True)
fig.savefig(args.upsetPlot, bbox_inches='tight', dpi=300)

#save metrics for future reference
metrics.to_csv(args.outmetrics, index=False)      


