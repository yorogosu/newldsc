import pandas as pd
import scipy as sp
import time
import sys

def get_bulik_sullivan_15_sids(sids_dir='/project/HeritPartition/faststorage/ld_scores_bulik_sullivan_2015/eur_w_ld_chr'):
    """
    Creates a by-chromosome dictionary of the ~1.2M SNP rsid's for which Bulik-Sullivan et al., (2015) reported LD scores.
    "Sids" stands for "SNP ID's".
    This function was originally created by Bjarni.
    """
    bulik_sullivan_snps = {}
    for chrom in range(1, 23):
        chrom_str = 'chr%d' % chrom
        filename = '%s/%d.l2.ldscore' % (sids_dir, chrom)
        ld_score_table = pd.read_table(filename)
        sids = sp.array(ld_score_table['SNP'])
        bulik_sullivan_snps[chrom_str] = sids
    return bulik_sullivan_snps


def calc_ld_table(norm_snps, cM_map, max_cM_dist=1, min_r2=0.2, verbose=True, normalize=True):
    """
    Calculate LD among all SNPs within a sliding window.
    
    This function only retains r^2 values above the given threshold.
    
    Written by Bjarni, modified by Yorgos to work with cM-based windows.
    """
    if min_r2 == 0.0:
        min_r2 = -1000 # to make sure that when we set r2 to 0 (no cutoff), even negative values will be included
    
    # Normalize SNPs (perhaps not necessary, but cheap)
    if normalize:
        norm_snps = norm_snps.T
        norm_snps = (norm_snps - sp.mean(norm_snps, 0)) / sp.std(norm_snps, 0, ddof=1)
        norm_snps = norm_snps.T
    if verbose:
        print 'Calculating LD table'
    t0 = time.time()
    num_snps, num_indivs = norm_snps.shape

    ld_table = {}
    for i in range(num_snps):
        ld_table[i] = {}

    for i in range(num_snps):
        start_i = norm_snps[abs(cM_map-cM_map[i])<=max_cM_dist].index.values[0]
        end_i = norm_snps[abs(cM_map-cM_map[i])<=max_cM_dist].index.values[-1]
        ld_vec = sp.dot(norm_snps[i:i+1], norm_snps[abs(cM_map-cM_map[i])<=max_cM_dist].T) / float(num_indivs-1)
        
        ld_vec = ld_vec**2 # square
        ld_vec = ld_vec - ((1-ld_vec)/(num_indivs-2)) # subtract error term
        
        ld_vec = sp.array(ld_vec).flatten()
        
        for k in range(start_i, end_i):
            ld_vec_i = k - start_i
            if ld_vec[ld_vec_i] > min_r2:
                ld_table[i][k] = ld_vec[ld_vec_i]
                ld_table[k][i] = ld_vec[ld_vec_i]

        if verbose:
            if i % 1000 == 0:
                sys.stdout.write('.')
                sys.stdout.flush()
    if verbose:
        sys.stdout.write('Done.\n')
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to calculate the LD table' % (t / 60, t % 60)
    del norm_snps
    return ld_table


def calc_ld_score(ld_table):
    """
    Calculates LD scores by summing over all r^2 values assigned to a specific tag SNP.
    
    Returns a list of LD scores
    """
    ld_scores = []
    
    for i in range(len(ld_table)):
        ld_score = 0
        for j in ld_table[i]:
            ld_score += ld_table[i][j]
        if i in ld_table[i]:
            ld_score = ld_score - ld_table[i][i] # ld_table[i][i] ~= 1
        ld_scores.append(ld_score)
    return ld_scores

