import scipy as sp
import pandas as pd
import numpy as np
import h5py
from pca_util import pca, calc_pca_score
from transform_util import transform_genotypes
from mapping_util import bp_to_cM
from ld_util import calc_ld_table, calc_ld_score

def get_genotype_data(in_h5f, chrom, maf_thres=0,
                      indiv_filter=None,
                      return_raw_snps=True, return_snps_info=True,
                      return_normalized_snps=True,
                      chrom_ok_snp_dict=None, hwe_filter=True, unscaled_grm=True):
    """
    Write something here.
    """
    chrom_str = 'chr%d' % chrom
    print 'Loading SNPs'
    snps = in_h5f[chrom_str]['snps'][...]
    
    if indiv_filter is not None:
        snps = snps[:, indiv_filter]
    
    if return_snps_info:
        positions = in_h5f[chrom_str]['positions'][...]
        snp_ids = in_h5f[chrom_str]['snp_ids'][...]
        nts = in_h5f[chrom_str]['nts'][...]
        hwe = in_h5f[chrom_str]['hwe'][...]
            
    if chrom_ok_snp_dict is not None:
        if not return_snps_info:
            snp_ids = in_h5f[chrom_str]['snp_ids'][...]
        ok_snp_filter = sp.in1d(snp_ids, chrom_ok_snp_dict[chrom_str])
        snps = snps[ok_snp_filter]
        if return_snps_info:
            positions = positions[ok_snp_filter]
            snp_ids = snp_ids[ok_snp_filter]
            nts = nts[ok_snp_filter]
            hwe = hwe[ok_snp_filter]
    
    snp_means = sp.mean(snps, 1)
    snp_means.shape = (len(snp_means), 1)
    snp_stds = sp.std(snps, 1)
    snp_stds.shape = (len(snp_stds), 1)
    
    ind_means = sp.mean(snps, 0)
    ind_means.shape = (1, len(ind_means))
    ind_stds = sp.std(snps, 0)
    ind_stds.shape = (1, len(ind_stds))   
    
    if maf_thres > 0:
        print 'Filtering SNPs with MAF <', maf_thres
        std_thres = sp.sqrt(2.0 * (1 - maf_thres) * (maf_thres))
        maf_filter = snp_stds.flatten() > std_thres
        snps = snps[maf_filter]
        snp_means = snp_means[maf_filter]
        snp_stds = snp_stds[maf_filter]
        if return_snps_info:
            positions = positions[maf_filter]
            snp_ids = snp_ids[maf_filter]
            nts = nts[maf_filter]
            hwe = hwe[maf_filter]
    
    if hwe_filter:
        print 'Filtering SNPs for departure from Hardy-Weinberg proportions with p-val<', 0.05/len(snps)
        hwe_thres = 0.05/len(snps)
        hwe_filter = hwe[()].flatten() > hwe_thres
        snps = snps[hwe_filter]
        snp_means = snp_means[hwe_filter]
        snp_stds = snp_stds[hwe_filter]
        if return_snps_info:
            positions = positions[hwe_filter]
            snp_ids = snp_ids[hwe_filter]
            nts = nts[hwe_filter]
            hwe = hwe[hwe_filter]  
    
    ret_dict = {'snp_means':snp_means, 'snp_stds':snp_stds}
    
    if return_normalized_snps:
        print 'Normalizing SNPs by row (marker)'
        norm_snps = sp.array((snps - snp_means) / snp_stds)
        ret_dict['norm_snps'] = norm_snps
        
        print 'Normalizing SNPs by column (individual)'
        norm_snps_col = sp.array((snps - ind_means) / ind_stds)
        ret_dict['norm_snps_col'] = norm_snps_col
    
    if return_raw_snps:
        ret_dict['snps'] = snps
    
    if return_snps_info:
        ret_dict['positions'] = positions
        ret_dict['snp_ids'] = snp_ids
        ret_dict['nts'] = nts
        ret_dict['hwe'] = hwe

    if unscaled_grm:
        grm = sp.dot(norm_snps.T, norm_snps)
        ret_dict['unscaled_grm'] = grm
        cov_mat = sp.dot(norm_snps_col.T, norm_snps_col)
        ret_dict['unscaled_cov_mat'] = cov_mat
    
    ret_dict['M'] = len(snps)
    
    return ret_dict


def add_leave_one_out_data(chromosome_dict, max_cM_dist, min_r2, n_pcs, calculate_original_ld_scores=False):
    """
    "Leave-one-out data refer to data that were calculated by omitting one chromosome at a time.
    
    These can only be calculated on top of a full chromosome dictionary (i.e. a dictionary that includes all 22 autosomes).
    """
    M, N = sp.shape(chromosome_dict['chr1']['snps'])
    all_chrom_strs = ['chr' + str(a) for a in range(1,23)]
    
    for chrom_str in all_chrom_strs:
        leave_one_out_chrom_strs = all_chrom_strs[:]
        leave_one_out_chrom_strs.remove(chrom_str)
        unscaled_grm_leave_one_out = sp.zeros(shape=(N, N))
        unscaled_cov_mat_leave_one_out = sp.zeros(shape=(N, N))
        total_M = 0
        
        for chrom_str_i in leave_one_out_chrom_strs:
            unscaled_grm_leave_one_out += chromosome_dict[chrom_str_i]['unscaled_grm']
            unscaled_cov_mat_leave_one_out += chromosome_dict[chrom_str_i]['unscaled_cov_mat']
            total_M += chromosome_dict[chrom_str_i]['M']
        
        grm_leave_one_out = unscaled_grm_leave_one_out / (total_M-1)
        quasi_grm_leave_one_out = ((1-0.999999)*unscaled_cov_mat_leave_one_out + 0.999999*unscaled_grm_leave_one_out) / (total_M-1)        
        del unscaled_grm_leave_one_out, unscaled_cov_mat_leave_one_out
        
        quasi_evals_leave_one_out, quasi_evecs_leave_one_out = pca(quasi_grm_leave_one_out)        
        transformed_genotypes_leave_one_out = transform_genotypes(chromosome_dict[chrom_str]['norm_snps'], quasi_grm_leave_one_out)
        transformed_genotypes_leave_one_out = pd.DataFrame(transformed_genotypes_leave_one_out)
        
        cM_map = bp_to_cM(chromosome_dict[chrom_str]['positions'], chrom_str)        
        big_ld_table = calc_ld_table(transformed_genotypes_leave_one_out, cM_map, max_cM_dist=max_cM_dist, min_r2=min_r2, verbose=True, normalize=True) # max_cM_dist and min_r2 are passed through the command line
        ld_scores_leave_one_out = calc_ld_score(big_ld_table)
        if calculate_original_ld_scores:
            big_ld_table_bulik_sullivan = calc_ld_table(pd.DataFrame(chromosome_dict[chrom_str]['norm_snps']), cM_map, max_cM_dist=max_cM_dist, min_r2=min_r2, verbose=True, normalize=True)
            ld_scores_bulik_sullivan = calc_ld_score(big_ld_table_bulik_sullivan)
        pc_scores_leave_one_out = calc_pca_score(chromosome_dict[chrom_str]['norm_snps'], quasi_evecs_leave_one_out, quasi_evals_leave_one_out, n_pcs=n_pcs) # n_pcs is passed through the command line
        
        #Dictionary additions:
        chromosome_dict[chrom_str]['ld_scores_leave_one_out'] = np.array(ld_scores_leave_one_out)
        del big_ld_table, ld_scores_leave_one_out
        if calculate_original_ld_scores:
            chromosome_dict[chrom_str]['ld_scores_bulik_sullivan'] = np.array(ld_scores_bulik_sullivan)
            del big_ld_table_bulik_sullivan, ld_scores_bulik_sullivan

        chromosome_dict[chrom_str]['cM'] = cM_map
        chromosome_dict[chrom_str]['transformed_norm_snps_leave_one_out'] = transformed_genotypes_leave_one_out
        del transformed_genotypes_leave_one_out
        
        chromosome_dict[chrom_str]['quasi_grm_leave_one_out'] = quasi_grm_leave_one_out
        chromosome_dict[chrom_str]['quasi_evals_leave_one_out'] = quasi_evals_leave_one_out
        chromosome_dict[chrom_str]['quasi_evecs_leave_one_out'] = quasi_evecs_leave_one_out
        chromosome_dict[chrom_str]['pc_scores_leave_one_out'] = pc_scores_leave_one_out
        del pc_scores_leave_one_out

