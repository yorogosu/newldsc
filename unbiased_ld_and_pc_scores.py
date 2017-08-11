import numpy as np
import h5py
import sys
from ld_util import get_bulik_sullivan_15_sids
from data_util import get_genotype_data, add_leave_one_out_data

maf_arg = sys.argv[1]
maf_arg = float(maf_arg)
pop_arg = sys.argv[2]
pop_arg = pop_arg.split("_")
n_pc_arg = sys.argv[3]
n_pc_arg = int(n_pc_arg)
cM_arg = sys.argv[4]
cM_arg = float(cM_arg)
r2_arg = sys.argv[5]
r2_arg = float(r2_arg)
bs_arg = sys.argv[6]
bs_arg = bool(bs_arg)

def run_code(KG_path='/home/yorgos/HeritPartition/faststorage/1Kgenomes_bjarni/phase3/',
                      out_path='/home/yorgos/HeritPartition/faststorage/pypcma/',
                      maf_thres=0.05, pop_list=['CEU', 'FIN', 'TSI', 'IBS', 'GBR'],
                      n_pcs=20, max_cM_dist=1, min_r2=0.2, calculate_original_ld_scores=False):
    """
    Write something.
    """
    print maf_thres, pop_list, n_pcs, max_cM_dist, min_r2
    kg_h5f = h5py.File(KG_path+'1k_genomes_unrelated.hdf5', 'r')
    chrom_ok_snp_dict = get_bulik_sullivan_15_sids() # a list of SNPs that were included in Bulik-Sullivan et al.
    
    n_indivs = len(kg_h5f['indivs']['ancestry'][...])
    indiv_index = np.zeros(n_indivs, dtype=bool)
    for pop in pop_list:
        indiv_index += kg_h5f['indivs']['ancestry'][...] == pop
    indiv_index = np.where(indiv_index==True)
    indiv_index = indiv_index[0] # to get rid of the tuple structure
    
    chromosome_dict = {}
    all_chroms = range(1,23)
    for chrom in all_chroms:
        genotype_data = get_genotype_data(kg_h5f, chrom, chrom_ok_snp_dict=chrom_ok_snp_dict, maf_thres=maf_thres, indiv_filter=indiv_index)
        chromosome_dict['chr'+str(chrom)] = genotype_data
    
    add_leave_one_out_data(chromosome_dict, max_cM_dist, min_r2, n_pcs, calculate_original_ld_scores)
        
    out_file = open(out_path+'scores_'+str(maf_thres)+'_'+'_'.join(pop_list)+'_'
                    +str(n_pcs)+'pcs_'+str(max_cM_dist)+'cM_'
                    +str(min_r2)+'_minr2.txt', "w") #creates a name for the output file based on the chosen options
    
    if calculate_original_ld_scores:
        out_file.write("rsid\tchr\tposition\tcM\tbs.ld.score\tnew.ld.score\tpca.score\n")
        for chrom in all_chroms:
            chrom_str = 'chr'+str(chrom)
            for i in range(len(chromosome_dict[chrom_str]['snp_ids'])):
                out_file.write(str(chromosome_dict[chrom_str]['snp_ids'][i])+'\t'+str(chrom)+'\t'
                               +str(chromosome_dict[chrom_str]['positions'][i])+'\t'
                               +str(chromosome_dict[chrom_str]['cM'][i])+'\t'
                               +str(chromosome_dict[chrom_str]['ld_scores_bulik_sullivan'][i])+'\t'
                               +str(chromosome_dict[chrom_str]['ld_scores_leave_one_out'][i])+'\t'
                               +str(chromosome_dict[chrom_str]['pc_scores_leave_one_out'][i])+'\n')
    else:
        out_file.write("rsid\tchr\tposition\tcM\tnew.ld.score\tpca.score\n")
        for chrom in all_chroms:
            chrom_str = 'chr'+str(chrom)
            for i in range(len(chromosome_dict[chrom_str]['snp_ids'])):
                out_file.write(str(chromosome_dict[chrom_str]['snp_ids'][i])+'\t'+str(chrom)+'\t'
                               +str(chromosome_dict[chrom_str]['positions'][i])+'\t'
                               +str(chromosome_dict[chrom_str]['cM'][i])+'\t'
                               +str(chromosome_dict[chrom_str]['ld_scores_leave_one_out'][i])+'\t'
                               +str(chromosome_dict[chrom_str]['pc_scores_leave_one_out'][i])+'\n')       
    out_file.close()
    
if __name__ == '__main__':
    run_code(maf_thres=maf_arg, pop_list=pop_arg, n_pcs=n_pc_arg, max_cM_dist=cM_arg, min_r2=r2_arg, calculate_original_ld_scores=bs_arg)
