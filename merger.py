import sys
import os
import scipy.stats

current_path = str(os.getcwd())+'/'
pc_ld_scores_name = sys.argv[1]
sum_stats_name = sys.argv[2]
out_file_name = sys.argv[3]

pc_ld_scores = open(current_path+pc_ld_scores_name, 'r')
sum_stats = open(current_path+sum_stats_name, 'r')
sum_stats2 = sum_stats
out_file = open(current_path+out_file_name, 'w')
pc_ld_scores_dict = {}
sum_stats_dict = {}



def qchisq_upper(p,df):
    """
    Calculates chi square statistics from p-values (upper tail).
    """
    return scipy.stats.chi2.ppf((1-p),df)


first_line = sum_stats.readline().split()

contains_pval = True

#find the snp column
if 'MarkerName' in first_line:
    s_index = first_line.index('MarkerName')
elif 'SNP' in first_line:
    s_index = first_line.index('SNP')
elif 'marker' in first_line:
    s_index = first_line.index('marker')
else:
    sys.exit('SNP column in sum stats file not found!')

#find the p-val/chi.sq column
if 'Pval' in first_line:
    p_index = first_line.index('Pval')
elif 'P' in first_line:
    p_index = first_line.index('P')
elif 'p' in first_line:
    p_index = first_line.index('p')
elif 'P_VALUE' in first_line:
    p_index = first_line.index('P_VALUE')
elif 'pvalue' in first_line:
    p_index = first_line.index('pvalue')
elif 'CHISQ' in first_line:
    p_index = first_line.index('CHISQ')
    contains_pval = False
else:
    sys.exit('P-value/chi.sq column in sum stats file not found!')

for line in sum_stats:
    try:
        sum_stats_dict[line.split()[s_index]] = line.split()[p_index] # parse p-val
    except IndexError:
        pass

for line in pc_ld_scores:
    pc_ld_scores_dict[line.split()[0]] = line.split()[1:]

for snp in pc_ld_scores_dict:
    if snp not in sum_stats_dict:
        continue
    if contains_pval:
        p = float(sum_stats_dict[snp])
        chisq = qchisq_upper(p, 1)
        out_file.write(snp+'\t'+'\t'.join(pc_ld_scores_dict[snp])+'\t'+sum_stats_dict[snp]+'\t'+str(chisq)+'\n')
    else:
        out_file.write(snp+'\t'+'\t'.join(pc_ld_scores_dict[snp])+'\t'+sum_stats_dict[snp]+'\n')

out_file.close()
