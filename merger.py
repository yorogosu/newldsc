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
    Calculates chi square statistics from p-values (upper tail) by use of the inverse survival function (isf).
    """
    return scipy.stats.chi2.isf(p,df)


first_line = sum_stats.readline().split()
print first_line
contains_pval = True

#find the snp column
if 'MarkerName' in first_line:
    s_index = first_line.index('MarkerName')
elif 'Markername' in first_line:
    s_index = first_line.index('Markername')
elif 'MARKERNAME' in first_line:
    s_index = first_line.index('MARKERNAME')
elif 'SNP' in first_line:
    s_index = first_line.index('SNP')
elif 'marker' in first_line:
    s_index = first_line.index('marker')
elif 'RSID' in first_line:
    s_index = first_line.index('RSID')
elif 'rs' in first_line:
    s_index = first_line.index('rs')
elif 'rs_number' in first_line:
    s_index = first_line.index('rs_number')
elif 'rsid' in first_line:
    s_index = first_line.index('rsid')
elif 'rsID' in first_line:
    s_index = first_line.index('rsID')
elif 'ID' in first_line:
    s_index = first_line.index('ID')
elif 'snp' in first_line:
    s_index = first_line.index('snp')
elif 'Snp' in first_line:
    s_index = first_line.index('Snp')
elif 'SNPID' in first_line:
    s_index = first_line.index('SNPID')
elif 'SNP_ID' in first_line:
    s_index = first_line.index('SNP_ID')
elif 'snpid' in first_line:
    s_index = first_line.index('snpid')
elif 'RSNUMBER' in first_line:
    s_index = first_line.index('RSNUMBER')
else:
    sys.exit('SNP column in sum stats file not found!')
print 'The s_index is: ', s_index
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
elif 'PVAL' in first_line:
    p_index = first_line.index('PVAL')
elif 'P_Value' in first_line:
    p_index = first_line.index('P_Value')
elif 'P_fix' in first_line:
    p_index = first_line.index('P_fix')
elif 'p_sanger' in first_line:
    p_index = first_line.index('p_sanger')
elif 'p-value' in first_line:
    p_index = first_line.index('p-value')
elif 'P-value' in first_line:
    p_index = first_line.index('P-value')
elif 'p_gc' in first_line:
    p_index = first_line.index('p_gc')
elif 'Pvalue' in first_line:
    p_index = first_line.index('Pvalue')
elif 'P_VALUE_GCADJ' in first_line:
    p_index = first_line.index('P_VALUE_GCADJ')
elif 'MainP' in first_line:
    p_index = first_line.index('MainP')
elif 'P-val' in first_line:
    p_index = first_line.index('P-val')
elif 'pval' in first_line:
    p_index = first_line.index('pval')
elif 'GWAS_P' in first_line:
    p_index = first_line.index('GWAS_P')
elif 'PVALUE' in first_line:
    p_index = first_line.index('PVALUE')
elif 'p.value' in first_line:
    p_index = first_line.index('p.value')
elif 'P.value' in first_line:
    p_index = first_line.index('P.value')
elif 'P_fathers_age_death' in first_line:
    p_index = first_line.index('P_fathers_age_death')
elif 'P_mothers_age_death' in first_line:
    p_index = first_line.index('P_mothers_age_death')
elif 'P_parents_age_death' in first_line:
    p_index = first_line.index('P_parents_age_death')
elif 'P_top_1_percent' in first_line:
    p_index = first_line.index('P_top_1_percent')
elif 'CHISQ' in first_line:
    p_index = first_line.index('CHISQ')
    contains_pval = False
else:
    sys.exit('P-value/chi.sq column in sum stats file not found!')
print 'The p_index is: ', p_index
for line in sum_stats:
    try:
        sum_stats_dict[line.split()[s_index]] = line.split()[p_index] # parse p-val
    except IndexError:
        pass

for line in pc_ld_scores:
    pc_ld_scores_dict[line.split()[0]] = line.split()[1:]

out_file.write("rsid\tchr\tposition\tcM\tbs.ld.score\tnew.ld.score\tpca.score\tchi.sq\n")
for snp in pc_ld_scores_dict:
    if snp not in sum_stats_dict:
        continue
    if contains_pval:
        p = float(sum_stats_dict[snp])
        chisq = qchisq_upper(p, 1)
        out_file.write(snp+'\t'+'\t'.join(pc_ld_scores_dict[snp])+'\t'+str(chisq)+'\n')
    else:
        out_file.write(snp+'\t'+'\t'.join(pc_ld_scores_dict[snp])+'\t'+sum_stats_dict[snp]+'\n')

out_file.close()
