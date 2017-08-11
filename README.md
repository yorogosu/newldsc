# newldsc
LD score estimates without long-range LD

The code is run with this command (example):

python unbiased_ld_and_pc_scores.py 0.01 CEU_FIN_TSI_IBS_GBR 20 1 0.2 True

There are six arguments:

  (i) MAF for SNP filtering [here 0.01]
 (ii) 1000GP population string (population strings connected with underscores) [here CEU_FIN_TSI_IBS_GBR]
(iii) Number of PC's (k) for PC score. Note to self: it will be more useful to report the PC weights before summing them to the final PC score so that users can sum over their k of choice [here 20]
 (iv) cM window for LD score [here 1]
  (v) whether original Bulik-Sullivan LD scores should be calculated (added for comparison purposes) [here True]
  
  
