# newldsc
This repository contains (A) code for LD score estimates *without long-range LD* & pc scores (unbiased_ld_and_pc_scores.py); (B) code for merging the resulting LD score files with a given set of summary statistics (merger.py); and (C) code for plotting LD score regression and outputting a table with intercepts (plotter.R)

A. LD score & PC score estimates

The code is run with this command (example):

python unbiased_ld_and_pc_scores.py 0.01 CEU_FIN_TSI_IBS_GBR 20 1 0.2 True

There are six arguments:

  (i) MAF for SNP filtering [here 0.01]
  
 (ii) 1000GP population string (population strings connected with underscores) [here CEU_FIN_TSI_IBS_GBR]
 
(iii) Number of PC's (k) for PC score. Note to self: it will be more useful to report the PC weights before summing them to the final PC score so that users can sum over their k of choice [here 20]

 (iv) cM window for LD score [here 1]
 
  (v) r-squared threshold for LD pruning withing the window (here 0.2)
 
  (vi) whether original Bulik-Sullivan LD scores should be calculated (added for comparison purposes) [here True]
  
  
B. Data merging

The code is run with this command (example):

python merger.py scores_0.01_CEU_TSI_IBS_GBR_20pcs_1.0cM_0.2_minr2.txt GIANT_Yang2012Nature_publicrelease_HapMapCeuFreq_BMI.txt outfile.txt

There are three arguments:

  (i) The LD and PC scores from a previous run [here scores_0.01_CEU_TSI_IBS_GBR_20pcs_1.0cM_0.2_minr2.txt]
  
 (ii) A file with summary statistics [here GIANT_Yang2012Nature_publicrelease_HapMapCeuFreq_BMI.txt]
 
(iii) An output filename of free choice [here outfile.txt]

C. Plotting

In order to make plot titles more understandable, I created a dictionary with a phenotype description for each sum stat file (dummy_dictionary.txt). This R code is a bit rough and needs refining. It runs like this:

R < plotter.R --no-save

but ideally we want something like this:

Rscript plotter.R args[1] args[2] args[3] etc ...
