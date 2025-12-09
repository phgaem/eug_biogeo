#Project: HISTORICAL BIOGEOGRAPHY AND NICHE EVOLUTION IN EUGENIINAE (MYRTACEAE)
#Author: Paulo Henrique Gaem
#University of Michigan, Ann Arbor

# PERFORMING ANCESTRAL RANGE RECONSTRUCTION WITH BIOGEOBEARS

# Setting working directory
setwd("~/Documents/Projetos/3_PhD/1_Project/Eugenia Biogeography")

#loading packages

#library(devtools)
#devtools::install_github(repo="nmatzke/BioGeoBEARS")

library(ape)
library(BioGeoBEARS)
library(cladoRcpp)
library(parallel)
library(dplyr)
library(tidyr)

#location of phylogenetic tree
trfn <- 'data/2_eugtree.tre'
moref(trfn) #checking the file

tr = read.tree(trfn)

#location of the pres-abs file
geogfn <- 'data/16_biogeobears_presabs.txt'
moref(geogfn) #checking the file

tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Check the maximum range size in your geographical data
max(rowSums(dfnums_to_numeric(tipranges@df)))

# Set the maximum range size based on your geographical data
max_range_size <- max(rowSums(dfnums_to_numeric(tipranges@df)))

# Check number of states which will influence the total analysis time (if more than 500-600 calculations will get really slow)
numstates_from_numareas(numareas=10, maxareas=3, include_null_range=T)

disp_multipl <- c('data/19_dispers_multipl_neutral.txt','data/17_dispers_multipl_east.txt',
                  'data/18_dispers_multipl_west.txt')


###################################################################
# Running BioGeoBEARS with the DEC model - unconstrained model + two hypotheses
# for details, check http://phylo.wikidot.com/biogeobears#script

rdatafiles.DEC <- c('results/1_Eugenia_DEC-Unconst.Rdata', 'results/2_Eugenia_DEC-East.Rdata',
                    'results/3_Eugenia_DEC-West.Rdata')

pdffiles.DEC <- c('results/1_Eug_DEC-Unconstr.pdf', 'results/2_Eug_DEC-East.pdf',
                  'results/3_Eug_DEC-West.pdf')

pdftitles.DEC <- c('BioGeoBEARS DEC (Unconstrained) on Eugenia',
                   'BioGeoBEARS DEC (Eastbound) on Eugenia', 'BioGeoBEARS DEC (Westbound) on Eugenia')

for (i in 1:length(disp_multipl)) {
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001
  BioGeoBEARS_run_object$include_null_range = TRUE
  
  # Use dispersal multiplier from the vector
  BioGeoBEARS_run_object$dispersal_multipliers_fn = disp_multipl[i]
  
  BioGeoBEARS_run_object$on_NaN_error = -1e50
  BioGeoBEARS_run_object$speedup = TRUE     
  BioGeoBEARS_run_object$use_optimx = "GenSA"   
  BioGeoBEARS_run_object$num_cores_to_use = 11 # Adjust according to computer specifications
  BioGeoBEARS_run_object$force_sparse = FALSE
  
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE
  
  BioGeoBEARS_run_object
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table
  check_BioGeoBEARS_run(BioGeoBEARS_run_object) # Check if run is ready -- should return "TRUE"
  
  # Use result file name from the vector
  resfn = rdatafiles.DEC[i]
  
  # Run DEC
  res = bears_optim_run(BioGeoBEARS_run_object)
  save(res, file=resfn) # Save results
  resDEC = res
  
  # Plot PDF of reconstruction
  pdffn = pdffiles.DEC[i]
  pdf(pdffn, width=10, height=40)
  
  # Use title from the vector
  analysis_titletxt = pdftitles.DEC[i]
  
  # Setup
  results_object = resDEC
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, plotwhat="text",
                           label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6,
                           titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir,
                           include_null_range=TRUE, tipranges=tipranges)
  
  # Pie chart
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, plotwhat="pie",
                           label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6,
                           titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir,
                           include_null_range=TRUE, tipranges=tipranges)
  
  dev.off()
}




################################################################
# Running BioGeoBEARS with the DEC+J model - unconstrained model + two hypotheses
# for details, check http://phylo.wikidot.com/biogeobears#script

load.DEC <- c('results/1_Eugenia_DEC-Unconst.Rdata', 'results/2_Eugenia_DEC-East.Rdata',
              'results/3_Eugenia_DEC-West.Rdata')

rdatafiles.DECJ <- c('results/4_Eugenia_DECJ-Unconst.Rdata', 'results/5_Eugenia_DECJ-East.Rdata',
                    'results/6_Eugenia_DECJ-West.Rdata')

pdffiles.DECJ <- c('results/4_Eug_DECJ-Unconstr.pdf', 'results/5_Eug_DECJ-East.pdf',
                  'results/6_Eug_DECJ-West.pdf')

pdftitles.DECJ <- c('BioGeoBEARS DEC+J (Unconstrained) on Eugenia',
                   'BioGeoBEARS DEC+J (Eastbound) on Eugenia', 'BioGeoBEARS DEC+J (Westbound) on Eugenia')

for (i in 1:length(disp_multipl)) {
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001
  BioGeoBEARS_run_object$include_null_range = TRUE
  
  BioGeoBEARS_run_object$dispersal_multipliers_fn = disp_multipl[i]
  
  BioGeoBEARS_run_object$on_NaN_error = -1e50    
  BioGeoBEARS_run_object$speedup = TRUE
  BioGeoBEARS_run_object$use_optimx = TRUE
  
  BioGeoBEARS_run_object$num_cores_to_use = 10
  BioGeoBEARS_run_object$force_sparse = FALSE
  
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
  
  # Set up DEC+J model
  # Get the ML parameter values from the 2-parameter nested model
  # (this will ensure that the 3-parameter model always does at least as good)
  
  load(paste(load.DEC[i]))
  
  dstart = res$outputs@params_table["d","est"]
  estart = res$outputs@params_table["e","est"]
  jstart = 0.0001
  
  # Input starting values for d, e
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  # Add j as a free parameter
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Name file in which to save the run
  resfn = rdatafiles.DECJ[i]
  
  # Run DEC+J
  res = bears_optim_run(BioGeoBEARS_run_object)
  save(res, file=resfn) # Save results
  resDEC = res
  
  # Plot PDF of reconstruction
  pdffn = pdffiles.DECJ[i]
  pdf(pdffn, width=10, height=40)
  
  # Plot ancestral states - DEC
  analysis_titletxt = pdftitles.DECJ[i]
  
  # Setup
  results_object = resDEC
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, plotwhat="text",
                           label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6,
                           titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir,
                           include_null_range=TRUE, tipranges=tipranges)
  
  # Pie chart
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, plotwhat="pie",
                           label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6,
                           titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir,
                           include_null_range=TRUE, tipranges=tipranges)
  
  dev.off()
}



###################################################################
# Running BioGeoBEARS with the DIVALIKE model - unconstrained model + two hypotheses
# for details, check http://phylo.wikidot.com/biogeobears#script

rdatafiles.DIVALIKE <- c('results/7_Eugenia_DIVALIKE-Unconst.Rdata',
                         'results/8_Eugenia_DIVALIKE-East.Rdata',
                         'results/9_Eugenia_DIVALIKE-West.Rdata')

pdffiles.DIVALIKE <- c('results/7_Eug_DIVALIKE-Unconstr.pdf', 'results/8_Eug_DIVALIKE-East.pdf',
                       'results/9_Eug_DIVALIKE-West.pdf')

pdftitles.DIVALIKE <- c('BioGeoBEARS DIVALIKE (Unconstrained) on Eugenia',
                        'BioGeoBEARS DIVALIKE (Eastbound) on Eugenia',
                        'BioGeoBEARS DIVALIKE (Westbound) on Eugenia')

for (i in 1:length(disp_multipl)) {
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001
  BioGeoBEARS_run_object$include_null_range = TRUE
  
  # Use dispersal multiplier from the vector
  BioGeoBEARS_run_object$dispersal_multipliers_fn = disp_multipl[i]
  
  BioGeoBEARS_run_object$on_NaN_error = -1e50
  BioGeoBEARS_run_object$speedup = TRUE     
  BioGeoBEARS_run_object$use_optimx = "GenSA"   
  BioGeoBEARS_run_object$num_cores_to_use = 11 # Adjust according to computer specifications
  BioGeoBEARS_run_object$force_sparse = FALSE
  
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE
  
  # Set up DIVALIKE model
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
  
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Use result file name from the vector
  resfn = rdatafiles.DIVALIKE[i]
  
  # Run DIVALIKE
  res = bears_optim_run(BioGeoBEARS_run_object)
  save(res, file=resfn) # Save results
  resDIVALIKE <- res
  
  # Plot PDF of reconstruction
  pdffn = pdffiles.DIVALIKE[i]
  pdf(pdffn, width=10, height=40)
  
  # Use title from the vector
  analysis_titletxt = pdftitles.DIVALIKE[i]
  
  # Setup
  results_object = resDIVALIKE
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt,
                                  addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7,
                                  statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir,
                                  include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  # Pie chart
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
                           label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                           cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  dev.off()
}



###################################################################
# Running BioGeoBEARS with the DIVALIKE+J model - unconstrained model + two hypotheses
# for details, check http://phylo.wikidot.com/biogeobears#script

load.DIVALIKE <- c('results/7_Eugenia_DIVALIKE-Unconst.Rdata', 'results/8_Eugenia_DIVALIKE-East.Rdata',
              'results/9_Eugenia_DIVALIKE-West.Rdata')

rdatafiles.DIVALIKEJ <- c('results/10_Eugenia_DIVALIKEJ-Unconst.Rdata', 'results/11_Eugenia_DIVALIKEJ-East.Rdata',
                     'results/12_Eugenia_DIVALIKEJ-West.Rdata')

pdffiles.DIVALIKEJ <- c('results/10_Eug_DIVALIKEJ-Unconstr.pdf', 'results/11_Eug_DIVALIKEJ-East.pdf',
                   'results/12_Eug_DIVALIKEJ-West.pdf')

pdftitles.DIVALIKEJ <- c('BioGeoBEARS DIVALIKE+J (Unconstrained) on Eugenia',
                    'BioGeoBEARS DIVALIKE+J (Eastbound) on Eugenia', 'BioGeoBEARS DIVALIKE+J (Westbound) on Eugenia')


for (i in 1:length(disp_multipl)) {
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001
  BioGeoBEARS_run_object$include_null_range = TRUE
  
  # Use dispersal multiplier from the vector
  BioGeoBEARS_run_object$dispersal_multipliers_fn = disp_multipl[i]
  
  BioGeoBEARS_run_object$on_NaN_error = -1e50
  BioGeoBEARS_run_object$speedup = TRUE     
  BioGeoBEARS_run_object$use_optimx = "GenSA"   
  BioGeoBEARS_run_object$num_cores_to_use = 11 # Adjust according to computer specifications
  BioGeoBEARS_run_object$force_sparse = FALSE
  
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE
  
  # Set up DIVALIKE+J model
  
  load(paste(load.DIVALIKE[i]))
  
  dstart = res$outputs@params_table["d","est"]
  estart = res$outputs@params_table["e","est"]
  jstart = 0.0001
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999
  
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object)
  
  # Use result file name from the vector
  resfn = rdatafiles.DIVALIKEJ[i]
  
  # Run DIVALIKE
  res = bears_optim_run(BioGeoBEARS_run_object)
  save(res, file=resfn) # Save results
  resDIVALIKEJ <- res
  
  # Plot PDF of reconstruction
  pdffn = pdffiles.DIVALIKEJ[i]
  pdf(pdffn, width=10, height=40)
  
  # Use title from the vector
  analysis_titletxt = pdftitles.DIVALIKEJ[i]
  
  # Setup
  results_object = resDIVALIKEJ
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt,
                                  addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7,
                                  statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir,
                                  include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  # Pie chart
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
                           label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                           cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  dev.off()
  
}




###################################################################
# Running BioGeoBEARS with the BAYAREALIKE model - unconstrained model + two hypotheses
# for details, check http://phylo.wikidot.com/biogeobears#script

rdatafiles.BAYAREALIKE <- c('results/13_Eugenia_BAYAREALIKE-Unconst.Rdata', 'results/14_Eugenia_BAYAREALIKE-East.Rdata',
                    'results/15_Eugenia_BAYAREALIKE-West.Rdata')

pdffiles.BAYAREALIKE <- c('results/13_Eug_BAYAREALIKE-Unconstr.pdf', 'results/14_Eug_BAYAREALIKE-East.pdf',
                  'results/15_Eug_BAYAREALIKE-West.pdf')

pdftitles.BAYAREALIKE <- c('BioGeoBEARS BAYAREALIKE (Unconstrained) on Eugenia',
                   'BioGeoBEARS BAYAREALIKE (Eastbound) on Eugenia', 'BioGeoBEARS BAYAREALIKE (Westbound) on Eugenia')

for (i in 1:length(disp_multipl)) {
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001
  BioGeoBEARS_run_object$include_null_range = TRUE
  
  # Use dispersal multiplier from the vector
  BioGeoBEARS_run_object$dispersal_multipliers_fn = disp_multipl[i]
  
  BioGeoBEARS_run_object$on_NaN_error = -1e50
  BioGeoBEARS_run_object$speedup = TRUE     
  BioGeoBEARS_run_object$use_optimx = "GenSA"   
  BioGeoBEARS_run_object$num_cores_to_use = 11 # Adjust according to computer specifications
  BioGeoBEARS_run_object$force_sparse = FALSE
  
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE
  
  # Set up BAYAREALIKE model

  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
  
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
  # BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object) # Check if run is ready -- should return "TRUE"
  
  # Use result file name from the vector
  resfn = rdatafiles.BAYAREALIKE[i]
  
  # Run DEC
  res = bears_optim_run(BioGeoBEARS_run_object)
  save(res, file=resfn) # Save results
  resBAYAREALIKE = res
  
  # Plot PDF of reconstruction
  pdffn = pdffiles.BAYAREALIKE[i]
  pdf(pdffn, width=10, height=40)
  
  # Use title from the vector
  analysis_titletxt = pdftitles.BAYAREALIKE[i]
  
  # Setup
  results_object = resBAYAREALIKE
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt,
        addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7,
        statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir,
        include_null_range=TRUE, tr=tr, tipranges=tipranges)
  
  # Pie chart
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"),
        plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6,
        titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE,
        tr=tr, tipranges=tipranges)
  
  dev.off()
}



###################################################################
# Running BioGeoBEARS with the BAYAREALIKE+J model - unconstrained model + two hypotheses
# for details, check http://phylo.wikidot.com/biogeobears#script

load.BAYAREALIKE <- c('results/13_Eugenia_BAYAREALIKE-Unconst.Rdata', 'results/14_Eugenia_BAYAREALIKE-East.Rdata',
                      'results/15_Eugenia_BAYAREALIKE-West.Rdata')

rdatafiles.BAYAREALIKEJ <- c('results/16_Eugenia_BAYAREALIKEJ-Unconst.Rdata', 'results/17_Eugenia_BAYAREALIKEJ-East.Rdata',
                            'results/18_Eugenia_BAYAREALIKEJ-West.Rdata')

pdffiles.BAYAREALIKEJ <- c('results/16_Eug_BAYAREALIKEJ-Unconstr.pdf', 'results/17_Eug_BAYAREALIKEJ-East.pdf',
                          'results/18_Eug_BAYAREALIKEJ-West.pdf')

pdftitles.BAYAREALIKEJ <- c('BioGeoBEARS BAYAREALIKE+J (Unconstrained) on Eugenia',
                           'BioGeoBEARS BAYAREALIKE+J (Eastbound) on Eugenia', 'BioGeoBEARS BAYAREALIKE+J (Westbound) on Eugenia')

for (i in 1:length(disp_multipl)) {
  BioGeoBEARS_run_object = define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn = trfn
  BioGeoBEARS_run_object$geogfn = geogfn
  BioGeoBEARS_run_object$max_range_size = max_range_size
  BioGeoBEARS_run_object$min_branchlength = 0.000001
  BioGeoBEARS_run_object$include_null_range = TRUE
  
  # Use dispersal multiplier from the vector
  BioGeoBEARS_run_object$dispersal_multipliers_fn = disp_multipl[i]
  
  BioGeoBEARS_run_object$on_NaN_error = -1e50
  BioGeoBEARS_run_object$speedup = TRUE     
  BioGeoBEARS_run_object$use_optimx = "GenSA"   
  BioGeoBEARS_run_object$num_cores_to_use = 11 # Adjust according to computer specifications
  BioGeoBEARS_run_object$force_sparse = FALSE
  
  BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
  BioGeoBEARS_run_object$return_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  BioGeoBEARS_run_object$calc_ancprobs = TRUE
  
  # Set up BAYAREALIKE+J model
  
  load(paste(load.BAYAREALIKE[i]))
  
  dstart = res$outputs@params_table["d","est"]
  estart = res$outputs@params_table["e","est"]
  jstart = 0.0001
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999
  
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
  BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999
  
  BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
  check_BioGeoBEARS_run(BioGeoBEARS_run_object) # Check if run is ready -- should return "TRUE"
  
  # Use result file name from the vector
  resfn = rdatafiles.BAYAREALIKEJ[i]
  
  # Run BAYAREALIKE+J
  res = bears_optim_run(BioGeoBEARS_run_object)
  save(res, file=resfn) # Save results
  resBAYAREALIKEJ = res
  
  # Plot PDF of reconstruction
  pdffn = pdffiles.BAYAREALIKEJ[i]
  pdf(pdffn, width=10, height=40)
  
  # Use title from the vector
  analysis_titletxt = pdftitles.BAYAREALIKEJ[i]
  
  # Setup
  results_object = resBAYAREALIKEJ
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  # States
  res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"),
                                  plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7,
                                  splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                                  cornercoords_loc=scriptdir, include_null_range=TRUE,
                                  tr=tr, tipranges=tipranges)
  
  # Pie chart
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
                           label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8,
                           plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr,
                           tipranges=tipranges)
  
  dev.off()
}

###################################################################
# Comparing results. First, testing the contribution of the J parameter for the three models and across hypotheses

# Set up empty tables to hold the statistical results
restable = NULL
teststable = NULL

# DEC vs. DEC+J

load.DEC <- c('results/1_Eugenia_DEC-Unconst.Rdata', 'results/2_Eugenia_DEC-East.Rdata',
              'results/3_Eugenia_DEC-West.Rdata')

load.DECJ <- c('results/4_Eugenia_DECJ-Unconst.Rdata', 'results/5_Eugenia_DECJ-East.Rdata',
              'results/6_Eugenia_DECJ-West.Rdata')


for (i in 1:length(load.DEC)) {
  
  load(paste(load.DEC[i]))
  resDEC <- res
  
  load(paste(load.DECJ[i]))
  resDECJ <- res
  
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECJ)
  
  numparams1 = 3
  numparams2 = 2
  
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDECJ, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  
  rbind(res2, res1)
  tmp_tests = conditional_format_table(stats)
  
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
}

# DIVALIKE vs. DIVALIKE+J

load.DIVALIKE <- c('results/7_Eugenia_DIVALIKE-Unconst.Rdata', 'results/8_Eugenia_DIVALIKE-East.Rdata',
              'results/9_Eugenia_DIVALIKE-West.Rdata')

load.DIVALIKEJ <- c('results/10_Eugenia_DIVALIKEJ-Unconst.Rdata', 'results/11_Eugenia_DIVALIKEJ-East.Rdata',
               'results/12_Eugenia_DIVALIKEJ-West.Rdata')


for (i in 1:length(load.DIVALIKE)) {
  
  load(paste(load.DIVALIKE[i]))
  resDIVALIKE <- res
  
  load(paste(load.DIVALIKEJ[i]))
  resDIVALIKEJ <- res
  
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEJ)
  
  numparams1 = 3
  numparams2 = 2
  
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEJ, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  
  rbind(res2, res1)
  tmp_tests = conditional_format_table(stats)
  
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
}

# DIVALIKE vs. DIVALIKE+J

load.BAYAREALIKE <- c('results/13_Eugenia_BAYAREALIKE-Unconst.Rdata', 'results/14_Eugenia_BAYAREALIKE-East.Rdata',
                   'results/15_Eugenia_BAYAREALIKE-West.Rdata')

load.BAYAREALIKEJ <- c('results/16_Eugenia_BAYAREALIKEJ-Unconst.Rdata', 'results/17_Eugenia_BAYAREALIKEJ-East.Rdata',
                    'results/18_Eugenia_BAYAREALIKEJ-West.Rdata')


for (i in 1:length(load.BAYAREALIKE)) {
  
  load(paste(load.BAYAREALIKE[i]))
  resBAYAREALIKE <- res
  
  load(paste(load.BAYAREALIKEJ[i]))
  resBAYAREALIKEJ <- res
  
  LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
  LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEJ)
  
  numparams1 = 3
  numparams2 = 2
  
  stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
  
  res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEJ, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  
  rbind(res2, res1)
  tmp_tests = conditional_format_table(stats)
  
  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)
}

teststable$alt = c("DEC+J Unconstr.", "DEC+J Eastbound", "DEC+J Westbound",
                   "DIVALIKE+J Unconstr.", "DIVALIKE+J Eastbound", "DIVALIKE+J Westbound", 
                   "BAYAREALIKE+J Unconstr.", "BAYAREALIKE+J Eastbound", "BAYAREALIKE+J Westbound")

teststable$null = c("DEC Unconstr.", "DEC Eastbound", "DEC Westbound",
                    "DIVALIKE Unconstr.", "DIVALIKE Eastbound", "DIVALIKE Westbound", 
                    "BAYAREALIKE Unconstr.", "BAYAREALIKE Eastbound", "BAYAREALIKE Westbound")

#saving testable
save(teststable, file='results/19_teststable.RData')

teststable <- teststable %>% unnest(cols = where(is.list))
write.csv(teststable, file='results/20_teststable.csv', row.names=F)

#prepping restable (Model weights of all models)
restable$model <- c("DEC Unconstr.", "DEC+J Unconstr.", "DEC Eastbound", "DEC+J Eastbound", "DEC Westbound","DEC+J Westbound",
                        "DIVALIKE Unconstr.", "DIVALIKE+J Unconstr.", "DIVALIKE Eastbound", "DIVALIKE+J Eastbound", "DIVALIKE Westbound","DIVALIKE+J Westbound",
                        "BAYAREALIKE Unconstr.", "BAYAREALIKE+J Unconstr.", "BAYAREALIKE Eastbound", "BAYAREALIKE+J Eastbound", "BAYAREALIKE Westbound","BAYAREALIKE+J Westbound")

restable <- restable %>% select(model, LnL, numparams, d, e, j)

# With AICs:
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
restable <- restable_AIC_rellike

# With AICcs -- factors in sample size
samplesize = length(tr$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable = cbind(restable, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike

restable$AICc_wt <- restable_AICc_rellike$AICc_wt

#saving restable
save(restable, file='results/21_restable.RData')
restable <- restable %>% unnest(cols = where(is.list))
write.csv(conditional_format_table(restable), file='results/22_restable.csv', row.names=F)