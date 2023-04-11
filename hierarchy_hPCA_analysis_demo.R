# Demo of hierarchical PCA analyses on simulated dataset,
#   using fns. in 'hierarchy_hPCA_analysis_fns.R', 'PCA_fns.R', & 'PCA_dataCompression_fns.R'
#


#------------------------------------
### Required functions for script ###
#------------------------------------
#   -PCA scripts from https://github.com/koreywylie/PCA.Rachakonda2016
#   -hPCA script from https://github.com/koreywylie/HierarchicalPrinCompAnalyses
#   -all scripts must be located in current dir.
script_fns = list()
script_fns[[length(script_fns)+1]] = 'PCA_fns.R'                  # PCA methods, optimized for big data & parallelizable
script_fns[[length(script_fns)+1]] = 'PCA_dataCompression_fns.R'  # ~GIFT toolbox data compression, subj. & group PCA, multi-model order ICA, etc.
script_fns[[length(script_fns)+1]] = 'hierarchy_hPCA_analysis_fns.R'  # treelet hPCA applied to temporally-concat voxels, back-reconstruction, etc.


#----------------------------------
### Load Exp.-specific params. ###
#----------------------------------
exp_dir = getwd()

exp_prefix = 'SimHier-'



#---------------------------
### Data Sources ###
#---------------------------
data_dir = 'Hierarchy'
# data_path = file.path(exp_dir, data_dir)
data_path = '/data/analysis/korey/SimTB_hierarchy/Hier-3Levels_3Degree'  # 4/11/2023 --kw-- testing...

data_wc = paste0('*_subject_.*_DATA.nii')  # wildcard for data generated with SimTB

hier_dir = data_dir
hier_path = file.path(exp_dir, hier_dir)

hier_wc = '_Hierarchy.RData'


#-------------------------
### Additional params. ###
#-------------------------
### Masks ###
mask_file = paste0(exp_prefix, '_MASK.nii')
mask_path = file.path(data_path, mask_file)    # path to SimTB mask volume
stopifnot(file.exists(mask_path))

pca_dir = hier_dir
pca_path = file.path(exp_dir, pca_dir)

pca_subj_wc = '_PCA1.RData'
pca_group_wc = '_groupPCA.RData'

### PCA parameters ###
pca_params = list('numOfPC1' = 300,    # to match application to HCP data
                  'numOfPC2' = 300,
                  'preproc_type' = 'variance.norm'
                  )


#---------------------------
### Save Info ###
#---------------------------
save_path = hier_path


#----------------------
### Preliminaries ###
#----------------------
### Load required functions ###
script_fns = file.path(exp_dir, script_fns)
if(!all(file.exists(script_fns))){  print(script_fns)   }
stopifnot(all(file.exists(script_fns)))
for (s in 1:length(script_fns)){ source(script_fns[[s]]) }

### Create Save dirs. ###
if (!all(dir.exists(save_path))){
  for (p in save_path){
    dir.create(p, recursive=T, showWarnings=F)
  }
}



#----------------------------------------------
### Get paths to subj. volumes ###
#----------------------------------------------
data_files = list.files(data_path, data_wc, recursive = T, full.names = T)

cat(paste0('\nAnalyzed subjs.:\n'))
print(data_files)
cat('\n\n')




######################################
###   Initial data compressions   ###
######################################
### Subj.-level PCA ###
subj_PCA(data_files,
         K1           = pca_params$numOfPC1,
         preproc.type = pca_params$preproc_type,
         discard.n    = NA,
         whitenPCs    = F,  # strongly recommend against whitening, results in large & small eigenvalues treated equally
         mask         = mask_path,
         prefix       = exp_prefix,
         save_path    = pca_path,
         verbose = T)
pca1_files = list.files(pca_path, pca_subj_wc, full.names=T)
stopifnot(length(pca1_files) > 0)
stopifnot(length(pca1_files) == length(data_files))


### Group-level PCA ###
group_PCA(pca1_files,
          K2           = pca_params$numOfPC2,
          prefix       = exp_prefix,
          whitenPCs    = F,  # strongly recommend against whitening, results in large & small eigenvalues treated equally
          save_path    = pca_path,
          verbose = T)
pca2_files = list.files(pca_path, pca_group_wc, full.names=T)
stopifnot(length(pca2_files) > 0)





##############################################
## Group-level hierarchical PCA on voxels ###
##############################################
hier_input_files = pca2_files

group_HierPCA(hier_input_files,
              prefix = exp_prefix,
              save_path = hier_path,
              verbose = T)
hier_file = list.files(hier_path, hier_wc, full.names=T)
stopifnot(length(hier_file) == 1)



#########################
### Set/Change stats. ###
#########################
change_stats_info(hier_file,
                  tests = c('SET',
                            't-test',
                            'FWE'),
                  gamma.eigs = 0.4,  # upper bound on 2nd eig., from calc. in manuscript
                  pvalue.eigs.unc = 0.001,
                  pvalue.MCP.base = 0.001,
                  verbose = T)



###############################################################
### Subj. & Group-level back-reconstruction of hPCA sources ###
###############################################################
### Subj.-level Source Back-reconstruction ###
group_HierBackReconstruction(hier_path,
                             prefix = exp_prefix,
                             separate_spatial_maps = T,  # recommended for leaf nodes in hier. level subset, but costly for final large levels
                             separate_time_series = F,  # recommended 2/2 accuracy
                             hier_stats_filter = T,
                             level_highest_redundant = T, # b.r. on last redundant merge level (last sig. by eigenvalue sum test)
                             levels_PC1 = T,  # options for more fine-grained control of returned spatial maps
                             levels_PC2 = F,
                             levels_show_mergers = F,
                             verbose = T)

### Create Mean & Var. vols ###
group_ScaleSources(hier_path,
                   scale='z-scores',
                   verbose=T)

group_summaryStats(hier_path,
                   prefix=exp_prefix,
                   create.csv=T,
                   create.nii=T,
                   verbose=T)



### Summary of b.r. output ###
csv_file = list.files(hier_path, '*Hier_SpatialMapsMeans.csv', full.names=T)[1]
stopifnot(!all(is.na(csv_file)) && file.exists(csv_file))
print(read.csv(csv_file, stringsAsFactors = F))  #NOTE:  appearance of NA, <NA> artifact of print(), not in .csv file



