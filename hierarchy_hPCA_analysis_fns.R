# Functions to perform a hierarchical analyses of ~fMRI data,
#   using hierarchical ("treelets") PCA,
#     w/ hierarchy constructed from voxel data,
#       following temporally concatenated & PCA-compressed group data
#
# Procedure:
#   1. Compress data using PCA applied individually to each subject,
#       w/ compression applied to temporal dimension (see 'subj_PCA' from PCA_dataCompression_fns.R).
#   2. Compress data using PCA applied to group data,
#       w/ subjs. concatenated in (PCA-reduced) temporal dimension (see 'group_PCA' from PCA_dataCompression_fns.R).
#   3. Estimate hierarchy as hPCA, using 'group_HierPCA' below.
#   4. Back-reconstruct spatial maps & time series for all sources specific to each subject,
#         using 'group_HierBackReconstruction' below.
#   5. Calculate descriptive statistics,
#         using 'group_ScaleSources' & 'group_summaryStats'.
#
# Requires:       (https://github.com/koreywylie/PCA.Rachakonda2016)
#   PCA_fns.R                 : optimized PCA method for fMRI data, simple PCA algorithms & MPOWIT
#   PCA_dataCompression_fns.R : fns. to load PCA data & decompress as part of b.r.
#
# References:
#   Lee, A.B., Nadler, B., and Wasserman, L. (2008).  Treelets--An adaptive
#     multi-scale basis for sparse unordered data.  The Annals of Applied
#     Statistics 2(2), p. 435-471.
#

##################################################################
group_HierPCA <- function (data, 
                           prefix = NA, 
                           save_path = NA,
                           return.data = F, 
                           verbose = F){
  ##################################################################
  # Perform hierarchical analysis, on group-level compressed & temporally-concatenated voxels,
  #   using hierarchical PCA ("treelets PCA").
  #
  # Input:
  #   data : file w/ group PCA space, output from group_PCA()
  #   prefix        : identifying prefix to attach to saved output
  #   save_path     : path to dir. for saved files
  # Output: 
  #   Saved file w/ suffix  *_hierarchy.RData
  #   Saved table of hierarchy levels as 'hierarchy_levels_aggregations.csv', w/ cols.:
  #     'k1'          : index of new sum comp, incorporated into subsequent levels of hieararchy 
  #     'k2'          : index of new diff. comp, stored & constant therafter
  #     'similarity'  : sim. between voxels/comps at prev. level, entered into level's PCA
  #     'Lambda1'     : 1st eigenvalue of level PCA, after centering & scaling vars. to mu=0, sd=1 
  #     'Lambda2'     : 2nd eigenvalue of level PCA, after centering & scaling vars. to mu=0, sd=1
  #     'c'           : constant c in diagonal of level Jacobi rotation matrix
  #     's'           : constant s in off-diagonal of level Jacobi rotation matrix
  # Requires:
  #   get_groupPCA_dat()     : load & format group-level PCA info, from 'PCA_dataCompression_fns.R'
  #   hierarchicalPCA()      : treelet PCA levels
  #   calc_Jacobi()          : numerically stable calculation of level Jacobi rotation matrices for above
  #  
  
  if (verbose){  cat('\nEstimating hierarchical PCA levels...')}
  
  ### Inputs & defaults ###
  s_vars = c()  # running char vector of vars. to save
  prefix = as.character(prefix)
  save_path = as.character(save_path)
  
  subj_pca_suffix  = '_PCA1.RData'       # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData'    # filesuffix for group-level temporal PCA from group_PCA()
  hier_suffix = 'Hierarchy.RData'   # filesuffix for hierarchical ICA/PCA from group_HierPCA()
  
  hier_wc = hier_suffix
  group_pca_wc = group_pca_suffix
  if (!all(is.na(prefix))){
    hier_wc = paste0(prefix, '.*_', hier_wc)
    group_pca_wc = paste0(prefix, '.*_', group_pca_wc)
  }
  s_vars = c(s_vars, 'prefix')
  
  
  ### Load data ###
  data_files = NA
  if (is.character(data) && all(file.exists(data))){
    data_files = data
    
    if (is.na(save_path)){
      save_path = unique(dirname(data_files))[1]
    }
    
      
    ### Initialize hierarchy level 0 w/ group-level PCA data ###
    group_pca_file = data_files
    group_pca = get_groupPCA_dat(group_pca_file, prefix=prefix,
                                 return.PCs=T,
                                 verbose=verbose)
    X_current = group_pca$X
    L1_current = NA
    mask = group_pca$mask
    
  }else if (is.numeric(data)){
    X_current = data
    data_files = NA
    mask = NA
    
  }else{
    stopifnot(is.numeric(data) || (is.character(data)) && file.exists(data))
  }
  s_vars = c(s_vars, 'data_files', 'mask')
  
  
  ### Hierarchy params. & data dims. ###
  L.max = dim(X_current)[1] - 1    # max. possible levels, as number of voxels - 1
  V = dim(X_current)[1]            # number of voxels/ROIs
  K2 = dim(X_current)[2]           # number of group-level PCs
  s_vars = c(s_vars, 'V', 'K2')
  
  
  #---------------------------------------------------
  ### Format stats. testing info w/ options:  ###
  #---------------------------------------------------
  # SET:  Smallest Eigenvalue Test to test for levels merging redundant vars.
  # T-test:  if similarity is a simple correlation, for merging unrelated vars.
  # FDR: False Discovery Rate
  # FWE: Family-Wise Error correction (~Bonferroni)
  #       (see fn. below for parameters for above)

  stats.info = format_stats_info(c('SET',
                                   't-test',
                                   'FWE'),
                                 'gamma.eigs' = 0.4,  # analytic calculation, based on CNR ~ U(0.65,2) in simtb code
                                 'pvalue.eigs.unc' = 0.001, # uncorrected p-value for ERT/SET
                                 'pvalue.sim.unc'  = 0.001, # uncorrected p-value for testing similarities w/ t-test
                                 'pvalue.MCP.base' = 0.001, # corrected p-value for FDR/FWE
                                 'K2' = K2)                 # d.f. for t-tests
  s_vars = c(s_vars, 'stats.info')
  #---------------------------------------------------
  
  
  #####################
  ###   Main fn.   ###
  #####################
  if (verbose){
    cat('\n...calculating hierarchical ("treelets") PCA algorithm of Lee et al. (2008)...')
  }
  hier_algorithm = 'hPCA'
  s_vars = c(s_vars, 'hier_algorithm')
  
  ### Calculated similarities between initial vars. ###
  if (verbose){  cat('\n\n......calculating similarities between initial vars...')  }
  similarity_table = calc_similarities(X_current,
                                       save_path = save_path,
                                       verbose = verbose)
  if (verbose){  cat('\n')  }
  
  
  ### Calc. hierarchy levels ###
  levels_table_save = file.path(save_path, 'hierarchy_levels_aggregations.csv')
  
  basis_L = hierarchicalPCA(X_current, 
                            sim.prev   = similarity_table,
                            levels.table.save = levels_table_save,
                            verbose = verbose)
  
  
  
  ### Get info to construct additional levels ###
  X_current = basis_L$X_l                # vars. in last level of hierarchy analyzed
  similarity_table = basis_L$similarity  # formatted sim. matrix, including path to (very long) .csv file(s)
  levels_table = basis_L$levels_table    # path to levels aggregation table, saved as (very long) .csv file(s)
  s_vars = c(s_vars, 'X_current', 'similarity_table', 'levels_table')
  
  
  ### Saving ###
  s_vars = unique(s_vars)
  save.fname = hier_suffix
  if (!is.na(prefix)){
    save.fname = paste0(prefix, '_', save.fname)
  }else{
    save.fhame = sub('^_', '', save.fname)
  }
  if (!is.na(save_path)){
    save.fname = file.path(save_path, save.fname)
  }
  if (verbose){  cat(paste0('\n...saving hierarchy info as:\n      ', save.fname, '\n\n'))  }
  save(list=s_vars, 
       file=save.fname)
  stopifnot(file.exists(save.fname))
  
  
  ##################################################################
  if (return.data){
    hier = hier_from_levels_table(levels_table, verbose=verbose)
    return(hier)
  }
} ##################################################################


##################################################################
calc_similarities <- function (X, save_path, verbose = F){
  ##################################################################
  # Calculates similarity matrix for rows (~voxels) of input data X,
  #   as correlation matrix formatted as list.
  #   
  # Input:
  #   X          : data matrix      [vars. x obs.]
  #   save_path  : location to save csv tables
  # Output:
  #
  
  if (verbose){  cat('\n......calculating simple correlations between vars...')  }
  
  p = dim(X)[1]
  n = dim(X)[2]
  R = cor(t(X))
  
  ix = c(1:p)
  i1 = i2 = c()
  for (i in 2:length(ix)){
    i1 = c(i1, ix[seq(1,(i-1))])
    i2 = c(i2, ix[rep(i,i-1)])
  }
  m = p * (p - 1) / 2
  sim_table1 = matrix(0, m, 3)
  sim_table1 = cbind(R[upper.tri(R)], i1, i2, deparse.level=0)
  sim_table1 = as.data.frame(sim_table1)
  names(sim_table1) = NULL
  
  if (nrow(sim_table1) > 0){
    save_table_path = file.path(save_path, 'similarity_table1.csv')  
    if (verbose){  cat(paste0('\n.........saving as:  ', save_table_path))  }
    write.csv(sim_table1, file=save_table_path, row.names=F)
    
    if (verbose && file.exists(save_table_path)){  cat('\n............done')  }
  }
  
  ##################################################################
  return(list(similarity_matrix = NA,   # sim. matrix formatted as [vars. by c(r, k1, k2)] table
              files          = save_table_path,   # filepaths to saved sim. matrix pieces
              files_vars     = 'sim_table1',      # var. names assigned to above sim. matrix pieces
              files_sims_max = max(sim_table1[,1])))  # max. similarities in each file
} ##################################################################



##################################################################
get_hier_dat <- function(hier_file, prefix = NA,
                         hier_levels = NA,
                         hier_active_all = F,
                         return.X_input = F,
                         return.X_current = F,
                         return.similarities = F,
                         return.level.vars = F,
                         return.init.redundant.vars = F,
                         return.premerge.vars = F,
                         return.level.rotations = F,
                         return.level.basis.fns = F,
                         apply.stats = F,
                         resave.hier_file = F,
                         verbose = F){
  ##################################################################
  # Loads & formats group-level hierarchical PCA
  #
  # Input:
  #   hier_file   : path to file w/ output from group_HierPCA()
  #   prefix      : prefix, used to verify check experiment info
  #   hier_levels                : integer vector, used to limit analysis to select levels
  #   hier_active_all            : Reconstruct all active vars. at requested hier_levels
  #   return.X_input             : data input to hierarchy
  #   return.X_current           : data from last level of hierarchy constructed
  #   return.similarities        : current similarity matrix, path to saved file(s)
  #   return.level.vars          : analyze sum & diff. vars. for levels
  #   return.init.redundant.vars : analyze initial levels which merge redundant vars.
  #                                 (~all active vars. in lowest non-redundant level pre-merge)
  #   return.premerge.vars       : analyze pre-merge state of vars. merged at select levels,
  #                                 used to show mergers at each level
  #   return.level.rotations     : analyze Jacobian rotation matrices at each level 
  #   return.level.basis.fns     : analyze scaling & detail fns. for returned vars.
  #   apply.stats : limit analysis to selected levels by apply statistical testing,
  #                   (e.g.: Smallest Eigenvalue Test to identify levels merging redundant vars.)
  #   resave.hier_file : save all hierarchy info in "hier_file" above
  #
  # Output:  
  #   hier : formatted list of level info, w/ elements depending on above options
  #
  # Requires:
  #   get_groupPCA_dat()  : from PCA_dataCompression_fns.R
  #   hier_from_levels_table()
  #   reconstruct_hier_levels()
  #   find_LevelEndSeq_byTrait()
  #   hier_apply_stats()
  #
  
  if (verbose){cat(paste0('\n...loading hierarchical analysis data from:\n      ', hier_file))}
  
  ### Prelims ###
  stopifnot(file.exists(hier_file))
  l.all_active_vars = c()
  l.all_inactive_vars = c()
  filter.by.stats = as.logical(apply.stats)  # returns all levels that either pass stats. criteria, or weren't tested
  
  Space.hier = new.env()
  load(hier_file, envir=Space.hier)
  stopifnot(all(c('levels_table', 'data_files', 'prefix') %in% names(Space.hier)))
  if (all(is.na(prefix))){
    prefix = get('prefix', envir=Space.hier)
  }else{
    stopifnot(prefix == get('prefix', envir=Space.hier))
  }
  
  levels_table = get('levels_table', envir=Space.hier)
  if (is.character(levels_table) && !is.na(levels_table)){
    stopifnot(file.exists(levels_table))
  }
  
  
  ### Format hierarchy info ###
  l.max = NA
  if (!all(is.na(hier_levels)) && is.numeric(hier_levels)){
    l.max = max(hier_levels, na.rm=T)
  }
  
  X.orig = get_groupPCA_dat(get('data_files', envir=Space.hier), prefix=prefix,
                            return.PCs=T, verbose=verbose)$X
  
  hier = hier_from_levels_table(levels_table, 
                                verbose = verbose)
  hier[[1]]$data.dims = dim(X.orig)  # original dimensions, needed for stats. testing
  
  
  if (filter.by.stats || return.init.redundant.vars){
    if (!any('stats.outcome' == names(hier[[1]]))){
      if ('stats.info' %in% names(Space.hier)){
        stats.info = get('stats.info', envir=Space.hier)
        stopifnot(any(c('p_value', 'p_value_cutoff_unc', 'p_value_cutoff_corrected') %in% names(stats.info[[1]])))
        
        hier = hier_apply_stats(hier, stats.info, verbose=verbose)
        
      }else if (verbose){
        cat(paste0('\n\nWARNING:  could not find formatted statistical testing info in:  ', hier_file))
      }
    }
  }
  
  ### Get info for saved levels of hierarchy ###
  h.levels.saved = unlist(lapply(hier, function(yy) return(yy$h.level)))
  
  ### Select hierarchy levels ###
  h.levels.requested = NA
  hier_levels.filter = NA
  if (!all(is.na(hier_levels))){
    if (is.character(hier_levels) && ((hier_levels == 'final') || (hier_levels == 'last'))){
      
      ### Return all sources at highest/final level of hierarchy ###
      h.levels.requested = hier[[length(hier)]]$h.level
      hier_levels.filter = length(hier)
      
      l.all_active_vars = c(l.all_active_vars,
                            h.levels.requested)
      
    }else if (is.character(hier_levels) && (hier_levels == 'all')){
      ### Return all sources at each/every level of hierarchy ###
      h.levels.requested = h.levels.saved
      hier_levels.filter = 1:length(hier)
      
    }else if (is.numeric(hier_levels) && hier_active_all){
      ### Return all sources at requested level(s) of hierarchy ###
      h.levels.requested = h.levels.saved[h.levels.saved %in% hier_levels]
      hier_levels.filter = which(h.levels.saved %in% hier_levels)
      
      l.all_active_vars = c(l.all_active_vars,
                            h.levels.requested)
    
    }else{
      ### Return sum & diff., pre-merge sources at select levels ###
      hier_levels = as.numeric(hier_levels)
      
      if (verbose && !all(hier_levels %in% h.levels.saved)){
        cat(paste0('\nWARNING: could not find requested hierarchy levels in saved data: ',
                   paste(hier_levels[!(hier_levels %in% h.levels.saved)], collapse=', ')))
      }
      
      h.levels.requested = h.levels.saved[h.levels.saved %in% hier_levels]
      hier_levels.filter = which(h.levels.saved %in% hier_levels)
    }
  }
  
  h.levels.stats = NA
  hier_stats.filter = NA
  if (filter.by.stats){
    if (verbose && !any('stats.outcomes' == names(hier[[1]]))){
      cat(paste0('\nWARNING:  could not filter hierarchy levels by stats, statistical testing not yet applied!'))
    }
    stopifnot(all(c('stats.tests', 'stats.outcomes') %in% names(hier[[1]])))
    
    h.levels.stats.pass = rep(T, length(h.levels.saved))
    
    if ('redundant.vars' %in% names(hier[[1]])){
      ### Skip returning levels merging redundant vars. ###
      h.levels.stats.redun = unlist(lapply(hier, function(yy) return(yy$redundant.vars)))
      h.levels.stats.pass[h.levels.stats.redun] = F
      if (verbose && any(h.levels.stats.redun)){
        cat(paste0('\n......redundant var. mergers found by Smallest Eigenvalue Test:   ',
                   sum(h.levels.stats.redun),'/',length(h.levels.stats.redun)))
        cat(paste0('\n.........excluding above levels from b.r.'))
      }else{
        cat(paste0('\n......_NO_ redundant var. mergers found by Smallest Eigenvalue Test:   ',
                   sum(h.levels.stats.redun),'/',length(h.levels.stats.redun)))
      }
    }
    
    if ('nonsimilar.vars' %in% names(hier[[1]])){
      ### Skip returning levels merging dissimilar, non-hierarchical vars. ###
      h.levels.stats.dissim = unlist(lapply(hier, function(yy) return(yy$nonsimilar.vars)))
      h.levels.stats.pass[h.levels.stats.dissim] = F
      if (verbose && any(h.levels.stats.dissim)){
        cat(paste0('\n......dissimilar, non-hierarchical var. mergers found:   ',
                   sum(h.levels.stats.dissim),'/',length(h.levels.stats.dissim)))
        cat(paste0('\n.........excluding above levels from b.r.'))
      }else{
        cat(paste0('\n......_NO_ dissimilar/non-hierarchical var. mergers found:   ',
                   sum(h.levels.stats.dissim),'/',length(h.levels.stats.dissim)))
      }
    }
    
    if (verbose && any(!h.levels.stats.pass)){
      cat(paste0('\n......total redundant/non-sig. levels:  ',sum(!h.levels.stats.pass),'/',length(h.levels.stats.pass)))
      cat(paste0('\n.........excluding above levels from further analyses'))
    }
    h.levels.stats = h.levels.saved[h.levels.stats.pass]
    hier_stats.filter = which(h.levels.stats.pass)
  }
  
  if (!all(is.na(h.levels.stats)) && !all(is.na(h.levels.requested))){
    h.levels.stats = h.levels.stats[h.levels.stats %in% h.levels.requested]
    h.levels.requested = h.levels.requested[h.levels.requested %in% h.levels.stats]
    hier_levels.filter = hier_levels.filter[hier_levels.filter %in% h.levels.requested]
    hier_stats.filter = hier_stats.filter[hier_stats.filter %in% h.levels.stats]
  }
  
  ### Return similarity matrix ###
  if (return.similarities){
    hier[[1]]$similarity_table = get('similarity_table', envir=Space.hier)
  }
  
  ### Find vars. at end of initial sequence of redundant levels ###
  hier_end_redun_seq.filter = NA
  if (return.init.redundant.vars){
    if (verbose && !any('stats.outcome' == names(hier[[1]])) && !('stats.info' %in% names(Space.hier))){
      cat(paste0('\nWARNING:  could not find sequence of initial redundant vars., missing info to test hier. level mergers for redundancy!'))
    }
    hier_end_redun_seq.filter = find_LevelEndSeq_byTrait(hier, 'redundant.vars', first=F)
    h.levels.end_redun_seq = hier[[hier_end_redun_seq.filter]]$h.level
    l.all_active_vars = c(l.all_active_vars, h.levels.end_redun_seq)
  }
  
  ### Check for missed/discarded levels ###
  l.missing = c()
  if (!all(is.na(h.levels.requested))){
    for (l in h.levels.requested){
      if (!any(l == h.levels.saved)){
        l.missing = c(l.missing, l)
      }else if (return.level.vars){
        l.ind = which(l == h.levels.saved)
        if (!('new.vars' %in% names(hier[[l.ind]])) ||
            !is.numeric(hier[[l.ind]]$new.vars)){
          l.missing = c(l.missing, l)
        }
      }
    }
  }
  if (!all(is.na(h.levels.stats))){
    for (l in h.levels.stats){
      if (!any(l == h.levels.saved)){
        l.missing = c(l.missing, l)
      }else if (return.level.vars){
        l.ind = which(l == h.levels.saved)
        if (!('new.vars' %in% names(hier[[l.ind]]))){
          l.missing = c(l.missing, l)
        }
      }
    }
  }
  l.missing = sort(unique(l.missing))
  if (length(l.missing) == 0){  l.missing = NA  }
  
  ### Select levels for pre-merge vars. ###
  l.premerge = NA
  if (return.premerge.vars){
    l.premerge = h.levels.saved
    
    l.premerge.filter = c()  # filter pre-merge levels by requested levels or by sig. stats, if indicated
    if (!all(is.na(hier_levels.filter))){  l.premerge.filter = c(l.premerge.filter, hier_levels.filter)}
    if (!all(is.na(hier_stats.filter))){  l.premerge.filter = c(l.premerge.filter, hier_stats.filter)}
    if (length(l.premerge.filter) > 0){
      l.premerge.filter = unique(l.premerge.filter)
      l.premerge = l.premerge[l.premerge.filter]
    }
  }
  
  ### Estimate missed/discarded levels ###
  if (length(l.all_active_vars) == 0){  l.all_active_vars = NA  }
  if (length(l.all_inactive_vars) == 0){  l.all_inactive_vars = NA  }
  if (return.level.basis.fns){
    l.basis = sort(unique(c(l.all_active_vars, l.missing, l.premerge)))
    l.basis = l.basis[!is.na(l.basis)]
  }else{  l.basis = NA  }
  if (!all(is.na(c(l.all_active_vars, l.missing, l.premerge)))){
    hier[[1]]$X = X.orig
    
    hier = reconstruct_hier_levels(hier,
                                   L.active = l.all_active_vars,
                                   L.inactive = l.all_inactive_vars,
                                   L.new_vars = l.missing,
                                   L.premerge = l.premerge,
                                   L.basis = l.basis,
                                   verbose = verbose)
  }
  
  ### Discard unused hierarchy levels ###
  hier.filter = c(hier_levels.filter,        # include hierarchy levels specified by input
                  hier_end_redun_seq.filter, # include last level in initial seq. of redundant var. mergers
                  hier_stats.filter)         # include stat.-tested levels/non-tested levels
  
  
  if (!all(is.na(hier.filter))){
    hier.filter = sort(unique(hier.filter[!is.na(hier.filter)]))
    
    if (!any(1 == hier.filter)){  # Preserve meta-info from 1st level of hierarchy
      meta_names = c('levels_table_file', 'data.dims',
                     'stats.tests', 'stats.MCP.method', 'stats.MCP.p.value.base', 
                     'X')
      meta_names = meta_names[meta_names %in% names(hier)]
      hier_meta = hier[[1]][meta_names]
      hier = hier[hier.filter]  # R vector indexing of list, keeps named vars. w/n list ind.
      hier[[1]] = c(hier[[1]], hier_meta)
    }else{
      hier = hier[hier.filter]  # R vector indexing of list, keeps named vars. w/n list ind.
    }
    
    if (length(hier) == 0){
      if (verbose){'\n......returning empty hierarchy after filtering out all non-requested levels!'}
      return(hier)
    }
  }
  
  ### Discard mixing matrices (~Jacobi rotation matrices), if indicated ###
  if (!return.level.rotations){
    for (h in 1:length(hier)){
      if ('J_l' %in% names(hier[[h]])){  hier[[h]]$J_l = NULL  }
    }
  }
  
  ### Add meta- info to 1st level ###
  hier[[1]]$prefix            = prefix
  hier[[1]]$group_pca_file    = get('data_files',        envir=Space.hier)
  hier[[1]]$mask              = get('mask',              envir=Space.hier)
  hier[[1]]$levels_table      = get('levels_table',      envir=Space.hier)
  if ('hier_algorithm' %in% names(Space.hier)){  
    hier[[1]]$hier_algorithm = get('hier_algorithm', envir=Space.hier)
  }
  
  ### Return data input to hierarchy ###
  if (return.X_input){
    if (!('X' %in% names(hier[[1]]))){
      if ('data' %in% names(Space.hier)){
        hier[[1]]$X = get('data', envir=Space.hier)  # get input to fn.
      }else{
        hier[[1]]$X = get_groupPCA_dat(hier[[1]]$group_pca_file, prefix = prefix,
                                       return.PCs = T,
                                       verbose = verbose)$X
      }
    }
  }else{  hier[[1]]$X = NULL  }
  
  
  ### Return current vars. in hierarchy ###
  if (return.X_current && ('X_current' %in% names(Space.hier))){
    hier[[1]]$X_current = get('X_current', envir=Space.hier)
  }else if (return.X_current){
    hier[[1]]$X_current = get('data', envir=Space.hier)
  }
  
  ### Save compiled & formatted hier info ###
  if (resave.hier_file){
    if (verbose){  cat('\n.........re-saving hierarchy info')  }
    assign('hier', hier, envir=Space.hier)
    save(list=ls(Space.hier), file=hier_file, envir=Space.hier)
  }
  
  ##################################################################
  return(hier)      # list w/ hierarchical PCA basis'
} ##################################################################


##################################################################
hier_from_levels_table <- function(levels_table, 
                                   verbose = F){
  ##################################################################
  # Creates hierarchy levels info from saved 'levels_table',
  #   formatting as if output from main hierarchy fn.
  #
  
  if (verbose){  cat('\n...constructing hierarchy from saved levels aggregation table...')}
  
  levels_table_file = NA
  if (is.character(levels_table) && file.exists(levels_table)){
    levels_table_file = levels_table
    levels_table = read.csv(levels_table, stringsAsFactors=F, header=T)
  }
  stopifnot(is.data.frame(levels_table) || is.numeric(levels_table))
  
  hier = list()
  for (l in 1:(nrow(levels_table))){
    hl.level = levels_table[l,]
    
    hier[[l]] = list('h.level'    = as.integer(l),       # saved hierarchy level
                     'k1'         = as.integer(hl.level[1]), # index of new comp, incorporated into subsequent levels of hieararchy 
                     'k2'         = as.integer(hl.level[2]), # index of new 2nd comp, stored & constant therafter
                     'similarity' = as.numeric(hl.level[3]), # sim. between vars. at prev. level, entered into level's PCA
                     'eigenvalues' = as.numeric(hl.level[4:5]))  # eigenvalues of PCA for level
    c = as.numeric(hl.level[6])
    s = as.numeric(hl.level[7])
    hier[[l]]$J_l = matrix(c(c,-s,
                             s, c), 2, 2, byrow=T)  # Jacobian rotation matrix from PCA
  }
  
  ### Add meta-info ###
  hier[[1]]$levels_table_file = levels_table_file
  
  ##################################################################
  return(hier)
} ##################################################################


##################################################################
format_stats_info <- function(tests, ..., verbose = F){
  ##################################################################
  # (Re-)Formats & returns statistical testing params. for hPCA,
  #   identifying levels merging redundant vars. w/ Smallest Eigenvalue Test,
  #   levels merging uncorrelated/unrelated vars. w/ t-test,
  #   & controlling for multiple comparisons using FWE or FDR
  #
  # Input:
  #   tests  : vector statistical tests to apply, accepted options:
  #     SET    : Smallest Eigenvalue Test, to identify redundant mergers
  #     t-test : t-test for correlations, to identify dissimilar mergers
  #     FDR    : False Discovery Rate
  #     FWE    : Family-Wise Error (~Bonferroni)
  #   ... : required parameters for above tests, w/ accepted names as input fn. arguments:
  #     stats.info      : output of prev. call to format_stats_info(), if changing stats. tests or params.
  #     gamma.eigs      : expected max. 2nd eigenvalue
  #     pvalue.eigs.unc : p-value for Smallest Eigenvalue Test or Eigenvalue Ratio Test, uncorrected
  #     pvalue.sim.unc  : p-value for merging levels based on sim. measure, uncorrected
  #     pvalue.MCP.base : corrected p-value for Multiple Comparison Procedure (ex: q-value for FDR)
  #     K2              : length of PCA-compressed, temp. concatenated time series
  #     h.levels.max    : total number of specified/possible levels in hierarchy (for FWE correction)
  # Output:  formatted list by test(s)
  #
  
  if (verbose){  cat('\n......formatting statistical inference w/ parameters:')  }
  
  stopifnot(is.character(tests))
  known.tests = c('SET', 'Smallest Eigenvalue Test',
                  'ttest', 't-test', 
                  'FWE', 'FDR')
  if (!all(tests %in% known.tests)){
    unknowns = tests[!(tests %in% known.tests)]
    if (verbose){
      cat(paste0('\n\nWARNING:  requested stats. tests not implemented:  ', paste(unknowns, collapse=', ')))
    }
    stopifnot(all(tests %in% known.tests))
  }
  
  ### Default params. ###
  gamma.eigs = 0.4         # max. 2nd eigenvalue
  pvalue.eigs.unc = 0.05   # p-value for Smallest Eigenvalue Test or Eigenvalue Ratio Test, uncorrected
  pvalue.sim.unc  = 0.05   # p-value for merging levels based on sim. measure, uncorrected
  pvalue.MCP.base = 0.05   # corrected p-value for Multiple Comparison Procedure (ex: q-value for FDR)
  K2 = NA                  # length of PCA-compressed, temp. concatenated time series
  h.levels.max = NA        # total number of specified/possible levels in hierarchy
  
  ### Input params. ###
  params = list(...)
  known.params = c('gamma.eigs', 'pvalue.eigs.unc', 'pvalue.sim.unc', 
                   'pvalue.MCP.base', 'K2', 'h.levels.max')
  param.names = names(params)
  
  if ('stats.info' %in% param.names){
    ### Preserve prev. input params, if unchanged ###
    p0 = which(param.names == 'stats.info')
    p1 = which(names(params[[p0]]) == 'input_params')
    prev.params = params[[p0]][[p1]]
    prev.params = prev.params[-which(names(prev.params) %in% param.names)]
    params = params[-which(names(params) == 'stats.info')]
    params = c(prev.params, params)
    param.names = names(params)
    
  }
  if (length(params) > 0){
    p.rm = c()
    for (p in 1:length(params)){
      ### Add new params. ###
      param.name = param.names[p]
      if (param.name %in% known.params){
        if (verbose){  cat(paste0('\n......   ',param.name,' = ', params[[param.name]]))  }
        assign(param.name, params[[param.name]])
      }else if (verbose  && !is.na(param.name) && (param.name != 'stats.info')){
        cat(paste0('\n\nWARNING:  could not understand input parameter:  ',param.name))
        p.rm = c(p.rm, p)
      }
    }
    if (length(p.rm) > 0){  # remove uninterpretable params from list
      params[p.rm] = NULL
      param.names = param.names[-p.rm]
    }
  }
  if (any(c('t-test', 'ttest') %in% tests)){  stopifnot('K2' %in% param.names)  }
  
  #-----------------------------------------
  ###     Stats. tests implemented     ###
  #-----------------------------------------
  SET = list('name' = 'Smallest Eigenvalue Test',
             'abbr' = paste0('SET(',gamma.eigs,')'),
             'param_name' = 'gamma.eigs',
             'parameter' = gamma.eigs,  # max. 2nd eigenvalue
             'p_value_cutoff_unc' = pvalue.eigs.unc)
  ttest = list('name' = 't-test',
               'abbr' = paste0('t(', K2 - 2,')'),
               'param_name' = 'df',
               'parameter' = K2 - 2,
               'p_value_cutoff_unc' = pvalue.sim.unc)
  MCP.FWE = list('name' = 'FWE',
                 'abbr' = 'FWE',
                 'param_name' = 'number of tests',
                 'parameter' = h.levels.max,  # if NA, defaults to number of p-values
                 'p_value_MCP_base' = pvalue.MCP.base)
  MCP.FDR = list('name' = 'FDR',
                 'abbr' = 'FDR',
                 'param_name' = NA,
                 'parameter' = NA,
                 'p_value_MCP_base' = pvalue.MCP.base)
  
  
  ### Format stats info ###
  stats.info = list()
  if (any(c('SET', 'Smallest Eigenvalue Test') %in% tests)){
    stats.info[[length(stats.info)+1]] = SET
  }
  if (any(c('t-test', 'ttest') %in% tests)){
    stats.info[[length(stats.info)+1]] = ttest
  }
  if ('FDR' %in% tests){
    stats.info[[length(stats.info)+1]] = MCP.FDR
  }else if ('FWE' %in% tests){
    stats.info[[length(stats.info)+1]] = MCP.FWE
  }
  ### Append input params ###
  stats.info$input_params = params
  
  ##################################################################
  return(stats.info)
} ##################################################################


##################################################################
change_stats_info <- function(hier_file, tests = NA, ..., verbose = F){
  ##################################################################
  # Change stats. info in 'hier_file' & re-save.
  #   See 'format_stats_info()' for params. & details.
  #
  
  if (verbose){  cat('\n......changing & re-formatting statistical testing info...')  }
  stopifnot(is.character(hier_file) && file.exists(hier_file))
  if (verbose){  cat(paste0('\n.........in file:  ', hier_file))  }
  
  params = list(...)
  
  Space.hier = new.env()
  load(hier_file, envir=Space.hier)
  if (!('stats.info' %in% names(Space.hier))){
    stats.info = format_stats_info(tests, ..., verbose = verbose)
    
  }else{
    prev.stats.info = get('stats.info', envir = Space.hier)
    prev.tests = unlist(lapply(prev.stats.info, function(yy) yy$name))
    tests = unique(c(prev.tests, tests))
    
    stats.info = format_stats_info(tests, stats.info = prev.stats.info, ..., verbose = verbose)
  }
  
  if (verbose){  cat('\n.........re-saving updated stats. info...')  }
  assign('stats.info', stats.info, envir=Space.hier)
  save(file = hier_file, list = ls(Space.hier), envir = Space.hier)
  if (verbose){  cat('done')  }
  
} ##################################################################


##################################################################
hier_apply_stats <- function(hier, stats.info, 
                             stats.table.save_path = NA,
                             verbose = F){
  ##################################################################
  # Applies statistical testing to hierarchy levels,
  #   using tests & parameters formatted w/ format_stats_info()
  #     & returns (thresholded/stats. tested) hierarchy
  #
  # Input:
  #   hier  : output from get_hier_dat() or hier_from_levels_table()
  #   stats.info : list of stats tests to apply, each formatted as
  #     name   : name of statistical test
  #     abbr   : abbreviation of test name
  #     param_name : stats. test param. name
  #     parameter  : stats. test param. value
  #     p_value_cutoff_unc : uncorrected p-value cutoff
  #     p_value_MCP_base   : MCP base (~q-value, optional)
  #     p_value_cutoff_corrected   : corrected, post-MCP p-value cutoff (optional)
  #   stats.table.save_path : if not NA, path to save location for .csv table showing stats for all levels
  #
  # Output: 
  #   hier : list as above, w/ stats. descriptive info in 1st level:
  #     stats.tests         : name of statistical test(s) applied
  #     stats.MCP.method    : name of MCP applied, if any
  #     stats.MCP.p.value.base : base corrected p-value/q-value for MCP test 
  #         ... & added fields for each level:
  #     stats.p.values      : p-value(s) for test(s), uncorrected
  #     stats.outcomes      : T/F if H0 rejected, pre- or post-MCP
  #     stats.MCP.corrected : T/F indicates if Multiple Comparison Procedure applied
  #     redundant.vars      : indicates if level mergers redundant vars.
  #     nonsimilar.vars     : indicates if level mergers unrelated (i.e., uncorrelated) vars.
  #
  
  if (verbose){  cat('\n...statistically testing hierarchy levels...')  }
  
  if (!is.list(hier) || (length(hier) == 0)){
    if (verbose){  cat('\n\nWARNING:  empty hierarchy input, cannot apply statistical testing')  }
    return(hier)
  }
  if (all(c('name', 'parameter', 'p_value_cutoff_unc') %in% names(stats.info))){
    stats.info = list(stats.info)
  }else if ((length(stats.info) < 1) ||
            !all(c('name', 'parameter', 'p_value_cutoff_unc') %in% names(stats.info[[1]]))){
    if (verbose){  cat(paste0('\n\nWARNING:  could not understand input stats test info:  ', stats.info))  }
    return(hier)
  }
  if (all(is.na(stats.table.save_path)) || !is.character(stats.table.save_path) || !file.exists(stats.table.save_path)){
    if ('levels_table_file' %in% names(hier[[1]])){
      stats.table.save_path = dirname(hier[[1]]$levels_table_file)
    }
  }
  stopifnot(is.character(stats.table.save_path))
  
  
  ### Main Statistical Testing ###
  for (t in 1:length(stats.info)){
    if (names(stats.info[t]) == 'input_params'){  next  }
    
    ### Overwrite prev. applied test with new parameters ###
    if (('stats.tests' %in% names(hier[[1]]))){
      prev.stats.tests = gsub('[0-9]', '', hier[[1]]$stats.tests)
      curr.stat.test = gsub('[0-9]', '', stats.info[[t]]$abbr)
      if (curr.stat.test %in% prev.stats.tests){
        if (verbose){
          i = which(curr.stat.test %in% prev.stats.tests)
          cat(paste0('\n......changing prev. applied test:  ',hier[[1]]$stats.tests[i],'  to new test:  ',stats.info[[t]]$abbr))
        }
        
        if (any(hier[[1]]$stats.tests %in% c('FDR', 'FWE'))){
          if (verbose){  cat(paste0('\n......resetting MCP control method...'))  }
          i = which(hier[[1]]$stats.tests %in% c('FDR', 'FWE'))
          hier[[1]]$stats.tests = hier[[1]]$stats.tests[-i]  # MCP method removed 1st to allow proper indexing of p-values & stats outcomes below
          hier[[1]]$stats.MCP.method = NULL
        }
        
        i = which(curr.stat.test == prev.stats.tests)
        for (h in 1:length(hier)){
          hier[[h]]$stats.outcomes = hier[[h]]$stats.outcomes[-i]
          hier[[h]]$stats.p.values = hier[[h]]$stats.p.values[-i]
          hier[[h]]$stats.MCP.corrected = hier[[h]]$stats.MCP.corrected[-i]
        }
        hier[[1]]$stats.tests = hier[[1]]$stats.tests[-i]
      }
    }
    hier[[1]]$stats.tests = c(hier[[1]]$stats.tests, stats.info[[t]]$abbr)
    
    p.cutoff = 0.05
    MCP.corrected = F
    if ('p_value_cutoff_corrected' %in% names(stats.info[[t]])){
      p.cutoff = stats.info[[t]]$p_value_cutoff_corrected
      MCP.corrected = T
    }else if ('p_value_cutoff_unc' %in% names(stats.info[[t]])){
      p.cutoff = stats.info[[t]]$p_value_cutoff_unc
    }
    
    if (stats.info[[t]]$name == 'Smallest Eigenvalue Test'){
      ### Apply Smallest Eigenvalue Test to test for redundancy in levels ###
      if (verbose){cat('\n......applying Smallest Eigenvalue Test to identify levels merging redundant vars....')}
      if (verbose){cat(paste0('\n.........   gamma = ', stats.info[[t]]$parameter,',  &  p < ', p.cutoff))}
      
      if (all(c('data.dims', 'eigenvalues') %in% names(hier[[1]])) && 
          all(is.numeric(hier[[1]]$eigenvalues))){
        vrbs1 = verbose
        for (h in 1:length(hier)){
          if (vrbs1){  cat(paste0('\n.........    example test for level = ',hier[[h]]$h.level,':'))  }
          
          test = smallestEigenval_test(eigs = hier[[h]]$eigenvalues,
                                       n = hier[[1]]$data.dims[2], # num. obs. as post-groupPCA dim. K2
                                       g = stats.info[[t]]$parameter,   
                                       pval = p.cutoff,
                                       return.pca = T, # return exact p-values & test outcome
                                       verbose = vrbs1)
          hier[[h]]$redundant.vars = test$lambda_sum.test  # H0: non-redundant vars, rejection indicates redundant vars.
          hier[[h]]$stats.outcomes  = c(hier[[h]]$stats.outcomes, test$lambda_sum.test)
          hier[[h]]$stats.p.values  = c(hier[[h]]$stats.p.values, test$lambda_sum.pvalue.exact)
          hier[[h]]$stats.MCP.corrected = c(hier[[h]]$stats.MCP.corrected, MCP.corrected)
          
          if (vrbs1){  cat(paste0('\n.........   & similarly for all other levels l = 2, 3, ..., ',length(hier)))  }
          vrbs1 = F
        }
      }else{
        if (verbose){cat(paste0('\n\nWARNING: could not find data dims., eigenvalues for Smallest Eigenvalue Test!'))}
        stopifnot(all(c('data.dims', 'eigenvalues') %in% names(hier[[1]])))
        stopifnot(all(is.numeric(hier[[1]]$eigenvalues)))
      }
      
    }else if (stats.info[[t]]$name == 't-test'){
      ### Apply t-test to correlations to test for levels merging associated vars. ###
      if (verbose){cat('\n......applying t-test to correlations to identify levels merging dissimilar vars....')  }
      
      for (h in 1:length(hier)){
        K2 = hier[[1]]$data.dims[2] # num. obs. as post-groupPCA dim. K2
        r = hier[[h]]$similarity
        tstat = r * sqrt(K2 - 2) / sqrt(1 - r^2)  # Anderson. (2003). Wiley, New Jersey. p. 121.
        test = list('statistic' = tstat,
                    'parameter' = K2 - 2,
                    'p.value'   = 1 - pt(tstat, K2 - 2), # one-sided t-test
                    'estimate'  = r,
                    'similarity'= r)
        hier[[h]]$nonsimilar.vars = test$p.value > p.cutoff  # H0: r = 0, failure to reject indicates invalid merger
        hier[[h]]$stats.outcomes   = c(hier[[h]]$stats.outcomes, test$p.value <= p.cutoff)
        hier[[h]]$stats.p.values   = c(hier[[h]]$stats.p.values, test$p.value)
        hier[[h]]$stats.MCP.corrected = c(hier[[h]]$stats.MCP.corrected, MCP.corrected)
      }
      
    }else if (!(stats.info[[t]]$name %in% c('Smallest Eigenvalue Test', 't-test', 'FWE', 'FDR'))){
      if (verbose){cat(paste0('\n\nWARNING: unknown/unimplemented stats. test  ', stats.info[[t]]$name))}
      stopifnot(stats.info[[t]]$name %in% c('Smallest Eigenvalue Test', 't-test', 'FWE', 'FDR'))
    }
  }
  
  
  ### Multiple Comparisons Procedure ###
  stopifnot('stats.p.values' %in% names(hier[[1]]))
  for (t in 1:length(stats.info)){
    if (names(stats.info[t]) == 'input_params'){  next  }
    pvals = unlist(lapply(hier, function(yy) return(yy$stats.p.values)))
    
    if (stats.info[[t]]$name == 'FWE'){
      if (!('parameter' %in% names(stats.info[[t]])) || !is.numeric(stats.info[[t]]$parameter)){
        stats.info[[t]]$parameter = length(pvals)
      }
      p.cutoff = stats.info[[t]]$p_value_MCP_base / stats.info[[t]]$parameter
      if (verbose){  cat(paste0('\nApplying FWE to ',stats.info[[t]]$parameter,' tests, base p = ',
                                stats.info[[t]]$p_value_MCP_base, ' (corrected)'))  }
      
    }else if (stats.info[[t]]$name == 'FDR'){
      p.cutoff = FDR(pvals, stats.info[[t]]$p_value_MCP_base, verbose = verbose)
      
    }else{  next  }
    
    if (verbose){  cat(paste0('\n...MCP-corrected p-value cutoff:  p < ', p.cutoff))  }
    
    for (h in 1:length(hier)){
      hier[[h]]$stats.outcomes = hier[[h]]$stats.p.values <= p.cutoff
      hier[[h]]$stats.MCP.corrected = rep(T, length(hier[[h]]$stats.p.values))
      
      if (any(grepl('^SET\\(.*\\)$', hier[[1]]$stats.tests))){  # correct SET check for merging redundant vars. for MCP
        i = which(grepl('^SET\\(.*\\)$', hier[[1]]$stats.tests))
        hier[[h]]$redundant.vars = hier[[h]]$stats.outcomes[i]  # H0 rejected, merged vars. likely redundant
      }else if (any(grepl('^t\\([0-9]*\\)$', hier[[1]]$stats.tests))){  # correct non-sim. t-test outcome for MCP
        i = which(grepl('^t\\([0-9]*\\)$', hier[[1]]$stats.tests))
        hier[[h]]$nonsimilar.vars = !hier[[h]]$stats.outcomes[i]  # H0 not rejected, merged vars. likely unrelated
      }
    }
    hier[[1]]$stats.MCP.method = stats.info[[t]]$name 
    hier[[1]]$stats.MCP.p.value.base = stats.info[[t]]$p_value_MCP_base
  }
  
  
  ### Save output as .csv table ###
  if (!all(is.na(stats.table.save_path))){
    create_hier_StatsTable(hier, stats.table.save_path, verbose)
  }
  
  ##################################################################
  return(hier)
} ##################################################################


##################################################################
find_LevelEndSeq_byTrait <- function(hier, 
                                     trait.name = 'redundant.vars',
                                     trait.present = T,
                                     sequential = T, first = T){
  ##################################################################
  # Finds first/last level in sequence of hierarchy levels output from get_hier_dat(),
  #   as determined by presence/absence of named trait.
  # Returns +/-Inf if no level of hierarchy has named trait.
  #
  
  stopifnot(length(hier) > 0)
  stopifnot(is.character(trait.name))
  stopifnot(trait.name %in% names(hier[[1]]))
  stopifnot(is.logical(get(trait.name, hier[[1]])))
  
  h.levels.trait = unname(unlist(lapply(hier, function(yy) yy[trait.name])))
  if (!trait.present){  h.levels.trait = !h.levels.trait  }
  if (sequential && (length(h.levels.trait) > 1)){
    # find largest run of True values in input, discounting isolated smaller runs
    h.levels.trait = with(rle(h.levels.trait),
                          rep(seq_along(values) == which.max(lengths * values), lengths))
  }
  
  h.levels.saved = unlist(lapply(hier, function(yy) return(yy$h.level)))
  h.levels.saved = h.levels.saved[h.levels.trait]
  
  if (first){  # find lowest hierarchy level with trait
    return(  min(h.levels.saved)  ) 
  }else{       # find highest hierarchy level with trait
    return(  max(h.levels.saved)  )
  }
  
} ##################################################################


##################################################################
reconstruct_hier_levels <- function (hier, 
                                     L.active = NA,
                                     L.inactive = NA,
                                     L.new_vars = NA,
                                     L.premerge = NA,
                                     L.basis = NA,
                                     verbose = F){
  ##################################################################
  # Reconstructs _all_ active vars. at level(s) 'L.active' of hierarchy,
  #   all inactive (i.e. diff) vars. at level(s) 'L.inactive' of hierarchy, 
  #   new "sum" & "diff" vars at levels specified by    'L.new_vars',
  #   pre-merged state of vars. merged at level         'L.premerge',
  #   & return basis fn. showing wts. for var. at level 'L.basis'
  #
  # 'hier' is a list output by 'group_HierPCA()' w/ elements:
  #   originally input data       'hier[[1]]$X' 
  #   & local rotation matrices   'hier[[h]]$J_l'
  # Called as part of 'get_hier_dat()'
  #
  
  if (verbose){  cat('\n......reconstructing hierarchy levels...')  }
  
  ### Sanity checks ###
  stopifnot(length(hier) > 0)
  stopifnot(all(c('X', 'k1', 'k2', 'J_l') %in% names(hier[[1]])))
  stopifnot(!all(is.na(L.active)) || !all(is.na(L.inactive)) ||
              !all(is.na(L.new_vars)) || !all(is.na(L.premerge)) ||
              !all(is.na(L.basis)))
  L = length(hier)
  hier.l = unlist(lapply(hier, function(yy) return(yy$h.level)))
  if (!all(is.na(L.active))){    stopifnot(max(L.active) <= max(hier.l))  }
  if (!all(is.na(L.inactive))){  stopifnot(max(L.inactive) <= max(hier.l))  }
  if (!all(is.na(L.new_vars))){  stopifnot(max(L.new_vars) <= max(hier.l))  }
  if (!all(is.na(L.premerge))){  stopifnot(max(L.premerge) <= max(hier.l))  }
  if (!all(is.na(L.basis))){     stopifnot(max(L.basis) <= max(hier.l))  }
  stopifnot(all(hier.l[2:L] - hier.l[1:(L-1)] == 1)) # check for sequential levels
  
  ### Prelims ###
  X = hier[[1]]$X
  active_inds = c(1:dim(X)[1])
  inactive_inds = c()
  if (!all(is.na(L.basis))){  B_l = diag(rep(1, nrow(X)))  }
  
  ### Level 0: store pre-merge vars. ###
  if (!all(is.na(L.premerge)) && (1 %in% L.premerge) && (1 == hier[[1]]$h.level)){
    k1k2 = c(hier[[1]]$k1, hier[[1]]$k2)
    hier[[1]]$premerge.vars = t(X[k1k2, ])
    hier[[1]]$premerge.vars.i = k1k2
  }
  
  ### Calc. rotations & merge vars. for each level ###
  for (h in 1:length(hier)) {
    l = hier[[h]]$h.level
    if (l > max(c(L.active, L.inactive, L.new_vars, L.premerge), na.rm=T)){  break  }
    k1 = hier[[h]]$k1
    k2 = hier[[h]]$k2
    
    a11 = var(X[k1,])
    a22 = var(X[k2,])
    a12 = cov(X[k1,], X[k2,])
    a12 = abs(a12)
    
    ### Merge vars. ###
    mXk12 = rowMeans(X[c(k1, k2),])
    s1s2 = sqrt(c(a11, a22))
    X_k1k2 = sweep(X[c(k1, k2),], 1, mXk12, FUN = "-")
    X_k1k2 = sweep(X[c(k1, k2),], 1, s1s2, FUN = "/")
    X[c(k1, k2),] = t(hier[[h]]$J_l) %*% X_k1k2
    mXk12 = mXk12 / s1s2
    
    ### Update matrix of basis functions ###
    if (!all(is.na(L.basis))){
      B_l[,c(k1,k2)] = B_l[,c(k1,k2)] %*% hier[[h]]$J_l
    }
    
    ### Add merged sum & diff. vars to info ###
    if (!all(is.na(L.new_vars)) && (l %in% L.new_vars)){
      hier[[h]]$new.vars = t(X[c(k1, k2), ])  # sum & diff. vars, respectively, in Lee et al. (2008) eq. (4)
      hier[[h]]$new.vars.i = c(k1, k2)
      
      ### Add basis fns. for above ###
      if (!all(is.na(L.basis)) && any(l == L.basis)){
        # Include scaling & detail fns. associated w/ sum & diff. vars, respectively, in Lee et al. (2008) eq. (4)
        hier[[h]]$basis.fns = B_l[,c(k1,k2)]
      }
    }
    
    ### Add pre-merged vars. merged in next level to info ###
    if (!all(is.na(L.premerge)) && ((l+1) %in% L.premerge)){
      k1k2 = c(hier[[l+1]]$k1, hier[[l+1]]$k2)
      hier[[l+1]]$premerge.vars = t(X[k1k2, ])
      hier[[l+1]]$premerge.vars.i = k1k2
    }
    
    if (!all(is.na(L.active)) || !all(is.na(L.inactive))){
      active_inds = active_inds[ -which(active_inds == k2) ]
      inactive_inds = c(inactive_inds, k2)
      
      if (!all(is.na(L.active)) && !all(is.na(L.inactive)) &&
          any(l == L.active) && any(l == L.inactive)){
        if (length(c(active_inds, inactive_inds)) > 1000){
          cat(paste0('\n   WARNING: greater than 1000 vars. added at level = ',l,
                     '\n      will be very computationally intensive to process!'))
        }
        if ('new.vars' %in% names(hier[[h]])){  hier[[h]]$new.vars = NULL  }
        if ('new.vars.i' %in% names(hier[[h]])){  hier[[h]]$new.vars.i = NULL  }
        
        all_inds = sort(c(active_inds, inactive_inds))
        hier[[h]]$new.vars = t(X[all_inds,, drop=F])
        hier[[h]]$new.vars.i = all_inds
        hier[[h]]$active_inds = active_inds
        hier[[h]]$inactive_inds = inactive_inds
        
        if (!all(is.na(L.basis)) && any(l == L.basis)){
          # include Scaling &/or detail fns. associated w/ above vars in Lee et al. (2008) eq. (4)
          hier[[h]]$basis.fns = B_l[,all_inds]
        }
        
      }else if (!all(is.na(L.active)) && any(l == L.active)){
        if (length(active_inds) > 1000){
          cat(paste0('\n   WARNING: greater than 1000 active vars. added at level = ',l,
                     '\n      will be very computationally intensive to process!'))
        }
        if ('new.vars' %in% names(hier[[h]])){  hier[[h]]$new.vars = NULL  }
        if ('new.vars.i' %in% names(hier[[h]])){  hier[[h]]$new.vars.i = NULL  }
        
        hier[[h]]$new.vars = t(X[active_inds,, drop=F])
        hier[[h]]$new.vars.i = active_inds
        hier[[h]]$active_inds = active_inds
        
        if (!all(is.na(L.basis)) && any(l == L.basis)){
          # Include scaling &/or detail fns. associated w/ above vars in Lee et al. (2008) eq. (4)
          hier[[h]]$basis.fns = B_l[,active_inds]
        }
        
      }else if (!all(is.na(L.inactive)) && any(l == L.inactive)){
        if (length(inactive_inds) > 1000){
          cat(paste0('\n   WARNING: greater than 1000 inactive vars. added at level = ',l,
                     '\n      will be very computationally intensive to process!'))
        }
        if ('new.vars' %in% names(hier[[h]])){  hier[[h]]$new.vars = NULL  }
        if ('new.vars.i' %in% names(hier[[h]])){  hier[[h]]$new.vars.i = NULL  }
        
        hier[[h]]$new.vars = t(X[inactive_inds,, drop=F])
        hier[[h]]$new.vars.i = inactive_inds
        hier[[h]]$inactive_inds = inactive_inds
        
        if (!all(is.na(L.basis)) && any(l == L.basis)){
          # Include scaling &/or detail fns. associated w/ above vars in Lee et al. (2008) eq. (4)
          hier[[h]]$basis.fns = B_l[,inactive_inds]
        }
      }
    }
  }
  
  ##################################################################
  return(  hier  )
} ##################################################################


##################################################################
smallestEigenval_test <- function(X = NA, eigs = NA, scale = F,
                                  n = NA, g = 0.1, m = 1, pval = 0.001, 
                                  return.pca = T, return.PCs = F, 
                                  verbose = F){
  ##################################################################
  # Smallest Eigenvalue Test,
  #   to test for sig. of sum of smallest m+1 eigenvalues:
  #
  # H_0: (l_m+1 + l_m+2 + ... l_p)   > g
  #      sum of smallest eigenvalues > g
  # H_A: sum of smallest eigenvalues < g
  #   first m principal components represent all measurements
  #
  #   Anderson, T.W. (2003). An introduction to multivariate 
  #     statistical analysis (3rd ed). Wiley, New Jersey. p. 480.
  #
  
  g = as.numeric(g)
  m = as.integer(m)
  pval = as.numeric(pval)
  
  if (is.numeric(X)){
    ### Calculate eigenvalues & test input data matrix ###
    X = as.matrix(X)
    pca = prcomp(X, 
                 retx=return.PCs, center = T, scale. = scale)
    n = nrow(X) - 1
    p = length(pca$sdev)
    l = pca$sdev[(m+1):p]^2  # smallest eigenvalues
    eigs = pca$sdev^2
    
  }else if (is.numeric(eigs)){
    ### test input eigenvalues & sample size ###
    stopifnot(length(eigs) >= m)
    stopifnot(!is.na(n) && is.numeric(n))  # need sample size
    n = n - 1
    
    pca = list()
    p = length(eigs)
    l = eigs[(m+1):p]
  }
  
  z = qnorm(1 - pval)            # upper critical point of N(0,1), w/ sig. level pval
  cr = sqrt(2 * sum(l^2) / n)
  critical_value = g - z * cr    # Anderson (2003) p. 479 eq. (3)
  stat = sum(l)
  test = stat < critical_value   # rejection region is less than critical point
  p.exact = pnorm((stat - g) / cr) # exact p-value, by solving for z in eq. (3)
  
  # Upper bound on Confidence Interval for sum of eigenvalues in population, with approximate confidence (1 - pval)
  CI.upper = sum(l) + z *cr      # Anderson (2003) p. 479 eq. (4)
  
  if (verbose){
    cat(paste0('\nSmallest Eigenvalue Test:'))
    cat(paste0('\n  Sum of smallest eigenvalues = ', signif(sum(l))))
    cat(paste0('\n  Sum of all eigenvalues = ', signif(sum(eigs))))
    cat(paste0('\n  C.I. upper bound for sum of smallest eigenvalues = ',
               signif(CI.upper)))
    cat(paste0('\n  Confidence level for upper bound = ', 1 - pval))
    cat(paste0('\n  Test statistic:  ', signif(stat)))
    cat(paste0('\n  Critical value: test stat.  < ', signif(critical_value)))
    cat(paste0('\n  Null hypothesis rejected at p<',pval,': ', test))
  }
  
  pca$lambda = eigs                     # eigenvalues of covariance matrix
  pca$lambda_sum = l                    # sum of smallest eigenvalues
  pca$lambda_sum.gamma = g              # input cutoff for sum
  pca$lambda_sum.m = m                  # number of largest eigenvalues to keep, if H_0 is rejected
  pca$lambda_sum.pvalue.cutoff = pval   # sig. level of test
  pca$lambda_sum.pvalue.exact = p.exact # exact p-value of test
  pca$lambda_sum.statistic = stat       # test statistic
  pca$lambda_sum.test   = test          # True if null hypothesis is rejected
  pca$lambda_sum.CI.upper = CI.upper    # upper bound of confidence interval for sum
  pca$lambda_sum.confidence = 1 - pval  # confidence level for above upper bound
  
  ##################################################################
  if (return.pca){  return(pca)  }else{  return(test)  }
} ##################################################################

################################################
FDR <- function(pvals, FDR.base, verbose = F){
  ################################################
  # Calculates p-value threshold for FDR,
  #  assumes positive dependence of variables.
  #
  # INPUT: vector of p-values & corrected false-discovery rate (q-value)
  # OUTPUT: p-value cuttoff threshold for given q-value
  #
  #		Benjamini & Hochberg (1995), Controlling
  #	the False Discovery Rate: a practical and powerful
  #	approach to multiple testing.  J.R. Statist. Soc. B.
  #	57(1), p. 289-300.
  
  if (verbose){  cat(paste0('\nApplying FDR to ',length(pvals),' tests, base q = ',FDR.base))  }
  
  q <- FDR.base
  m <- length(pvals)
  cutoff <- 0
  
  p <- sort(pvals)
  for (ii in 1:m){
    if (p[ii] <= q * ii/m){
      cutoff <- p[ii]
    }
  }
  if (verbose){  cat(paste0('\n...   FDR cutoff:  p < ', signif(cutoff,3)))  }
  return(cutoff)
} ################################################


##################################################################
calc_Jacobi <- function(a12, a11 = NA, a22 = NA){
  ##################################################################
  # Calculates numerically-stable Jacobi matrix,
  #   given input covariance a12 & variances a11, a22.
  #   Anderson, E. (2000). Discontinuous plane rotations and the symmetric eigenvalue problem.
  #
  
  a12 = as.numeric(a12)
  if (isTRUE(all.equal(a12, 0))){
    c = 1  # return identity matrix since angle = 0
    s = 0
  }else if (isTRUE(all.equal(a11, a22))){
    c = cos(pi / 4)  # max angle, since |angle| <= pi /4
    s = sin(pi / 4)
    if (a12 < 0){  s = -s  }
  }else if (is.numeric(a11) && is.numeric(a22)){
    tau = (a22 - a11) / (2 * a12)
    u = sqrt(1 + tau^2)
    t = min(abs(c(-tau + u, -tau - u)))
    if (a22 > a11){  # preserve ordering of input variances
      s = 1 / sqrt(1 + t^2)
      c = s*t
    }else{
      c = 1 / sqrt(1 + t^2)
      s = c*t
    }
    if (a12 < 0){  s = -s  }
  }else{  return(NA)  }
  
  return(matrix(c(c, -s,
                  s, c), 2, 2, byrow=T))
} ##################################################################



##################################################################
hierarchicalPCA <- function (X, 
                             sim.prev,
                             levels.table.save = NA,
                             verbose = F) {
  ##################################################################
  #  Implements the hierarchical Treelets PCA algorithm of Lee et al. (2008),
  #     & applies stats. on eigenvalues to flag levels merging redundant vars.
  #
  # Input:
  #   X           : data matrix      [vars. x obs.]
  #   sim.prev    : initial similarity matrix between vars., used to construct additional levels of hierarchy,
  #                   output from calc_similarities().
  #   levels.table.save : save path or filename for saved hierarchy aggregation by level as .csv table,
  #                         with cols. formatted as:  k1, k2, similarity.
  #
  # Output:  list w/ elements:
  #   X_l          : coordinate vars. for final level, created by:
  #                       X_l = ... %*% t(J_2) %*% t(J_1) %*% X0
  #   similarity   : current similarity matrix (or path to saved .csv), used to contruct additional levels
  #   levels_table : path to saved .csv table, or matrix with merging info, structured as:
  #                   1st col.:  indice of 1st comp., used in higher levels
  #                   2nd col.:  indice of 2nd comp., stored & unchanged
  #                   3rd col.:  similarity value (~corr.)
  #                   4th col.:  1st/leading eigenvalue for level PCA
  #                   5th col.:  2nd eigenvalue for level PCA
  #                   6th col.:  constant c on diagonal of level rotation matrix J_l
  #                   7th col.:  constant s on diagonal of level rotation matrix J_l
  #
  # Requires:
  #   calc_similarities()  : calculates similarity matrix by blocks, from calc_Similarities_BigData_fns.R
  #   smallestEigenval_test() : tests relevance of PCs, if used as stopping or pausing criteria
  #   calc_Jacobi()           : numerically stable calculation of level Jacobi rotation matrices
  #
  
  ### Prelims ###
  p = dim(X)[1]
  n = dim(X)[2]
  levels.max = p - 1
  
  if (!all(is.na(sim.prev))){
    stopifnot(is.list(sim.prev))
    stopifnot(all(c('similarity_matrix', 'files', 'files_vars') %in% names(sim.prev)))
    if (verbose){
      cat('\n...similarity matrix detected as saved file:  ')
      cat(paste0('\n      ', sim.prev$files))
      cat('\n   ......starting hierarchical PCA construction using data in above...')
    }
    sim_hpca = sim.prev
  }
  
  if (all(c('files', 'files_vars') %in% names(sim_hpca)) && !all(is.na(sim_hpca$files))){
    stopifnot(length(sim_hpca$files) == length(sim_hpca$files_vars))
    sim_files      = sim_hpca$files
    sim_files_vars = sim_hpca$files_vars
    sim_files_maxs = sim_hpca$files_sims_max
    
    if (verbose){  cat('\n...reading similarity table files:')  }
    for (f in 1:length(sim_files)){
      if (verbose){  cat(paste0('\n      ', sim_files[f]))  }
      if (!all(is.na(sim_files_vars[f]))){
        assign(sim_files_vars[f], read.csv(sim_files[f], header=F, stringsAsFactors=F))
        stopifnot(exists(sim_files_vars[f]))
        sim_files_maxs[f] = max(get(sim_files_vars[f])[,1], na.rm=T)
      }
    }
  }
  if (!all(is.na(sim_files_vars))){
    for (v in sim_files_vars){
      if (!is.na(v) && !exists(v)){
        cat(paste0('\n\n...WARNING: cannot find loaded sim. table piece in R environment:  ', v,'\n\n'))
        stopifnot(exists(v))
      }
    }
  }
  
  ### Initiatize & pre-allocate new levels table for requested levels ###
  levels_table = as.data.frame(matrix(0, levels.max, 7))
  names(levels_table) = c('Index1', 'Index2', 'Similarity', 'Lambda1', 'Lambda2', 'c', 's')
  mode(levels_table$Index1) = 'integer'  # ind. for k1, becomes ind. for PC1
  mode(levels_table$Index2) = 'integer'  # ind. for k2, becomes ind. for PC2
  
  ### Pre-allocate for main similarity table(s) ###
  # NOTE: main 'sim_table' stores similarities for all new pairs of vars., separate from old pieces stored in 'sim_files'
  k2.prev = c()  # indices for inactive "diff" vars. (PC2), already incorporated into hierarchy
  n.active.inds = p
  n.pairs.new = sum(c((n.active.inds - 2):0)[1:levels.max])  # number of possible new pairs
  sim_table = matrix(NA, n.pairs.new, 3)
  
  
  #####################################################
  ### Calc. rotations & merge vars. for each level  ###
  #####################################################
  h.final = F
  if (!exists('l.start')){  l.start = 0  }  # starting level of hierarchy
  old.k1k2 = c()   # indices of outdated k1 & k2 indices in sim. table files
  for (h in 1:levels.max){
    if (verbose){  cat(paste0("\nLevel ", h,"  /  ", levels.max,"  :"))  }

    ### Get max. similarity & var. indices ###
    if ((!all(is.na(sim_table)) || is.ff(sim_table)) && (nrow(sim_table) > 0) &&
        (max(sim_table[,1], na.rm=T) >= max(sim_files_maxs, na.rm=T))){
      r.max.i = which.max(sim_table[,1])   # find max. sim. from main sim_table, stored in memory
      sk1k2 = sim_table[r.max.i, ]
      s.max = as.numeric(sk1k2[1])
    }else{  # find max. similarity from stored pieces of sim. table, stored in ff files on disk
      f = which.max(sim_files_maxs)
      s.max = as.numeric(sim_files_maxs[f])
      r.max.i = which(get(sim_files_vars[f])[,1] == s.max)
      sk1k2 = get(sim_files_vars[f])[r.max.i, ]
    }
    k1 = as.integer(unname(unlist(sk1k2[2])))
    k2 = as.integer(unname(unlist(sk1k2[3])))
    old.k1k2 = sort(unique(c(old.k1k2, k1, k2)))
    
    if (verbose){
      cat(paste0("      k1 = ",k1,",  k2 = ",k2,",  max. similarity = ",signif(s.max, 3)))
    }
    
    
    #----------------
    ### Local PCA ###
    #----------------
    a11 = var(X[k1,])
    a22 = var(X[k2,])
    a12 = cov(X[k1,], X[k2,])
    
    mXk12 = rowMeans(X[c(k1, k2),])
    s1s2 = sqrt(c(a11, a22))
    X_k1k2 = sweep(X[c(k1, k2),], 1, mXk12, FUN = "-")
    X_k1k2 = sweep(X[c(k1, k2),], 1, s1s2, FUN = "/")
    
    ### Get PCA & Eigenvalue stats ###
    pca = prcomp(t(X_k1k2), retx = F, center = T, scale. = T)
    
    ### Get PCA rotation matrix & update data ###
    #     Theory: t(X) * J = PCs in cols., where X is subset X[c(k1,k2),]
    #           t(PCs) = t(J) * X to obtain new vars. as PCs, by merging rows of X.
    #       Thus, t(J) is the mixing matrix for the rows of X.
    J = calc_Jacobi(a12, 1, 1)
    X[c(k1, k2),] = t(J) %*% X[c(k1, k2),]
    mXk12 = mXk12 / s1s2
    
    ### Update levels aggregation table ###
    levels_table$Index1[h] = k1
    levels_table$Index2[h] = k2
    levels_table$Similarity[h] = s.max
    levels_table$Lambda1[h] = pca$sdev[1]^2
    levels_table$Lambda2[h] = pca$sdev[2]^2
    levels_table$c[h] = J[1,1]
    levels_table$s[h] = J[2,1]
    
    
    ### Remove old indices from similarity table(s) of vars. created in current run ###
    old.inds = which((sim_table[,2] %in% c(k1,k2)) | (sim_table[,3] %in% c(k1,k2)))
    sim_table[old.inds,1] = NA
    
    if (!all(is.na(sim_files_vars))){
      new_maxs = rep(NA, length(sim_files_vars))
      for (f in 1:length(sim_files_vars)){
        r.max = -Inf
        old.inds = (get(sim_files_vars[f])[,2] %in% c(old.k1k2)) | (get(sim_files_vars[f])[,3] %in% c(old.k1k2))
        r.max = max(-Inf, get(sim_files_vars[f])[which(!old.inds), 1], na.rm=T)
        new_maxs[f] = r.max
      }
      for (f in 1:length(sim_files_vars)){
        sim_files_maxs[f] = new_maxs[f]  # will be used to index similarities, w/o removing prev. entries
      }
    }
    
    ### Update similarity table w/ new var. ###
    #             NOTE:             # k1 always contains PC1 "sum" after above sorting, regardless of ranking & order of inds
    k2.prev = sort(c(k2.prev, k2))  # k2 always contains PC2 "diff" var., not incorporated in higher levels of hierarchy
    k.update = c(1:p)               # indices of vars. to update similarities w/ k.select...
    k.update = k.update[ -which(k.update %in% c(k1, k2.prev)) ]  # ...skip calc. sim. w/ self & w/ diff. vars. from prev. levels
    if (length(k.update) == 0){
      if (verbose){  cat('\n\n...final level complete, no remaining vars. in input!')  }
      h.final = T
      break
    }else{
      # if (verbose){  cat('\n...calculating similarities w/ new var...')  }
      r = cor(t(X[k.update,, drop=F]), X[k1,])

      if (all(is.na(sim_table[,2]))){
        r.start.i = 1
        r.end.i = length(k.update)
      }else{
        r.start.i = max(which(!is.na(sim_table[,2]))) + 1
        r.end.i = r.start.i - 1 + length(k.update)
      }
      sim_table[r.start.i:r.end.i, 1] = c(r)
      sim_table[r.start.i:r.end.i, 2] = k.update
      sim_table[r.start.i:r.end.i, 3] = rep(k1, length(k.update))
    }
  } #------------  (end of loop for hierarchy level)  -------------#
  
  
  #------------------------
  ### Saving & clean-up ###
  #------------------------
  if (verbose){  cat('\n.....final level reached, cleaning up similarity table files(s)...')  }
  sim_table = NA
  
  if (!all(is.na(sim_files_vars))){
    for (f in 1:length(sim_files_vars)){
      if (is.na(sim_files[f]) || is.na(sim_files_vars[f])){
        next
      }else if (is.character(sim_files[f])){
        if (file.exists(sim_files[f])){
          file.remove(sim_files[f])
        }else if (verbose){
          cat(paste0('\n  WARNING: could not clean-up sim. table file:  ', sim_files[f]))
          next
        }
        sim_files[f] = NA
        sim_files_vars[f] = NA
        sim_files_maxs[f] = -Inf
      }
    }
  }
  
  if (any(levels_table$Index1 == 0)){
    empty_rows = levels_table$Index1 == 0
    inds = which(empty_rows)
    
    if (all(diff(inds) == 1) && (max(inds) == nrow(levels_table))){
      if (verbose){  cat('\n......trimming unused pre-allocated rows from current levels table...')}
      
      ### Trim unused rows from levels table ###
      levels_table_slimmed = as.data.frame(matrix(0, sum(!empty_rows), 7))
      names(levels_table_slimmed) = names(levels_table)
      mode(levels_table_slimmed$Index1) = 'integer'
      mode(levels_table_slimmed$Index2) = 'integer'
      levels_table_slimmed$Index1 = as.integer(levels_table$Index1[!empty_rows])
      levels_table_slimmed$Index2 = as.integer(levels_table$Index2[!empty_rows])
      levels_table_slimmed$Similarity = levels_table$Similarity[!empty_rows]
      levels_table_slimmed$Lambda1 = levels_table$Lambda1[!empty_rows]
      levels_table_slimmed$Lambda2 = levels_table$Lambda2[!empty_rows]
      levels_table_slimmed$c = levels_table$c[!empty_rows]
      levels_table_slimmed$s = levels_table$s[!empty_rows]
      
      levels_table = levels_table_slimmed
    }
  }
  if (!all(is.na(levels.table.save)) && is.character(levels.table.save)){
    if (!grepl('.csv$', levels.table.save)){  levels.table.save = paste0(levels.table.save, '.csv')  }
    if (!dir.exists(dirname(levels.table.save))){  dir.create(dirname(levels.table.save), recursive=T)  }
    if (verbose){  cat(paste0('\n...saving & appending hierarchy levels table:\n      ', levels.table.save))  }
    
    write.csv(levels_table, levels.table.save, row.names=F)
    if (file.exists(levels.table.save)){
      levels_table = levels.table.save  # return path to saved aggregation table to min. mem. use
    }else{
      print(paste0('WARNING: levels table .csv not saved, could not find  ', levels.table.save))
    }
  }
  
  ##################################################################
  return(list('X_l' = X,
              'similarity' = list('similarity_matrix' = sim_table,
                                  'files'  = sim_files,
                                  'files_vars' = sim_files_vars,
                                  'files_sims_max' = sim_files_maxs),
              'levels_table' = levels_table))
} ##################################################################



##################################################################
group_HierBackReconstruction <- function(groupHier_path, prefix = NA, 
                                         align_sources = T,
                                         separate_spatial_maps = T,
                                         separate_time_series = F,
                                         hier_levels = NA,
                                         hier_stats_filter = F,
                                         level_highest_redundant = F,
                                         levels_active_all = F,
                                         levels_PC1 = T,
                                         levels_PC2 = F,
                                         levels_show_mergers = F,
                                         verbose = F,
                                         ...){
  ##################################################################
  # Subj.-specific source back-reconstruction for all subjs. in group
  #   & selected levels of a hierarchical analysis
  # Estimates subj-level spatial maps & time series
  #
  # Input:
  #   groupHier_path : Path to dir. w/ saved outputs of 
  #                     subj_PCA(), group_PCA(), group_HierPCA().
  #   prefix         : Identifying prefix to attach to saved output.
  #   separate_spatial_maps : if T, separate hPCA comp. spatial maps that are independent w/n hierarchy,
  #                             resulting in maximally distinct & non-overlapping spatial maps.
  #   separate_time_series  : if T, separate hPCA comp. time series that are independent w/n hierarchy,
  #                             resulting in maximally orthogonal time series,
  #                         (i.e., for each requested comp. & given subset of hierarchy, 
  #                           remove effects of returned leaf nodes that are not ancestors of comp.),
  #                       if F, hPCA comp. spatial maps are constructed from hPCA time series as in seed-FC.
  #   align_sources     : Align all hPCA comps. signs to a positive 3rd moment (skew).
  #   hier_levels         : Reconstruct comps from select level of hierarchy. Options
  #                         "final"            : reconstructs _all_ comps from highest level only
  #                         vector of integers : reconstructs select comps from specific levels.
  #   hier_stats_filter   : if T, filter reconstruct levels by applying statistical inference testing,
  #                           else back-reconstruct levels as above, regardless of any stats testing applied.
  #   level_highest_redundant : if T, reconstruct all active vars. in from last redundant level by stats testing,
  #                               essentially summarizing previous redundant mergers in early hierarchy levels.
  #   levels_active_all   : Reconstruct all active vars. at requested levels.
  #   levels_PC1          : Reconstruct primary comp. (sum var.) at each level.
  #   levels_PC2          : Reconstruct secondary comp. (difference var.) at each level.
  #   levels_show_mergers : Shows merger created by each level, by displaying (in order, for each level):
  #                         1. 1st in pairs of vars. merged, in pre-merged state from prev. level
  #                         2. 2nd in pairs of vars. merged, in pre-merged state from prev. level
  #                         3. "sum" var., created at level by combining variances of above pair of vars. (PC1)
  #                         4. "diff" var., created by subtracting variances of above pair of vars. (PC2).
  # Output: 
  #   saved files w/ *_ICAsources.RData added as suffix & vars.:
  #     S_i  : estimated subj.-specific source spatial maps as [voxels x sources] matrix
  #     R_i  : estimated subj.-specific source time series  as [time x sources] matrix
  #     hier_inds             : indices of hierarchical ICA/PCA sources w/n above matrices
  #     Hierarchy_level       : vector of levels of the hierarchy for b.r. sources (0=initial ICA data)
  #     hier_k1_inds          : index of primary new comp., created at each level
  #     hier_k2_inds          : index of secondary new comp., created at each level & stored thereafter
  #     hier_ns_inds          : indices of initial/redundant/n.s. sources, summarizing early mergers
  #     hier_algorithm        : algorithm used to construct hierarchy, hPCA = hierarchical PCA
  # Requires:
  #   load.data()          : data munging from PCA_fns.R
  #   zeroMean_Yi()        :  "      "      "     "
  #   flatten_img()        :  "      "      "     "
  #   detrend()            : detrending for pre-processing, from PCA_dataCompression_fns.R
  #   subj_preprocData()   : pre-processes data used for subj.-level PCA, from PCA_dataCompression_fns.R
  #   get_groupPCA_dat()   : retrieves & formats group-level PCA info, from PCA_dataCompression_fns.R
  #   get_subjPCA_dat()    : retrieves & formats subj.-level PCA info, from PCA_dataCompression_fns.R
  #   subj_HierBackReconstruct() : calculates & formats sources' from levels of hPCA
  #   hier_find_dependent_levels()   : finds dependent levels w/n filtered levels of hierarchy
  #   hier_find_leaves()             : finds leaf nodes in hierarchy
  #   hier_find_independent_levels() : finds hierarchically independent vars., from different branches
  #   hier_separate_sources()        : separates sources using orthogonal projection
  #
  
  if (verbose){  cat('\nStarting subj.-level source back-reconstruction...')  }
  
  ### Optional additional inputs ###
  #     Reconstruct scaling & detail basis fns. associated w/ analyzed vars.,
  #       saved in 'hier' w/n file output from group_HierPCA().
  #     Can be computationally intensive & may result in a large file, 
  #       if hierarchy is deep and many levels are requested
  reconstruct_basis = F  # recommended for hierarchy constructed from voxels, avoids large [voxels x voxels] matrix
  # reconstruct_basis = T  # potentially useful for smaller hierarchies
  
  #     Reconstructing all vars. at lower levels can be excessive,
  #       input vector of integers to arg below can be used to reconstruct & return select subset of requested vars.
  levels_subset_PCs = F  # if F, back-reconstruct all requested vars. at each level
  # levels_subset_PCs = c(20:30)  # ...ex.: reconstruct a subset of vars. 20:30 at each level
  
  args = list(...)
  if ('reconstruct_basis' %in% names(args)){
    reconstruct_basis = as.logical(args$reconstruct_basis)
  }
  if ('levels_subset_PCs' %in% names(args)){
    levels_subset_PCs = as.integer(args$levels_subset_PCs)
    if (!as.logical(level_highest_redundant) || !as.logical(levels_subset_PCs)){
      levels_subset_PCs = F
    }
  }
  
  ### Inputs ###
  groupHier_path = as.character(groupHier_path) 
  prefix = as.character(prefix) 
  separate_spatial_maps = as.logical(separate_spatial_maps)
  separate_time_series = as.logical(separate_time_series)
  separate_sources = separate_spatial_maps || separate_time_series
  align_sources = as.logical(align_sources)
  if (!all(is.na(hier_levels))){
    if (all(hier_levels == F)){
      hier_levels = NA
    }else if (is.character(hier_levels)){
      stopifnot((hier_levels %in% c('final', 'last', 'all')))
    }else if (is.numeric(hier_levels)){
      hier_levels = as.integer(hier_levels)
    }
  }
  hier_stats_filter = as.logical(hier_stats_filter)
  level_highest_redundant = as.logical(level_highest_redundant)
  levels_active_all = as.logical(levels_active_all)
  levels_PC1 = as.logical(levels_PC1)
  levels_PC2 = as.logical(levels_PC2)
  levels_show_mergers = as.logical(levels_show_mergers)
  verbose = as.logical(verbose)
  
  if (verbose){
    cat('\n...for hierarchical sources...')
    if (level_highest_redundant){
      cat('\n......including all sources in initial non-sig. levels of hierarchy, after merging redundant vars....')
    }
    if (is.character(hier_levels) && ((hier_levels == 'final') || (hier_levels == 'last'))){
      cat('\n......including all sources at last (highest) level of hierarchy...')
    }else if (is.character(hier_levels) && (hier_levels == 'all')){
      cat('\n......including hierarchical sources at every level...')
      cat('\n.........(Warning: may be computationally intensive!)')
    }else if (is.numeric(hier_levels) && (length(hier_levels) > 0)){
      cat('\n......including sources at select level(s):  ')
      cat(paste(hier_levels, collapse = ','))
    }
    if (levels_active_all){
      cat('\n.........including all active vars. at above levels...')
      cat('\n.........(Warning: may be computationally intensive!)')
    }
    if (levels_show_mergers){
      cat('\n......including 4 sources per level of hierarchy, showing merger created at each level, in order:')
      cat('\n.........1. 1st in pair of pre-merged vars., from previous level')
      cat('\n.........2. 2nd in pair of pre-merged vars., from previous level')
      if (levels_PC1){
        cat('\n.........3. merged "sum" var., created at hierarchy level')
      }
      if (levels_PC2 && levels_PC1){
        cat('\n.........4. merged "difference" var., created at hierarchy level')
      }else if (levels_PC2 && !levels_PC1){
        cat('\n.........3. merged "difference" var., created at hierarchy level')
        cat('\n             (note: "sum" var., showing similarities of merged vars., not requested!')
      }
    }else if (levels_PC1 && levels_PC2){
      cat('\n......including 2 sources created at each sig. level of hierarchy, in order:')
      cat('\n.........1. merged "sum" var., showing shared variance (PC1)')
      cat('\n.........2. merged "difference" var., showing difference in variances (PC2)')
    }else if (levels_PC1){
      cat('\n......including new primary source (i.e., sum of variances for merged vars.) at each hPCA level...')
    }else if (levels_PC2){
      cat('\n......including new secondary source (i.e., difference of variances for merged vars.) at each hPCA level...')
    }
    if (hier_stats_filter){
      cat('\n......after filtering out levels merging redundant/unrelated vars. with statistical testing...')
    }else{
      cat('\n......including non-sig. levels of hierarchy, if stats. testing was applied')
    }
    if (any(as.logical(levels_subset_PCs))){
      cat('\n......while limiting each level to a subset of all requested vars.')
    }
  }
  if (separate_sources){
    if (separate_time_series){
      cat('\n......while separating all returned hPCA component time series...')
    }
    if (separate_spatial_maps){
      cat('\n......while separating all returned hPCA component spatial map...')
    }
    cat('\n         by orthogonally projecting hPCA components onto all hierarchically independent leaf nodes')
  }
  cat('\n')
  
  ### Defaults ###
  subj_pca_suffix  = '_PCA1.RData'    # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData' # filesuffix for group-level temporal PCA from group_PCA()
  hier_suffix = 'Hierarchy.RData'   # filesuffix for hierarchical PCA from group_HierPCA()
  br_suffix   = '_HierComps.RData'  # filesuffix for back-reconstructed subj. sources  
  
  hier_wc = hier_suffix
  group_pca_wc = group_pca_suffix
  if (!all(is.na(prefix))){
    hier_wc = paste0(prefix, '.*_', hier_wc)
    group_pca_wc = paste0(prefix, '.*_', group_pca_wc)
  }
  
  save_path = groupHier_path         # save to same dir.
  
  ### Sanity checks ###
  stopifnot(is.character(groupHier_path) && dir.exists(groupHier_path))
  stopifnot(length(list.files(groupHier_path, hier_wc)) > 0)
  if (length(list.files(groupHier_path, hier_wc)) > 1){
    cat(paste0('WARNING: multiple hierarchical ICA runs detected in ', groupHier_path,
               ',\n only the 1st will be processed'))
  }
  
  ### Prelims ###
  apply.stats = hier_stats_filter | level_highest_redundant
  if (file.exists(groupHier_path) && !dir.exists(groupHier_path)){
    hier_file = groupHier_path
  }else{
    hier_file  = list.files(groupHier_path, hier_wc, full.names=T)[1]
  }
  stopifnot(file.exists(hier_file))
  
  hier = get_hier_dat(hier_file, 
                      prefix = prefix, 
                      hier_levels = hier_levels,
                      hier_active_all = levels_active_all,
                      return.level.vars = T,
                      return.init.redundant.vars = level_highest_redundant,
                      return.premerge.vars = levels_show_mergers,
                      return.level.basis.fns = reconstruct_basis,
                      apply.stats = apply.stats,
                      verbose = verbose)
  
  if (length(hier) > 0){
    group_pca_file = hier[[1]]$group_pca_file
    mask           = hier[[1]]$mask
    levels_table   = hier[[1]]$levels_table
    
    if (verbose && !any(names(hier[[1]]) == 'new.vars')){
      cat('\n\n.........Stopping back-reconstruction, no levels/vars. selected in hierarchy based on input parameters!\n\n')
    }
    stopifnot('new.vars' %in% names(hier[[1]]))
    stopifnot('new.vars.i' %in% names(hier[[1]]))
    stopifnot(!all(is.na(group_pca_file)) && file.exists(group_pca_file))
  }else{
    if (verbose){
      cat('\n.........Stopping back-reconstruction, no vars. in hierarchy selected based on input parameters!')
    }
    stopifnot(length(hier) > 0)
  }
  
  
  PCA_g = get_groupPCA_dat(group_pca_file, prefix,
                           return.PCs = F,
                           verbose = verbose)
  subj_pca_files = PCA_g$subj_pca_files
  stopifnot(is.character(subj_pca_files) 
            && all(file.exists(subj_pca_files)))
  S = length(subj_pca_files)
  
  
  ### Additional selection & filtering of un-needed vars. in hierarchy ###
  if (!levels_PC1 || !levels_PC2){
    k.rm = which(c(!levels_PC1, !levels_PC2))
    
    for (h in 1:length(hier)){
      if (!any(c('active_inds', 'inactive_inds') %in% names(hier[[h]]))){
        hier[[h]]$new.vars = hier[[h]]$new.vars[,-k.rm, drop=F]
        hier[[h]]$new.vars.i = hier[[h]]$new.vars.i[-k.rm] 
      }
    }
  }
  if (any(as.logical(levels_subset_PCs))){
    for (h in 1:length(hier)){
      if ((length(hier[[h]]$new.vars.i) >= min(levels_subset_PCs)) &&
          (length(hier[[h]]$new.vars.i) <= max(levels_subset_PCs))){
        if (verbose){
          cat(paste0('\n...level = ',h, ', total vars. = ', length(hier[[h]]$new.vars.i),
                     '\n...... selecting subset of ',length(levels_subset_PCs), 
                     ' vars.:', paste(levels_subset_PCs, collapse=', ')))
          k.rm = c(1:length(new.vars.i))
          k.rm = k.rm[-levels_subset_PCs]
          hier[[h]]$new.vars = hier[[h]]$new.vars[,-k.rm, drop=F]
          hier[[h]]$new.vars.i = hier[[h]]$new.vars.i[-k.rm]
          if (any(c('active_inds', 'inactive_inds') %in% names(hier[[h]]))){
            select.inds = c(hier[[h]]$active_inds, hier[[h]]$inactive_inds)
            select.inds = select.inds[-levels_subset_PCs]
            if ('active_inds' %in% names(hier[[h]])){
              select.act.inds = hier[[h]]$active_inds %in% kept.inds
              hier[[h]]$active_inds = hier[[h]]$active_inds[select.act.inds]
            }
            if ('inactive_inds' %in% names(hier[[h]])){
              select.inact.inds = hier[[h]]$inactive_inds %in% kept.inds
              hier[[h]]$inactive_inds = hier[[h]]$inactive_inds[select.inact.inds]
            }
          }
        }
      }
    }
  }
  
  ### Remove empty/non-selected levels ###
  h.rm = c()
  for (h in 1:length(hier)){
    if (ncol(hier[[h]]$new.vars) < 1){  h.rm = c(h.rm, h)  }
  }
  if (length(h.rm) > 0){  hier = hier[-h.rm]  }
  
  
  ### Find ancestor/descendent nodes w/n selected levels of hierarchy ###
  hier[[1]]$separate_spatial_maps = separate_spatial_maps
  hier[[1]]$separate_time_series = separate_time_series
  if (separate_sources){
    levels_table_file = hier[[1]]$levels_table
    h.levels          = unlist(lapply(hier, function(yy) yy$h.level))
    new.vars.i        = lapply(hier, function(yy) yy$new.vars.i)
    
    h.dependent = hier_find_dependent_levels(levels_table_file, h.levels, new.vars.i, 
                                             verbose = verbose)
    h.leaves = hier_find_leaves(levels_table_file, h.levels, new.vars.i, 
                                verbose = verbose)
    h.independent = hier_find_independent_levels(h.dependent, h.levels, new.vars.i, h.leaves, 
                                                 verbose = verbose)
    hier = hier_separate_sources(hier, h.independent, 
                                 verbose = verbose)
  }
  
  
  ### Group hPCA Back-reconstruction ###
  vrbs1 = verbose    # limit printed output to 1st subj. only
  for (s in 1:S){
    if (verbose){
      cat(paste0('\n...Back-reconstructing sources for:  ', basename(subj_pca_files[s])))
    }
    
    ### Get subject's PCA ###
    PCA_s = get_subjPCA_dat(subj_pca_files[s], prefix,
                            return.PCs = F, return.orig.data = T)
    
    ### Back-reconstruction(s) ###
    br = list('S_i' = NULL, 'R_i'= NULL)
    
    br = subj_HierBackReconstruct(hier,
                                  s, PCA_s, PCA_g,
                                  align_sources,
                                  br, verbose = vrbs1)
    vrbs1 = F
    
    ### Append meta info ###
    br$mask = mask
    br$levels_table = levels_table
    
    ### Append info for excluded time points ###
    if (('discard.n' %in% names(PCA_s)) && 
        is.numeric(PCA_s$discard.n) && (PCA_s$discard.n > 0)){
      br$discard.n = PCA_s$discard.n
    }
    if (('good.timepoints' %in% names(PCA_s)) && 
        is.numeric(PCA_s$good.timepoints) && (length(PCA_s$good.timepoints) > 0)){
      br$good.timepoints = PCA_s$good.timepoints
    }
    if (('total.timepoints' %in% names(PCA_s)) && 
        is.numeric(PCA_s$total.timepoints) && (PCA_s$total.timepoints > 0)){
      br$total.timepoints = PCA_s$total.timepoints
    }
    
    ### Saving ###
    save.fname = basename(subj_pca_files[s])
    if (!is.na(prefix) && !grepl(prefix, save.fname)){
      save.fname = paste0(prefix, save.fname)
    }
    save.fname = sub(subj_pca_suffix, br_suffix, save.fname)
    save.fname = file.path(save_path, save.fname)
    save(file=save.fname, list=c('br'))
  }
  
  if (verbose){  cat(paste0('\n...saved as:  ', save_path, '/*', br_suffix, '\n\n'))  }
} ##################################################################



##################################################################
hier_find_dependent_levels <- function(levels_table_file, h.levels, new.vars.i, verbose = F){
  ##################################################################
  # Finds levels that are dependent within hierarchy,
  #   used to create distinct & non-overlapping spatial maps during b.r.
  #
  # Input:  elements of list "hier" output from get_hier_dat():
  #   levels_table_file : path to saved .csv of merged vars. at each level of hierarchy
  #   h.levels          : vector of select hierarchy levels included in hier
  #   new.vars.i        : list of indices vectors for new vars. included in each level of hier
  # Output:
  #   h.dependent : list of dependencies for each var. w/n each level of hier,  formatted as
  #                   if var. v1 w/n hier list element h1 is dependent on var. v2 level l2, then
  #     h.dependent[[h1]]$levels[[v1]][l2] : indexes dependency on var. w/n level l2
  #     h.dependent[[h1]]$voxels[[v1]][v2] : indexes dependent var. w/n level l2
  #
  if (verbose){
    cat('\n...finding dependent vars. based on sequence of mergers within hierarchy...')
  }
  
  levels_table_file = as.character(levels_table_file)
  h.levels = as.integer(h.levels)
  new.vars.i = as.list(new.vars.i)
  stopifnot(file.exists(levels_table_file))
  H = length(h.levels)
  stopifnot(length(new.vars.i) == H)
  
  h.dependent = vector("list", H)  # list of dependencies, in terms of hier levels & var. indices, for each level w/n input "hier"
  v.current = vector("list", H)  # list of most recent var. inds. in dependency chain for each level w/n input "hier"
  for (h1 in 1:H){
    h.dependent[[h1]] = list('levels' = vector("list", length(new.vars.i[[h1]])),
                             'voxels' = vector("list", length(new.vars.i[[h1]])))
    v.current[[h1]] = rep(NA, length(new.vars.i[[h1]]))
  }
  
  if (H == 1){  # nothing to do, all vars. w/n level are independent w/n hierarchy
    return(h.dependent)  
  }
  
  ### Load table of vars. merged at each level ###
  levels_table = read.csv(levels_table_file, stringsAsFactors=F, header=T)
  stopifnot(is.data.frame(levels_table) || is.numeric(levels_table))
  
  
  ### Find dependencies by collapsing hierarchy & tracking indices ###
  for (l in min(h.levels):max(h.levels)){
    ### Update current var. indices in dependency chains ###
    if (any(levels_table[l,2] %in% unlist(v.current))){  # update 2nd var. index only, 1st var. index unchanged
      for (h1 in 1:H){
        if (any(levels_table[l,2] %in% v.current[[h1]])){
          v2 = which(levels_table[l,2] == v.current[[h1]])
          v.current[[h1]][v2] = levels_table[l,1]  # sum var. for level l indexed by 1st index at subsequent levels
        }
      }
    }
    
    if (l %in% h.levels){
      ### Initialize index vars. if initial appearance of level ###
      h1 = which(l == h.levels)  # current level index
      if (all(is.na(v.current[[h1]]))){
        v.current[[h1]] = new.vars.i[[h1]]
      }
      
      ### Check for dependencies on prev. levels ###
      if (any(v.current[[h1]] %in% unlist(v.current[-h1]))){
        for (h2 in (1:H)[-h1]){
          if (any(levels_table[l,c(1,2)] %in% v.current[[h2]])){
            ### Append dependent level & var. index to appropriate dependency chains ###
            v2 = c(which(levels_table[l,1] == v.current[[h2]]),
                   which(levels_table[l,2] == v.current[[h2]]))
            if (length(v2) > 0){
              for (v1 in 1:length(v.current[[h1]])){
                for (n in 1:length(v2)){
                  h.dependent[[h1]]$levels[[v1]] = c(h.dependent[[h1]]$levels[[v1]], 
                                                     h.levels[h2])  # append dep. level to current node v1
                  h.dependent[[h1]]$voxels[[v1]] = c(h.dependent[[h1]]$voxels[[v1]], 
                                                     v2[n])  # append var. index of dep. node on current node
                  h.dependent[[h2]]$levels[[v2[n]]] = c(h.dependent[[h2]]$levels[[v2[n]]], 
                                                        h.levels[h1]) # append current level to dep. lower-level node
                  h.dependent[[h2]]$voxels[[v2[n]]] = c(h.dependent[[h2]]$voxels[[v2[n]]], 
                                                        v1) # append var. index of current node to dep. lower level node
                }
              }
            }
          }
        }
      }
    }
  }
  
  ##################################################################
  return(h.dependent)
} ##################################################################

##################################################################
hier_find_leaves <- function(levels_table_file, h.levels, new.vars.i, verbose = F){
  ##################################################################
  # Find leaf variables in subset of hierarchy,
  #   based on initial appearance among subset of levels.
  # Used for orthogonal projection of minimal independent subset of vars.
  #
  # Input:  elements of list "hier" output from get_hier_dat():
  #   levels_table_file : path to saved .csv of merged vars. at each level of hierarchy
  #   h.levels          : vector of select hierarchy levels included in hier
  #   new.vars.i        : list of indices vectors for new vars. included in each level of hier
  # Output:
  #   h.leaves : list of vars. acting leaf nodes w/n each level of hier., 
  #               formatted as vector of indices
  # 
  if (verbose){
    cat('\n...finding leaf nodes within subset of hierarchy...')
  }
  
  levels_table_file = as.character(levels_table_file)
  h.levels = as.integer(h.levels)
  new.vars.i = as.list(new.vars.i)
  stopifnot(file.exists(levels_table_file))
  H = length(h.levels)
  stopifnot(length(new.vars.i) == H)
  
  h.leaves = vector("list", H)  # list of leaf vars., in terms of hier levels & var. indices, for each level w/n input "hier"
  v.prev = c()  # vector of previously used var. indices
  
  if (H == 1){  # nothing to do, all vars. are leaf nodes w/n subset of hierarchy
    h.leaves[[H]] = c(1:length(new.vars.i[[H]]))
    return(h.leaves)  
  }
  
  ### Load table of vars. merged at each level ###
  levels_table = read.csv(levels_table_file, stringsAsFactors=F, header=T)
  stopifnot(is.data.frame(levels_table) || is.numeric(levels_table))
  
  
  ### Find dependencies by collapsing hierarchy & tracking indices ###
  for (l in min(h.levels):max(h.levels)){
    if (l %in% h.levels){
      ### Find new leaves in hier. level ###
      h1 = which(l == h.levels)  # current level index
      v.new = which(!(new.vars.i[[h1]] %in% v.prev))  # indices of vars. in level h1 not in prev. hier levels
      if (length(v.new) > 0){
        h.leaves[[h1]] = v.new
      }
      
      ### Append new var. indices ###
      v.prev = unique(c(v.prev, new.vars.i[[h1]]))
      
    }else{
      ### Update current dependent var. indices ###
      if (levels_table[l,2] %in% v.prev){  # update 2nd var. index only, 1st var. index unchanged
        v2 = which(levels_table[l,2] == v.prev)
        v.prev = v.prev[-v2]
        v.prev = unique(v.prev, levels_table[l,1])  # sum var. for level l indexed by 1st index at subsequent levels
      }
    }
  }
  
  ##################################################################
  return(h.leaves)
} ##################################################################


##################################################################
hier_find_independent_levels <- function(h.dependent, h.levels, new.vars.i, 
                                         h.leaves = NA, new.k1.vars.i = NA, new.k2.vars.i = NA, 
                                         verbose = F){
  ##################################################################
  # Finds levels that are independent within hierarchy,
  #   using output from get_hier_dat() & hier_find_dependent_levels().
  #
  # Input:
  #   h.dependent   : list of var. & level dependencies, from hier_find_dependent_levels().
  #   h.levels      : vector of select hierarchy levels included in "hier" output from get_hier_dat().
  #   new.vars.i    : list of indices vectors for new vars. included in each level of hier.
  #   h.leaves      : (optional) list of var. indices by level for leaf nodes in subset of hier levels,
  #                     output from hier_find_leaves().
  #                     If included, will limit orthogonal projection to these leaf vars.
  #   new.k1.vars.i : (optional) vector of indices for hPCA sum vars. at each level of hierarchy.
  #   new.k2.vars.i : (optional) vector of inds. for hPCA diff. vars.,
  #                     not used for orthogonal projection since PCA comps. are already orthogonal.
  # Output:
  #   h.independent : list of independent vars. to each var. w/n hier list in other levels, indexed as
  #       h.independent[[h1]][[v1]][[h2]][v2]  where:
  #         h1 : index of hier list level
  #         v1 : index of var. w/n hier list level h1
  #         h2 : hier list index of independent vars. for v1, w/n other hier list levels
  #         v2 : index of independent var. w/n hier level h2, for var. v1 w/n hier list level h1
  #
  
  if (verbose){
    cat('\n......finding independent vars. in separate branches of hierarchy...')
  }
  
  h.dependent = as.list(h.dependent)
  h.levels = as.integer(h.levels)
  new.vars.i = as.list(new.vars.i)
  H = length(h.levels)
  stopifnot(length(new.vars.i) == H)
  stopifnot(length(h.dependent) == H)
  if (!all(is.na(h.leaves))){
    stopifnot(length(h.leaves) == H)
  }
  if (!all(is.na(new.k1.vars.i))){
    stopifnot(length(new.k1.vars.i) == H)
  }
  if (!all(is.na(new.k2.vars.i))){
    stopifnot(length(new.k2.vars.i) == H)
  }
  
  h.independent = vector("list", H)
  for (h1 in 1:H){
    V1 = length(new.vars.i[[h1]])
    h.independent[[h1]] = vector("list", V1)
    for (v1 in 1:V1){  # Create list of independent vars. in other levels
      h.independent[[h1]][[v1]] = vector("list", H)
    }
  }
  
  if (H == 1){  # nothing to do, all vars. w/n level are independent w/n hierarchy
    for (v1 in 1:V1){
      h.independent[[1]][[v1]][[1]] = c(1:V1)[-v1]
    }
    return(h.independent)
  }
  
  
  for (h1 in 1:H){
    V1 = length(new.vars.i[[h1]])
    
    for (v1 in 1:V1){
      for (h2 in 1:H){
        v.ind = c(1:length(new.vars.i[[h2]]))  # by default, assume all vars. ind. w/n hierarchy
        if (h1 == h2){
          v.ind = v.ind[-v1]  # non-zero vars. are not self-orthogonal
          if (!all(is.na(new.k1.vars.i)) || !all(is.na(new.k2.vars.i))){
            # hPCA sum & diff vars. w/n same level are already orthogonal, skip projection
            if ((new.vars.i[[h1]][v1] == new.k1.vars.i[h1]) && 
                any(new.vars.i[[h1]][v.ind] == new.k2.vars.i[h1])){
              v.ind = v.ind[-which(new.vars.i[[h1]][v.ind] == new.k2.vars.i[h1])]
            }
            if ((new.vars.i[[h1]][v1] == new.k2.vars.i[h1]) && 
                any(new.vars.i[[h1]][v.ind] == new.k1.vars.i[h1])){
              v.ind = v.ind[-which(new.vars.i[[h1]][v.ind] == new.k1.vars.i[h1])]
            }
          }
        }else{
          h.dependent.l = h.dependent[[h1]]$levels[[v1]]
          h.dependent.v = h.dependent[[h1]]$voxels[[v1]]
          if (length(h.dependent.l) > 0){
            for (v2 in 1:length(h.dependent.l)){
              if (length(h.dependent.l[[v2]]) > 0){
                for (h3 in 1:length(h.dependent.l[[v2]])){
                  if (h.levels[h2] == h.dependent.l[[v2]][h3]){
                    v.ind = v.ind[-which(v.ind == h.dependent.v[[v2]][h3])]
                  }
                }
              }
            }
          }
        }
        h.independent[[h1]][[v1]][[h2]] = v.ind
      }
    }
  }
  
  if (!all(is.na(h.leaves))){
    ### Find minimal independent set from leaf nodes ###
    if (verbose){
      cat('\n.........finding minimal independent set of vars. as leaf nodes of hierarchy...')
    }
    for (h1 in 1:H){
      for (v1 in 1:length(h.independent[[h1]])){
        for (h2 in 1:H){
          v.ind = h.independent[[h1]][[v1]][[h2]]
          v.leaf = v.ind[which(v.ind %in% h.leaves[[h2]])]
          h.independent[[h1]][[v1]][[h2]] = v.leaf
        }
      }
    }
  }
  
  ##################################################################
  return(h.independent)
} ##################################################################



##################################################################
hier_separate_sources <- function(hier, h.independent, verbose = F){
  ##################################################################
  # Separates source signals while preserving hierarchical dependencies,
  #   using pseudo-inverse to project onto orthogonal subspace.
  # Used to increase specificity of hierarchical spatial maps.
  #
  # Input:
  #   hier  : list w/ hierarchy info by level, from get_hier_dat().
  #   h.independent : list of independent vars. w/n hier, from hier_find_independent_levels().
  # Output:
  #   hier  : formatted list as above, with added elements for each level:
  #     $premerge.vars.sep  : compressed group t.s. for pre-merge vars. (if included in input),
  #                             after removing effects of all other independent vars. in hier.
  #     $new.vars.sep       : compressed group t.s. for vars. in level,
  #                             after removing effects of all other independent vars. in hier.
  #
  
  if (verbose){
    cat('\n......projecting hPCA t.s. onto orthogonal subspace of independent vars...')
  }
  
  stopifnot(length(hier) > 0)
  stopifnot(all(c('new.vars', 'h.level') %in% names(hier[[1]])))
  h.levels = unlist(lapply(hier, function(yy) return(yy$h.level)))
  H = length(h.levels)
  K2 = nrow(hier[[1]]$new.vars)
  
  for (h1 in 1:H){
    ### Get hierarchical independent vars. for pre-merge vars. ###
    if (all(c('premerge.vars', 'premerge.vars.i') %in% names(hier[[h1]]))){
      V1 = ncol(hier[[h1]]$premerge.vars)
      premerge.vars.sep = matrix(NA, K2, V1)
      for (v1 in 1:V1){
        if (hier[[h1]]$premerge.vars.i[v1] %in% unlist(hier[[h1]][c('k1', 'k2')])){
          V2 = length(unlist(h.independent[[h1]][[1]]))  # if var. is merged at level, use dependencies for sum var. from level...
        }else{
          premerge.vars.sep[,v1] = hier[[h1]]$premerge.vars[,v1]  # skip orthog. projection if premerge.var is leaf node after filtering
          next
        } 
        if (V2 == 0){
          premerge.vars.sep[,v1] = hier[[h1]]$premerge.vars[,v1]  # skip orthog. projection if premerge.var is leaf node after filtering
        }else{
          R.indep = matrix(NA, K2, V2)
          v.start = 0
          v.end = 0
          for (h2 in 1:H){
            v.inds = h.independent[[h1]][[1]][[h2]]  # use main var. dependencies for level
            if (length(v.inds) > 0){
              v.start = v.end + 1
              v.end = v.end + length(v.inds)
              R.indep[, v.start:v.end] = hier[[h2]]$new.vars[, v.inds]
            }
          }
          
          ### Project pre-merge vars. onto independent subspace ###
          R_pinv = MASS::ginv(R.indep)
          P = diag(K2) - R.indep %*% R_pinv  # projection operator onto orthogonal subspace
          premerge.vars.sep[,v1] = P %*% hier[[h1]]$premerge.vars[,v1]
        }
      }
      hier[[h1]]$premerge.vars.sep = premerge.vars.sep
    }
    
    ### Get hierarchical independent vars. for main vars. ###
    V1 = ncol(hier[[h1]]$new.vars)
    new.vars.sep = matrix(NA, K2, V1)
    for (v1 in 1:V1){
      V2 = length(unlist(h.independent[[h1]][[v1]]))
      if (V2 == 0){
        new.vars.sep[,v1] = hier[[h1]]$new.vars[,v1]  # skip orthog. projection for leaf nodes after filtering
        
      }else{
        R.indep = matrix(NA, K2, V2)
        v.start = 0
        v.end = 0
        for (h2 in 1:H){
          v.inds = h.independent[[h1]][[v1]][[h2]]
          if (length(v.inds) > 0){
            v.start = v.end + 1
            v.end = v.end + length(v.inds)
            R.indep[, v.start:v.end] = hier[[h2]]$new.vars[, v.inds]
          }
        }
        
        ### Project onto independent subspace ###
        R_pinv = MASS::ginv(R.indep)
        P = diag(K2) - R.indep %*% R_pinv  # projection operator onto orthogonal subspace
        new.vars.sep[,v1] = P %*% hier[[h1]]$new.vars[,v1]
      }
    }
    hier[[h1]]$new.vars.sep = new.vars.sep
  }
  
  ##################################################################
  return(hier)
} ##################################################################


##################################################################
subj_HierBackReconstruct <- function(hier,
                                     s, PCA_s, PCA_g, 
                                     align_sources = T,
                                     br = list(), verbose = F){
  ##################################################################
  # Subj.-specific source back-reconstruction
  #   for all sig. levels of group hierarchical PCA
  #     called as part of group_HierBackReconstruction().
  #
  # Input:
  #   hier     : list w/ hierarchy info by level, from get_hier_dat()
  #   s     : subject index
  #   PCA_s : subj.-specific PCA data, from get_subjPCA_dat()
  #   PCA_g : group PCA data, from get_groupPCA_dat()
  #   align_sources : align sources to skew
  #   br    : previous back-reconstructed components for subj.
  # Output:     list w/ elements
  #   S_i              : matrix of source spatial maps [voxels x comps.]
  #   R_i              : matrix of source time series  [time x comps.]
  #   Hierarchy_level  : hierarchical ICA/PCA level for return comps.
  #   hier_prev1_inds  : indices of _pre-merge_ var1's for each level of hierarchy
  #   hier_prev2_inds  : indices of _pre-merge_ var2's for each level of hierarchy
  #   hier_k1_inds     : indices of post-merge sum comps. for each level of hierarchy
  #   hier_k2_inds     : indices of post-merge diff. comp. for each level of hierarchy
  #   hier_inds        : indices for all sig. comps. from hICA / hPCA
  # Requires:
  #   zeroMean_Yi()  : datamunging fn. from PCA_fns.R
  #   flatten_img()  :  "      "      "     "
  #   align_toSkew() : flips signs of comps. to ensure pos. skew
  #
  
  ### Inputs / defaults ###
  stopifnot(length(hier) > 0)
  stopifnot('new.vars' %in% names(hier[[1]]))
  stopifnot('new.vars.i' %in% names(hier[[1]]))
  s = as.integer(s)
  
  br = as.list(br)
  if (length(br) == 0){  br = list('S_i' = NULL, 'R_i'= NULL)  }
  if (!('R_i' %in% names(br))){  br$R_i = NULL  }
  if (!('S_i' %in% names(br))){  br$S_i = NULL  }
  
  stopifnot(all(c('U_i','L_i', 'Y', 'whitenedPCs') %in% names(PCA_s)))
  U_i = PCA_s$U_i
  K1 = ncol(U_i)
  if (PCA_s$whitenedPCs){
    L_i = PCA_s$L_i
  }else{
    L_i = rep(1, K1)
  }
  if (('Y_vx' %in% names(PCA_s)) && is.numeric(PCA_s$Y_vx)){  
    Y_i = PCA_s$Y_vx  # voxel-level data required for voxel-level wholebrain spatial maps, may be different from (~ROI) t.s. in subj.-level PCA
  }else{  Y_i = PCA_s$Y  }
  Y_i = zeroMean_Yi(Y_i)
  
  stopifnot(all(c('G', 'L', 'whitenedPCs') %in% names(PCA_g)))
  if (!is.na(s)){ # find subj.specific column space w/n group PCA space
    i.start = K1 * (s-1) + 1
    i.end =   K1 * s
    G_i = PCA_g$G[i.start:i.end,] # subj.-specific slice of group reduction matrix
  }else{
    G_i = PCA_g$G
  }
  K2 = ncol(PCA_g$G)
  if (PCA_g$whitenedPCs){
    L_gr = PCA_g$L     # re-scale PCs to reverse whitening
  }else{
    L_gr = rep(1, K2)
  }
  
  if (!('S_i' %in% names(br)) || is.null(br$S_i)){
    hier_i0 = 1  # start index of hierarchical ICA/PCA sources in output
  }else{
    hier_i0 = ncol(br$S_i) + 1  # start index of hierarchical ICA/PCA sources in output
    if (!('Hierarchy_level' %in% names(br))){  br$Hierarchy_level = rep(NA, ncol(br$S_i))  }
    if (!('hier_algorithm' %in% names(br))){  br$hier_algorithm = rep(NA, ncol(br$S_i))  }
    if (!('orig_var_inds' %in% names(br))){  br$orig_var_inds = rep(NA, K.prev)  }
  }
  
  ### Sanity checks ###
  stopifnot(is.numeric(L_i))
  stopifnot(is.numeric(L_gr))
  stopifnot(is.matrix(Y_i) && is.numeric(Y_i))
  stopifnot(ncol(G_i) == nrow(hier[[1]]$new.vars))
  if ('premerge.vars' %in% names(hier[[1]])){
    stopifnot(ncol(G_i) == nrow(hier[[1]]$premerge.vars))
  }
  stopifnot(ncol(U_i) == nrow(G_i))
  stopifnot(ncol(Y_i) == nrow(U_i))
  
  
  ### Main Fn. ###
  h_levels = c()
  vrbs1 = verbose
  for (h in 1:length(hier)){ # index hierarchical sources starting from bottom of hierarchy
    l = hier[[h]]$h.level
    
    #----------------------------------------------------------------------
    ### Back-reconstruct new vars. in pre-merged state from prev. level ###
    #----------------------------------------------------------------------
    if ('premerge.vars' %in% names(hier[[h]])){
      k.inds = hier[[h]]$premerge.vars.i
      h_levels = c(h_levels, rep(l, length(k.inds)))
      
      if (!('S_i' %in% names(br)) || is.null(br$S_i)){
        br$hier_prev1_inds = c(br$hier_prev1_inds, which(k.inds == hier[[h]]$k1))  # pre-merge var1 ind. always post-merge k1 ind.
        br$hier_prev2_inds = c(br$hier_prev2_inds, which(k.inds == hier[[h]]$k2))  # pre-merge var2 ind. always post-merge k2 ind.
      }else{
        if (hier[[h]]$k1 %in% k.inds){
          br$hier_prev1_inds = c(br$hier_prev1_inds, ncol(br$S_i) + which(k.inds == hier[[h]]$k1))
        }
        if (hier[[h]]$k2 %in% k.inds){
          br$hier_prev2_inds = c(br$hier_prev2_inds, ncol(br$S_i) + which(k.inds == hier[[h]]$k2))
        }
      }
      
      ########################################################
      ### Un-compress time dimension & create spatial maps ###
      ########################################################
      if (hier[[1]]$separate_time_series || hier[[1]]$separate_spatial_maps){
        R_i_orth = U_i %*% diag(sqrt(L_i)) %*% G_i %*% diag(sqrt(L_gr)) %*% hier[[h]]$premerge.vars.sep
      }
      if (!hier[[1]]$separate_time_series || !hier[[1]]$separate_spatial_maps){
        R_i_orig = U_i %*% diag(sqrt(L_i)) %*% G_i %*% diag(sqrt(L_gr)) %*% hier[[h]]$premerge.vars
      }
      
      if (hier[[1]]$separate_time_series){
        R_i = R_i_orth
      }else{
        R_i = R_i_orig
      }
      
      if (hier[[1]]$separate_spatial_maps){ 
        S_i = cor(t(Y_i), R_i_orth)
      }else{
        S_i = cor(t(Y_i), R_i_orig)
      }
      
      
      ### Align spatial maps & calibrate t.s. ###
      if (align_sources){
        S_i = align_toSkew(S_i, verbose=F)
        
        sgn = rep(1,ncol(S_i))
        for (k in 1:ncol(S_i)){
          s.max = max(S_i[,k])  # calibrate comp. based on highest magnetude voxel
          if (max(abs(S_i[,k])) > s.max){  s.max = -max(abs(S_i[,k]))  }
          v = which(S_i[,k] == s.max)
          if (cor(c(Y_i[v,]), R_i[,k]) < -.Machine$double.eps){
            sgn[k] = -1
          }
        }
        if (any(sgn < 0)){
          R_i = R_i %*% diag(sgn)
        }
      }
      
      ### Format output & meta-info ###
      br$R_i = cbind(br$R_i, R_i, deparse.level=0)
      br$S_i = cbind(br$S_i, S_i, deparse.level=0)
      if ('hier_algorithm' %in% names(hier[[h]])){
        br$hier_algorithm = c(br$hier_algorithm, rep(hier[[h]]$hier_algorithm, ncol(S_i)))
      }else if ('hier_algorithm' %in% names(hier[[1]])){
        br$hier_algorithm = c(br$hier_algorithm, rep(hier[[1]]$hier_algorithm, ncol(S_i)))
      }
      br$orig_var_inds = c(br$orig_var_inds, k.inds)
    }
    
    
    #--------------------------------------------
    ### Back-reconstruct new vars. for level ###
    #--------------------------------------------
    k.inds = hier[[h]]$new.vars.i
    if (!('S_i' %in% names(br)) || is.null(br$S_i)){
      br$hier_k1_inds = c(br$hier_k1_inds, which(k.inds == hier[[h]]$k1))
      br$hier_k2_inds = c(br$hier_k2_inds, which(k.inds == hier[[h]]$k2))
    }else{
      if (hier[[h]]$k1 %in% k.inds){
        br$hier_k1_inds = c(br$hier_k1_inds, ncol(br$S_i) + which(k.inds == hier[[h]]$k1))
      }
      if (hier[[h]]$k2 %in% k.inds){
        br$hier_k2_inds = c(br$hier_k2_inds, ncol(br$S_i) + which(k.inds == hier[[h]]$k2))
      }
    }
    h_levels = c(h_levels, rep(l, length(k.inds)))
    
    
    ########################################################
    ### Un-compress time dimension & create spatial maps ###
    ########################################################
    if (hier[[1]]$separate_time_series || hier[[1]]$separate_spatial_maps){
      R_i_orth = U_i %*% diag(sqrt(L_i)) %*% G_i %*% diag(sqrt(L_gr)) %*% hier[[h]]$new.vars.sep
    }
    if (!hier[[1]]$separate_time_series || !hier[[1]]$separate_spatial_maps){
      R_i_orig = U_i %*% diag(sqrt(L_i)) %*% G_i %*% diag(sqrt(L_gr)) %*% hier[[h]]$new.vars
    }
    
    if (hier[[1]]$separate_time_series){
      R_i = R_i_orth
    }else{
      R_i = R_i_orig
    }
    
    if (hier[[1]]$separate_spatial_maps){
      S_i = cor(t(Y_i), R_i_orth)
    }else{
      S_i = cor(t(Y_i), R_i_orig)
    }
    
    ### Align spatial maps & calibrate t.s. ###
    if (align_sources){
      S_i = align_toSkew(S_i, verbose=vrbs1)
      
      sgn = rep(1,ncol(S_i))
      for (k in 1:ncol(S_i)){
        s.max = max(S_i[,k])  # calibrate comp. based on highest magnetude voxel
        if (max(abs(S_i[,k])) > s.max){  s.max = -max(abs(S_i[,k]))  }
        v = which(S_i[,k] == s.max)
        
        if (cor(c(Y_i[v,]), R_i[,k]) < -.Machine$double.eps){
          sgn[k] = -1
        }
      }
      
      if (any(sgn < 0)){
        if (length(sgn) == 1){  sgn = as.matrix(sgn)  }
        R_i = R_i %*% diag(sgn)
      }
    }
    vrbs1 = F
    
    
    ### Format output & meta-info ###
    br$R_i = cbind(br$R_i, R_i, deparse.level=0)
    br$S_i = cbind(br$S_i, S_i, deparse.level=0)
    if ('hier_algorithm' %in% names(hier[[h]])){
      br$hier_algorithm = c(br$hier_algorithm, rep(hier[[h]]$hier_algorithm, ncol(S_i)))
    }else if ('hier_algorithm' %in% names(hier[[1]])){
      br$hier_algorithm = c(br$hier_algorithm, rep(hier[[1]]$hier_algorithm, ncol(S_i)))
    }
    br$orig_var_inds   = c(br$orig_var_inds, k.inds)
  }
  br$Hierarchy_level = c(br$Hierarchy_level, h_levels)
  br$hier_inds = c(br$hier_inds, c(hier_i0:ncol(br$R_i)))
  
  ##################################################################
  return(br)
} ##################################################################


##################################################################
align_toSkew <- function(data, k.indices = NA,
                         return.list = F,
                         verbose = F){
  ##################################################################
  # Tweaks output from hierarchical Principal Components Analysis (hPCA)
  #   aligning all hPCA comp. signs to a positive 3rd moment (skew),
  #   a rough work around for PCA sign indeterminacy,
  #   that ensures +/- PC loadings have some sort of meaning in the data
  #
  # Input:
  #   data : matrix formatted as [voxels x sources]
  #   k.indices : list of integer vectors, 
  #                 limits alignment & estimation to subset of sources
  #   return.list : return original skew signs & aligned data as list
  # Output:  either aligned input matrix data, 
  #           or (if return.list) list w/ elements:
  #     data : input matrix, w/ all cols. +/- flipped to skew pos.
  #     skew_sign : sign of original skew
  #
  
  if (verbose){ cat('\n......aligning +/- sources to skew positively...') }
  
  ##################
  skew <- function(dat){
    dat = as.vector(dat)
    n = length(dat)
    m = mean(dat)
    sigma = sd(dat)
    return( (1/n) * sum((dat - m)^3) / (sigma^3) )
  } ################
  
  stopifnot(is.numeric(data))
  stopifnot(length(dim(data))==2)
  stopifnot(prod(dim(data)) > 0)
  if (all(is.na(k.indices))){
    k.inds = c(1:ncol(data))
  }else{
    k.inds = sort(as.integer(k.indices))  # preserve ordering of sources in mixing matrix & spatial maps
  }
  return.list = as.logical(return.list)
  
  sk = rep(1, ncol(data))
  sk[k.inds] = apply(data, 2, skew)
  if (any(sk < -.Machine$double.eps)){
    sk = sign(sk)
    if (length(sk) == 1){  sk = as.matrix(sk)  }
    data = data %*% diag(sk)  # right-multiplication flips signs of cols.
  }
  
  #########################################
  if (return.list){
    return(list('data'       = data,
                'skew_signs' = sk,
                'k.indices'  = k.inds))
  }else{  return(data)  }
} #########################################


##################################################################
group_ScaleSources <- function(group_info_path,
                               scale = 'z-scores',
                               insert.NAs = T,
                               verbose = F){
  ##################################################################
  # Calibrates & scales subj.-specific sources signals
  #   following group hierchical PCA & back-reconstruction
  #
  # Input:
  #   group_info_path : Path to dir. w/ saved outputs of 
  #                       subj_PCA(), group_PCA(), group_HierPCA().
  #   scale           : Scaling method for subj. source signals, options:
  #                       z-scores : subract mean & scale to standard deviation
  #                       L2       : L2 norm, scale to square root of sum of squared weights
  #                       Sup      : supremum norm, scale to max. abs. value
  #   insert.NAs      : Insert NA values for discarded/scrubbed time points,
  #                       requires vector of 'good.timepoints' & scalar 'total.timepoints' saved in above files.
  # Output: 
  #   -saved matrix as 2D nifti vol. w/ time series for each subj.
  #
  
  if (verbose){  cat('\nScaling subj.-level ICA/hICA/hPCA sources...')  }
  
  ### Inputs & defaults ###
  subj_pca_suffix  = '_PCA1.RData'    # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData' # filesuffix for group-level temporal PCA from group_PCA()
  hier_suffix = 'Hierarchy.RData'   # filesuffix for hierarchical PCA from group_HierPCA()
  br_suffix   = '_HierComps.RData'  # filesuffix for back-reconstructed subj. sources  
  spatial_maps_suffix = '_HierSpatialMaps.nii' # filesuffix for 4D nii vol. w/ back-reconstructed spatial maps
  time_series_suffix  = '_HierTimeSeries.nii'  # filesuffix for 2D nii vol. w/ back-reconstructed t.s.
  
  group_info_path = as.character(group_info_path)
  scale = as.character(scale)
  if (length(scale) > 1){  scale = scale[1]  }
  
  ### Prelims ###
  if (is.character(scale)){
    stopifnot(scale %in% c('z-scores', 
                           'L2', 
                           'Sup', 'sup', 'supremum'))
  }else if (all(is.na(scale))){
    scale = ''  # skip scaling ICs
  }else{
    stopifnot(is.character(scale))
  }
  br_wc = paste0('*', br_suffix)
  
  ### Get file names & paths ###
  br_files = list.files(group_info_path, br_wc, full.names=T)
  S = length(br_files)
  
  ### Scale back-reconstructed components ###
  for (s in 1:S){
    if (verbose){  cat(paste0('\n...scaling subj.:  ', basename(br_files[s])))  }
    if (s == 1){  vrbs = verbose  }else{  vrbs = F  }  # shorten output after 1st subj.
    load(br_files[s])
    
    if (all(is.na(br$S_i))){
      if (verbose){  cat('\n......Warning:  no back-reconstructed to scale for this subject')}
      next
    }
    
    if (scale == 'z-scores'){
      if (vrbs){  cat('\n......subtracting mean & scaling to s.d. for all sources  ~  z-scores')  }
      
      br$S_i_means = colMeans(br$S_i, na.rm=T)
      br$S_i = sweep(br$S_i, 2, colMeans(br$S_i, na.rm=T))
      br$S_i_sd = apply(br$S_i, 2, sd, na.rm=T)
      br$S_i = apply(br$S_i, 2, function(yy) return(yy / sd(yy, na.rm=T)))
      
      br$R_i_means = colMeans(br$R_i, na.rm=T)
      br$R_i = sweep(br$R_i, 2, colMeans(br$R_i, na.rm=T))
      br$R_i_sd = apply(br$R_i, 2, sd, na.rm=T)
      br$R_i = apply(br$R_i, 2, function(yy) return(yy / sd(yy, na.rm=T)))
      
    }else if (scale == 'L2'){
      if (vrbs){  cat('\n......scaling all sources to L2 norm')  }
      
      L2.norm <- function(Lambda){  sqrt(sum(Lambda^2, na.rm=T))  }
      
      br$S_i_L2.norm = apply(br$S_i, 2, L2.norm)
      br$S_i = apply(br$S_i, 2, function(yy) return(yy / L2.norm(yy)))
      br$R_i_L2.norm = apply(br$R_i, 2, L2.norm)
      br$R_i = apply(br$R_i, 2, function(yy) return(yy / L2.norm(yy)))
      
    }else if (scale %in% c('Sup', 'sup', 'supremum')){
      if (vrbs){  cat('\n......scaling all sources to supremum norm')  }
      
      Sup.norm <- function(Lambda){  max(abs(Lambda), na.rm=T)    } # Suprenum norm for finite-dim. vectors
      
      br$S_i_Sup.norm = apply(br$S_i, 2, Sup.norm)
      br$S_i = apply(br$S_i, 2, function(yy) return(yy / Sup.norm(yy)))
      br$R_i_Sup.norm = apply(br$R_i, 2, Sup.norm)
      br$R_i = apply(br$R_i, 2, function(yy) return(yy / Sup.norm(yy)))
    }
    br$scale.method = scale
    
    if (insert.NAs && all(c('good.timepoints', 'total.timepoints') %in% names(br))){
      if (vrbs){  cat('\n......buffering time series w/ NA values in discarded/scrubbed times:')  }
      stopifnot(max(br$good.timepoints) <= br$total.timepoints)
      stopifnot(nrow(br$R_i) <= br$total.timepoints)
      if (br$total.timepoints != nrow(br$R_i)){
        stopifnot(length(br$good.timepoints) == nrow(br$R_i))
      }
      
      if (nrow(br$R_i) < br$total.timepoints){
        n.NA = br$total.timepoints - length(br$good.timepoints)
        n = br$total.timepoints
        K = ncol(br$R_i)
        if (vrbs){  cat(paste0('\n        inserting NAs in ',n.NA,'/',n,' total original time points'))  }
        
        R_i = matrix(NA, n, K)
        R_i[br$good.timepoints,] = br$R_i
        br$R_i = R_i
      }
    }
    
    save.fname = br_files[s]
    if (vrbs){  cat(paste0('\n......& saving as:  ', save.fname))  }
    save(list='br', file=save.fname)
  }
  
  if (verbose){  cat(paste0('\n.........saved in dir.:  ', unique(dirname(br_files)), '\n\n'))  }
} ##################################################################


##################################################################
group_summaryStats <- function(group_info_path, prefix = NA, 
                               calc.var = F,
                               create.csv = T,
                               create.nii = T,
                               verbose = F){
  ##################################################################
  # Calculates summary statistics
  #   & saves as Nifti volume
  #   for group-level hierarchical PCA or spatial ICA output
  #
  # Input:
  #   group_info_path : Path to dir. w/ saved outputs of 
  #                       subj_PCA(), group_PCA(), group_HierPCA()
  #   prefix   : identifying prefix to attach to saved output
  #   calc.var : calculate standard deviation for each component
  #   create.csv : creates .csv table w/ component info
  #   create.nii : save output as nifti volume
  # Output: 
  #   -saved 4D nifti vol. w/ mean of subj. spatial maps for each IC
  #   -if indicated, saved 4D nifti vol. w/ standard deviation
  #     of subj. spatial maps for each IC
  #   -saved matrix as 2D nifti vol. w/ mean of subj. time series
  #   -csv table w/ columns:
  #     Component               : filename & 4d-nifti index for component
  #     Hierarchy Level         : Level of component in hierarchy, NA for non-hPCA/hICA comps.
  #     Hierarchy Variable Type : explanation:
  #                               'initial sig.' : comp. is from lowest (significant) level of hierarchy
  #                               'sum' : comp. represents merged shared variance in hPCA (~PC1)
  #                               'diff': comp. represents difference of merged vars. in hPCA (~PC2)
  #     Variable Index          : Index for comp. w/n ordered set of original variables,
  #                                 used for locating component within data
  #     Variable Source Level   : Previous level appearance of variable within hierarchy,
  #                                 used to compare merging of components across levels
  #     Algorithm               : Hierarchical PCA or ICA algorithm
  #

  suppressMessages(library(RNifti))
  
  if (verbose){ cat('\nCalculating group-level summary statistics (mean, etc.) of sources...')}
  
  ### Inputs & defaults ###
  subj_pca_suffix  = '_PCA1.RData'    # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData' # filesuffix for group-level temporal PCA from group_PCA()
  hier_suffix = 'Hierarchy.RData'   # filesuffix for hierarchical PCA from group_HierPCA()
  br_suffix   = '_HierComps.RData'  # filesuffix for back-reconstructed subj. sources  
  spatial_maps_suffix = '_HierSpatialMaps.nii' # filesuffix for 4D nii vol. w/ back-reconstructed spatial maps
  time_series_suffix  = '_HierTimeSeries.nii'  # filesuffix for 2D nii vol. w/ back-reconstructed t.s.
  mean_spatial_map_suffix = 'Hier_SpatialMapsMeans.nii' # filesuffix for mean of spatial maps
  var_spatial_map_suffix  = 'Hier_SpatialMapsVars.nii'  # filesuffix for variance of spatial maps
  mean_time_series_suffix = 'Hier_TimeSeriesMeans.nii'  # filesuffix for mean of time series
  var_time_series_suffix  = 'Hier_TimeSeriesVars.nii'  # filesuffix for variance of time series
  
  group_info_path = as.character(group_info_path)
  prefix = as.character(prefix)
  create.csv = as.logical(create.csv)
  create.nii = as.logical(create.nii)
  
  save_path = group_info_path         # save to same dir.
  
  ### Prelims ###
  br_wc           = paste0('*', br_suffix)
  spatial_maps_wc = paste0('*', spatial_maps_suffix)
  time_series_wc  = paste0('*', time_series_suffix)
  if (!all(is.na(prefix))){
    br_wc           = paste0(prefix,'.', br_wc)
    spatial_maps_wc = paste0(prefix,'.', spatial_maps_wc)
    time_series_wc  = paste0(prefix,'.', time_series_wc)
  }
  
  br_files = list.files(group_info_path, br_wc, full.names=T)
  S = length(br_files)
  if (S==0 && verbose){
    cat(paste0('\n...Warning:  Could not find any back-reconstructed files in:  ',group_info_path,'\n'))
  }
  stopifnot(S > 0)
  
  ###################
  get_Time_minMax <- function(ts_files, nii=T){
    # Finds limiting number of timepoints shared among subjs.
    if (!nii){  Space.br = new.env()  }
    S = length(ts_files)
    T.max = Inf
    for (s in 1:S){
      if (nii){
        T.max = min(T.max, niftiHeader(ts_files[s])$dim[2])
      }else{
        load(ts_files[s], envir=Space.br)
        T.max = min(T.max, nrow(get('br',envir=Space.br)$R_i))
      }
    }
    return(T.max)
  } #################
  
  
  ### Main fn. ###
  T.limit = get_Time_minMax(br_files, nii=F)
  
  for (s in 1:S){
    load(br_files[s])
    na.flags = is.na(br$R_i[1:T.limit,1])
    non.NA.inds = as.numeric(!na.flags)
    if (any(na.flags)){  br$R_i[is.na(br$R_i)] = 0  }
    
    if (s==1){
      timespace_comps = rbind(cbind(br$R_i[1:T.limit,], non.NA.inds),
                              cbind(br$S_i, rep(NA,nrow(br$S_i))))
    }else{
      timespace_comps = timespace_comps + rbind(cbind(br$R_i[1:T.limit,], non.NA.inds),
                                                cbind(br$S_i, rep(NA,nrow(br$S_i))))
    }
  }
  
  load(br_files[1])
  br_mean = br
  k.NA = ncol(timespace_comps)  # last col. is count of number of subjs. w/o NA values for each time point
  br_mean$R_i = timespace_comps[1:T.limit, -k.NA] / timespace_comps[1:T.limit, k.NA]  # av. over non-NA subjs. by time point
  br_mean$R_i[is.infinite(br_mean$R_i)] = 0  # set time points w/o data from any subj. to 0
  br_mean$R_i[is.na(br_mean$R_i)] = 0
  br_mean$S_i = timespace_comps[-(1:T.limit), -k.NA] / S

  if (calc.var){
    for (s in 1:S){
      load(br_files[s])
      na.flags = is.na(br$R_i[1:T.limit,1])
      non.NA.inds = as.numeric(!na.flags)
      if (any(na.flags)){  br$R_i[is.na(br$R_i)] = 0  }
      
      if (s==1){
        timespace_comps = rbind(cbind((br$R_i[1:T.limit,] - br_mean$R_i)^2, non.NA.inds),
                                cbind((br$S_i - br_mean$S_i)^2, rep(NA,nrow(br$S_i))))
      }else{
        timespace_comps = timespace_comps + rbind(cbind((br$R_i[1:T.limit,] - br_mean$R_i)^2, non.NA.inds),
                                                  cbind((br$S_i - br_mean$S_i)^2, rep(NA,nrow(br$S_i))))
      }
    }
    timespace_comps[1:T.limit, k.NA] = pmax(0, timespace_comps[1:T.limit, k.NA] - 1)
    
    br_var = br
    k.NA = ncol(timespace_comps)  # last col. is count of number of subjs. w/o NA values for each time point
    br_var$R_i = timespace_comps[1:T.limit, -k.NA] / timespace_comps[1:T.limit, k.NA]  # av. over non-NA subjs. by time point
    br_var$R_i[is.infinite(br_var$R_i)] = 0  # set time points w/o data from any subj. to 0
    br_var$R_i[is.na(br_var$R_i)] = 0
    br_var$S_i = timespace_comps[-(1:T.limit), -k.NA] / (S-1)
  }
  
  
  ### Saving ###
  save.fname1 = mean_spatial_map_suffix
  save.fname2 = mean_time_series_suffix
  save.fname3 = var_spatial_map_suffix
  save.fname4 = var_time_series_suffix
  if (!all(is.na(prefix))){
    save.fname1 = paste0(prefix, save.fname1)
    save.fname2 = paste0(prefix, save.fname2)
    save.fname3 = paste0(prefix, save.fname3)
    save.fname4 = paste0(prefix, save.fname4)
  }
  save.fname1 = file.path(save_path, save.fname1)
  save.fname2 = file.path(save_path, save.fname2)
  save.fname3 = file.path(save_path, save.fname3)
  save.fname4 = file.path(save_path, save.fname4)
  
  if (create.nii){
    if (verbose){  cat(paste0('\n...saving vol. of comp. mean spatial maps as:  ', basename(save.fname1)))  }
    if (verbose){  cat(paste0('\n...saving ~2D vol. of comp. mean time series as:  ', basename(save.fname2)))  }
    group_createNiftiVols(group_info_path, 
                          br_mean, 
                          prefix=prefix,
                          spatial_maps_save=save.fname1, 
                          time_series_save=save.fname2,
                          verbose=verbose)
    
    if (calc.var){
      if (verbose){  cat(paste0('\n...saving vol. of comp. variance spatial maps as:  ', basename(save.fname3)))  }
      if (verbose){  cat(paste0('\n...saving ~2D vol. of comp. variance time series as:  ', basename(save.fname4)))  }
      group_createNiftiVols(group_info_path, 
                            br_var, 
                            prefix=prefix,
                            spatial_maps_save=save.fname3, 
                            time_series_save=save.fname4,
                            verbose=verbose)
    }
  }
  
  
  if (create.csv){
    save.fname5 = save.fname1
    save.fname5 = sub('\\.nii$', '', save.fname5)
    save.fname5 = sub('\\.RData$', '', save.fname5)
    save.fname5 = paste0(save.fname5, '.csv')
    group_create_b.r.sourceTable(br_mean, save.fname5, verbose)
    
    if (calc.var){
      save.fname6 = save.fname3
      save.fname6 = sub('\\.nii$', '', save.fname6)
      save.fname6 = sub('\\.RData$', '', save.fname6)
      save.fname6 = paste0(save.fname6, '.csv')
      group_create_b.r.sourceTable(br_var, save.fname6, verbose)
    }
  }
  
  if (verbose){  cat(paste0('\n......saved in dir.:  ', save_path, '\n\n'))  }
} ##################################################################


##################################################################
group_createNiftiVols <- function(group_info_path, 
                                  br_files = NA, prefix = NA,
                                  spatial_maps_save = NA, 
                                  time_series_save = NA,
                                  verbose = F){
  ##################################################################
  # Creates nifti volumes for subj.-specific sources signals
  #   following group hierarchical PCA & back-reconstruction
  #
  # Input:
  #   group_info_path : Path to dir. w/ saved outputs of 
  #                       subj_PCA(), group_PCA(), group_HierPCA()
  #   br_files : data or files to project into brain mask (~MNI space),
  #               defaults to back-reconstructed data output from group_HierBackReconstruction,
  #               w/ or w/o scaling by group_ScaleSources().
  #              Data must be formatted as list w/ elements 'S_i' & 'R_i'
  #   prefix : identifying prefix to attach to saved output
  #   spatial_maps_save : file name & path for saved subjects' spatial maps nii vols.
  #   time_series_save  : file name & path for saved subjects' time series 2D-nii vol.
  # Output: 
  #   -saved 4D nifti vol. spatial maps for each subj.
  #   -saved matrix as 2D nifti vol. w/ time series for each subj.
  # Requies:
  #   unflatten_img() to load spatial maps into 4D nifti vol., from PCA_fns.R
  #
  
  suppressMessages(library(RNifti))
  
  if (verbose){ cat('\nProjecting subj.-level ICA/hICA/hPCA sources onto brain/ROI mask...')}
  
  ### Inputs & defaults ###
  subj_pca_suffix  = '_PCA1.RData'    # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData' # filesuffix for group-level temporal PCA from group_PCA()
  hier_suffix = 'Hierarchy.RData'   # filesuffix for hierarchical PCA from group_HierPCA()
  br_suffix   = '_HierComps.RData'  # filesuffix for back-reconstructed subj. sources  
  spatial_maps_suffix = '_HierSpatialMaps.nii' # filesuffix for 4D nii vol. w/ back-reconstructed spatial maps
  time_series_suffix  = '_HierTimeSeries.nii'  # filesuffix for 2D nii vol. w/ back-reconstructed t.s.
  mean_spatial_map_suffix = 'Hier_SpatialMapsMeans.nii' # filesuffix for mean of spatial maps
  var_spatial_map_suffix  = 'Hier_SpatialMapsVars.nii'  # filesuffix for variance of spatial maps
  mean_time_series_suffix = 'Hier_TimeSeriesMeans.nii'  # filesuffix for mean of time series
  var_time_series_suffix  = 'Hier_TimeSeriesVars.nii'  # filesuffix for variance of time series
  
  group_info_path = as.character(group_info_path)
  prefix = as.character(prefix)
  spatial_maps_save = as.character(spatial_maps_save)
  time_series_save = as.character(time_series_save)
  
  ### Prelims ###
  if (all(is.na(prefix))){
    group_pca_wc = group_pca_suffix
    hier_wc      = hier_suffix
    br_wc        = paste0('*', br_suffix)
  }else if (!all(is.na(prefix))){
    group_pca_wc = paste0(prefix, '.*_', group_pca_suffix)
    hier_wc      = paste0(prefix, '.*_', hier_suffix)
    br_wc        = paste0(prefix, '.*_', br_suffix)
  }
  
  ### Get file names & paths ###
  if (all(is.na(br_files))){
    br_files = list.files(group_info_path, br_wc, full.names=T)
    S = length(br_files)
  }else if (is.character(br_files)){
    stopifnot(all(file.exists(br_files)))
    S = length(br_files)
  }else{
    stopifnot(all(c('S_i','R_i') %in% names(br_files)))
    stopifnot(!all(is.na(spatial_maps_save)))
    S = 1
  }
  
  
  ### Main fn. ###
  for (s in 1:S){
    if (is.character(br_files)){
      if (verbose){  cat(paste0('\n...projecting subj.:  ', br_files[s]))}
      load(br_files[s])
      stopifnot(exists('br'))
    }else{
      br = br_files
    }
    if (all(is.na(br$S_i))){
      if (verbose){  cat(paste0('\n......Warning:  no back-reconstructed sources for ', basename(br_files[s])))  }
      next
    }
    
    spatialmaps.nii = unflatten_img(br$S_i, br$mask)
    timeseries.nii = asNifti(br$R_i)
    
    if (!all(is.na(spatial_maps_save))){
      save.fname1 = spatial_maps_save
    }else{
      save.fname1 = sub(br_suffix, spatial_maps_suffix, br_files[s])
    }
    if (!is.na(prefix) && !grepl(prefix, save.fname1)){
      save.fname = paste0(prefix, save.fname1)
    }
    if (!all(is.na(time_series_save))){
      save.fname2 = time_series_save
    }else{
      save.fname2 = sub(br_suffix, time_series_suffix, br_files[s])
    }
    if (!is.na(prefix) && !grepl(prefix, save.fname2)){
      save.fname = paste0(prefix, save.fname2)
    }
    if (verbose){  cat(paste0('\n......& saving as:  ', save.fname1))  }
    writeNifti(spatialmaps.nii, save.fname1)
    if (verbose){  cat(paste0('\n......& saving as:  ', save.fname2))  }
    writeNifti(timeseries.nii, save.fname2)
  }
  
  if (verbose){  cat(paste0('\n...saved in dir.:  ', unique(dirname(save.fname1)), '\n\n'))  }
} ##################################################################


##################################################################
group_create_b.r.sourceTable <- function(br, save_path, verbose = F){
  ##################################################################
  # Creates table of info for back-reconstructed spatial maps
  #   for group-level hierarchical PCA or spatial ICA output
  #
  # Input:
  #   br        : back-reconstructed data output of group_HierBackReconstruction,
  #                 w/ or w/o scaling by group_ScaleSources()
  #   save_path : file name & path for saved .csv table w/ component info
  # Output:    .csv table w/ columns:
  #     Component               : filename & 4d-nifti index for component
  #     Hierarchy Level         : Level of component in hierarchy, NA for non-hPCA comps.
  #     Hierarchy Variable Type : Options:
  #                               'initial sig.' : comp. is from lowest (significant) level of hierarchy
  #                               'pre-merge_var1' : comp. is 1st of merged vars. at hierarchy level, pre-merged state
  #                               'pre-merge_var2' : comp. is 2nd of merged vars. at hierarchy level, pre-merged state
  #                               'sum' : comp. represents merged shared variance in hPCA (~PC1)
  #                               'diff': comp. represents difference of merged vars. in hPCA (~PC2)
  #     Variable Index          : Index for comp. w/n ordered set of original variables,
  #                                 used for locating component within data
  #     Variable Source Level   : Previous level appearance of variable within hierarchy,
  #                                 used to compare merging of components across levels
  #     Algorithm               : Hierarchical PCA algorithm
  #
  
  ### Inputs/Prelims ###
  br = as.list(br)
  save_path = as.character(save_path)
  
  stopifnot('S_i' %in% names(br))
  K = ncol(br$S_i)
  
  csv.header = c()
  csv.body = matrix(NA,K,0)
  
  ### Add "filename,k" index for all comps. w/n 4d nii vol. ###
  k_inds = paste0(basename(save_path))
  k_inds = sub('\\.nii$', '', k_inds)
  k_inds = sub('\\.RData$', '', k_inds)
  k_inds = sub('\\.csv$', '', k_inds)
  k_inds = paste0(k_inds,',',as.integer(1:K))
  stopifnot(length(k_inds) == nrow(csv.body))
  csv.header = c(csv.header, 'Component')
  csv.body = cbind(csv.body, k_inds, deparse.level=0)
  
  ### Add Hierarchy Level ###
  if ('Hierarchy_level' %in% names(br)){
    Hierarchy_level = br$Hierarchy_level
  }else if (exists('Space.br') && ('Hierarchy_level' %in% names(Space.br))){
    Hierarchy_level = get('Hierarchy_level', envir=Space.br)
  }else{
    Hierarchy_level = rep(NA, K)
  }
  if (!all(is.na(Hierarchy_level))){
    stopifnot(length(Hierarchy_level) == nrow(csv.body))
    csv.header = c(csv.header, 'Hierarchy Level')
    csv.body = cbind(csv.body, Hierarchy_level, deparse.level=0)
  }
  
  ### Add Hierarchy algorithm ###
  if ('hier_algorithm' %in% names(br)){
    varType = rep(NA,K)
    init0 = pc1 = pc2 = ic1 = ic2 = F
    
    hpca = which(br$hier_algorithm == 'hPCA')
    if (length(hpca) > 0){
      if ('hier_k1_inds' %in% names(br)){
        if (!all(is.na(br$hier_k1_inds))){  pc1 = br$hier_k1_inds[br$hier_k1_inds %in% hpca]  }
        if (any(as.logical(pc1))){
          varType[pc1] = 'sum'  # shared variance of merged vars. at each level of hierarchy (~PC1 for hPCA)
        }
      }
      if ('hier_k2_inds' %in% names(br)){
        if (!all(is.na(br$hier_k2_inds))){  pc2 = br$hier_k2_inds[br$hier_k2_inds %in% hpca]  }
        if (any(as.logical(pc2))){
          varType[pc2] = 'diff' # differences between merged vars. at each levelof hierarchy (~PC2 for hPCA)
        }
      }
    }
    
    if ('hier_prev1_inds' %in% names(br)){
      if (!all(is.na(br$hier_prev1_inds))){
        varType[br$hier_prev1_inds] = 'pre-merge_var1'  # 1st merged var. at each level of hierarchy
      }
    }
    if ('hier_prev2_inds' %in% names(br)){
      if (!all(is.na(br$hier_prev2_inds))){
        varType[br$hier_prev2_inds] = 'pre-merge_var2'  # 2nd merged var. at each level of hierarchy
      }
    }
    if (('hier_ns_inds' %in% names(br)) && !all(is.na(br$hier_ns_inds))){
      init0 = br$hier_ns_inds[br$hier_ns_inds %in% c(hpca, hica)]
      varType[init0] = 'initial sig.'  # indicates vars. from lowest significant level of hierarchy
    }
    
  }
  if (exists('varType') && !all(is.na(varType))){
    stopifnot(length(varType) == nrow(csv.body))
    csv.header = c(csv.header, 'Hierarchy Variable Type')
    csv.body = cbind(csv.body, varType, deparse.level=0)
  }
  
  ### Add prev. level of vars. in hierarchy ###
  prev_level_inds = rep(NA, K)
  if ((c('levels_table', 'orig_var_inds') %in% names(br)) && exists('Hierarchy_level')){
    levels_table = br$levels_table
    if (is.character(levels_table) && file.exists(levels_table) && grepl('.csv$', levels_table)){
      levels_table = read.csv(levels_table, stringsAsFactors=F, header=T)
    }else{  stopifnot(grepl('.csv$', levels_table))  }
    for (k in 1:K){
      if (grepl('pre-merge_var', varType[k]) && is.numeric(br$orig_var_inds[k]) && (Hierarchy_level[k] > 1)){
        l = Hierarchy_level[k]
        prev_level_inds[k] = max(which(levels_table$Index1[1:(l-1)] == br$orig_var_inds[k]), 0)
      }
    }
  }
  if (!all(is.na(prev_level_inds))){
    stopifnot(length(prev_level_inds) == nrow(csv.body))
    csv.header = c(csv.header, 'Prev. Hierarchy Source Level')
    csv.body = cbind(csv.body, prev_level_inds, deparse.level=0)
  }
  
  ### Add indices of vars. in saved Hierarchy files ###
  orig_var_inds = rep(NA, K)
  if ('orig_var_inds' %in% names(br)){
    orig_var_inds = br$orig_var_inds
  }
  if (!all(is.na(orig_var_inds))){
    stopifnot(length(orig_var_inds) == nrow(csv.body))
    csv.header = c(csv.header, 'Variable Index')
    csv.body = cbind(csv.body, orig_var_inds, deparse.level=0)
  }
  
  ### Add Hierarchy algorithm ###
  if ('hier_algorithm' %in% names(br)){
    hier_algorithm = br$hier_algorithm
  }else if (exists('Space.br') && ('hier_algorithm' %in% names(Space.br))){
    hier_algorithm = get('hier_algorithm', envir=Space.br)
  }else{     # default to hierarchical PCA
    hier_algorithm = rep('hPCA', K)
  }
  if (!all(is.na(hier_algorithm))){
    stopifnot(length(hier_algorithm) == nrow(csv.body))
    csv.header = c(csv.header, 'Algorithm')
    csv.body = cbind(csv.body, hier_algorithm, deparse.level=0)
  }
  
  ### Save table ###
  csv.body = rbind(csv.header, csv.body, deparse.level=0)
  csv.fname = save_path
  if (!grepl('\\.csv$', csv.fname)){
    csv.fname = paste0(csv.fname, '.csv')
  }
  write.table(csv.body, file=csv.fname, 
              sep=",", qmethod="double",
              row.names=F, col.names=F)
  if (verbose){  cat(paste0('\n...saving table of spatial map info as:  ', basename(csv.fname)))  }
  
} ##################################################################


##################################################################
create_hier_StatsTable <- function(hier, save_path, verbose = F){
  ##################################################################
  # Create table w/ statistical testing info,
  #   uncorrected p-values by hierarchy level
  #
  # Input:
  #   hier      : output from get_hier_dat() & hier_apply_stats()
  #   save_path : path to dir. for saved files
  # Output:    .csv table w/ hierarchy level in 1st col.
  #               & p-values for each test in subsequent cols.
  #
  
  stopifnot(is.list(hier) && (length(hier) > 0))
  stopifnot(all(c('stats.tests',
                  'stats.p.values',
                  'stats.MCP.corrected',
                  'stats.MCP.method',
                  'stats.MCP.p.value.base') %in% names(hier[[1]])))
  if (!dir.exists(save_path)){  dir.create(save_path)  }
  
  n = length(hier[[1]]$stats.p.values)
  levels_stats = as.data.frame(matrix(NA, length(hier), n+2))
  MCP.colname = paste0(hier[[1]]$stats.MCP.method,'_', 
                       hier[[1]]$stats.MCP.p.value.base)
  stats.colnames = hier[[1]]$stats.tests
  stats.colnames = stats.colnames[-which(stats.colnames %in% c('FDR', 'FWE'))]
  names(levels_stats) = c('Level', stats.colnames, MCP.colname)
  mode(levels_stats$Level) = 'integer'
  mode(levels_stats[[MCP.colname]]) = 'logical'
  
  for (h in 1:length(hier)){
    levels_stats[h,1]          = hier[[h]]$h.level
    levels_stats[h,c(2:(n+1))] = hier[[h]]$stats.p.values
    levels_stats[h,n+2]        = all(hier[[h]]$stats.MCP.corrected)
  }
  
  save.fname = file.path(save_path, 'hierarchy_levels_stats.csv')
  if (verbose){  cat(paste0('\n......saving table of stats. applied to hierarchy levels:  \n         ', save.fname))  }
  write.csv(levels_stats, file=save.fname, row.names=F)
  
} ##################################################################
