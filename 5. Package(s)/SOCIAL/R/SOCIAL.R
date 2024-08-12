#IRIS package
#Sahil Sahni & Kun Wang 07/31/24

## SOCIAL

## ------- Functions for SOCIAL -------

## Step 1 Functions:

#' SOCIAL step 1
#'
#' A function to extract ligand receptor pairs that are found within our sc data.
#' @param db a data frame of known ligand receptor gene pairs, column V1 is the ligand, column V2 is the receptor (complex with up-to 3 genes seperated by ";", i.e. A;B;C).
#' @param expr expression matrix (normalized units), rows are genes, columns are cells.
#'
#' @return pairs: a data frame of ligand receptor gene pairs, column V1 is the ligand, column V2 is the receptor
#'
#' @export
SOCIAL.query_LRdb <- function(db, expr){

  sc_gene = toupper(rownames(expr))

  pairs = toupper(as.matrix(db))
  pairs = unique(pairs)

  pairs_check = pairs %>% as.data.frame(.) %>% dplyr::mutate(receptor=str_split_fixed(pairs[,2], ";", 3)) %>% dplyr::mutate_all(na_if, "") %>% as.matrix(.)

  ind1 = (pairs_check[,1] %in% sc_gene);
  ind2 = (pairs_check[,2] %in% sc_gene);
  ind3 = (pairs_check[,3] %in% sc_gene)|(pairs_check[,3] %in% c(""));
  ind4 = (pairs_check[,4] %in% sc_gene)|(pairs_check[,4] %in% c(""));

  pairs = pairs[ind1&ind2&ind3&ind4,]

  return(pairs)
}

## Step 2 & 3 Functions:

#'SOCIAL step 2 & 3
#'
#' A function to infer cell-2-cell interactions from sc data.
#' @param expr expression matrix (normalized units), rows are genes, columns are cells. Set rownames to gene symbols and column names to cell ids.
#' @param ct_map a dataframe describing the mapping of each cell to individual and cell type.
#' @param pairs a data frame of ligand receptor gene pairs, column V1 is the ligand, column V2 is the receptor
#' @param n_iterations number of bootstrap iterations to infer the null distribution
#'
#' @return list including a interaction score, empirical p-value, mean ligand expression, and mean receptor expression matrices, rows are cell-cell interactions, columns are samples.
#'
#' @export
SOCIAL.calculate_interaction_score <- function(expr,ct_map, pairs, n_iterations = 100){ #for spatial data

  LR_interest = c(pairs$ligand, pairs$receptor) %>% strsplit(., ";") %>% unlist(.) %>% unique(.) # LR of interest to speed of data aggregation

  ctypes = unique(ct_map$cell_type)
  ctype_groupings = t(combn(ctypes,2))
  ctype_groupings = rbind(ctype_groupings,ctype_groupings[,c(2,1)])

  #compute real scores (Step 2)
  summary_expr = stats::aggregate.data.frame(as.data.frame(t(expr[LR_interest,])), by = list(ct_map$cell_type), FUN = mean, na.rm = T)
  rownames(summary_expr) = summary_expr$Group.1
  summary_expr = summary_expr[,-1]
  summary_expr = t(summary_expr)

  summary_expr2 = stats::aggregate.data.frame(as.data.frame(t(expr[LR_interest,]) > 0), by = list(ct_map$cell_type), FUN = mean, na.rm = T)
  rownames(summary_expr2) = summary_expr2$Group.1
  summary_expr2 = summary_expr2[,-1]
  summary_expr2 = t(summary_expr2)

  real_scores = sapply(1:nrow(ctype_groupings), function(j) { # compute interaction scores
    ct1 = ctype_groupings[j,1]
    ct2 = ctype_groupings[j,2]
    scores = suppressMessages(t(apply(pairs,1,function(y){
      g1 = strsplit(y[1],split = ";")[[1]]
      g2 = strsplit(y[2],split = ";")[[1]]
      score = c()
      ligand_exp = NULL
      receptor_exp = NULL
      if(length(g1) > 1)
        ligand_exp = min(as.numeric(summary_expr[g1,ct1]))
      else
        ligand_exp = as.numeric(summary_expr[g1,ct1])
      if(length(g2) > 1)
        receptor_exp = min(as.numeric(summary_expr[g2,ct1]))
      else
        receptor_exp = as.numeric(summary_expr[g2,ct1])

      #interaction score
      sc = ligand_exp*receptor_exp
      return(sc)
    })))
    return(scores)
  })

  ligand_scores = sapply(1:nrow(ctype_groupings), function(j) { # compute mean ligand expression
    ct1 = ctype_groupings[j,1]
    ct2 = ctype_groupings[j,2]
    scores = suppressMessages(t(apply(pairs,1,function(y){
      g1 = strsplit(y[1],split = ";")[[1]]
      score = c()
      ligand_exp = NULL
      if(length(g1) > 1)
        ligand_exp = min(as.numeric(summary_expr[g1,ct1]))
      else
        ligand_exp = as.numeric(summary_expr[g1,ct1])

      #ligand expression
      sc = ligand_exp
      return(sc)
    })))
    return(scores)
  })

  receptor_scores = sapply(1:nrow(ctype_groupings), function(j) { # compute mean receptor expression
    ct1 = ctype_groupings[j,1]
    ct2 = ctype_groupings[j,2]
    scores = suppressMessages(t(apply(pairs,1,function(y){
      g2 = strsplit(y[2],split = ";")[[1]]
      score = c()
      receptor_exp = NULL
      if(length(g2) > 1)
        receptor_exp = min(as.numeric(summary_expr[g2,ct1]))
      else
        receptor_exp = as.numeric(summary_expr[g2,ct1])

      #receptor expression
      sc = receptor_exp
      return(sc)
    })))
    return(scores)
  })

  rownames(real_scores) = sapply(1:nrow(pairs), function(i) paste(pairs[i,1],"_",pairs[i,2],sep = ""))
  colnames(real_scores) = sapply(1:nrow(ctype_groupings), function(i) paste(ctype_groupings[i,1],"_",ctype_groupings[i,2],sep = ""))

  rownames(ligand_scores) = sapply(1:nrow(pairs), function(i) paste(pairs[i,1],"_",pairs[i,2],sep = ""))
  colnames(ligand_scores) = sapply(1:nrow(ctype_groupings), function(i) paste(ctype_groupings[i,1],"_",ctype_groupings[i,2],sep = ""))

  rownames(receptor_scores) = sapply(1:nrow(pairs), function(i) paste(pairs[i,1],"_",pairs[i,2],sep = ""))
  colnames(receptor_scores) = sapply(1:nrow(ctype_groupings), function(i) paste(ctype_groupings[i,1],"_",ctype_groupings[i,2],sep = ""))

  #compute null scores (Step 3)
  null_scores = lapply(1:n_iterations, function(i) {
    start = Sys.time()
    usam = unique(ct_map$samples)
    null_labels = ct_map$cell_type
    for(samp in usam){
      null_labels[ct_map$samples == samp] = sample(ct_map$cell_type[ct_map$samples == samp]) #shuffle cell labels
    }
    null_map = data.frame(cell = ct_map$Cell, cell_type = sample(null_labels))
    summary_expr = stats::aggregate.data.frame(as.data.frame(t(expr[LR_interest,])), by = list(null_map$cell_type), FUN = mean, na.rm = T)
    rownames(summary_expr) = summary_expr$Group.1
    summary_expr = summary_expr[,-1]
    summary_expr = t(summary_expr)

    summary_expr2 = stats::aggregate.data.frame(as.data.frame(t(expr[LR_interest,]) > 0), by = list(null_map$cell_type), FUN = mean, na.rm = T)
    rownames(summary_expr2) = summary_expr2$Group.1
    summary_expr2 = summary_expr2[,-1]
    summary_expr2 = t(summary_expr2)
    null_score = sapply(1:nrow(ctype_groupings), function(j) {
      ct1 = ctype_groupings[j,1]
      ct2 = ctype_groupings[j,2]
      scores = suppressMessages(t(apply(pairs,1,function(y){
        g1 = strsplit(y[1],split = ";")[[1]]
        g2 = strsplit(y[2],split = ";")[[1]]
        score = c()
        ligand_exp = NULL
        receptor_exp = NULL
        if(length(g1) > 1)
          ligand_exp = min(as.numeric(summary_expr[g1,ct1]))
        else
          ligand_exp = as.numeric(summary_expr[g1,ct1])
        if(length(g2) > 1)
          receptor_exp = min(as.numeric(summary_expr[g2,ct1]))
        else
          receptor_exp = as.numeric(summary_expr[g2,ct1])

        #null score
        sc = ligand_exp*receptor_exp
        return(sc)
      })))
      return(scores)
    })
    print(Sys.time() - start)
    return(null_score)
  })

  #compute empirical p-value (Step 3)
  pvals = matrix(0,nrow = nrow(real_scores), ncol = ncol(real_scores))
  for(i in 1:n_iterations){
    pvals = pvals + (null_scores[[i]] >= real_scores) #assumption is >= so that 0 >= 0 if considered true
  }
  pvals = pvals/n_iterations
  rownames(pvals) = sapply(1:nrow(pairs), function(i) paste(pairs[i,1],"_",pairs[i,2],sep = ""))
  colnames(pvals) = sapply(1:nrow(ctype_groupings), function(i) paste(ctype_groupings[i,1],"_",ctype_groupings[i,2],sep = ""))

  return(list(interaction_scores = real_scores,p.value = pvals, ligand_exp=ligand_scores, receptor_exp=receptor_scores))
}

## Step 4 Functions:

#' SOCIAL step 4
#'
#' A function to infer ligand receptor interaction activity from SOCIAL object.
#' @param input SOCIAL object outputted from SOCIAL.cis_rslurm() or SOCIAL.calculate_interaction_score.
#' @param p empirical p-value cutoff (what empirical p value is considered significantly active)
#' @param median whether activity should implement ligand and receptor median expression (above median = active)
#'
#' @return binary matrix (0 cci is inactive; 1 cci is active), rows are cell-cell interactions, columns are samples.
#'
#' @export
SOCIAL.infer_activity <- function(input, p_cutoff=0.05, median_cutoff=T){

  # empirical p value (active if p < 0.05)
  input_p = input$empiricalp %>% magrittr::set_rownames(.$int) %>% .[,-1] %>% as.matrix(.)
  input_p = ifelse(input_p < p_cutoff, 1, 0) %>% ifelse(is.na(.), 0, .)

  if(median_cutoff==T){
    # ligand median expression (active if mean ligand expr > median)
    input_l = input$ligandexp %>% magrittr::set_rownames(.$int) %>% .[,-1] %>% as.matrix(.)
    input_l = apply(input_l, 1, function(x){ifelse(x > median(x, na.rm=T), 1,0) }) %>% t(.) %>% ifelse(is.na(.), 0, .)

    # receptor median expression (active if mean receptor expr > median)
    input_r = input$receptorexp %>% magrittr::set_rownames(.$int) %>% .[,-1] %>% as.matrix(.)
    input_r = apply(input_r, 1, function(x){ifelse(x > median(x, na.rm=T), 1,0) }) %>% t(.) %>% ifelse(is.na(.), 0, .)

    # create a 3D array
    exp_array <- abind::abind(input_p, input_l, input_r, along=3)
    dimnames(exp_array)[[3]] = list('cci','ligand','receptor')
    dim(exp_array)

    # sum across 3 array cci, ligand, receptor (1+1+1 = 3 if active)
    sum_array <- apply(exp_array, c(1, 2), sum)

    # spot score where sum == 3 is 1 (active); others is 0 (inactive)
    output_3 = ifelse(sum_array == 3, 1, 0) %>% t(.)
    return(output_3)
  }

  output_1 = input_p
  return(output_1)
}

## Additional Functions:

#' SOCIAL step 2 & 3 rslurm
#'
#' A function to infer cell-2-cell interactions (SOCIAL steps 2-3) from sc data (for each sample) using rslurm (ideal for when multiple samples are available).
#' @param expr expression matrix (normalized units), rows are genes, columns are cells. Set rownames to gene symbols and column names to cell ids.
#' @param ct_map a dataframe describing the mapping of each cell to individual and cell type.
#' @param pairs a data frame of ligand receptor gene pairs, column V1 is the ligand, column V2 is the receptor
#' @param n_iterations number of bootstrap iterations to infer the null distribution
#'
#' @return list including a interaction score, empirical p-value, mean ligand expression, and mean receptor expression matrices
#'
#' @export
SOCIAL.cis_rslurm <- function(expr,ct_map, pairs, n_iterations = 100){

  ## create jobname
  jobname = paste('SOCIAL', Sys.time())
  print(jobname)

  ## create directory to save intermediate input
  dir = getwd()
  dir.create(paste(dir,'/',jobname,'/',sep=""))
  temp_path_cis = paste(dir,'/',jobname,'/',sep="")

  ## save sc expr files as sparse matrix
  sparse_expr=Matrix::Matrix(expr,sparse = T)
  sparse_path = paste(temp_path_cis,'sparse_expr.RDS',sep='')
  saveRDS(sparse_expr, file=sparse_path)

  ## generate and save intermediate ct_map files
  ct_map_input = generate_ctmap_per_sample(ct_map, expr)
  ct_map_path = paste(temp_path_cis,'intermediate_ctmap.RDS',sep='')
  saveRDS(ct_map_input, file=ct_map_path)

  ## save pairs files
  pairs_path = paste(temp_path_cis,'LR_pairs_prefiltered.csv',sep='')
  write.csv(pairs, pairs_path)

  ## generate rslurm data.frame
  run = 1:length(unique(ct_map$samples))
  samples = unique(ct_map$samples)
  expr_path = sparse_path
  ct_map_path = ct_map_path
  pairs_path = pairs_path
  n_iterations= n_iterations

  rslurm_df = data.frame(run=run, samples=samples, expr_path=expr_path, ct_map_path=ct_map_path, pairs_path=pairs_path, n_iterations=n_iterations)

  ## run rslurm
  sjob = rslurm::slurm_apply(run_cis_rslurm, rslurm_df, jobname=jobname, nodes=400, cpus_per_node=1,submit=TRUE, slurm_options=list(time='24:00:00', mem='20g')) #CHANGE HYPERPARAMETER AS NEEDED
  profile = rslurm::get_slurm_out(sjob, outtype = 'raw', wait = TRUE) #SOCIAL results
  names(profile) = rslurm_df$samples
  rslurm::cleanup_files(sjob)

  ## aggregate results
  pval_list = lapply(names(profile), function(x){
    m_pv = reshape2::melt(profile[[x]]$p.value)
    colnames(m_pv) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_pv = m_pv %>% dplyr::mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)]
    return(m_pv)
  })

  score_list = lapply(names(profile), function(x){
    m_sc = reshape2::melt(profile[[x]]$interaction_scores)
    colnames(m_sc) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_sc = m_sc %>% dplyr::mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)]
    return(m_sc)
  })

  ligand_list = lapply(names(profile), function(x){
    m_pv = reshape2::melt(profile[[x]]$ligand_exp)
    colnames(m_pv) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_pv = m_pv %>% dplyr::mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)]
    return(m_pv)
  })

  receptor_list = lapply(names(profile), function(x){
    m_sc = reshape2::melt(profile[[x]]$receptor_exp)
    colnames(m_sc) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_sc = m_sc %>% dplyr::mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)]
    return(m_sc)
  })

  pval_df = pval_list %>% purrr::reduce(full_join, by='int')
  scr_df = score_list %>% purrr::reduce(full_join, by='int')
  ligand_df = ligand_list %>% purrr::reduce(full_join, by='int')
  receptor_df = receptor_list %>% purrr::reduce(full_join, by='int')

  unlink(temp_path_cis, recursive = TRUE) #remove temporary path

  social_profile=list(empiricalp=pval_df, interactionscore=scr_df, ligandexp=ligand_df, receptorexp=receptor_df)
  return(social_profile)
}

## Internal Functions:

#generate_ctmap_per_sample(): A internal function to generates ct_map for each sample to load faster for generating interaction score using cis_rslurm().

#ct_map: a dataframe describing the mapping of each cell to individual and cell type. it is of the following format:
#Cell cell_type samples
#         cy79_p4_CD45_neg_PDL1_neg_E11_S1115_comb       EVT   Mel79
#  cy79_p1_CD45_neg_PDL1_pos_AS_C1_R1_F07_S67_comb       EVT   Mel79
#  cy79_p1_CD45_neg_PDL1_pos_AS_C1_R1_G01_S73_comb       EVT   Mel79
# cy79_p1_CD45_neg_PDL1_neg_AS_C4_R1_D09_S141_comb       EVT   Mel79
# cy79_p1_CD45_neg_PDL1_neg_AS_C4_R1_D06_S138_comb       EVT   Mel79
#                    cy53_1_CD45_neg_C06_S318_comb       EVT   Mel53

#expr: expression matrix (normalized units), rows are genes, columns are cells. columns is non-subsetted by samples!

generate_ctmap_per_sample <- function(ct_map, expr){
  samples = unique(ct_map$samples)

  ct_map_list = parallel::mclapply(1:length(samples), function(i){ #generates ct_map files for each sample
    sample_i = samples[i]
    ct_map_sample = subset(ct_map, (Cell %in% colnames(expr)) & (samples == sample_i)) #ct_map
    return(ct_map_sample)
  })

  names(ct_map_list) = samples
  return(ct_map_list)
}

#run_cis_rslurm(): A internal function to execute calculate interaction score on rslurm efficiently

run_cis_rslurm <- function(run, samples, expr_path, ct_map_path, pairs_path, n_iterations){

  #CODE START
  set.seed(run+1)

  #load sc data
  expr = readRDS(expr_path) %>% as.matrix(.)
  print('expr loaded')

  #load filtered LR pairs
  pairs = read.csv(pairs_path) %>% dplyr::select(ligand,receptor)
  print('pairs loaded')

  #load ct_map
  load_ctmap = readRDS(ct_map_path)
  ct_map = load_ctmap[[samples]] %>% dplyr::select(Cell, cell_type, samples)
  print('ctmap loaded')

  #filter expr (assuming duplicate cell IDs):
  ct_map = ct_map %>% dplyr::mutate(uID = paste(Cell,1:nrow(ct_map), sep="nci")) #append character to help prevent duplicates
  id_df = data.frame(Cell=ct_map$Cell, uID=ct_map$uID)
  expr_filter = merge(id_df, t(expr), by.x=1, by.y=0) %>% magrittr::set_rownames(.$uID) %>% .[,-1:-2] %>% as.matrix(.) %>% t(.)
  ct_map = ct_map %>% dplyr::mutate(Cell = uID) %>% dplyr::select(-uID)
  expr_filter = expr_filter[,ct_map$Cell] #reorder
  colnames(expr_filter)

  SOCIAL_res = SOCIAL.calculate_interaction_score(expr_filter,ct_map, pairs, n_iterations)
  print('SOCIAL Step 2&3 Complete')

  return(SOCIAL_res)
}

#' Generate pseudopatients ct_map
#'
#' A function used to randomly down sample x% of cells from each cell type for each group (i.e. post-ICI or naive) to generate n number of pseudopatients for each group.
#'@param n number of pseudopatients for each group
#'@param prop fraction (in decimal) of cells to be downsampled from each cell type for each group
#'@param group_map data.frame with columns Cell, cell_type, samples, **and** group
#'
#'@return ct_map: a SOCIAL input
#'
#'@export
generate_pseudopatient <- function(n, ct_info, prop=0.4){
  output = ct_info %>% group_by(group, cell_type) %>% slice_sample(prop=prop) %>% ungroup(.) %>% mutate(samples=paste(group,n,sep="")) %>% dplyr::select(Cell, cell_type, samples)
  return(output)
}
