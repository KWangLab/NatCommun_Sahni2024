#IRIS package
#Sahil Sahni & Kun Wang 07/31/24

## SPECIAL

## ------- Functions for SPECIAL -------

## Main Function:

#' SPECIAL
#'
#' A function to infer cell-2-cell interactions (SPECIAL steps 1-2) from a single spatial sample using rslurm.
#' @param expr expression matrix (normalized units), rows are genes, columns are cells. Set rownames to gene symbols and column names to cell ids
#' @param pairs a data frame of ligand receptor gene pairs, column V1 is the ligand, column V2 is the receptor (complex with up-to 3 genes). This data.frame does not have to be prefiltered, and could be a database of plausible ligand-receptor pairs (as in SOCIAL.query_LRdb).
#' @param loc a dataframe describing the mapping of each cell to an individual location and cell type. Direct output of CytoSPACE
#' @param name name of tumor slide or sample
#' @param n_iterations number of bootstrap iterations to infer the null distribution
#' @param distance distance of individual spatial region (for bulk it is the radius; for SlideSeqV2 it is the diameter)
#' @param platform 'slideseqv2' or 'bulk'
#' @param puck_diameter diameter (in unit "x") of the SlideSeqV2 puck (for SlideSeqV2 ONLY)
#'
#' @return list including empirical p-value, mean ligand expression, and mean receptor expression matrices, rows are cell-cell interactions, columns are samples. In addition, cell fraction and cell count of each cell type in a specific region, as well as regions removed due to only a single cell type being present (thus paracrine interactions cannot be calculated)
#'
#' @export
SPECIAL <- function(expr, pairs, loc, name, n_iterations, distance, platform, puck_diameter=NA){
  print('running SPECIAL')

  ## source SPECIAL

  ## initialize seed
  set.seed(20892)

  ## create directory to save intermediate input
  DIR = getwd()
  dir.create(paste(DIR,'/temporary/',sep=""))
  path_1 = paste(DIR,'/temporary/',sep="")
  dir.create(paste(path_1,name,'/',sep=""))
  temp_path_cis = paste(path_1, name, '/', sep='') #path to temporary folder

  ## run step 1 of SPECIAL ##

  # if legacy or visium
  if (toupper(platform) == 'BULK'){
    print('input bulk')

    region = SPECIAL.sliding_window(loc, distance)
  }

  # if slideseqv2
  if (toupper(platform) == 'SLIDESEQV2'){
    print('input SLIDESEQV2')

    region = SPECIAL.cluster_region(loc, distance, puck_diameter)
  }

  ## calculate cell fraction

  cellcount = region %>%
    dplyr::group_by(`samples`, cell_type) %>%
    dplyr::summarise(count=n()) %>%
    tidyr::pivot_wider(., names_from=cell_type, values_from=count) %>% as.data.frame(.) %>%
    magrittr::set_rownames(.$samples) %>% .[,-1] %>% replace(., is.na(.), 0)
  cellfraction = cellcount/rowSums(cellcount)

  ## remove regions ("samples" ~ to seamlessly apply to SOCIAL) and cellfraction with only one cell type
  removed_regions <- names(which(apply(cellfraction == 1, 1, any)))

  cellfraction = cellfraction[!(rownames(cellfraction) %in% removed_regions),]
  region = subset(region, !(samples %in% removed_regions))

  ## run step 1 of SOCIAL ##
  pairs = SOCIAL.query_LRdb(pairs,expr)

  ## save intermediate files for SOCIAL step 2 ##

  ## generate intermediate ct_map per region ("samples")
  ct_map_input = generate_ctmap_per_sample(region, expr)
  ct_map_path = paste(temp_path_cis,'intermediate_ctmap.RDS',sep='')
  saveRDS(ct_map_input, file=ct_map_path)

  ## save sc expr files as sparse matrix
  sparse_expr=Matrix::Matrix(expr,sparse = T)
  sparse_path = paste(temp_path_cis,'sparse_expr.RDS',sep='')
  saveRDS(sparse_expr, file=sparse_path)

  ## save pairs files
  pairs_path = paste(temp_path_cis,'LR_pairs_prefiltered.csv',sep='')
  write.csv(pairs, pairs_path)

  ## generate rslurm data.frame
  run = 1:length(unique(region$samples))
  samples = unique(region$samples)
  expr_path = sparse_path
  ct_map_path = ct_map_path
  pairs_path = pairs_path
  n_iterations = n_iterations

  rslurm_df = data.frame(run=run, samples=samples, expr_path=expr_path, ct_map_path=ct_map_path, pairs_path=pairs_path, n_iterations=n_iterations)

  ## run rslurm (SOCIAL Step 2-3) ##
  sjob = rslurm::slurm_apply(run_cis_rslurm, rslurm_df, jobname=name, nodes=nrow(rslurm_df), cpus_per_node=8,submit=TRUE, slurm_options=list(time='72:00:00', mem='280g'), preschedule_cores=FALSE)
  profile = rslurm::get_slurm_out(sjob, outtype = 'raw', wait = TRUE)
  rslurm::cleanup_files(sjob) #remove slurm folder
  names(profile) = rslurm_df$samples

  ## creates data.frame of all possible ligand_receptor pairs x samples; values are empirical p
  pval_list = lapply(names(profile), function(x){
    m_pv = reshape2::melt(profile[[x]]$p.value)
    colnames(m_pv) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_pv = m_pv %>% dplyr::mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)]
    return(m_pv)
  })

  pval_df = pval_list %>% reduce(full_join, by='int')

  ## creates data.frame of all possible ligand_receptor pairs x samples; values are ligand exp
  ligand_list = lapply(names(profile), function(x){
    m_pv = reshape2::melt(profile[[x]]$ligand_exp)
    colnames(m_pv) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_pv = m_pv %>% dplyr::mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)]
    return(m_pv)
  })

  ligand_df = ligand_list %>% reduce(full_join, by='int')

  ## creates data.frame of all possible ligand_receptor pairs x samples; values are receptor exp
  receptor_list = lapply(names(profile), function(x){
    m_pv = reshape2::melt(profile[[x]]$receptor_exp)
    colnames(m_pv) = c('ligand_receptor', "Lcell_Rcell", names(profile[x]))
    m_pv = m_pv %>% dplyr::mutate(int=paste(Lcell_Rcell, ligand_receptor, sep='_')) %>% .[,c(4,3)]
    return(m_pv)
  })

  receptor_df = receptor_list %>% purrr::reduce(full_join, by='int')

  ## delete created directory
  unlink(temp_path_cis, recursive = TRUE)

  ## return objects
  return(list(name = name, cellfraction=cellfraction, empiricalp=pval_df, removed_regions=removed_regions, ligandexp=ligand_df, receptorexp=receptor_df, cellcount=cellcount))
}

## Step 1 Functions:

#SPECIAL.cluster_region(): function to partiion spatial transcriptomics slide to physiological relevant ligand-receptor interaction regions
#loc: a dataframe describing the mapping of each cell to spatial coordinate and cell type (output of CytoSPACE. it is of the following format:
#   UniqueCID: unique sc barcode (generated by CytoSPACE)
#   OriginalCID:  original sc barcode from study
#   CellType: cell type of individual single-cell
#   SpotID: barcode or coordinate of the spatial location of the aligned single-cell
#   row: y coordinate of aligned single-cell
#   col: x coordinate of aligned single-cell
#cluster_diameter: diameter (in unit "x") of the individual regions
#puck_diameter: diameter (in unit "x") of the SlideSeqV2 puck

SPECIAL.cluster_region <- function(loc, cluster_diameter, puck_diameter){

  # create puck coordinate map
  map = as.data.frame(loc) %>% dplyr::mutate(x=row, y=col)
  map = map %>% mutate(x = (x-min(x))/max(x)*puck_diameter, y = (y-min(y))/max(y)*puck_diameter) #rescale to puck dimension
  coord = as.matrix(map[,c('x','y')]) #create matrix

  # initialize
  puck_radius = puck_diameter/2
  cluster_radius = cluster_diameter/2
  centers = as.integer((pi*puck_radius^2)/(pi*cluster_radius^2)) #number of spatial spots that can fit in the puck

  # calculate kmean clustering with ideal centers
  set.seed(20892)
  kmeans_result = stats::kmeans(coord, centers=centers)
  map$cluster = kmeans_result$cluster

  # cluster center data.frame
  cluster_coord = data.frame(cluster = 1:centers, kmeans_result$centers) %>% dplyr::mutate(regions=paste(round(x,4),round(y,4),sep='x'))

  # merge data.frame and filter
  res = merge(map, cluster_coord, by.x='cluster', by.y='cluster') %>% dplyr::select(Cell=OriginalCID, cell_type=CellType, samples=regions)
  return(res)
}

#SPECIAL.sliding_window(): function to partiion spatial transcriptomics slide to physiological relevant ligand-receptor interaction regions
#loc: a dataframe describing the mapping of each cell to spatial coordinate and cell type (output of CytoSPACE. it is of the following format:
#   UniqueCID: unique sc barcode (generated by CytoSPACE)
#   OriginalCID:  original sc barcode from study
#   CellType: cell type of individual single-cell
#   SpotID: barcode or coordinate of the spatial location of the aligned single-cell
#   row: y coordinate of aligned single-cell
#   col: x coordinate of aligned single-cell
#radius: radius (in unit "x") of the individual region

SPECIAL.sliding_window <- function(loc, radius=1){
  map = data.frame(x=loc$row, y=loc$col) %>% distinct(.)

  windows = lapply(1:nrow(map), function(o){
    origin = map[o,]

    window = subset(loc, sqrt((row - origin$x)^2 + (col - origin$y)^2) <= radius) #euclidean distance
    window$window = paste(origin$x,'x',origin$y,sep='')

    window = window %>% dplyr::select(Cell=OriginalCID, cell_type=CellType, samples=window)
    return(window)
  })
  windows = do.call(rbind,windows)
  return(windows)
}

## Additional Functions:

#' SPECIAL for Cohorts
#'
#' A function to infer cell-2-cell interactions (SPECIAL steps 1-2) from a spatial cohort (for each spatial sample) using rslurm (ideal for when multiple samples are available).
#' @param path_input SPECIAL organized input
#' @param n_iterations number of bootstrap iterations to infer the null distribution
#' @param distance distance of individual spatial regions (for bulk it is the radius; for SlideSeqV2 it is the diameter)
#' @param platform 'slideseqv2' or 'bulk'
#' @param puck_diameter diameter (in unit "x") of the SlideSeqV2 puck (for SlideSeqV2 ONLY)
#'
#' @return list including empirical p-value, mean ligand expression, and mean receptor expression matrices, rows are cell-cell interactions, columns are samples. In addition, cell fraction and cell count of each cell type in a specific region, as well as regions removed due to only a single cell type being present (thus paracrine interactions cannot be calculated).
#'
#' @export
SPECIAL.cohort <- function(path_input, n_iterations, distance, platform, puck_diameter=NA){

  print('running SPECIAL.cohort')

  #load SPECIAL and its input
  special_input = readRDS(path_input)

  #create data.frame to run rslurm across different slides/samples
  rslurm_df = data.frame(index=1:length(special_input), path_input=path_input, n_iterations=n_iterations, distance=distance, platform=platform, puck_diameter=puck_diameter)
  jobname = paste('per','index',Sys.time(), sep='_')

  #run rslurm
  sjob_out = rslurm::slurm_apply(SPECIAL.per_sample, rslurm_df, jobname=jobname, nodes=nrow(rslurm_df), cpus_per_node=8,submit=TRUE, slurm_options=list(time='72:00:00', mem='10g'), preschedule_cores=FALSE)
  output = rslurm::get_slurm_out(sjob_out, outtype = 'raw', wait = TRUE)
  rslurm::cleanup_files(sjob_out)

  return(output)
}

## Internal Functions:

#SPECIAL.per_sample(): A internal function to execute SPECIAL on each sample. Can be expanded in the future to add Null distribution.

SPECIAL.per_sample <- function(index, path_input, n_iterations, distance, platform, puck_diameter=NA){
  print('running SPECIAL.per_sample')

  #load SPECIAL

  #data.frame for rslurm
  rslurm_df = data.frame(index=index, path_input=path_input, path_special=path_special, n_iterations=n_iterations, distance=distance, platform=platform, puck_diameter=puck_diameter)

  #run rslurm per slide
  sjob_in = rslurm::slurm_apply(run_special_rslurm, rslurm_df, jobname=paste('index',index,sep="_"), nodes=nrow(rslurm_df), cpus_per_node=10,submit=TRUE, slurm_options=list(time='72:00:00', mem='40g'), preschedule_cores=FALSE) #2 before
  output = rslurm::get_slurm_out(sjob_in, outtype = 'raw', wait = TRUE)
  rslurm::cleanup_files(sjob_in)

  #outputs
  return(output)
}

#run_special_rslurm(): A internal function to execute SPECIAL on rslurm

run_special_rslurm <- function(index, path_input, name, n_iterations, distance, platform, puck_diameter=NA){
  print('running run_special_rslurm')

  #load special functions

  #load input path
  input = readRDS(path_input)

  #load indexed sample
  indexed_sample = input[[index]]
  expr = indexed_sample$expr
  loc = indexed_sample$loc
  name = indexed_sample$name
  pairs = indexed_sample$pairs # could be prefiltered or a unfiltered (ligand-receptor data.base as in SOCIAL.query_LRdb)
  distance = distance
  platform = platform
  puck_diameter = puck_diameter

  #run SPECIAL step 1 and 2
  results = SPECIAL(expr, pairs, loc, name, n_iterations, distance, platform, puck_diameter)
  return(results)
}
