# Data availability
The CODEFACS, LIRICS, SOCIAL, CytoSPACE, and SPECIAL data (and relevant data) generated in this study have been deposited to https://zenodo.org/records/13172848.

# Data Directory:
## a. benchmark
results regarding applying existing transcriptomics-based biomarkers to predict ICB response.

### i. signatures
1. **```impres_auc_supp.RDS```**  AUCs for IMPRES on bulk RNA-seq cohort (mono./comb. therapy separated)
2. **```sig_auc.RDS```**          AUCs for signatures (including IMPRES) on bulk RNA-seq cohort
3. **```sig_input_thrane.RDS```** raw signature scores on Thrane et al. spatial cohort
4. **```sig_input.RDS```**        raw signature scores on bulk RNA-seq cohort

### ii. TIDE
1. **```tide_auc_supp.RDS```**    AUCs for TIDE on bulk RNA-seq cohort (mono./comb. therapy seperated)
2. **```tide_auc.RDS```**         AUCs for TIDE on bulk RNA-seq cohort
3. **```tide_input.RDS```**       raw signature scores for TIDE on bulk RNA-seq cohort

## b. CODEFACES
deconvolved expression data genereated for this study. 

### i. output
1. **```GE_imp_auslander.rds```** CODEFACS' deconvolved expression for Auslander et al.
2. **```GE_imp_puch.rds```**      CODEFACS' deconvolved expression for puch et al. 

## c. CytoSPACE
CytoSPACE output for Biermann et al. and Thrane et al. For Biermann et al., CytoSPACE was ran on matched SlideSeqV2 and snRNA-seq data. The "firstctassigned" and "secondctassigned" corresponds to the location of the first and second cell type inferred respectively. For Thrane et al., CytoSPACE ran on unmatched Legacy ST and scRNA-seq (Jerby-Arnon et al. Cell 2018).

### i. results
1. **Biermann et al.**
2. **Thrane et al.**

## c. IRIS
results regarding applying IRIS to identify resistance relevant ligand-receptor interactions.

### i. results
1. **```IRIS_RDI_output_mc.data```** inferred RDIs' for testing on mono/comb ICB therapy cohorts seperately
2. **```IRIS_RDI_output.data```** inferred RDIs' for testing on ICB cohorts
3. **```IRIS_RUI_output_mc.data```** inferred RUIs' for testing on mono/comb ICB therapy cohorts seperately
4. **```IRIS_RUI_output.data```** inferred RUIs' for testing on ICB cohorts

## d. LIRICS
LIRICS data genereated for this study. 

### i. LIRICS database
1. **```lirics_skcm_database_interactions.RDS```** list of all possible literature-curated cell-type-specific ligand-receptor interactions within the Melanoma TME across the 10 cell types deconvolved.

### ii. organized data/ICB
1. **```IRIS_LIRICS_final_ICB_cohort_input_OS.Rdata```** Organized LIRICS data for ICB cohorts with overall survival timeline information
2. **```IRIS_LIRICS_final_ICB_cohort_input_PFS.Rdata```** Organized LIRICS data for ICB cohorts with progression free survival timeline information
3. **```IRIS_LIRICS_final_ICB_cohort_input_response.Rdata```** Organized LIRICS data for ICB cohorts with ICB response information (input for IRIS)

### iii. organized data/TCGA
1. **```TCGA_LIRICS_SURV.Rdata```** Organized LIRICS data for TCGA-SKCM with prorgession free and overall survival timeline information
2. **```TCGA_LIRICS.Rdata```** Organized LIRICS data for TCGA-SKCM with additional meta information (deconvolved cell fraction, hot/cold TME, etc.)
   * *note that "tmb", "aneuoploidy", and "purity" were **NOT** used in this study and are from: [Taylor et al. Cell 2018](https://pubmed.ncbi.nlm.nih.gov/29622463/)*

## e. miscellaneous
1. **```coding_genes.txt```** list of protein coding genes

## f. sc ICB data
### i. Jerby-Arnon et al. 2018 Cell
1. **```Livnat_updated_celltypes_table.rds```** Reclustered cell type meta informatoin used for Jerby-Arnon et al.'s single-cell cohort (i.e. for SOCIAL)

## g. SOCIAL
SOCIAL data genereated for this study. 

### i. output/Jerby-Arnon et al. Cell 2018
1. **```livnat2018_TPM_400pseudopat_40ds_SOCIAL.rds```** SOCIAL output

## h. SPECIAL
SPECIAL data genereated for this study. 

### i. Biermann et al. 2022/organized input
1. **```biermann2022_SPECIAL_input.RDS```** SPECIAL input

### ii. Biermann et al. 2022/output
1. **```biermann2022_SPECIAL_input.RDS```** SPECIAL output

### iii. Thrane et al. 2018/organized input
1. **```thrane2018_SPECIAL_input.RDS```** SPECIAL input

### iv. Thrane et al. 2018/output
1. **```thrane2018_SPECIAL_input.RDS```** SPECIAL output

## i. Spatial Seurat
Cell fraction data generated for input to CytoSPACE.

### i. Thrane et al./cell fraction
Spatial Seurat deconvolved cell fraction at each spatial coordinate

### i. Thrane et al./weights
Spatial Seurat deconvolved cell fraction across the entire spatial slide (used for CytoSPACE)

# Citation
If using IRIS, SOCIAL, SPECIAL, or any of the pertaining results and data, please cite:

Sahni et al. "A machine learning model reveals expansive downregulation of ligand-receptor interactions enhancing lymphocyte infiltration in melanoma with acquired resistance to Immune Checkpoint Blockade. *N C* **X**, XXXX (XXXX). https://doi.org

# Corresponding Author(s)
Kun Wang (kun.wang@nih.gov) and Eytan Ruppin (eytan.ruppin@nih.gov)
