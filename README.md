# A machine learning model reveals expansive downregulation of ligand-receptor interactions enhancing lymphocyte infiltration in melanoma with acquired resistance to Immune Checkpoint Blockade
**Last Updated: 08/14/2024**

Code repository to reproduce the results and findings that are published in: 

Sahni et al. "A machine learning model reveals expansive downregulation of ligand-receptor interactions enhancing lymphocyte infiltration in melanoma with acquired resistance to Immune Checkpoint Blockade. *N C* **X**, XXXX (XXXX). https://doi.org

### IRIS: **I**mmunotherapy **R**esistance *cell-cell* **I**nteraction **S**canner
<img src="https://github.com/sahil-sahni/IRIS/blob/main/4.%20Figure/images/png/IRIS%20figure%201%20Final%20Version%20%5Bnc%20acc%5D.png" alt="grouping">

We developed **I**mmunotherapy **R**esistance cell-cell **I**nteraction **S**canner (IRIS), a computational method specifically designed to identify immune checkpoint blockade (ICB) resistance relevant ligand-receptor interactions in the tumor microenvironment (TME), given a patients cohort including tumor bulk expression data and ICB treatment response data. The gene expression data is deconvolved using [**CODEFACS**](https://pubmed.ncbi.nlm.nih.gov/34983745/) such that the input to IRIS in a given patients cohort is comprised of two components: 1. Literature-curated cell-type-specific ligand-receptor interaction activity profiles (denoting either activation: 1 or inactivation: 0) in each tumor sample, which is inferred using [**LIRICS**](https://pubmed.ncbi.nlm.nih.gov/34983745/) from the deconvolved expression – an interaction is considered as activated if the (deconvolved) expression of both its ligand and receptor genes is above their median expression values across the cohort samples, and inactivated otherwise;  2. The corresponding ICB response outcome for each patient. 

IRIS consists of two steps: Step I uses a Fisher’s test to identify differentially activated ligand-receptor interactions in the pre-treatment and non-responder post-treatment samples. These interactions are categorized as either resistant downregulated interactions (RDI) or resistant upregulated interactions (RUI) based on their differential activity state in the post-treatment vs. the pre-treatment state; that is, RDIs are downregulated in post-treatment resistant patients and vice versa for RUIs. Step II employs a hill climbing aggregative feature selection algorithm to choose the optimal set of RDIs or RUIs for classifying responders and non-responders in pre-treatment samples. The final output of IRIS is a selected set of RDIs and RUIs hypothesized to facilitate in ICB resistance, that can be used to predict ICB therapy response in a new ICB cohort.

See **Tutorial** & **Package** here: COMING SOON

#### Installation of IRIS
```r
COMING SOON
```

### SOCIAL: **S**ingle-cell transcrip**O**mics **C**ell-cell **I**nteraction **AL**gorithm
<img src="https://github.com/sahil-sahni/IRIS/blob/main/4.%20Figure/images/png/SOCIAL%20%5Bnc%20acc%5D.png" alt="grouping">

We developed an R method, **SOCIAL** (**S**ingle-cell transcript**O**mics **C**ell-cell **I**nteraction **AL**gorithm), to identify significant ligand-receptor interactions between two specific cell types, drawing upon insights from Kumar et al.'s (https://pubmed.ncbi.nlm.nih.gov/30404002/), Vento-Tormo et al.'s (https://pubmed.ncbi.nlm.nih.gov/30429548/), and our own [**LIRICS**](https://pubmed.ncbi.nlm.nih.gov/34983745/) framework. Our decision to create our own code stemmed from four primary motivations: 1. Leveraging the strengths of previous methods: By combining aspects of the three approaches, we aimed to maximize the accuracy and robustness of our ligand-receptor interaction predictions. 2. Implementing an R-based solution: While the first method lacked publicly accessible code and the second was in Python, we sought to create an R-based solution for accessibility and ease of use. 3. Incorporating our comprehensive database: Our ligand-receptor interaction database (LIRICS) provided rich and informative annotations, enhancing the depth of our analysis. 4. Accommodating variations in ligand-receptor interaction activity observed across patients.

SOCIAL comprises three main steps: 1. Querying the LIRICS database: Initially, we queried the LIRICS database to identify plausible ligand-receptor interactions; 2. Computing interaction scores: Next, we computed the ligand-receptor interaction score by multiplying the average expression levels of the ligand and receptor complexes for each interaction pair and cell type. 3. Permutation testing: Following that, we performed permutation tests (utilizing 100 iterations in our study) by randomly shuffling cell type labels. This allowed us to derive empirical p-values by calculating the fraction of permutation tests resulting in a higher interaction score than the foreground score determined in step 2. A lower p-value suggests a higher likelihood of the interaction occurring. 4. Optionally, ligand-receptor interactions can be further denoted as significantly activated if the average expression level of both the ligand and receptor genes is greater than the median across all samples.

See **Tutorial** & **Package** here: XXXX

### SPECIAL: **SP**atial c**E**ll-**C**ell **I**nteraction **AL**gorithm
<img src="https://github.com/sahil-sahni/IRIS/blob/main/4.%20Figure/images/png/SPECIAL%20%5Bnc%20acc%5D.png" alt="grouping">

To quantify the activity of cell-type-specific ligand-receptor interactions within each spatial transcriptomics slide, we further developed our in-house single-cell ligand-receptor inference tool called SOCIAL, into SPECIAL (**SP**atial c**E**ll-**C**ell **I**nteraction **AL**gorithm). This novel iteration is customized specifically for spatial transcriptomics with aligned single-cell transcriptomes, the direct output of [**CytoSPACE**](https://github.com/digitalcytometry/cytospace/tree/main). It consists of three major steps: Step I utilizes either a sliding window or k-means clustering approach on bulk (i.e. Visium 10X and Legacy) and SlideSeqV2 spatial transcriptomics, respectively, to divide spatial slides into “regions” of approximately 250 μm in diameter. Step II employs SOCIAL steps 1 through 3 to infer cell-type-specific interaction activity within each ~250 μm region. Step III, ligand-receptor interactions are further denoted as significantly activated if the average expression levels of both the ligand and receptor genes within the respective cell type is greater than the median across all regions. The final output of SPECIAL is a cell-type-specific ligand-receptor interaction activity profile across all regions in a spatial transcriptomics slide.

See **Tutorial** & **Package** here: COMING SOON!!

#### Installation of SOCIAL & SPECIAL
```r
COMING SOON
```

## Data availability 
The CODEFACS, LIRICS, SOCIAL, CytoSPACE, and SPECIAL data (and relevant data) generated in this study have been deposited to XX.

## Folders
0. **Data**: details regarding data files from zenodo
1. **IRIS**: code to reproduce results for **I**mmunotherapy **R**esistance *cell-cell* **I**nteraction **S**canner
2. **Benchmark**: code to reproduce results for previously established transcriptomics biomarkers' of ICB
3. **SOCIAL & SPECIAL**: code to reproduce results for **S**ingle-cell transcrip**O**mics **C**ell-cell **I**nteraction **AL**gorithm and **SP**atial c**E**ll-**C**ell **I**nteraction **AL**gorithm
4. **Figure**: code to reproduce figures and results
5. **Package(s)** **[coming soon]**: R packages for user implementation of **IRIS**, **SOCIAL**, and **SPECIAL** in their own work
    1. **IRIS** R package with vignettes
    2. **SOCIAL** R package (includes SPECIAL) with vignettes

## System requirements
IRIS, SOCIAL, and SPECIAL were developed on R (v4.4.1) using R packages: dplyr (v1.1.4), magrittr (v2.0.3), parallel (v4.4.1), pROC (v1.18.5), rBayesianOptimization (v1.2.1), tidyr (v1.3.1), abind (v1.4-5), Matrix (v1.7-0),  urr (v1.0.2), reshape2 (1.4.4), rslurm (v0.6.2), and stats (v4.4.1). All statistical analyses were done on R (v4.4.1).

## Citation
If using IRIS, SOCIAL, SPECIAL, or any of the pertaining results, please cite:

Sahni et al. "A machine learning model reveals expansive downregulation of ligand-receptor interactions enhancing lymphocyte infiltration in melanoma with acquired resistance to Immune Checkpoint Blockade. *N C* **X**, XXXX (XXXX). https://doi.org

## Contact
### Corresponding Author(s)
1. **Kun Wang** (kun.wang@nih.gov)
2. **Eytan Ruppin** (eytan.ruppin@nih.gov)

### Developer(s)
* IRIS was developed by Sahil Sahni ([sahil-sahni](https://github.com/sahil-sahni)) and Kun Wang ([kwangcb](https://github.com/kwangcb))
* SOCIAL was developed by Kun Wang ([kwangcb](https://github.com/kwangcb)), Sahil Sahni ([sahil-sahni](https://github.com/sahil-sahni)), and Sushant Patkar ([spatkar94](https://github.com/spatkar94))
* SPECIAL was developed by Sahil Sahni ([sahil-sahni](https://github.com/sahil-sahni)) and Kun Wang ([kwangcb](https://github.com/kwangcb))


## Acknowledgement(s)
IRIS, SOCIAL, and SPECIAL figures were created with BioRender.com. IRIS figure was developed by Di Wu, PhD.
