<div align = "center">
  <h1>IRIS R Package</h1>
</div>
  
<img src="https://github.com/sahil-sahni/IRIS/blob/main/4.%20Figure/images/png/IRIS%20figure%201%20Final%20Version%20%5Bnc%20acc%5D.png" alt="grouping">

## IRIS: **I**mmunotherapy **R**esistance *cell-cell* **I**nteraction **S**canner
#### Last Updated: 08/11/24

We developed **I**mmunotherapy **R**esistance cell-cell **I**nteraction **S**canner (IRIS), a computational method specifically designed to identify immune checkpoint blockade (ICB) resistance relevant ligand-receptor interactions in the tumor microenvironment (TME), given a patients cohort including tumor bulk expression data and ICB treatment response data. The gene expression data is deconvolved using [**CODEFACS**](https://pubmed.ncbi.nlm.nih.gov/34983745/) such that the input to IRIS in a given patients cohort is comprised of two components: 1. Literature-curated cell-type-specific ligand-receptor interaction activity profiles (denoting either activation: 1 or inactivation: 0) in each tumor sample, which is inferred using [**LIRICS**](https://pubmed.ncbi.nlm.nih.gov/34983745/) from the deconvolved expression – an interaction is considered as activated if the (deconvolved) expression of both its ligand and receptor genes is above their median expression values across the cohort samples, and inactivated otherwise;  2. The corresponding ICB response outcome for each patient. 

IRIS consists of two steps: Step I uses a Fisher’s test to identify differentially activated ligand-receptor interactions in the pre-treatment and non-responder post-treatment samples. These interactions are categorized as either resistant downregulated interactions (RDI) or resistant upregulated interactions (RUI) based on their differential activity state in the post-treatment vs. the pre-treatment state; that is, RDIs are downregulated in post-treatment resistant patients and vice versa for RUIs. Step II employs a hill climbing aggregative feature selection algorithm to choose the optimal set of RDIs or RUIs for classifying responders and non-responders in pre-treatment samples. The final output of IRIS is a selected set of RDIs and RUIs hypothesized to facilitate in ICB resistance, that can be used to predict ICB therapy response in a new ICB cohort.

See **Tutorial** & **Package** here: XXXX

#### Installation of IRIS
```r
require(devtools)
devtools::install_github("sahil-sahni/5. Package(s)/IRIS")
```
