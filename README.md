
<!-- README.md is generated from README.Rmd. Please edit that file -->

# unCTC: Characterising single circulating tumor cell transcriptomes

## Introduction

To identify Circulating Tumor Cells (CTCs) from White Blood Cells, unCTC
uses pathway-based unsupervised clustering of single-cell RNA-Seq data
(WBCs). It accepts a list of raw Countdata/TPM expression single-cell
RNA-Seq matrices. In the expression matrix, genes must be arranged in
rows and must put cells in columns. unCTC combines all matrices based on
commongenes, filters out low-expression genes and cells, removes the
batch effect, and normalizes the integrated matrix. For unsupervised
clustering of Circulating Tumor Cells (CTCs) and White Blood Cells
(WBCs), a normalized matrix is transformed to pathway space, and deep
dictionary learning utilizing k means clustering is implemented. unCTC
also calculates Copy Amount Variations (CNVs) in chromosomes across
clusters, revealing the number of CNVs and the p/q arm variation
location. With the use of Stouffer’s Z-score, unCTC allows for the
determination of various canonical markers identifying
malignant/epithelial/immune origins (Stouffer et al., 1949). The
expression of other canonical markers confirms the lineage of
Circulating Tumour Cells (CTCs).

## Installation

You can install the released version of unCTC from ….

``` r

#First, you need to install the devtools package. 
install.packages("devtools")

#Load the devtools package.
library(devtools)

#Install the unCTC package
install_github("SaritaPoonia/unCTC")

#Load unCTC
library(unCTC)
```

## Software Requirements

  - R (tested in R version 4.0.3 (2020-10-10))  
  - R libraries required include : PCAtools, viridis, magrittr, stringr,
    limma, rworldmap, data.table, GSVA, pheatmap, qusage, umap, ggplot2,
    SingleCellExperiment, Linnorm, readxl, readtext, infercnv, ggpubr,
    reticulate, GenomicRanges, readtext, S4Vectors, cowplot, edgeR,
    D3GB, dplyr, SummarizedExperiment, IRanges, Matrix, utils, stats
  - Python 3 with installed modules: h5py, scipy, pandas, numpy,
    sklearn, matplotlib.pyplot, scipy.io

## Data Requirements

unCTC requires:  
\* List of expression data matrices.  
\* List of expression data matrices name in the same order.  
\* Gene list: List of specific genes  
\* Genesets: list of pathways \* A gene/chromosome positions file

## Usage and workflow

Two data matrices Poonia\_et\_al.\_TPMData and
Ding\_et\_al.\_WBC1\_TPMData, and their meta data are given with this
package.

### Load Data and meta data from package

``` r
Poonia_et_al._TPMData = unCTC::Poonia_et_al._TPMData
Ding_et_al._WBC1_TPMData = unCTC::Ding_et_al._WBC1_TPMData
Poonia_et_al._metaData = unCTC::Poonia_et_al._metaData
Ding_et_al._WBC1_metaData = unCTC::Ding_et_al._WBC1_metaData
```

Here we are using two other dataset Ding\_et\_al.\_WBC2\_TPMData and
Ebright\_et\_al.\_TPMData.

  - Download all datasets from below
    link:  
    <https://drive.google.com/file/d/1GGnDZFw40ULjQ7_RGQGW1m-hEkcGKcBD/view?usp=sharing>  
  - download unCTC\_data.zip folder.
  - Unzip unCTC\_data.zip..
  - It contains all the data used in the creating and validating unCTC
    package.
  - Here we are using 4 expression data and 4 corresponding meta data
    files. Out of which Poonia\_et\_al.\_TPMData,
    Poonia\_et\_al.\_metaData, Ding\_et\_al.\_WBC1\_TPMData and
    Ding\_et\_al.\_WBC1\_metaData are given with the package. So we can
    direct load from unCTC package.
  - Load Ding\_et\_al.\_WBC2\_TPMData, Ding\_et\_al.\_WBC2\_metaData,
    Ebright\_et\_al.\_TPMData and Ebright\_et\_al.\_metaData from
    downloaded folder.

<!-- end list -->

``` r
load("/path/of/downloaded/folder/Ebright_et_al._TPMData.RData")
Ebright_et_al._TPMData = Ebright_et_al._TPMData

load("/path/of/downloaded/folder/Ding_et_al._WBC2_TPMData.RData")
Ding_et_al._WBC2_TPMData = Ding_et_al._WBC2_TPMData

load("/path/of/downloaded/folder/Ebright_et_al._metaData.RData")
Ebright_et_al._metaData = Ebright_et_al._metaData

load("/path/of/downloaded/folder/Ding_et_al._WBC2_metaData.RData")
Ding_et_al._WBC2_metaData = Ding_et_al._WBC2_metaData
```

### Load geneset

This package includes one geneset, which is taken from molecular
signature database.

``` r
c2.all.v7.2.symbols = unCTC::c2.all.v7.2.symbols
```

``` r
#Create Expression data list
dataList = list(Poonia_et_al._TPMData,Ebright_et_al._TPMData,
                Ding_et_al._WBC1_TPMData,Ding_et_al._WBC2_TPMData)

#Create Data Id's list
dataId = list("Poonia_et_al._TPMData","Ebright_et_al._TPMData",
              "Ding_et_al._WBC1_TPMData","Ding_et_al._WBC2_TPMData")

#Create Meta data list
MetaData = list(Poonia_et_al._metaData, Ebright_et_al._metaData, 
                Ding_et_al._WBC1_metaData, Ding_et_al._WBC2_metaData )

#Genesets given with this package
genesets = c2.all.v7.2.symbols
```

## Calculate pathway enrichment score

PathwayEnrichmentScore uses the following steps:

Integrate data passed in the list based on common genes.

  - Filter out low expressed genes and cells.  
  - Use Linnorm.Norm() for normalization and batch effect correction.  
  - Create singleCellObject instance.  
  - Calculate pathway enrichment score.  
  - Calculate metadata.  
  - Return pathway enrichment score and pathway metadata.

PathwayEnrichmentScore requires following inputs:  
\* data\_list: List of expression data matrices  
\* data\_id: List of expression data matrices name in the same order.  
\* Genesets: List of pathways \* min.size: Minimum size of genes in
pathways/Genesets,Default is 10  
\* max.size: Maximum size of genes in pathways/Genesets,Default is 500  
\* min\_Sample: filter out genes which are not expressedin at least
min\_Sample cells, Default is 5.  
\* min\_Gene: Filter out those cells which do not express at least
min\_Gene genes, Default is 1500. \* Parallel\_threads : Number of
threads in parallel to execute process

``` r
#Call PathwayEnrichmentScore
Pathway_score = unCTC::PathwayEnrichmentScore(data_list =dataList,
                                        data_id = dataId,
                                        Genesets = genesets,
                                        min.size=10,
                                        max.size = 500,
                                        min_Sample = 5,
                                        min_Gene = 1500,
                                       Parallel_threads=8L)
```

### Calculate the optimal number of clusters for pathway enrichment score matrix

For the above pathway enrichment score matrix, we calculate the number
of clusters using the Elbow method.

``` r
library(factoextra)
library(NbClust)
fviz_nbclust(Pathway_score$Pathway_score, kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Elbow method")
```

<img src="man/figures/K_clustNo.png" width="80%" style="display: block; margin: auto;" />

## DDLK Clusteing

DDLk\_Clust need the following inputs

  - PathwayScore: Pathway enrichment score matrix. We get this from
    PathwayEnrichmentScore\_output$Pathway\_score  
  - PathwayMetaData: Pathways metadata. we get this from
    PathwayEnrichmentScore\_output$Pathway\_metadata  
  - n: Number of clusters for K-means clustering. We can determine the
    optimal number of clusters for pathway enrichment score matrix using
    the Elbow method, Average silhouette method, Gap statistic method.  
  - out.dir: Directory to save enrichment score. This input is
    mandatory. Default is your current directory.  
  - MetaData: list of metadata for all expression matrices. This input
    is optional. Only give if you have the same columns in all metadata
    matrices.

<!-- end list -->

``` r
# DDLK_Clust use python environment so set python environment before running DDLK_Clust()
# Set Python3 Path: Example: Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python3")

Sys.setenv(RETICULATE_PYTHON = "/path/to/python3") 
library(reticulate)

#Retrive information about the version of python being used by reticulate
reticulate::py_config()
#If version is different from the the given path then restart session and
#give path again can change path

DDLK_Clusters = unCTC::DDLK_Clust(PathwayScore = Pathway_score$Pathway_score,
                           PathwayMetaData = Pathway_score$Pathway_metadata,
                           n = 4,
                           out.dir = getwd(),
                           MetaData = MetaData
                           )
```

## unCTC plots:

unCTC\_plots Plots principal components of pathway enrichment score.

  - Return list of three plots:group\_by\_Class,
    group\_by\_Cluster,p1\_p5\_Pairssplot

Required input for unCTC\_plots method is:

  - PathwayScore: Pathway enrichment score matrix. We get this from
    PathwayEnrichmentScore\_output$Pathway\_score.  
  - PathwayMetaData: Pathways metadata. we get this from
    PathwayEnrichmentScore\_output$Pathway\_metadata.  
  - colorby: Any column name from PathwayMetaData, default is
    “Data\_id”.  
  - Color\_cluster: Any column name from PathwayMetaData, default is
    “Clusters”.  
  - pairsplotLegend: Legend Position in pairsplot. Choose one from
    “left”,“right” and “none”

<!-- end list -->

``` r
plots = unCTC::unCTC_plots(Pathway_score = DDLK_Clusters$Pathway_score,
                    Pathway_metadata = DDLK_Clusters$PathwayDDLK_clust,
                    colorby = "Data_id",
                    Color_cluster = "Clusters",
                    pairsplotLegend = "none")
```

### Color by datasets id

  - Color by any column from metadata, default is
“Data\_id”.

<!-- end list -->

``` r
plots$group_by_Class
```

<img src="man/figures/Group-by-class-1.png" width="70%" style="display: block; margin: auto;" />

### Color by clusters

``` r
plots$group_by_Cluster
```

<img src="man/figures/Group-by-Cluster-1.png" width="70%" style="display: block; margin: auto;" />

stacked bar plot shows the count of CTCs (red) and WBCs (blue), while
the x-axis shows clusters.

``` r
library(ggplot2)
ggplot(DDLK_Clusters$PathwayDDLK_clust, aes(x=Clusters, fill = Cell_type))+
       theme_classic()+
       geom_bar(stat="count")+
       scale_color_manual()+
       scale_fill_manual(
       values = c("deepskyblue3","darkred",
               "darkgreen","dark turquoise"))
```

<img src="man/figures/Counts-per-cluster-1.png" width="80%" style="display: block; margin: auto;" />

## Differential genes

Provide differential genes between given groups.

Differential\_genes require following inputs:

  - data\_list: List of expression data matrices  
  - min\_Sample: filter out genes which are not expressedin at least
    min\_Sample cells, Default is 5  
  - min\_Gene: Filter out those cells which do not express at least
    min\_Gene genes, Default is 1500.  
  - DDLK\_Clusters: output of DDLK\_Clust()  
  - data\_id: List of expression data matrices name in same order.  
  - data\_type: choose from given vector,c(“Normalised”,“Raw”).
  - DifferentiateBy: Default is Clusters. We can choose any column name
    from DDLK\_Clusters$PathwayDDLK\_clust  
  - p\_val = Threshold p-value, Default is 0.05
  - lfc = Threshold log fold change value, Default is 0
  - up\_gene\_number = Select number of upregulated genes you want to
    obtain. Default is 10. If you get an error when computing the number
    of upregulated genes, relax the p val parameter.

<!-- end list -->

``` r
Diff_matrix = unCTC::Differential_genes(data_list=dataList,
                                 min_Sample = 5,
                                 min_Gene = 1500,
                                 DDLK_Clusters,
                                 Genesets = genesets,
                                 data_id=dataId,
                                 data_type = "Normalised",
                                 DifferentiateBy = "Clusters",
                                 up_gene_number = 100)
```

### Heatmap showing the top 100 upregulated genes in the 4 clusters.

``` r
library(pheatmap)
library(viridis)
annotation = Diff_matrix$annotations
annotation$Data_id <- NULL
annotation$GroupID <- NULL
annotation$Cell_type <- NULL
ann = annotation[,c("HormoneStatus","Class","Clusters")] 

pheatmap(t(scale(t(Diff_matrix$DiffMat))),cluster_cols = FALSE,
         show_colnames = FALSE,cluster_rows = FALSE, show_rownames = FALSE,
         color = viridis(1000),annotation = ann)
```

<img src="man/figures/Differential-gene-matrix-pheatmap-1.png" width="80%" style="display: block; margin: auto;" />

## Differential Pathways

Provide differential pathways between given groups.

Differential\_pathways require following inputs:

  - Pathway\_score: Output of PathwayEnrichmentScore.R method  
  - DDLK\_Clusters: output of DDLK\_Clust.R method  
  - DifferentiateBy: Default is Clusters. We can choose any column name
    from DDLK\_Clusters$PathwayDDLK\_clust  
  - p\_val = Threshold p-value, Default is 0.05
  - lfc = Threshold log fold change value, Default is 0
  - up\_pathways\_number = Select number of upregulated pathways you
    want to obtain. Default is 10. If you get an error when computing
    the number of upregulated pathways, relax the p val parameter.

<!-- end list -->

``` r
Diff_path = unCTC::Differential_pathways(Pathway_score,
                      DDLK_Clusters = DDLK_Clusters,
                      DifferentiateBy = "Clusters",
                      up_pathways_number = 10
                      )
```

### Top 10 upregulated pathways in each cluster

``` r
annotation = Diff_path$annotations
annotation$Data_id <- NULL
annotation$GroupID <- NULL
annotation$Cell_type <- NULL
ann = annotation[,c("HormoneStatus","Class","Clusters")] 

pheatmap(t(scale(t(Diff_path$DiffMatpathway))),cluster_cols = FALSE,
         show_colnames = FALSE,cluster_rows = FALSE, show_rownames = TRUE,
         color = viridis(1000),annotation = ann)
```

<img src="man/figures/Top10-Diff-Pathways-1.png" width="100%" style="display: block; margin: auto;" />

## Calcuate Stouffers score

Stouffer\_score method uses the following steps:

  - Calculate Z score of log normalized data  
  - Calculate Stouffer score for the specific gene list  
  - Assign Stouffer score to different groups of data and results are
    shown in boxplot

The followings input are required to calculate Stouffer’s score:

  - data\_list: List of expression data matrices.  
  - min\_Sample: filter out genes which are not expressedin at least
    min\_Sample cells, Default is 5.  
  - min\_Gene: Filter out those cells which do not express at least
    min\_Gene genes, Default is 1500.
  - data\_id: List of expression data matrices name in the same order.  
  - gene\_list: Specific genes in vector.  
  - MetaData: Optional, List of metadata of expression matrices. If
    given, then columns of all metadata in the list must be the same.  
  - metaColPos: default = 1,Only require if metadata is given. Position
    of the column in metadata for which we want to calculate and compare
    Stouffer score for given gene list.  
  - metaColName: default = “Class”, Only require if metadata is given.  
    Name of the column in metadata for which we want to calculate and
    compare Stouffer score for given gene list.

With this package, we have given two types of gene list:

### Load genelists

``` r
Breast_elevated_genes = unCTC::Breast_elevated_genes
Immune_signature_genes = unCTC::Blood_specific_gene
```

### Stouffer’s Score:

``` r
#Calculate Stouffer's score for Blood gene
S_WBC = unCTC::Stouffer_score(data_list = dataList,
                         min_Sample = 5,
                         min_Gene = 1500,
                         gene_list =Immune_signature_genes,
                         data_id = dataId,
                         Groupby = "Clusters",
                         DDLKCluster_data = DDLK_Clusters)


#Calculate Stouffer's score for Breast elevated genes
S_Breast = unCTC::Stouffer_score(data_list = dataList,
                           min_Sample = 5,
                           min_Gene = 1500,
                           gene_list = Breast_elevated_genes,
                           data_id = dataId,
                           Groupby = "Clusters",
                           DDLKCluster_data = DDLK_Clusters)
```

## Stouffer’s Score Plot

For better colour visualization we are using following color key:

``` r
library(ggplot2)
library(ggpubr)
ColorKey = c("darkred","deepskyblue3","darkolivegreen4",
             "dark turquoise","pale violet red",
             "steelblue","forestgreen","gray2",
             "gray50","hotpink","lightslateblue",
             "tan4","yellow3","sienna4","orchid4")
```

### For Immune genes:

``` r
ggplot(S_WBC$Stouffer_score,aes(x=Clusters,y= Stouffer_score,fill=Clusters))+
geom_boxplot(outlier.shape = NA) + theme_classic() +
scale_fill_manual(values = ColorKey) +
ggtitle("Immune gene signature")+
stat_compare_means(comparisons = S_WBC$comparisons,
label = "p.signif", method = "t.test",ref.group = ".all.")
```

<img src="man/figures/Immune-genes-Stouffer-1.png" width="80%" style="display: block; margin: auto;" />

### For Breast elevated genes:

``` r
ggplot(S_Breast$Stouffer_score,aes(x=Clusters,y= Stouffer_score,fill=Clusters))+
geom_boxplot(outlier.shape = NA) + theme_classic() +
scale_fill_manual(values = ColorKey)+
ggtitle("Breast elevated gene signature")+
stat_compare_means(comparisons = S_Breast$comparisons,
label = "p.signif", method = "t.test",ref.group = ".all.")
```

<img src="man/figures/Breast-genes-Stouffer-1.png" width="80%" style="display: block; margin: auto;" />

## Copy Number Variation Analysis:

inferCNV R package is used for analysing copy number variation for raw
Count/TPM data. Along with all analysis of inferCNV,
unCTC::CNV\_alterations calculate addition and deletion position in p
and q arms in test/cancerous/ diseased data as compared to
reference/normal/healthy. To calculate p and q arm location from
inferCNV events, we used GRCh37 cytoband information.

CNV\_alterations require the following inputs:

  - data\_list: List of expression data matrices.  
  - data\_id: List of expression data matrices name in the same order.  
  - path= path to save output files.
  - min\_Sample: filter out genes which are not expressedin at least
    min\_Sample cells, Default is 5.  
  - min\_Gene: Filter out those cells which do not express at least
    min\_Gene genes, Default is 1500.
  - GenePositionFile: Gene/chromosomal order file. This package includes
    genecode hg19 gene order file. Either you can use this or download
    it from outer sources.  
  - threads\_no: Thread number for parallel processes, default is 8.  
  - MetaData: Optional, List of metadata of expression matrices. If
    given, then columns of all metadata in the list must be the same.  
  - Groupby: Any column name from MetaData, which we want to use as an
    annotation file. Only applicable if MetaData is given.  
  - Reference\_name: Any cell type from the data\_id list or any cell
    type from the column assign to Groupby.  
  - obs.title: Title of test/observation matrix. Default is
    “Observations”.  
  - ref.title: Title of reference matrix. Default is “References”.  
  - out.Filename: Output filename-prefix. Default is “inferCNV”.

### Load gene order file

``` r
gencode_v19_gene_pos =unCTC::gencode_v19_gene_pos
```

### inferCNV between clusters by taking cluster\_0 as reference

``` r
CNV_Alterations = unCTC::CNV_alterations(
                     data_list= dataList,
                     data_id= dataId,
                     min_Sample = 5,
                     min_Gene = 1500,
                     path= "/path/to store/output",      
                     GenePositionFile= gencode_v19_gene_pos,  
                     threads_no=8, 
                     MetaData = list(DDLK_Clusters$PathwayDDLK_clust),
                     Groupby = "Clusters",
                     Reference_name = "Cluster_0",
                     obs.title ="Observations", 
                     ref.title = "References",
                     out.Filename = "inferCNV_Cluster" 
                     )
```

### CNV States

``` r
print("HMM state (1 = 2 copies loss, 2 = 1 copy loss, 3 = neutral,4 = 1 copy gain, 5 = 2 copies gain, 6 = 3+ copies gain),")

head(CNV_Alterations)
```

<img src="man/figures/p_q_arm_clusters.png" width="100%" style="display: block; margin: auto;" />
<img src="man/figures/infercnv_clusters.png" width="100%" style="display: block; margin: auto;" />

### inferCNV by taking WBC cells as reference and CTCs as observation

``` r
CNV_Alterations_dataId = unCTC::CNV_alterations(
                     data_list= dataList,
                     data_id= dataId,
                     min_Sample = 5,
                     min_Gene = 1500,
                     path= "/path/to store/output",      
                     GenePositionFile= gencode_v19_gene_pos,  
                     threads_no=8, 
                     MetaData = list(DDLK_Clusters$PathwayDDLK_clust),
                     Groupby = "Data_id",
                     Reference_name = c("Ding_et_al._WBC1_TPMData",
                      "Ding_et_al._WBC2_TPMData"), # WBC data as reference
                     obs.title ="Observations", 
                     ref.title = "References",
                     out.Filename = "inferCNV" 
                     )
```

### CNV States

``` r
print("HMM state (1 = 2 copies loss, 2 = 1 copy loss, 3 = neutral,4 = 1 copy gain, 5 = 2 copies gain, 6 = 3+ copies gain),")

head(CNV_Alterations_dataId)
```

<img src="man/figures/p_q_arm_Dataid.png" width="100%" style="display: block; margin: auto;" />
<img src="man/figures/infercnv_DataId.png" width="100%" style="display: block; margin: auto;" />

## Gene\_Violin\_plots

Give violin plot for a given Canonical marker expression.

Gene\_Violin\_plots require input:  
\* data\_list: List of expression data matrices  
\* data\_id: List of expression data matrices name in the same order.  
\* min\_Sample: filter out genes which are not expressedin at least
min\_Sample cells, Default is 5.  
\* min\_Gene: Filter out those cells which do not express at least
min\_Gene genes, Default is 1500. \* gene\_symbol: Specific gene for
which we want to see expression.  
\* MetaData: Optional, list of metadata of expression matrices. If given
then columns of all metadata in the list must be identical.  
\* Groupby: Any column name from MetaData, which we want to use to see
differential expression of the gene. Default is “data\_id”.

``` r
# Gene Violin plot
PTPRC = Gene_Violin_plots(data_list =dataList,
                  data_id = dataId,
                  min_Sample = 5,
                  min_Gene = 1500,
                  gene_symbol = "PTPRC",
                  DDLKCluster_data = DDLK_Clusters$PathwayDDLK_clust,
                  Groupby = "Clusters")
NKG7 = Gene_Violin_plots(data_list =dataList,
                  data_id = dataId,
                  min_Sample = 5,
                  min_Gene = 1500,
                  gene_symbol = "NKG7",
                  DDLKCluster_data = DDLK_Clusters$PathwayDDLK_clust,
                  Groupby = "Clusters")

EPCAM = Gene_Violin_plots(data_list =dataList,
                  data_id = dataId,
                  min_Sample = 5,
                  min_Gene = 1500,
                  gene_symbol = "EPCAM",
                  DDLKCluster_data = DDLK_Clusters$PathwayDDLK_clust,
                  Groupby = "Clusters")
KRT18 = Gene_Violin_plots(data_list =dataList,
                  data_id = dataId,
                  min_Sample = 5,
                  min_Gene = 1500,
                  gene_symbol = "KRT18",
                  DDLKCluster_data = DDLK_Clusters$PathwayDDLK_clust,
                  Groupby = "Clusters")
```

### Canonical marker expression

``` r
library(cowplot)
plot_grid(PTPRC,NKG7,EPCAM,KRT18)
```

<img src="man/figures/Canonical-markers-expression-1.png" width="100%" style="display: block; margin: auto;" />
