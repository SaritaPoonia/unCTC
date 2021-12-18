#------------------------------------------------------------------------
#-------------------- Get Differential genes  ---------------------------
#------------------------------------------------------------------------


#' Claculate differential gene between two type of groups.
#' @description Calculate the top ten most differentially elevated genes
#' in group 1 that have a significant p-value difference.
#' For differentiation, limma voom is utilised.
#'
#' @param data_mat Raw or Normalized data matrix in which genes in the row
#'  and cells in columns.
#' @param group Different groups in a vector of size equals to the
#' sample size of data_mat
#' @param data_type whether data is Normalised or Raw (without normalization)
#' @param Normalization_method Only used if expresion matrix is not
#' normalized. All normalization methods are explained in calcNormFactor
#' function of edgeR package
#' @param p_val Threshold p value, default is 0.05
#' @param lfc Threshold log fold change value, default is 0
#' @param up_gene_number select number of upregulated genes in each group
#'
#' @importFrom  edgeR cpm
#' @importFrom edgeR calcNormFactors
#' @importFrom edgeR DGEList
#' @import limma
#' @importFrom stats model.matrix
#' @importFrom stats coef
#'
#' @return Up_gene_mat return expression matrix of upregulated genes
#' in each group
#'
#' @examples
#' data = unCTC::Poonia_et_al._TPMData
#' groups = c(rep("TNBC",11),rep("NonTNBC",61))
#' output = GroupsDiffGenes(data_mat=data,
#'                           group=groups,
#'                           data_type="Raw",
#'                            Normalization_method="TMM")
#' @export GroupsDiffGenes
#'

GroupsDiffGenes = function(data_mat, #Expression Matrix
                           group,    #Vector of 2 groups
                           data_type = c("Normalised","Raw"),
                           Normalization_method = c("TMM","TMMwsp","RLE",
                                                    "upperquartile","none"),
                           p_val = 0.05,
                           lfc = 0,
                           up_gene_number = 10)

{
  #Create DGEList object
  d0 <- DGEList(data_mat)

  #Calculate normalization factors
  if (data_type== "Normalised")
    d0_calcnorm <- calcNormFactors(d0,method = "none")
  else
    d0_calcnorm <- calcNormFactors(d0,method = Normalization_method)

  #Specify the model to be fitted
  mm <- model.matrix(~0 + group)
  colnames(mm) = c("group1","group2")
  #Voom : A linear model is fitted to the log2 CPM for each gene
  voom_norm  =  voom(d0_calcnorm,mm,plot = FALSE,d0_calcnorm$samples$lib.size)

  #lmfits: fits a linear model using weighted least squares for each gene
  fit1 <- lmFit(voom_norm, mm)

  #makeContrasts :Specify which groups to compare:
  fit_contrast <- makeContrasts(group1 - group2,
                                levels = colnames(coef(fit1)))

  #contrasts.fit : Estimate contrast for each gene
  tmp_contrast <- contrasts.fit(fit1, fit_contrast)

  #eBayes: Empirical Bayes smoothing of standard errors
  error_smooth <- eBayes(tmp_contrast)

  #topTable: genes which are most differentially expressed
  top.table1 <- topTable(error_smooth, sort.by = "P", n = Inf)


  #Select genes which have adjusted P value < given threshold p value
  top.table2 = top.table1[top.table1$adj.P.Val < p_val,]

  #Select genes which have log fold change > given threshold lfc
  top.table3 = top.table2[top.table2$logFC > lfc,]
  table3_sort = top.table3[order(top.table3$logFC,decreasing = TRUE),]
  Upregulated_genes = head(rownames(table3_sort),up_gene_number)


  #Matrix
  Up_gene_mat = voom_norm$E[Upregulated_genes,]

  ####
  return(Up_gene_mat)
}




#------------------------------------------------------------------------
#-----------  Differential genes between clusters  ----------------------
#------------------------------------------------------------------------


#' Calculate differential genes between cells/clusters.
#' @description Provide differential genes between given groups.
#'
#' @param data_list List of expression matricies
#' @param min_Gene cell filter, filter out those cells which do not
#' express at least min_Gene genes
#' @param min_Sample gene filter, filter out genes which are not expressed
#' in at least min_Sample cells
#' @param DDLK_Clusters Output of DDLK_Clust.R method
#' @param data_id List of names/ids of expression matrix
#' @param data_type list of expression data passed in data.
#' Valid inputs are either raw or normalised.
#' @param Genesets list of genesets/pathways used to calculate
#' PathwayEnrichmentScore.
#' @param DifferentiateBy Any column name from
#' DDLK_Clusters$PathwayDDLK_clust, default is "Clusters"
#' @param p_val Threshold p value, default is 0.05
#' @param lfc Threshold log fold change value, default is 0
#' @param up_gene_number select number of upregulated genes in each group
#'
#' @return Diff_mat list of three objects
#' 1. DifferentialMatrix = Data matrix with top selected up_gene_number
#' in each group
#' 2. Diffup_genes = list of differential gene in each group
#' with selected up_gene_number
#' 3. annotations = Cell wise annotation of DifferentialMatrix

#' @examples
#' data1 = unCTC::Poonia_et_al._TPMData
#' data2 = unCTC::Ding_et_al._WBC1_TPMData
#' Data_list = list(data1,data2)
#' Data_Id = list("data1","data2")
#' Genesets = unCTC::c2.all.v7.2.symbols
#' Pathway_score = PathwayEnrichmentScore(data_list=Data_list,
#'                                         data_id= Data_Id,
#'                                          min_Sample = 5,
#'                                          min_Gene = 1500,
#'                                         Genesets=Genesets,
#'                                         min.size=70,
#'                                         max.size=100)
#'
#' DDLK_Clusters = DDLK_Clust(PathwayScore = Pathway_score$Pathway_score,
#'                            PathwayMetaData=Pathway_score$Pathway_metadata,
#'                             n=3,
#'                             out.dir = getwd())
#'
#' Output = Differential_genes(data_list=Data_list,
#'                               min_Sample = 5,
#'                               min_Gene = 1500,
#'                               DDLK_Clusters=DDLK_Clusters,
#'                               data_id = Data_Id,
#'                               Genesets=Genesets,
#'                               data_type="Normalised",
#'                               DifferentiateBy = "Clusters")
#
#' @export Differential_genes
Differential_genes = function(
                       data_list=list(), # list of expression matricies
                       min_Sample = 5,
                       min_Gene = 1500,
                       DDLK_Clusters, #Output of DDLK_Clust.R method
                       data_id, #list of expression matrices
                       #name in same order
                       Genesets,
                       data_type = c("Normalised","Raw"),
                       DifferentiateBy = "Clusters",
                       p_val = 0.05,
                       lfc = 0,
                       up_gene_number = 10)
{
  # Create SingleCellObject
  if(data_type=="Raw"){
  sce_obj = CreateSingleCellObject(data_list,min_sample = min_Sample,
                                   min_gene=min_Gene)
  sce_data = sce_obj$sce_norm

  #Select normalised expression matrix
  data_norm = as.matrix(sce_data@assays@data$Normalized_Data)

  # Assign colnames to expression matrix
  colnames(data_norm) = sce_data@colData@listData$Samples

  # Assign rownames to expression matrix
  rownames(data_norm) = rowData(sce_data)[,1]
  }else
  {
    sce_obj = CreateSingleCellObject(data_list,min_sample = min_Sample,
                                     min_gene=min_Gene)
    sce_data = sce_obj$sce_raw
    #Select normalised expression matrix
    data_norm = as.matrix(sce_data@assays@data$Raw_filtered_data)

    # Assign colnames to expression matrix
    colnames(data_norm) = sce_data@colData@listData$FilterDataSample

    # Assign rownames to expression matrix
    rownames(data_norm) =  rowData(sce_data)[,1]
  }
  #Clusters metadata
  metadata = DDLK_Clusters$PathwayDDLK_clust


  #Genes which are used to calculate PathwayEnrichmentScore
  GeneList = list()
  for (i in seq_along(Genesets)) {
    gene2 = intersect((rownames(data_norm)),
                      (unlist(Genesets[i])))
    GeneList[[i]] = gene2
  }
  names(GeneList) = names(Genesets)
  GeneList_genes = unique(unlist(GeneList))

  #Set column position according to groups/clusters
  data_norm1 = data_norm[GeneList_genes,rownames(metadata)]


  up_genes = list()
  up_geneList = list()
  for (i in seq_along(unique(metadata[,DifferentiateBy])))
  {
    a = rownames(metadata[metadata[,DifferentiateBy]==unique
                          (metadata[,DifferentiateBy])[i],])
    b=  rownames(metadata[metadata[,DifferentiateBy]!=unique
                          (metadata[,DifferentiateBy])[i],])
    data_a = data_norm1[,a]
    data_b = data_norm1[,b]
    data_ab = cbind(data_a,data_b)
    group1 = c(rep("group1",length(a)))
    group2 = c(rep("group2",length(b)))
    group = c(group1,group2)
    mat = GroupsDiffGenes(data_ab,group,
                          data_type="Normalised",
                          Normalization_method="none",
                          p_val = p_val,
                          lfc = lfc,
                          up_gene_number = up_gene_number
                          )
    up_genes[[i]] = mat[,rownames(metadata)]
    up_geneList[[i]] = rownames(mat)
    #print(rownames(mat))


    ##
  }
  Differential_mat =do.call("rbind", up_genes)
  Differential_mat_1 = Differential_mat[,rownames(metadata)]

  Diff_mat = list(DiffMat=Differential_mat_1,
                  Diffup_genes = up_geneList,
                  annotations = metadata
                  )
  return(Diff_mat)
}


