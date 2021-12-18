#------------------------------------------------------------------------
#-------------------- Get Differential pathways  ---------------------------
#------------------------------------------------------------------------


#' Claculate differential pathways between two type of groups.
#' @description Calculate top most differenciated upregulated pathways
#' with significant p-value between in group1, given in the form of vector.
#' Limma is used for differentiaition
#'
#' @param Pathway_mat Pathway enrichment score matrix in which cells in the
#' columns and pathways in the rows.
#' @param group Different groups in a vector of size equals to the
#' sample size of Pathway_mat
#' @param p_val Threshold p value, default is 0.05
#' @param lfc Threshold log fold change value, default is 0
#' @param up_pathways_number select number of upregulated pathways
#' in each group

#'
#' @importFrom  edgeR cpm
#' @importFrom edgeR calcNormFactors
#' @importFrom edgeR DGEList
#' @import limma
#' @importFrom stats model.matrix
#' @importFrom stats coef
#'
#' @return Up_pathway_mat return expression matrix of 10 and 100 most
#' upregulated pathways respectively
#'
#' @examples
#' data = unCTC::Poonia_et_al._TPMData
#' genesets = unCTC::c2.all.v7.2.symbols
#' Pathways_mat = PathwayEnrichmentScore(data_list=list(data),
#'                                        data_id = list("data"),
#'                                        Genesets = genesets,
#'                                        min.size = 50,
#'                                        max.size = 100)
#' groups = c(rep("TNBC",11),rep("NonTNBC",61))
#' output = GroupsDiffPathways(Pathway_mat=Pathways_mat$Pathway_score,
#'                           group=groups)
#'
#' @export GroupsDiffPathways
GroupsDiffPathways = function(
                            Pathway_mat, #Pathway enrichment score Matrix
                            group, #Vector of 2 groups, size = ncol(data_mat)
                            p_val = 0.05,
                            lfc = 0,
                            up_pathways_number = 10)
{
   #Specify the model to be fitted
  mm <- model.matrix(~0 + group)
  colnames(mm) = c("group1","group2")

  #lmfits: fits a linear model using weighted least squares for each pathways
  fit1 <- lmFit(Pathway_mat, mm)

  #makeContrasts :Specify which groups to compare:
  fit_contrast <- makeContrasts(group1 - group2,
                                levels = colnames(coef(fit1)))

  #contrasts.fit : Estimate contrast for each pathway
  tmp_contrast <- contrasts.fit(fit1, fit_contrast)

  #eBayes: Empirical Bayes smoothing of standard errors
  error_smooth <- eBayes(tmp_contrast)

  #topTable: pathwayss which are most differentially expressed
  top.table1 <- topTable(error_smooth, sort.by = "P", n = Inf)


  #Select pathwayss which have adjusted P value < given threshold p value
  top.table2 = top.table1[top.table1$adj.P.Val< p_val,]

  #Select pathways which have log fold change > given threshold lfc
  top.table3 = top.table2[top.table2$logFC > lfc,]
  table3_sort = top.table3[order(top.table3$logFC,decreasing = TRUE),]

  #Select top upregulated pathways
  Upregulated_pathways = head(rownames(table3_sort),up_pathways_number)


  #Matrix
  Up_pathway_mat = Pathway_mat[Upregulated_pathways,]
  return(Up_pathway_mat)
}


#------------------------------------------------------------------------
#-----------  Differential Pathways between clusters  ----------------------
#------------------------------------------------------------------------


#' Calculate differential pathways between cells/clusters.
#' @description Provide differential pathways between given groups.
#'
#' @param Pathway_score Output of PathwayEnrichmentScore.R method
#' @param DDLK_Clusters Output of DDLK_Clust.R method
#' @param DifferentiateBy Any column name from
#' DDLK_Clusters$PathwayDDLK_clust, default is "Clusters"
#' @param p_val Threshold p value, default is 0.05
#' @param lfc Threshold log fold change value, default is 0
#' @param up_pathways_number select number of upregulated pathways
#' in each group
#'
#' @return Diff_Pathways list of three objects
#' 1. DiffMatpathway = Pathway matrix with top most
#' up_pathways_number differential pathways in each group
#' 2. Diffup_pathways = top most up_pathways_number differential
#'  pathways in each group
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
#' Output = Differential_pathways(Pathway_score=Pathway_score,
#'                               DDLK_Clusters=DDLK_Clusters,
#'                               DifferentiateBy = "Clusters")
#
#' @export Differential_pathways
Differential_pathways = function(
                       Pathway_score,#output of PathwayEnrichmentScore.R method
                       DDLK_Clusters, #Output of DDLK_Clust.R method
                       DifferentiateBy = "Clusters",
                       p_val = 0.05,
                       lfc = 0,
                       up_pathways_number = 10)
{

  #Clusters metadata
  metadata = DDLK_Clusters$PathwayDDLK_clust


  #PathwayEnrichmentScore
  Pathway_score = Pathway_score$Pathway_score

  #Set column position according to groups/clusters
  Pathway_score1 = Pathway_score[,rownames(metadata)]


  up_pathways = list()
  up_pathwaysList = list()
  for (i in seq_along(unique(metadata[,DifferentiateBy])))
  {
    a = rownames(metadata[metadata[,DifferentiateBy]==unique
                          (metadata[,DifferentiateBy])[i],])
    b=  rownames(metadata[metadata[,DifferentiateBy]!=unique
                          (metadata[,DifferentiateBy])[i],])
    data_a = Pathway_score1[,a]
    data_b = Pathway_score1[,b]
    data_ab = cbind(data_a,data_b)
    group1 = c(rep("group1",length(a)))
    group2 = c(rep("group2",length(b)))
    group = c(group1,group2)
    mat = GroupsDiffPathways(data_ab,
                             group,
                             p_val = p_val,
                             lfc = lfc,
                             up_pathways_number = up_pathways_number)

    up_pathways[[i]] = mat[,rownames(metadata)]
    up_pathwaysList[[i]] = rownames(mat)
    #print(rownames(mat))
  }
  Differential_mat =do.call("rbind", up_pathways)
  Differential_mat_1 = Differential_mat[,rownames(metadata)]

  Diff_Pathways = list(DiffMatpathway=Differential_mat_1,
                  Diffup_pathways = up_pathwaysList,
                  annotations = metadata)
  return(Diff_Pathways)
}

