#------------------------------------------------------------------------
#------------------ Z score computation ---------------------------------
#------------------------------------------------------------------------

#' Claculate z-score of log transformed expression matrix
#' @description Claculate Z score of matrix
#'
#' @param data_mat log transformed Expresion matrix,
#' genes in row and cells in column
#'
#' @import Linnorm
#'
#' @return data_z zscore matrix
#' @export get_Zscore
#' @examples
#' data = unCTC::Poonia_et_al._TPMData
#' Filtered_data = data_filtering(data_mat= data,
#'                                min_sample =5,
#'                                min_gene=1500)
#' #Filtere data normalizationand log transformtion
#' library(Linnorm)
#' data_norm = log1p(Linnorm.Norm(Filtered_data))
#' data_zscore = get_Zscore(data_norm)
#'
get_Zscore <- function(data_mat){
  data_mat[data_mat == 0]<-1
  mean=colMeans(data_mat)
  std=apply(data_mat,2,sd)
  data_z = scale(data_mat,scale=std,center=mean)
  return(data_z)
}

#------------------------------------------------------------------------
#------------- Stouffer's score compuation ------------------------------
#------------------------------------------------------------------------

#' Calculate Stouffer's score
#' @description Stouffer's score computation
#'
#' @param data1 z score matrix
#' @param axis 1 = applies over rows, 2 = applies over columns
#'
#' @return Stouffer's score
#'
#' @import Linnorm
#'
#' @export get_stouffer
#'
#' @examples
#' data = unCTC::Poonia_et_al._TPMData
#' Filtered_data = data_filtering(data_mat= data,
#'                                min_sample =5,
#'                                min_gene=1500)
#' #Filtere data normalizationand log transformtion
#' library(Linnorm)
#' data_norm = log1p(Linnorm.Norm(Filtered_data))
#' data_zscore = get_Zscore(data_norm)
#' StoufferScore = get_stouffer(data_zscore,axis=2)
#'
get_stouffer <- function(data1,axis=2)
{
  Stouffers_score = apply(data1,axis,function(x){sum(x)/sqrt(length(x))})
  return(Stouffers_score)
}

#------------------------------------------------------------------------
#-------------------- Stouffer's Score Calculation ----------------------
#------------------------------------------------------------------------


#' Calculate Stouffer's Score for specific list of genes and plot
#' stouffer's score between different type of cells.
#' @description Calculate Stouffer,s score of genes passed as gene_list in
#' Stouffer_score()
#' @param data_list List of expression matrix
#' @param min_Gene cell filter, filter out those cells which do not
#' express at least min_Gene genes
#' @param min_Sample gene filter, filter out genes which are not expressed
#' in at least min_Sample cells
#' @param gene_list Vector of specific genes for which we
#' want to calculate Stouffer's score
#' @param data_id List of names/ids of expression matrix
#' @param MetaData Optional, List of metadata of expression matricies in same
#' order in which expression matricies in data_list, Column number and names
#' of all the MetaData in the list must be same
#' @param metaColPos column position in MetaData according
#' to which we want to group data, default is 1
#' @param metaColName Name of the column in MetaData according
#' to which we want to group data, default is "Class"
#' @param Groupby column name in MetaData according
#' to which data is grouped.
#' @param DDLKCluster_data output of DDLK_Clust.R method
#'
#' @importFrom ggplot2 ggplot
#'
#' @return genes_Stouffer_score1 Cell wise stouffer score with Class info
#' @export Stouffer_score
#'
#' @examples
#' data1 = unCTC::Poonia_et_al._TPMData
#' data2 = unCTC::Ding_et_al._WBC1_TPMData
#' Data_list = list(data1,data2)
#' Data_Id = list("data1","data2")
#' #LincrNA specific gene list
#' gene_list = unCTC::Breast_elevated_genes
#' S = Stouffer_score(data_list=Data_list,
#'                     min_Sample = 5,
#'                     min_Gene = 1500,
#'                     gene_list = gene_list,
#'                     Groupby = "data_id",
#'                     data_id= Data_Id)
#'
Stouffer_score = function(data_list=list(),
                          min_Sample = 5,
                          min_Gene = 1500,
                          gene_list=list(),
                          data_id=list(),
                          MetaData = list(),
                          metaColPos = 1,
                          metaColName = "Class",
                          Groupby = "Clusters",
                          DDLKCluster_data)
{
  #Colours for better visualization
  ColorKey = c("darkred","deepskyblue3","darkolivegreen4",
               "dark turquoise","pale violet red",
               "steelblue","forestgreen","gray2",
               "gray50","hotpink","lightslateblue",
               "tan4","yellow3","sienna4","orchid4")

  #Create comparison list
  if(Groupby=="Clusters")
  {clust = unique(DDLKCluster_data$PathwayDDLK_clust$Clusters)
  my_comparisons=list()
  for(i in seq_along(clust))
  {
    for(j in seq_along(clust))
    {
      if (i<j)
      {
        my_comparisons[[paste(i,j)]]=c(clust[[i]],clust[[j]])
      }
    }
  }
  }else if(missing(MetaData)){
    #Create metadata
    my_comparisons=list()
    for(i in seq_along(data_id))
    {
      for(j in seq_along(data_id))
      {
        if (i<j)
        {
          my_comparisons[[paste(i,j)]]=c(data_id[[i]],data_id[[j]])
        }
      }
    }

  }else
  {
    if (Groupby=="data_id")
    {
      my_comparisons=list()
      for(i in seq_along(data_id))
      {
        for(j in seq_along(data_id))
        {
          if (i<j)
          {
            my_comparisons[[paste(i,j)]]=c(data_id[[i]],data_id[[j]])
          }
        }
      }
    } else{
      my_comparisons=list()
      MetaData1 = as.data.frame(do.call("rbind", MetaData))
      Group = unique(MetaData1[,Groupby])
      for(i in seq_along(Group))
      {
        for(j in seq_along(Group))
        {
          if (i<j)
          {
            my_comparisons[[paste(i,j)]]=c(Group[[i]],Group[[j]])
          }
        }
      }
    }
  }

  #Create SingleCellObject of all expression data
  sce_obj = CreateSingleCellObject(data_list,min_sample = min_Sample,
                                   min_gene=min_Gene)
  sce_data = sce_obj$sce_norm

  #log transformation of expression matrix
  data_mat = log1p(sce_data@assays@data@listData$Normalized_Data)
  colnames(data_mat) = sce_data@colData@listData$Samples
  rownames(data_mat) = rowData(sce_data)[,1]
  #Covert gene names ito upper case
  row.names(data_mat) = make.names(toupper(rownames(data_mat)),unique = TRUE)

  if(missing(MetaData)){
  #Create metadata
  meta_data1 = Meta_data(data_list,data_id)
  metadata = meta_data1[colnames(data_mat), ,drop=FALSE]
  }else
  {
    MetaData1 = as.data.frame(do.call("rbind", MetaData))
    if(ncol(MetaData1)!=1){
      metadata = MetaData1[colnames(data_mat),]
    } else
    metadata = MetaData1[colnames(data_mat), ,drop=FALSE]
  }
  #Convert all genes in gene_list into uppercase to match with genes in data
  gene_list = toupper(gene_list)

  #Find common genes between gene_list and data
  ins_genes = intersect(gene_list,rownames(data_mat))

  #Calculate z score of complete matrix
  data_zscore = get_Zscore(data_mat)
  colnames(data_zscore) = colnames(data_mat)

  #Subset z score matrix with common genes
  gene_list_specificZscore = data_zscore[ins_genes,]

  #Calculate sample/Cell wise stouffer score
  genes_Stouffer_score = as.data.frame(get_stouffer(
                                       gene_list_specificZscore,2))
  genes_Stouffer_score$Class = metadata[,metaColPos]
  colnames(genes_Stouffer_score) = c("Stouffer_score",metaColName)
  if(missing(DDLKCluster_data)){
    genes_Stouffer_score1 = genes_Stouffer_score
  }else
  {
    genes_Stouffer_score1 = genes_Stouffer_score[
      rownames(DDLKCluster_data$PathwayDDLK_clust),]
    genes_Stouffer_score1$Clusters =
      DDLKCluster_data$PathwayDDLK_clust$Clusters

  }
  #Return dataframe with Stouffer Score and Class column
  Stouffer_score_and_pvalue = list(Stouffer_score = genes_Stouffer_score1,
                                   comparisons=my_comparisons)
  return(Stouffer_score_and_pvalue)
}
