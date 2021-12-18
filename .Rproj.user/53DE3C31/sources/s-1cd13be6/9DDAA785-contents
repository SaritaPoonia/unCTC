#------------------------------------------------------------------------
#-------------------- DDLK clusters  ------------------------------------
#------------------------------------------------------------------------

#' Cluster pathway mapped expression mtrix using Deep Dictionary Learning
#' with K-means clustering cost
#' @description clustering datasets with thousands of samples, and
#' incorporating the K-means clustering cost into the deep dictionary
#' learning (DDL) framework
#'
#' @param PathwayScore Pathway_score matrix form PathwayEnrichmentScore output
#' @param PathwayMetaData Pathway_metadata form PathwayEnrichmentScore output
#' @param n integer number of clusters for k means.
#' @param out.dir = Output directory to write Pathwayscore.
#' Default is current directory.
#' @param MetaData Optional, List of metadata of expression matricies in same
#' order in which expression matricies in data_list, Column number and names
#' of all the MetaData in the list must be same
#'
#' @importFrom utils write.csv
#' @importFrom utils read.delim
#'
#' @return Pathways_score_cluster Return:
#' 1. Pathway_score = Cell wise pathway enrichment score matrix,
#'    Pathways in row and cells/samples in column.
#' 2. PathwayDDLK_clust = sample/Cell wise DDLK cluster information.
#'
#' @examples
#' data1 = unCTC::Poonia_et_al._TPMData
#' data2 = unCTC::Ding_et_al._WBC1_TPMData
#' Data_list = list(data1,data2)
#' Data_Id = list("data1","data2")
#' Genesets = unCTC::c2.all.v7.2.symbols
#' Pathway_score = PathwayEnrichmentScore(data_list=Data_list,
#'                                         data_id= Data_Id,
#'                                         Genesets=Genesets,
#'                                         min.size=70,
#'                                         max.size=100)
#'
#' cluster_output = DDLK_Clust(PathwayScore = Pathway_score$Pathway_score,
#'                            PathwayMetaData=Pathway_score$Pathway_metadata,
#'                             n=3,
#'                             out.dir = getwd())
#'
#' @export DDLK_Clust
DDLK_Clust =  function(PathwayScore,
                       PathwayMetaData,
                       n,      #Number of cluster for k means (int)
                       out.dir = getwd(),  #Out directory to save file.
                       MetaData = list()
                       )
{
  #Call PathwayEnrichmentScore()
  Pathway_score = PathwayScore

  #Write Pathway score in a directory
  Pathway_score_csv = paste0(out.dir,"/Pathway_score.csv")
  write.csv(Pathway_score,file =Pathway_score_csv)

  #Pathway metadata
  Pathway_metadata = PathwayMetaData

  #Write metadata in out dir
  Pathway_metadata_csv = paste0(out.dir,"/Pathway_metadata.csv")
  write.csv(Pathway_metadata,Pathway_metadata_csv)


  #As we can not call 2 environment together so we can not use R Pathway score
  #in python script directly.
  #So you need to write Pathway enrichment scores in a same directory
  path <- paste(system.file(package="unCTC"), "DDLK_V2.py", sep="/")
  # Call DDLK python script
  system(paste('python3',path,'--n_clusters',n,' --out_dir',out.dir,sep=" "),
         wait=TRUE)
  #Read DDLK cluter info

  #Read DDLK cluster output
  DDLK_Cluster = read.delim(paste0(out.dir,"/DDLKlabels_nclusters_",n,".csv"),
                            header=FALSE)
  rownames(DDLK_Cluster) = colnames(Pathway_score)

  #Add cluster info to Pathway metadata
  Pathway_metadata$Clusters = paste0("Cluster_",DDLK_Cluster$V1)
  colnames(Pathway_metadata) = c("Data_id","Clusters")


  #Add cluster info in Pathway metadata
  Pathway_metadata1 = Pathway_metadata[order(Pathway_metadata$Clusters,
                                             decreasing = FALSE),]
  Pathway_score1 = Pathway_score[,rownames(Pathway_metadata1)]

  #Add given Metadata to the Pathway_metadata1
  if(missing(MetaData)){
    Pathway_metadata2 = Pathway_metadata1
  }else
  {
    MetaData1 = as.data.frame(do.call("rbind", MetaData))
    if(ncol(MetaData1)!=1){
      metadata = MetaData1[rownames(Pathway_metadata1),]
    } else {
      metadata = MetaData1[rownames(Pathway_metadata1), ,drop=FALSE]}
    Pathway_metadata2 = cbind(Pathway_metadata1,metadata)
  }

  #list Pathway_score and Pathway_metadata1 in a list
  Pathways_score_cluster = list(Pathway_score = Pathway_score1,
                               PathwayDDLK_clust = Pathway_metadata2)

  #Return Pathways_score_cluster list
  return(Pathways_score_cluster)
}

