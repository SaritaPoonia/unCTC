#------------------------------------------------------------------------
#-------------------- Pathway Enrichments score  ------------------------
#------------------------------------------------------------------------

#' Calculate Pathway enrichment score of all cells for all pathways
#' @description Calculate gene set enrichment across samples/cells using
#' the GSVA package in R, a non-parametric method and unsupervised
#' software programme.
#'
#' @param data_list List of gene expression data matricies. Genes/Features
#' should be in rows and cells/ samples in columns.
#' @param data_id List of names/ids of expression matrix
#' @param Genesets list of genesets/pathways.
#' @param min.size Minimum size of the resulting gene sets.
#' @param max.size Maximum size of the resulting gene sets.
#' @param min_Gene cell filter, filter out those cells which do not
#' express at least min_Gene genes
#' @param min_Sample gene filter, filter out genes which are not expressed
#' in at least min_Sample cells
#' @param Parallel_threads = Number of threads in parallel to execute process.
#'
#' @importFrom qusage read.gmt
#' @import GSVA
#' @import SingleCellExperiment
#'
#' @return PathwayData list of athway enrichment score and pathway metadata.
#' @export PathwayEnrichmentScore
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
#'                                         max.size=100,
#'                                         min_Sample = 5,
#'                                         min_Gene = 1500
#'                                         )
PathwayEnrichmentScore = function(data_list=list(),
                                  data_id=list(),
                                  Genesets,
                                  min.size=10,
                                  max.size=500,
                                  min_Sample = 5,
                                  min_Gene = 1500,
                                  Parallel_threads = 4L)
{
  #Create SingleCellObject of all expression data
  sce_obj = CreateSingleCellObject(data_list,min_sample =min_Sample,
                                   min_gene=min_Gene)
  sce_data = sce_obj$sce_norm

  #Fetch normalised data from SingleCellObject
  data_norm = as.matrix(sce_data@assays@data$Normalized_Data)
  colnames(data_norm) = sce_data@colData@listData$Samples
  rownames(data_norm) = rowData(sce_data)[,1]

  #Log transformation of normalised data matrix
  data_norm = log1p(data_norm)


  #Filter out those genes from genesets which are not present in data
  cat("Filter out those genes from genesets which are not present in data")
  cat("Filtering..")
  cat(" ")
  #row.names(data_norm) = toupper(rownames(data_norm))

  GeneList = list()
  for (i in seq_along(Genesets)) {
    gene2 = intersect((rownames(data_norm)),
                      (unlist(Genesets[i])))
    GeneList[[i]] = gene2
  }
  names(GeneList) = names(Genesets)

  cat("---------- Calculate pathway score with gsva()-----------")
  cat("Time to caculate GSVA score is directly depends on")
  cat("sample size of data and genesets size")
  cat(" ")
  Pathway_score = gsva(as.matrix(data_norm),
                       GeneList,
                       method= "gsva",
                       kcdf = "Gaussian",
                       min.sz=min.size,
                       max.sz=max.size,
                       parallel.sz=Parallel_threads,
                       mx.diff=FALSE,
                       verbose=TRUE)

  #Calcualte Pathway metadata

    meta_data = as.matrix(Meta_data(data_list,data_id))
    Pathway_metadata = as.data.frame(meta_data[colnames(Pathway_score),])

    PathwayData = list(Pathway_score = Pathway_score,
                     Pathway_metadata=Pathway_metadata )
  return(PathwayData)
}
