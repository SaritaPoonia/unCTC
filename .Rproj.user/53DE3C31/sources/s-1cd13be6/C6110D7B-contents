#------------------------------------------------------------------------
#-------------  Data filtering function ---------------------------------
#------------------------------------------------------------------------

#' Filter out poor quality data
#' @description Filter out poor quality genes and cells
#'
#' @param data_mat Expression matrix, genes in row and cells in column
#' @param min_gene cell filter, filter out those cells which do not
#' express at least min_gene genes
#' @param min_sample gene filter, filter out genes which are not expressed
#' in at least min_sample cells
#'
#' @import viridis
#' @import magrittr
#' @import umap
#' @import reticulate
#'
#' @return Filtered_data return filtered data
#' @examples
#' data = unCTC::Poonia_et_al._TPMData
#' Filtered_data = data_filtering(data_mat= data,
#'                                min_sample =5,
#'                                min_gene=1500)
#' @export
data_filtering = function(data_mat,
                          min_sample =5,
                          min_gene=1500)
{
  #gene filtering
  good_genes = apply(
    data_mat, 1,
    function(x) sum(x>0)) >= min_sample

  #subsetting the data
  data1= as.matrix(data_mat[good_genes,])

  ##cell filtering
  good_cells = apply(data1,2,
    function(x) sum(x>0)) >= min_gene

  ##subsetting the data
  Filtered_data = data1[,good_cells]

  return(Filtered_data)
}

#------------------------------------------------------------------------
#------------- Meta Data computation ------------------------------------
#------------------------------------------------------------------------

#' Meta Data computation
#' @description data frame with one column 'data_id'.
#' data_id is a vector of length equals to total number
#' of cells in all matricies.
#'
#' @param data_list List of expression matricies
#' @param data_id List of names/ids of expression matricies
#'
#' @return batch_id1 return dataframe with data_id column
#' @import Linnorm
#' @import SingleCellExperiment
#' @export Meta_data
#'
#' @examples
#' data1 = unCTC::Poonia_et_al._TPMData
#' data2 = unCTC::Ding_et_al._WBC1_TPMData
#' Data_list = list(data1,data2)
#' Data_Id = list("data1","data2")
#' metadata = Meta_data(data_list=Data_list,
#'                       data_id = Data_Id)
#'
Meta_data =  function(data_list=list(),
                      data_id = list()
                      )
{
  Batch_id = list()
  for (i in seq_along(data_list))
  {
    b = rep(paste0(data_id[i]),ncol(data_list[[i]]))
    Batch_info = as.data.frame(b)
    rownames(Batch_info) = colnames(data_list[[i]])
    Batch_id[[i]] = Batch_info

  }
  batch_id1 = as.data.frame(do.call("rbind", Batch_id))
  colnames(batch_id1) = "data_id"
  return(batch_id1)
}

#------------------------------------------------------------------------
#----------------------- Data integration  ------------------------------
#------------------------------------------------------------------------

#' Integrate single cell expression matricies in the list
#' @description Integrates all the data matrix present
#'  in the list on the basis of common genes.
#'
#' @param data_list List of data matrix. Genes/Features
#'  should be in rows and cells/ samples in columns.
#' @return integrated_data Return integrated data matrix
#'  in which genes in rows and cells in column.
#' @export Data_integration
#'
#' @examples
#' data1 <- unCTC::Poonia_et_al._TPMData
#' data2 = unCTC::Ding_et_al._WBC1_TPMData
#' Data_list = list(data1,data2)
#' Integrated_data = Data_integration(data_list=Data_list)

Data_integration = function(data_list)
{
  #If only one expression data matrix in list
  common_genes=rownames(data_list[[1]])
  integrated_data=data_list[[1]]

  #If more than one expression data matricies in list
  if (length(data_list)>1)
  {
    for(i in 2:length(data_list))
    {
      temp= rownames(data_list[[i]])
      common_genes=intersect(common_genes,temp)
    }
    integrated_data=integrated_data[common_genes,]
    for(i in 2:length(data_list))
    {
      integrated_data=cbind(integrated_data,
                    data_list[[i]][common_genes,])
    }
  }
  return(integrated_data)
}

#------------------------------------------------------------------------
#-------------- Creating SingleCellExperiment instances -----------------
#------------------------------------------------------------------------

#' Create SingleCellExperiment instance
#' @description Creating SingleCellExperiment
#'  instance of the data. Steps in creating
#'  SingleCellExperiment instance are like this.
#' 1. Integrate all the data matrix passed in the
#' list on the basois of common genes.
#' 2. Filter out low quality data.
#' 3. Use Linnorm.Norm function to remove batch
#' effect between data taken from different
#' experiments/bathes
#' 4. Create SingleCellExperiment instances. In
#' this instance two type of data format present.
#'
#' Raw_filtered_data = Low quality gene and cell filtered
#' out matrix without batch effect correction.
#'
#' @param data_list List of expression matrix.
#' Genes should be in rows and cells should be
#' in columns in each data in the list.
#' @param min_gene cell filter, filter out those cells which do not
#' express at least min_gene genes
#' @param min_sample gene filter, filter out genes which are not expressed
#' in at least min_sample cells
#'
#' @return sce SingleCellExperimemnt instance od data
#'  passed in a list.
#'
#' @import Linnorm
#' @import SingleCellExperiment
#' @import S4Vectors
#' @import SummarizedExperiment
#'
#' @export CreateSingleCellObject
#'
#' @examples
#' data1 = unCTC::Poonia_et_al._TPMData
#' data2 = unCTC::Ding_et_al._WBC1_TPMData
#' Data_list = list(data1,data2)
#' sce_obj = CreateSingleCellObject(data_list=Data_list,
#'                                   min_sample =5,
#'                                   min_gene=1500)
#'
CreateSingleCellObject = function(data_list=list(),
                                  min_sample =5,
                                  min_gene=1500
                                  )
{
  int_data= as.data.frame(Data_integration(data_list))

  #------- data filtering ---------
  filter_data = data_filtering(int_data,min_sample = min_sample,
                               min_gene=min_gene)
  #-------- Normalise data --------
  data_norm = Linnorm.Norm(filter_data)

  #-------- Normalised data --------
  Normalized_Data = as.matrix(data_norm)
  Features = rownames(data_norm)
  Samples = colnames(data_norm)
  colnames(Normalized_Data) = NULL
  rownames(Normalized_Data) = NULL
  norm_d =  Normalized_Data
  #-------- Filtered data --------
  Samples_raw = colnames(filter_data)
  Features_raw = rownames(filter_data)
  colnames(filter_data) = NULL
  rownames(filter_data) = NULL

  sce_norm <- SingleCellExperiment(
    list(Normalized_Data = norm_d),
    colData = DataFrame(Samples=Samples),
    rowData = DataFrame(Features=Features)
  )

  sce_raw <- SingleCellExperiment(
    list(Raw_filtered_data = filter_data),
    colData = DataFrame(FilterDataSample = Samples_raw ),
    rowData = DataFrame(FilterDataFeatures = Features_raw))
  return(list(sce_norm=sce_norm,sce_raw = sce_raw))
}



