#------------------------------------------------------------------------
#-------------------- Gene Violin plot  --------------------------------
#------------------------------------------------------------------------

#' Gene Violin plot
#' @description Create Violin plot of specific gene
#' and compare expression in all data in data list.
#'
#' @param data_list List of expression matricies
#' @param data_id List of names/ids of expression matrix
#' @param min_Gene cell filter, filter out those cells which do not
#' express at least min_Gene genes
#' @param min_Sample gene filter, filter out genes which are not expressed
#' in at least min_Sample cells
#' @param gene_symbol Gene should present in all expression data of data_list
#' @param MetaData Optional, List of metadata of expression matricies in same
#' order in which expression matricies in data_list, Column number and names
#' of all the MetaData in the list must be same
#' @param Groupby column name in MetaData according
#' to which data is grouped.
#' @param DDLKCluster_data PathwayDDLK_clust from DDLK_Clust.R method
#'
#' @import ggpubr
#' @import ggplot2
#'
#' @return Violin_plot violin plot and plot object of given gene.
#'
#' @export Gene_Violin_plots
#'
#' @examples
#' data1 = unCTC::Poonia_et_al._TPMData
#' data2 = unCTC::Ding_et_al._WBC1_TPMData
#' Data_list = list(data1,data2)
#' Data_Id = list("data1","data2")
#' Violin_plots = Gene_Violin_plots(data_list=Data_list,
#'                                  data_id= Data_Id,
#'                                  min_Sample = 5,
#'                                  min_Gene = 1500,
#'                                  Groupby = "data_id",
#'                                   gene_symbol ="TSPAN6" )

Gene_Violin_plots = function(data_list =list(),
                             data_id = list(),
                             min_Sample = 5,
                             min_Gene = 1500,
                             gene_symbol = "",
                             MetaData = list(),
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
  {clust = unique(DDLKCluster_data$Clusters)
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

  #Normalised data from data expression list
  sce_obj = CreateSingleCellObject(data_list,min_sample = min_Sample,
                                   min_gene=min_Gene)
  sce_data = sce_obj$sce_norm

  data_norm = as.matrix(sce_data@assays@data$Normalized_Data)
  colnames(data_norm) = sce_data@colData@listData$Samples
  rownames(data_norm) = rowData(sce_data)[,1]
  data_norm = log1p(data_norm)

  #Transpose data matrix
  data_norm_t = as.data.frame(t(data_norm))

  #If MetaData is not given then calculate metadata with Meta_data()
  if(missing(MetaData)){
  #Create meta data
  meta_data1 = Meta_data(data_list,data_id)
  metadata = meta_data1[colnames(data_norm), ,drop=FALSE]
  data_norm_t1 = cbind(data_norm_t,metadata)
  }else
  {
    MetaData1 = as.data.frame(do.call("rbind", MetaData))
    if(ncol(MetaData1)!=1){
      MetaData2 = MetaData1[colnames(data_norm),]
    } else
    {
      MetaData2 = MetaData1[colnames(data_norm), ,drop=FALSE]
    }
    meta_data1 = Meta_data(data_list,data_id)
    metadata = meta_data1[colnames(data_norm), ,drop=FALSE]
    data_norm_t1 = cbind(data_norm_t,metadata,MetaData2)
  }

  if(missing(DDLKCluster_data)){
    data_norm_t2 = data_norm_t1
  }else
  {
    data_norm_t2 = data_norm_t1[rownames(DDLKCluster_data),]
    data_norm_t2$Clusters = DDLKCluster_data$Clusters

  }

  #Vioin plot

  Violin_plot = ggviolin(data_norm_t2, x=Groupby, y=gene_symbol, fill = Groupby,
           legend = "top",palette =ColorKey) +
    stat_compare_means(comparisons = my_comparisons,label = "p.signif",
                       method = "t.test",ref.group = ".all.")

  Violin_plot = list(Violin_plot_obj = data_norm_t2,Violin_plot = Violin_plot)
  return(Violin_plot)
}

