#------------------------------------------------------------------------
#-------- Copy number variation in different arms of chromosome ---------
#------------------------------------------------------------------------

#' CNV in p and q arms of chromosome
#' @description Use single cell RNA-Seq expression data to identify copy
#' number variation at chromosomal level such as deletions or gains of
#' entire chromosome or large segments of chromosome
#'
#' @param data_list List of raw expression matrix. Genes should be in
#' rows and cells should be in columns in each data in the list.
#' @param data_id List of names/ids of expression matrix
#' @param min_Gene cell filter, filter out those cells which do not
#' express at least min_Gene genes
#' @param min_Sample gene filter, filter out genes which are not expressed
#' in at least min_Sample cells
#' @param path Path of output directory to save results
#' @param GenePositionFile A gene/chromosome positions file with chromosome
#' name, start, end position. "genecode hg19" positional file is given with
#' this package. Either you can use same using unCTC::gencode_v19_gene_pos
#' or can download from other sources.
#' @param threads_no (int) number of threads for parallel steps (default: 8)
#' @param obs.title Title of test/observation matrix.
#' Default is "Observations"
#' @param ref.title Title of reference matrix. Default is "References".
#' @param out.Filename Store results with out.Filename prefix.
#' Default is "inferCNV".
#' @param MetaData Optional, List of metadata of expression matricies in same
#' order in which expression matricies in data_list, Column number and names
#' of all the MetaData in the list must be same
#' @param Groupby Any column name from MetaData,which we want to use as
#' annotation file. Only applicable if MetaData is included.
#' @param cutoff The minimum average read counts per gene among
#' reference cells. (The default value is 1)
#' @param Reference_name Any one type of cell from data_id list or
#' any one cell type from column assign to Groupby.
#'
#' @import readtext
#' @import infercnv
#' @import D3GB
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#' @importFrom utils head
#' @importFrom utils write.csv
#' @importFrom utils read.delim
#'
#' @return p_and_q_arm_CNV
#'
#' @examples
#' data1 = unCTC::Poonia_et_al._TPMData
#' metadata = list(unCTC::Poonia_et_al._metaData)
#' Data_list = list(data1)
#' Data_Id = list("data1")
#' GenePoFile = unCTC::gencode_v19_gene_pos
#' GroupID = "GroupID"
#' ref = "P1"
#' path = getwd()
#' CNV_alterations(data_list=Data_list,
#'                 data_id= Data_Id,
#'                 min_Sample = 5,
#'                 min_Gene = 1500,
#'                 path=path,
#'                 GenePositionFile=GenePoFile,
#'                 threads_no=1,
#'                 MetaData=metadata,
#'                 Groupby=GroupID,
#'                 Reference_name=ref,
#'                 obs.title ="Observations",
#'                 ref.title = "References",
#'                 out.Filename = "inferCNV"
#'                 )
#'
#' @export CNV_alterations
CNV_alterations = function(
                     data_list=list(), # list of expression matrix
                     data_id=list(), # Expression matrix's name list
                     min_Sample = 5,
                     min_Gene = 1500,
                     path=" ",  # Path to save outputs
                     GenePositionFile=" ", # Gene order file
                     threads_no=8, # Number of threads
                     MetaData = list(),
                     Groupby = " ",
                     cutoff=1,
                     Reference_name=" ", # Reference Cells
                     obs.title ="Observations", # Observations matrix title
                     ref.title = "References",  # Reference matrix title
                     out.Filename = "inferCNV"  # prefix to save output
                     )
  {
  #Create SingleCellObject from all expression data in data_list
  sce_obj = unCTC::CreateSingleCellObject(data_list,min_sample = min_Sample,
                                   min_gene=min_Gene)
  sce_data = sce_obj$sce_raw


  data_mat = sce_data@assays@data@listData$Raw_filtered_data
  colnames(data_mat) = sce_data@colData@listData$FilterDataSample
  rownames(data_mat) = rowData(sce_data)[,1]

  #Assign Data id to each sample/cell
  if(missing(MetaData)){
  meta_data = as.matrix(Meta_data(data_list,data_id))
  metadata = as.data.frame(meta_data[colnames(data_mat),])
  }else
  {if(length(MetaData)==1){
    MetaData1 = as.data.frame(MetaData[[1]])
    MetaData1 = MetaData1[colnames(data_mat), ,drop=FALSE]
  }else
    MetaData1 = as.data.frame(do.call("rbind", MetaData))
    if(ncol(MetaData1)==1){
      metadata = MetaData1[colnames(data_mat), ,drop=FALSE]
    } else
      metadata = MetaData1[colnames(data_mat),Groupby ,drop=FALSE]
  }
  colnames(metadata) = NULL

  #Gene order file,colnames must be NULL
  pos_file = GenePositionFile

  #Creation of an infercnv object
  infercnv_obj1.2 = CreateInfercnvObject(raw_counts_matrix=data_mat,
                                         annotations_file = metadata,
                                         delim="\t",
                                         gene_order_file=pos_file,
                                         ref_group_names = Reference_name)

  infercnv_obj2_h.1 = infercnv::run(infercnv_obj1.2,
                                    cutoff=cutoff,
                                    # cutoff=1 works well for Smart-seq2,
                                    #and cutoff=0.1 works well for 10x Genomics
                                    out_dir=path,
                                    cluster_by_groups=TRUE,
                                    denoise=TRUE,
                                    HMM=TRUE,
                                    num_threads=threads_no)

  plot_cnv(infercnv_obj2_h.1,
           out_dir= path,
           obs_title=obs.title,
           ref_title=ref.title,
           cluster_by_groups=TRUE,
           x.center=1,
           x.range="auto",
           hclust_method='ward.D',
           color_safe_pal=FALSE,
           output_filename = out.Filename,
           output_format="pdf",
           png_res=300,
           dynamic_resize=0
  )

  # Find out arms which calculate P and q arm deletion
  #Read CNV region file from saved output directory
  cnv_region_file = readtext(
                     paste0(path,
                     "/HMM_CNV_predictions.*pred_cnv_regions.dat")
                     )

  cnv_region = read.delim(paste0(path,"/",cnv_region_file$doc_id))

  #Read CNV gene file from saved output directory
  cnv_gene_file = readtext(
                   paste0(path,
                  "/HMM_CNV_predictions.*pred_cnv_genes.dat")
                  )

  cnv_gene  =read.delim(paste0(path,"/",cnv_gene_file$doc_id))

  ins = intersect(rownames(cnv_region),rownames(cnv_gene))

  GRCh37_bands = GRCh37.bands # Chromosome name with p and q arm positions

  #Assign unique row name
  row.names(GRCh37_bands) = make.names(
                               paste0("chr",
                               GRCh37_bands$start,"_",GRCh37_bands$end),
                               unique = TRUE)
  GRCh37_bands$chr = paste0("chr",GRCh37_bands$chr)

  #Find genomic locations and their associated annotations.
  gr0 = with(GRCh37_bands, GRanges(chr, IRanges(start=start, end=end)))
  gr1 = with(cnv_region, GRanges(chr, IRanges(start=start, end=end)))

  #Finding/counting overlaps between objects containing genomic ranges
  hits = findOverlaps(gr1, gr0)
  ranges(gr1)[queryHits(hits)] = ranges(gr0)[subjectHits(hits)]
  class(gr1)

  ranges = as.data.frame(gr1@ranges)
  cnv_region$start = ranges$start

  cnv_region$end = ranges$end

  ranges2 = as.data.frame(gr0@ranges)
  GRCh37_bands$start = ranges2$start
  GRCh37_bands$end = ranges2$end

  row.names(GRCh37_bands) = make.names(
                                      paste0(GRCh37_bands$chr,"_",
                                      GRCh37_bands$start,"_",
                                      GRCh37_bands$end),unique = TRUE
                                      )


  row.names(cnv_region) = make.names(
                                  paste0(cnv_region$chr,"_",
                                         cnv_region$start,"_",
                                         cnv_region$end),unique = TRUE
                                    )
  ins2 = intersect(rownames(cnv_region),rownames(GRCh37_bands))

  GRCh37_bands_2 = GRCh37_bands[ins2,]
  cnv_region_2 = cnv_region[ins2,]
  cnv_region_2$chr_arm = paste0(cnv_region_2$chr,"_",GRCh37_bands_2$name)
  #View(cnv_region_2)
  cnv_region_3 = as.data.frame(cbind(cnv_region_2$cell_group_name,
                                     cnv_region_2$state,cnv_region_2$chr_arm))
  colnames(cnv_region_3) = c("cell_group","state","chr_arm")
  print("HMM state (1 = 2 copies loss, 2 = 1 copy loss, 3 = neutral,
        4 = 1 copy gain, 5 = 2 copies gain, 6 = 3+ copies gain),")
  print(head(cnv_region_3))
  write.csv(cnv_region_3,paste0(path,"/cnv_region_with_armInfo.csv"))
  return(p_and_q_arm_CNV = cnv_region_3 )
}





