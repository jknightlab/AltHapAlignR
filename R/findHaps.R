#' To read multiple bam files
#' 
#' @param bamFiles bam files 
#' @param bamNames sample names of each bam file
#' @param mapper mapping program used, only tophat format available 
#' @export
#' @import Rsamtools
#' @import GenomicAlignments
#' @examples
#' output <- readBamFiles(bamfiles, name_of_bamfiles, mapper="tophat")
#' 
#' # Output data
#' # mapping statistics
#' mapping_stat <- output$mapping_stat
#' # read ids
#' uniqueFragmentID <- output$uniqueFragmentID
#' # paired reads
#' list_pairedReads <- output$list_pairedReads


readBamFiles <- function(bamFiles, bamNames, mapper="tophat"){
  n_bams <- length(bamFiles)
  message( n_bams, " BAM files ", bamFiles, " (mapped by ", mapper, ")" )
  
  what= c("qname", "flag", "cigar")
  if(mapper=="tophat"){
    tags <- c("NM", "NH", "HI")
  }
  if(mapper=="bwa"){
    tags <- c("NM", "X0", "X1")
  }
  
  flags <- scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE)
  param=ScanBamParam(what=what, tag=tags, flag=flags)
  
  pairedReads <- list()
  dumped <- list()
  fid <- c()
  hap_stat <- data.frame()
  for(i in 1:n_bams){
    message("Load paired reads from ", bamNames[i], " sample ")
    pairedReads[i] <- readGAlignmentPairs(bamFiles[i], use.names=FALSE, param=param)
    dumped[i] <- getDumpedAlignments()
    
    fid <- c(fid, pairedReads[[i]]@first@elementMetadata$qname)
    tmp <- data.frame(hap=bamNames[i],
                      NH=pairedReads[[i]]@first@elementMetadata$NH,
                      NM=(pairedReads[[i]]@first@elementMetadata$NM + pairedReads[[i]]@last@elementMetadata$NM) )
    hap_stat <- rbind(hap_stat, tmp)
  }
  fid <- data.frame(FragmentID=unique(fid) )
  hap_stat <- data.frame(table(hap_stat))
  names(pairedReads) <- bamNames
  
  output <- list( mapping_stat=hap_stat,
                  uniqueFragmentID=fid,
                  list_pairedReads=pairedReads)
  return(output)
}




#' To filter gene feature from a GTF file
#' 
#' @param gtf GTF file
#' @param type selecting gene feature
#' @param max_level annotation level
#' @param status gene status
#' @export
#' @import GenomicAlignments
#' @examples
#' output <- filterGenes(gtf, type="protein_coding", max_level=2, status="KNOWN")
#' 
#' # Output data
#' # summary of gene information, format: data.frame


filterGenes <- function(gtf, type="protein_coding", max_level=2, status="KNOWN"){

       anno <- gtf[gtf@elementMetadata@listData$type=="gene"] 
      anno.chr <- seqnames(anno)
       anno.start <- anno@ranges@start
       anno.end <- anno@ranges@start + anno@ranges@width

       anno.info <- gsub("\"|;", "", values(anno)$group)
       anno.info <- strsplit(gsub("  ", " ", anno.info), split=" ")
       anno.info <- lapply(anno.info, function(x) {
                       data <- x[seq(2, length(x), 2)]
                       names(data) <- x[seq(1, length(x), 2)]
                       data })
       attribs <- names(anno.info[[1]])
       anno.info <- do.call(data.frame, c(lapply(attribs, function(x) sapply(anno.info, "[", x)), stringsAsFactors = FALSE))
       colnames(anno.info) <- attribs
       anno.info <- data.frame(chr=anno.chr, start=anno.start, end=anno.end, anno.info )

       anno.info <- subset(anno.info, gene_type==as.character(type) & gene_status==as.character(status) & max_level <= as.numeric(level) )
	

	return(anno.info)
}




#' To filter exon feature from a GTF file
#' 
#' @param gtf GTF file
#' @param type selecting exon feature
#' @param max_level annotation level
#' @param status exon status
#' @export
#' @import GenomicAlignments
#' @examples
#' output <- filterExons(gtf, type="protein_coding", max_level=2, status="KNOWN")
#' 
#' # Output data
#' # summary of exon information, format: data.frame

filterExons <- function(gtf, type="protein_coding", max_level=2, status="KNOWN"){

  anno <- gtf[gtf@elementMetadata@listData$type=="exon"]
  anno.chr <- seqnames(anno)
  anno.start <- anno@ranges@start
  anno.end <- anno@ranges@start + anno@ranges@width
  
  anno.info <- gsub("\"|;", "", values(anno)$group)
  anno.info <- strsplit(gsub("  ", " ", anno.info), split=" ")
  anno.info <- lapply(anno.info, function(x) {
    data <- x[seq(2, length(x), 2)]
    names(data) <- x[seq(1, length(x), 2)]
    data })
  attribs <- names(anno.info[[1]])
  anno.info <- do.call(data.frame, c(lapply(attribs, function(x) sapply(anno.info, "[", x)), stringsAsFactors = FALSE))
  colnames(anno.info) <- attribs
  exon.info <- data.frame(chr=anno.chr, start=anno.start, end=anno.end, anno.info )
  
  exon.info <- subset(exon.info, as.character(transcript_type)==as.character(type) & as.character(transcript_status)==as.character(status) )
  exon.info <- subset(exon.info, as.numeric(level) <= as.numeric(max_level) )
  
  return(exon.info)
}




#' To get overlapping exons in different genes
#' 
#' @param overlappingTrx output from getOverlappingTrx
#' @param anno_exon output from filterExons
#' @export
#' @import GenomicAlignments
#' @examples
#' output <- getOverlappingExons(overlappingTrx, anno_exon)
#' 
#' # Output data
#' # overlapping exons in different genes

getOverlappingExons <- function(overlappingTrx, anno_exon){
 
      overlappingExons <- c()
      for(ch in unique(as.character(overlappingTrx$chr))){
            message(paste0("filtering overlapping exons in different genes in ", ch))
            s_trx <- subset(overlappingTrx, chr==ch)
            n=0
            
            for(i in 1:max(s_trx$overlapping)){
                  g1 <- s_trx[s_trx$overlapping==i,]$gene_name[1]
                  g2 <- s_trx[s_trx$overlapping==i,]$gene_name[2]
                  s_exon <- subset(anno_exon, chr==ch & ( gene_name==g1 | gene_name==g2) )
                  s_exon <- s_exon[,c("chr", "start", "end", "gene_name", "exon_id")]
                  s_exon <- s_exon[order(as.numeric(s_exon$start)), ]
                  s_exon <- unique(s_exon)
                  
                  pos_start <- s_exon$start[1]
                  pos_end <- s_exon$end[1]
                  
                  for(i in 2:nrow(s_exon)){
                        if(s_exon$start[i] >= pos_start & s_exon$start[i] <= pos_end){
                            if(as.character(s_exon$gene_name[i])!=as.character(s_exon$gene_name[i-1])){
                                  n = n +1
                                  overlappingExons <- rbind(overlappingExons, 
                                              data.frame(s_exon[i-1,], overlapping=n ), 
                                              data.frame(s_exon[i,], overlapping=n ) )
                            }
                            else{
                              if(i==2){
                                  overlappingExons <- rbind(overlappingExons,
                                                          data.frame(s_exon[1,], overlapping=0 ) )
                              }
                                  overlappingExons <- rbind(overlappingExons,
                                                        data.frame(s_exon[i,], overlapping=0 ) )
                            }
                            
                        }
                        else{
                            if(i==2){
                              overlappingExons <- rbind(overlappingExons,
                                                        data.frame(s_exon[1,], overlapping=0 ) )
                            }
                            overlappingExons <- rbind(overlappingExons,
                                                    data.frame(s_exon[i,], overlapping=0 ) )
                        }
                        pos_start <- s_exon$start[i]
                        pos_end <- s_exon$end[i]
                  }
                  overlappingExons <- unique(overlappingExons)
                  
            }
      }
      
      return(overlappingExons)
}



#' To get overlapping transcripts in genes
#' 
#' @param gene_cordi data.table format (chr, gene_name, start, end)
#' @export
#' @import GenomicAlignments
#' @examples
#' output <- getOverlappingTrx(gene_cordi)
#' 
#' # Output data
#' # overlapping transcripts in genes

getOverlappingTrx <- function(gene_cordi){

        overlappingTrx <- c()
        for(ch in unique(as.character(gene_cordi$chr))){
              print(ch)
              s <- subset(gene_cordi, chr==ch)
              s <- s[order(as.numeric(s$start)),]
              
              pos_start <- s$start[1]
              pos_end <- s$end[1]
              n=0
              for(i in 2:nrow(s)){
                    if(s$start[i] >= pos_start & s$start[i] <= pos_end){
                        n = n +1
                        overlappingTrx <- rbind(overlappingTrx, 
                                     data.frame(s[i-1,], overlapping=n ), 
                                     data.frame(s[i,], overlapping=n ) )
                    }
                    pos_start <- s$start[i]
                    pos_end <- s$end[i]
              }
              overlappingTrx <- unique(overlappingTrx)
        }
        return(overlappingTrx)
}




#' predictRecombHapByGenes: To predict if a sample is  heterogeneous or homogeneous AND to predic haplotypes if the sample is heterogeneous
#' 
#' @param bamFiles bam files
#' @param bamChrs haplotype names
#' @param bamNames bam file names
#' @param sample_name output file name
#' @param gtf_MHC gtf file
#' @param anno_type feature of gene (ex, "All", protein_coding", "processed_transcript", etc.)
#' @param anno_level feature level (ex, 1, 2, or 3)
#' @param anno_status annotation status of gene feature (ex, "KNOWN", "NOVEL", or "PUTATIVE")
#' @param minDP_per_gene minimum read depth in genes
#' @param numWorkers the number of nodes
#' @param uniq_mapping to get only uniquely mapped read 
#' @param editing_distance number of mismatch
#' @param read_length read length in bam files
#' @param save_mapping_info 1 for saving all mapping data, it will produce big files. 
#' @param mapper output from only tophat availble 
#' 
#' @export
#' @import Rsamtools
#' @import GenomicAlignments
#' @import reshape2
#' @import foreach
#' @import doParallel
#' @import plyr
#' @importFrom data.table data.table
#' @examples
#' # Input data: bam files from mapping to two haplotype references  
#' output <- predictRecombHapByGenes(bamFiles, bamChrs, bamNames, sample_name="sample_name", gtf_MHC, anno_type="protein_coding", anno_level=2, anno_status="KNOWN", minDP_per_gene=30, numWorkers=2, uniq_mapping=1, editing_distance=0, read_length=50, save_mapping_info=0, mapper="tophat")
#' 
#' # Output data
#' # *mhc_heatmap.pdf : heatmap of selected haplotypes
#' # *max_comb.txt : mapping rates and gene counts in all pairs of haplotypes
#' # *selected_genes.txt : Mapping ranks in all pairs of haplotypes
#' # *mapping_stats_to_genes.txt : gene counts in haployptes
#' # *rates_of_weighed_reads_by_genes.txt :unique  mapping rates in all pairs of haplotypes




predictRecombHapByGenes <- function(bamFiles, bamChrs, bamNames, sample_name="sample_name", gtf_MHC, anno_type="protein_coding", anno_level=2, anno_status="KNOWN", minDP_per_gene=30, numWorkers=2, uniq_mapping=1, editing_distance=0, read_length=50, save_mapping_info=1, mapper="tophat"){


  ed_dis <- as.numeric(editing_distance)
  output_name <- as.character(sample_name)
  read_len <- as.numeric(read_length)
  
	if(anno_type!="All"){
		anno_exon <- filterExons(gtf, type=as.character(anno_type), max_level=as.numeric(anno_level), status=as.character(anno_status))
		anno_trx <- filterTrx(gtf, type=as.character(anno_type), max_level=as.numeric(anno_level), status=as.character(anno_status))
	}

	if(anno_type=="All"){
		anno_exon <- filterExons(gtf, max_level=as.numeric(anno_level) )
                anno_trx <- filterTrx(gtf, max_level=as.numeric(anno_level) )
	}
  
  overlappingTrx <- getOverlappingTrx(anno_trx)
  overlappingExons <- getOverlappingExons(overlappingTrx, anno_exon)
 


  #####################################
  ### Aliases and manual movement
  #####################################

  # 1) VARSL -> VARS2, C6orf205 -> MUC21
 
  anno_exon$gene_name <- gsub("VARSL", "VARS2", anno_exon$gene_name)
  anno_trx$gene_name <- gsub( "VARSL", "VARS2", anno_trx$gene_name)
  overlappingTrx$gene_name <- gsub("VARSL", "VARS2", overlappingTrx$gene_name)
  overlappingExons$gene_name <- gsub("VARSL", "VARS2", overlappingExons$gene_name)
 
  anno_exon$gene_name <- gsub("C6orf205", "MUC21", anno_exon$gene_name)
  anno_trx$gene_name <- gsub( "C6orf205", "MUC21", anno_trx$gene_name)
  overlappingTrx$gene_name <- gsub("C6orf205", "MUC21", overlappingTrx$gene_name)
  overlappingExons$gene_name <- gsub("C6orf205", "MUC21", overlappingExons$gene_name)
  
  
  # 2) OR10C1 reads in other haps goes to OR11A1 in pgf.
  ## if reads in OR11A1 in pgf mapped to OR10C1 in other haplotypes, descard reads in OR11A1 in pgf
  


  # exclude overlapping trx to different genes
  useTrx <- anno_trx[!(paste0(anno_trx$chr, "_", anno_trx$gene_name) %in% paste0(overlappingTrx$chr, "_", overlappingTrx$gene_name)), ]
  
  
  
#  1) exclude overlapping exons to different genes
#  useExon <- subset(overlappingExons, overlapping==0)


# 2) use all trx and exons
  useExon <- overlappingExons 
  


  
  # gtf: GRanges format
  n_bams <- length(bamFiles) 
  output <- readBamFiles(bamFiles, bamNames, mapper="tophat")
  message( paste0( n_bams, " BAM files loaded: ", Sys.time() )  )




  
  mapping_stat <- output$mapping_stat
  uniq_fid <- output$uniqueFragmentID
  list_pairedReads <- output$list_pairedReads
 
   
	# getting mapping information  
  get_mapped_reads <- function(paired_reads, bam_name){
      tmp <- data.frame( FragmentID=paired_reads@first@elementMetadata$qname,
                        hap_first_start=paired_reads@first@start,
                        hap_last_start=paired_reads@last@start,
                        hap=(paired_reads@first@elementMetadata$NM+ paired_reads@last@elementMetadata$NM),
                        NH=paired_reads@first@elementMetadata$NH,
                        hap_first_cigar=paired_reads@first@cigar,
                        hap_last_cigar=paired_reads@last@cigar,
			insert_size=(paired_reads@last@start-paired_reads@first@start) )
 
     if(uniq_mapping==1){
                tmp <- subset(tmp, as.numeric(as.character(NH))==1)
                #tmp <- subset(tmp, NH <= 3)
      }
      tmp <- subset(tmp, as.numeric(hap)<=as.numeric(editing_distance))
      colnames(tmp)[4] <- as.character(bam_name)
      return(tmp)
  }





  message("Get data.frame of mapped reads...")
  base.paired <- list()
  for(i in 1:n_bams){
      base.paired[[i]] <- get_mapped_reads(list_pairedReads[[i]], bamNames[i])
      base.paired[[i]]$hap_first_end=gsub("M|S|N|I|D", " ", base.paired[[i]]$hap_first_cigar )
      base.paired[[i]]$hap_last_end=gsub("M|S|N|I|D", " ", base.paired[[i]]$hap_last_cigar )
      
      a <- apply(base.paired[[i]], 1, function(x){
        sum( as.numeric(strsplit(x[9], " ")[[1]] ) )
      } )
      base.paired[[i]]$hap_first_end <- a + base.paired[[i]]$hap_first_start -1
      
      a <- apply(base.paired[[i]], 1, function(x){
              sum( as.numeric(strsplit(x[10], " ")[[1]] ) )
      } )
      base.paired[[i]]$hap_last_end <- a + base.paired[[i]]$hap_last_start -1
  }    
  names(base.paired) <- bamChrs
  




    ###########################################################
    #  
    #  to exclude reads
    #  1) start position is futher from trx start position
    #  2) end position is further from trx end position. 
    #
    ###########################################################

    message("Filtering reads mapped to the outside of gene regions...")

    ft.base.paired <- list()
    for(b in bamChrs){     
            message(b)
            tmp <- base.paired[[b]]
            tmp_useTrx <- subset(anno_trx, as.character(chr)==b)
            tmp_useTrx <- tmp_useTrx[order(tmp_useTrx$start),]
            tmp_useTrx$start <- tmp_useTrx$start-1000
            tmp_useTrx$end <- tmp_useTrx$end+1000
            
            ft_list <- apply(tmp_useTrx, 1, function(x){
                    a <- tmp[tmp$hap_first_start > as.numeric(x[3]) & tmp$hap_first_end < as.numeric(x[4]) & tmp$hap_last_start > as.numeric(x[3]) & tmp$hap_last_end < as.numeric(x[4]) ,]
                    return(a)
            })
            ft_list <- do.call(rbind.data.frame, ft_list)
            

            print(dim(tmp))
            print(dim(ft_list))
            ft_list <- unique(ft_list)
            print(dim(ft_list))
            
            ft.base.paired[[b]] <- ft_list
    }
        base.paired <-  ft.base.paired





    ###################################################
    #
    #  to exclude exons overlapping in different genes
    #
    ###################################################

      message("Filtering reads mapped to exons overlapping in different genes...")

    anno_trx_and_exon_ftrd <- rbind(data.frame(useTrx, fea="trx"),
                           	data.frame(useExon[,c("chr", "gene_name", "start", "end")], fea="exon" ) )

    tmp <- anno_exon[,c("chr", "gene_name", "start", "end")]
    anno_trx_and_exon_all <- rbind(data.frame(anno_trx, fea="trx"),
                           data.frame(unique(tmp), fea="exon" ) )

    excludeOverlappingExons <- subset(overlappingExons, overlapping >0)
    excludeOverlappingExons <- excludeOverlappingExons[, c("chr", "gene_name", "start", "end", "exon_id", "overlapping")]
    

    
  message( paste0( "Reads to genes (overlapping exons to different genes): ", Sys.time() )  )

  cl<-makeCluster(numWorkers)
  registerDoParallel(cl)

  reads_to_MultiGenes <- list()
  reads_to_MultiGenes <- foreach(i=1:n_bams) %dopar% {
    
    base_hap <- base.paired[[i]]
    sub_anno <- subset(excludeOverlappingExons, as.character(chr)==names(base.paired)[i])
    sub_anno <- unique(sub_anno[,1:5])
	
    base_hap_to_gene <- apply(sub_anno, 1, function(x){

        message(as.character(x[2]))
        tmp <- subset(base_hap, as.numeric(hap_first_start) >= as.numeric(x[3]) & as.numeric(hap_first_start) <= as.numeric(x[4]) )
        tmp <- subset(tmp, as.numeric(hap_first_end) >= as.numeric(x[3]) & as.numeric(hap_first_end) <= as.numeric(x[4]) )
        tmp <- subset(tmp, as.numeric(hap_last_start) >= as.numeric(x[3]) & as.numeric(hap_last_start) <= as.numeric(x[4]) )
        tmp <- subset(tmp, as.numeric(hap_last_end) >= as.numeric(x[3]) & as.numeric(hap_last_end) <= as.numeric(x[4]) )
        
        
        tmp_sub <- subset(sub_anno, gene_name==as.character(x[2]) )       
        
        if( nrow(tmp) >= minDP_per_gene){
          
            ft_anno <- c()
            for(e in 1:nrow(tmp_sub)){
                print(e)
                tmp2 <- subset(tmp, as.numeric(hap_first_start) >= as.numeric(tmp_sub$start)[e] & hap_first_start <= as.numeric(tmp_sub$end)[e] )
              
                if(nrow(tmp2) >=5 ){
              
                  tmp2 <- data.frame(tmp2, gene_name=as.character(tmp_sub$gene_name[1])  )
                  tmp2 <- unique(tmp2)
                  ft_anno <- rbind(ft_anno, tmp2)
                }
            }
            ft_anno <- unique(ft_anno) 

            return(ft_anno)
        }   

    })

      if(!is.null(base_hap_to_gene)) {
          base_hap_to_gene <- do.call(rbind.data.frame, base_hap_to_gene)
          base_hap_to_gene <- unique(base_hap_to_gene)
    
          reads_to_MultiGenes[[i]] <- base_hap_to_gene
      }
  }

  names(reads_to_MultiGenes) <- bamChrs
  message( paste0( "Done: ", Sys.time() )  )






  message( paste0( "Combine haplotypes by read: ", Sys.time() )  )  
  comb_mapping <- c()
  for(i in 1:n_bams){
    comb_mapping <- c(comb_mapping, as.character( reads_to_MultiGenes[[i]]$FragmentID)  )
  }
  comb_mapping <- unique(comb_mapping)
  comb_mapping <- data.frame(FragmentID=comb_mapping)

  if(nrow(comb_mapping) >0){
        haps.reads_to_MultiGenes <- list()
        haps.reads_to_MultiGenes <- foreach(i=1:n_bams) %dopar% {
              b <- reads_to_MultiGenes[[i]]
              b <- b[,c(1,4,5,11)]
              b <- unique(b)
              tmp <- merge(comb_mapping, b, by="FragmentID", all=TRUE)
              tmp <- data.frame(tmp, bamNames[i])
              names(tmp) <- c("FragmentID", "NM", "NH", "gene_name", "Hap")
    
              haps.reads_to_MultiGenes[[i]] <- tmp
        }
        names(haps.reads_to_MultiGenes) <- bamChrs
        haps.reads_to_MultiGenes <- do.call(rbind.data.frame, haps.reads_to_MultiGenes)
        message( paste0("Done: ", Sys.time() ) ) 

  }

	ft_reads <- data.frame(table(haps.reads_to_MultiGenes[,c(1,4)]))
	ft_reads <- subset(ft_reads, Freq >0)

	ft_reads <- data.frame(table(ft_reads[,1]))
	ft_reads <- subset(ft_reads, Freq >1)

	ft_reads <- haps.reads_to_MultiGenes[haps.reads_to_MultiGenes$FragmentID %in% ft_reads$Var1, ]
	imsi <- data.frame(table(ft_reads$gene_name) )
  write.table(imsi, paste0(sample_name, ".overlappingGenes.txt"), sep="\t", quote=F)





  message( paste0( "Reads to genes: ", Sys.time() )  )

  #cl<-makeCluster(numWorkers)
  #registerDoParallel(cl)

  reads_to_genes <- list()


for(i in 1:n_bams ){
	   
    base_hap <- base.paired[[i]]
    base_hap <- base_hap[!(as.character(base_hap$FragmentID) %in% ft_reads$FragmentID) , ]



#     # using trx anno
      base_anno_trx <- subset(anno_trx, as.character(chr)==names(base.paired)[i])

    # filtered 
    sub_anno <- subset(anno_trx_and_exon_ftrd, as.character(chr)==names(base.paired)[i] )
    # all
    #sub_anno <- subset(anno_trx_and_exon_all, as.character(chr)==names(base.paired)[i] )
    sub_anno <- unique(sub_anno)
  
 
    base_hap_to_gene <- apply(base_anno_trx, 1, function(x){
        message( paste0( bamChrs[i], ": ", as.character(x[2]) ) )

        tmp <- rbind( subset(base_hap, as.numeric(hap_first_start) >= as.numeric(x[3]) & as.numeric(hap_first_start) <= as.numeric(x[4]) ),
                      subset(base_hap, as.numeric(hap_first_end) >= as.numeric(x[3]) & as.numeric(hap_first_end) <= as.numeric(x[4]) ),
                      subset(base_hap, as.numeric(hap_last_start) >= as.numeric(x[3]) & as.numeric(hap_last_start) <= as.numeric(x[4]) ),
                      subset(base_hap, as.numeric(hap_last_end) >= as.numeric(x[3]) & as.numeric(hap_last_end) <= as.numeric(x[4]) )  )
        tmp <- unique(tmp)


        tmp_sub <- subset(sub_anno, gene_name==as.character(x[2]) )
	      tmp_exclude <- subset(excludeOverlappingExons, chr==as.character(x[1]) & gene_name==as.character(x[2]) )
	      tmp_exclude <- unique(tmp_exclude[,1:4])


	      if(nrow(tmp_sub) >0 & nrow(tmp) >= minDP_per_gene ) {
		            ft_anno <- c()
          
            		for(e in 1:nrow(tmp_sub)){
                		tmp2 <- rbind( subset(tmp, as.numeric(hap_first_start) >= as.numeric(tmp_sub$start)[e] & hap_first_start <= as.numeric(tmp_sub$end)[e] ),
                               			subset(tmp, as.numeric(hap_first_end) >= as.numeric(tmp_sub$start)[e] & hap_first_end <= as.numeric(tmp_sub$end)[e] ),
                               			subset(tmp, as.numeric(hap_last_start) >= as.numeric(tmp_sub$start)[e] & hap_last_start <= as.numeric(tmp_sub$end)[e] ),
                               			subset(tmp, as.numeric(hap_last_end) >= as.numeric(tmp_sub$start)[e] & hap_last_end <= as.numeric(tmp_sub$end)[e] ) )
                		tmp2 <- unique(tmp2)
				            tmp2 <- subset(tmp2, insert_size > 10 | insert_size < (-10) )                

                		if(nrow(tmp2) >0 ){
                     			a1 <- subset(tmp2, as.numeric(NH)==1)
                     			a2 <- subset(tmp2, as.numeric(NH)>1)
 					
                     			tmp2 <- a1
                     			if(nrow(a2) > 0){
						#a2 <- ddply(a2, c("FragmentID"), function(x2)x2[which.min(x2$hap_first_start),] )	# too slow
						a3 <- group_by(a2, as.character(FragmentID), package = "dplyr" )
						a3$insert_size <- abs(a3$insert_size)
						a3 <- filter(a3, insert_size==min(insert_size))
  						a3 <- filter(a3, hap_first_start==min(hap_first_start))
						a3 <- filter(a3, hap_first_end==min(hap_first_end))
						a3 <- filter(a3, hap_last_end==min(hap_last_end))
						a3 <- as.data.frame(a3)

                           			tmp2 <- rbind(a1, a3[,1:10])   
               	      			}
               
                  			tmp2 <- data.frame(tmp2, gene_name=as.character(tmp_sub$gene_name[1])  )
                  			tmp2 <- unique(tmp2)
                  			ft_anno <- rbind(ft_anno, tmp2)
               	 		}
            		}
            		ft_anno <- unique(ft_anno)
     
                # exclude overmapping exons
             		if( nrow(tmp_exclude)>0  & !is.null (ft_anno)  ){	#LY6G5B
                       		#message( paste0( "Excluding exons overlapping in different genes: ", as.character(x[2]), " in ", as.character(x[1]), ", ", Sys.time() )  )
                       		for(ec in 1:nrow(tmp_exclude)){
                               		ft_anno <- subset(ft_anno, !( hap_first_start >= as.numeric(tmp_exclude$start[ec]) & hap_first_start <= as.numeric(tmp_exclude$end[ec] ) ) )
                       		}
             		}

	            return(ft_anno)	
        }
    })
    
      base_hap_to_gene <- do.call(rbind.data.frame, base_hap_to_gene)
      base_hap_to_gene <- unique(base_hap_to_gene)


#       # filtering reads mapped to different genes: this filteres all reads mapped to different exons in the same gene and makes reads depth too low... 
       a <- data.frame(table(base_hap_to_gene$FragmentID ))
       a <- subset(a, Freq==1)
       base_hap_to_gene <- base_hap_to_gene[as.character(base_hap_to_gene$FragmentID) %in% as.character(a$Var1), ]

    
      reads_to_genes[[i]] <- base_hap_to_gene
  }

  names(reads_to_genes) <- bamChrs
  message( paste0( "Done: ", Sys.time() )  )

 
  




  message( paste0( "Combine haplotypes by read: ", Sys.time() )  )
  
  comb_mapping <- c()
  for(i in 1:n_bams){
    comb_mapping <- c(comb_mapping, as.character( reads_to_genes[[i]]$FragmentID)  )
  }
  comb_mapping <- unique(comb_mapping)
  comb_mapping <- data.frame(FragmentID=comb_mapping)
  
  
 
  
  haps.combine_all_reads <- list()
  haps.combine_all_reads <- foreach(i=1:n_bams) %dopar% {
    b <- reads_to_genes[[i]]
    b <- b[,c(1,4,5,11)]
    b <- unique(b)
    tmp <- merge(comb_mapping, b, by="FragmentID", all=TRUE)
    tmp <- data.frame(tmp, bamNames[i])
    names(tmp) <- c("FragmentID", "NM", "NH", "gene_name", "Hap")
    
    haps.combine_all_reads[[i]] <- tmp
  }
  names(haps.combine_all_reads) <- bamChrs
  
  haps.combine_all_reads <- do.call(rbind.data.frame, haps.combine_all_reads)

  message( paste0("Done: ", Sys.time() ) ) 

   #write.table(haps.combine_all_reads, paste0(output_name, ".haps.combine_all_reads.tmp.txt"), sep="\t", quote=F, row.names=F )




  message( paste0("Filtering: ", Sys.time() ) ) 


  ########################################################################### 
  #
  # Filtering reads 
  # 1) Choose reads with minimum ED across haplotypes
  # 2) Clean HLA-DRB*
  # 3) OR11A1 OR10C1 reads in other haps goes to OR11A1 in pgf.
  #    # if reads in OR11A1 in pgf mapped to OR10C1 in other haplotypes, descard reads in OR11A1 in pgf
  #
  #
  ############################################################################

 
  # 1) choose reads with minimum ED across haplotypes
  ft_1 <- haps.combine_all_reads
  ft_1 <- subset(ft_1, NH==1)		# only uniquely mapped reads
  ft_1$NM[is.na(ft_1$NM)] <- 100
  ft_1 <- group_by(ft_1, as.character(FragmentID) )
  ft_1 <- filter(ft_1, NM==min(NM))
  ft_1 <- data.frame(ft_1)
  ft_1 <- ft_1[,1:5]


  

   test <- data.frame(table(ft_1[,c(1,4)]))
   test <- subset(test, Freq >0)

   test2 <- data.frame(table(test[,1]))
   dim(test2)
   test2 <- subset(test2, Freq>1)
   dim(test2)

   #test3 <- ft_1[ft_1$FragmentID %in% test2$Var1, ]
   #data.frame(table(test3$gene_name) )
   ft_1 <- ft_1[ !( ft_1$FragmentID %in% test2$Var1 ), ]


  #write.table(ft_1, paste0("testing_", output_name, "_ft_1.txt"), sep="\t", quote=F, row.names=F)  ### <- first need to get this list.. make an another script for this..





    ##########################################
    ######## get mapping rates from pairs
    ##########################################

    ed_all <- ft_1
    #write.table(ed_all, paste0(sample_name, ".ed_all.txt"), sep="\t", quote=F, row.names=F)


	    message(paste0("Mapping rates from pairs of haplotypes: ", Sys.time() ) )

        output <- getMappingRatesFromPairs(ed_all, editing_distance=ed_dis, minDP_per_gene=30, anno="hg38", sample_name=output_name, read_len=read_length)
        write.table(output$selected_genes, paste0(sample_name, ".selected_genes.txt"), sep="\t", quote=F, row.names=F)
        write.table(output$max_comb, paste0(sample_name, ".max_comb.txt"), sep="\t", quote=F, row.names=F)
  


	    message(paste0("Weighting reads:", Sys.time() ) )


      output2 <- getWeightReads(ed_all, minDP_per_gene=30)
      write.table(output2$ed_summary, paste0(sample_name, ".ed_summary.txt"), sep="\t", quote=F, row.names=F)
      write.table(output2$mapping_stats_to_genes, paste0(sample_name, ".mapping_stats_to_genes.txt"), sep="\t", quote=F, row.names=F)
      write.table(output2$rates_of_weighed_reads_by_genes, paste0(sample_name, ".rates_of_weighed_reads_by_genes.txt"), sep="\t", quote=F, row.names=F) 
      
      
      message("Permutation test...")
      
      overMax <- c()
      for( i in 1:nrow( output2$rates_of_weighed_reads_by_genes )){
                print(as.character(output2$rates_of_weighed_reads_by_genes$gene_name)[i])
                output3 <- permutation_test(bamNames, output2$rates_of_weighed_reads_by_genes[i,], p_number=10000)
                #output3 <- output3[order(as.character(output3$hap) ),]
                
                overMax <- rbind(overMax,
                                 data.frame(gene_name=as.character(output2$rates_of_weighed_reads_by_genes$gene_name)[i],  
                                            output3) )
      }
      write.table(overMax, paste0(sample_name, ".permutation_test.txt"), sep="\t", quote=F, row.names=F) 
      
      
      
      

      if(save_mapping_info==1){
              mapping_rst <- list(mapping_stat=mapping_stat,
                                #editingDistances=ed,     ### only for testing, need to drop it.. it's too big 
                                editingDistances_all=ed_all,
				                        editingDistances_multiHits=haps.reads_to_MultiGenes,
				                        overlappingExons=excludeOverlappingExons)
      }
      if(save_mapping_info!=1){
            mapping_rst <- list(mapping_stat=mapping_stat)
      }


      return(mapping_rst)
      message(paste0("Finished:", Sys.time() ) )

}




#' To read multiple bam files
#' 
#' @param hapNames haploytpe names
#' @param weighted_gene output of 'predictRecombHapByGenes' (*rates_of_weighed_reads_by_genes.txt)
#' @param p_number permutation number
#' @export
#' @importFrom data.table data.table
#' @import Rsamtools
#' @import GenomicAlignments
#' @examples
#' output <- permutation_test(hapNames, weighted_gene, p_number=100000)
#' 
#' # Output data
#' # P-values of mapping results in all pairs of haplotypes by permutaion

permutation_test <- function(hapNames, weighted_gene, p_number=100000){
  
  
      colMax <- function(data) sapply(data, max, na.rm = TRUE)
  
      hapNames_ori <- hapNames
      weighted_gene_ori <- weighted_gene[,hapNames_ori]
      n_haps_ori <- length(hapNames_ori)
      
      weighted_gene <- weighted_gene_ori[which(weighted_gene_ori!=0)]
      hapNames <- colnames(weighted_gene)
      n_haps <- length(hapNames)
      
      
      ##### making matrix
      if(n_haps>=2){
            tmp <- t(combn(hapNames,2))
            b <- apply(tmp, 1, function(x){
                    x <- (weighted_gene[grep(paste0("^", x[1], "$"), names(weighted_gene))] ) + ( weighted_gene[grep(paste0("^", x[2], "$"), names(weighted_gene))] )
            })
            b <- unlist(b)
            tmp <- cbind(b, tmp)
            rownames(tmp) <- NULL
            tmp[tmp[,1]==max(b), ]
      
      
            # matrix structure
            vals<-sort(hapNames)
            nm<-matrix(NA, nrow=length(vals), ncol=length(vals), dimnames=list(vals, vals))
      }
      
      
      
      if(n_haps==2){
              nm[1,1] <- as.numeric(weighted_gene[,hapNames])[1]
              nm[2,2] <- as.numeric(weighted_gene[,hapNames])[2]
              nm[1,2] <- sum(as.numeric(weighted_gene[,hapNames]))
              nm[2,1] <- sum(as.numeric(weighted_gene[,hapNames]))
              #colMax(data.frame(nm))
              
              ##### sampling for 2 haplotypes
              shuffled_max_mean <- c()
              for(i in 1:p_number){
                #print(i)
                sf<-matrix(NA, nrow=length(vals), ncol=length(vals), dimnames=list(vals, vals))
                shuffled_vals <- sample(b, (n_haps * n_haps - n_haps)/2 , replace=TRUE)
                shuffled_diag <- as.numeric( sample(weighted_gene, n_haps, replace=TRUE ) )
                sf[1,1] <- as.numeric(shuffled_vals)
                sf[2,2] <- as.numeric(shuffled_vals)
                sf[1,2] <- shuffled_diag[1]
                sf[2,1] <- shuffled_diag[2]
                shuffled_max_mean <- c(shuffled_max_mean, max(colMeans(sf) ) )
              }
      }
      
      if(n_haps >2) {
            nm[as.matrix(tmp[, 2:3])] <- as.numeric(as.character(tmp[,1]))          # fill
            nm[as.matrix(tmp[, 3:2])] <- as.numeric(as.character(tmp[,1]))          #symmetric
            diag(nm)<- as.numeric(weighted_gene[,hapNames])
            nm_colMeans <- colMeans(nm)
            #max(colMeans(nm) )
            #colMax(data.frame(nm))

            ##### sampling for over 3 haplotypes
            shuffled_max_mean <- c()
            for(i in 1:p_number){
                  #print(i)
                  sf<-matrix(NA, nrow=length(vals), ncol=length(vals), dimnames=list(vals, vals))
                  shuffled_vals <- sample(b, (n_haps * n_haps - n_haps)/2 , replace=TRUE)
                  shuffled_diag <- as.numeric( sample(weighted_gene, n_haps, replace=TRUE ) )
                  sf[as.matrix(tmp[, 2:3])] <- as.numeric(shuffled_vals)
                  sf[as.matrix(tmp[, 3:2])] <- as.numeric(shuffled_vals)
                  diag(sf)<-shuffled_diag
                  shuffled_max_mean <- c(shuffled_max_mean, max(colMeans(sf) ) )
            }
            max(shuffled_max_mean)
      }
      
      
      hap_wi_no_mapping <- hapNames_ori[!(hapNames_ori %in% hapNames)]
      if(n_haps >=2){
            overMax <- c()
            for(i in 1:n_haps){
                  x <- as.data.frame(table(shuffled_max_mean > colMeans(nm)[i] ) ) 
                  t1=ifelse(length(x[x$Var1=="TRUE",c("Freq")]), x[x$Var1=="TRUE",c("Freq")], 0)
                  t2=ifelse(length(x[x$Var1=="FALSE",c("Freq")]), x[x$Var1=="FALSE",c("Freq")], 0 )
                  #print(t1)
                  #print(t2)
                  overMax <- rbind(overMax, data.frame(bigger=t1, smaller=t2 ) )
            }
            overMax$hap <- hapNames
            overMax$pvalues <- overMax$bigger / p_number

            if( length(hap_wi_no_mapping)>0 ){
                  overMax <- rbind(overMax,
                                  data.frame(bigger=p_number, smaller=0, hap=hap_wi_no_mapping, pvalues=1))
            }
      }
      
      if(n_haps==1){
            overMax <- rbind(data.frame(bigger=0, smaller=p_number, hap=hapNames, pvalues=0),
                             data.frame(bigger=p_number, smaller=0, hap=hap_wi_no_mapping, pvalues=1))
      }
      
      overMax <- overMax[order(as.character(overMax$hap) ),]
      
      return(overMax)
}




#' Weighting reads by editing distance to reference sequences
#' 
#' @param bamFiles bam files 
#' @param read_to_gene output from 'predictRecombHapByGenes' (*ed_all.txt)
#' @param minDP_per_gene minimum read depth required in gene
#' @export
#' @examples
#' output <- getWeightReads(read_to_gene, minDP_per_gene=30)
#' 
#' # Output data
#' # mapping rates in each haploypte
#' mapping_stat <- output$ed_summary
#' # gene counts in haploytpes
#' uniqueFragmentID <- output$mapping_stats_to_genes
#' # unique mapping rates in all pairs of haploytpes
#' list_pairedReads <- output$rates_of_weighed_reads_by_genes


getWeightReads <- function(read_to_gene, minDP_per_gene=30){
  
    weightED <- function(ed, i){
          X <- unique(ed)
          N <- length(X)
          i <- i-1
          ord <- which(ed == sort(X, partial = N-i)[N-i])
          return(ord)
    }
  

  
  all_genes <- unique(as.character(read_to_gene$gene_name) )
  all_genes <- all_genes[!is.na(all_genes)]
  
  all_genes_ed <- c()
  mapping_stats_to_genes <- c()      # add statistics... : total number of reads mapped to each gene...  
  for(g in all_genes){
    message(g)
    
    #t1 <- subset(tmp.combine_all_reads, as.character(gene_name)==g & as.numeric(NH)==1)
    t1 <- subset(read_to_gene, as.character(gene_name)==g )
    # a <- t1
    # a$ori <- gsub(".*rt\\d+_\\d+_.*_(.*)_\\d+_\\d+_\\d+_", "\\1", a$FragmentID)
    
    #subset(haps.combine_all_reads, FragmentID=="ENST00000307137_qbl_qbl_rt1_1_qbl_HLA-DRB3_100_275_128_")
    #GL000255.2.7032 ENST00000307137_qbl_qbl_rt1_1_qbl_HLA-DRB3_100_275_128_  0  1  HLA-DRB1  qbl
    #GL000255.2.7033 ENST00000307137_qbl_qbl_rt1_1_qbl_HLA-DRB3_100_275_128_  0  1  HLA-DRB3  qbl
    
    #gene_name       apd        cox        dbb      mann mcf      pgf        qbl     ssto
    #1  HLA-DRB1 0.3858184 0.07294606 0.08680229 0.3924232   0 0.202395 0.05607165 0.204884
    
    
    if(nrow(t1) > minDP_per_gene){
      f_freq <- length(unique(as.character(t1$FragmentID )))
      t1 <- as.data.frame(table(t1[,c("NM", "Hap")]), stringsAsFactors=FALSE )
      
      mapping_stats_to_genes <- rbind(mapping_stats_to_genes, 
                                      data.frame(gene_name=g, NM=NA, Hap="all", Freq=f_freq), 
                                      data.frame(gene_name=g, t1) )
      
      t1$nm2 <- 1
      i_num=length(unique(as.numeric(t1$NM )))
      
      if(i_num >1){
          # m_rate=1/i_num
          m_rate=1/max(as.numeric(t1$NM ) )
        
          for(i in i_num:1){
              ord <- weightED(ed=as.numeric(t1$NM), i)
              t1[ord,]$nm2 <- as.numeric(t1[ord,]$nm2) - (as.numeric(t1[ord,]$NM) * m_rate)
          }
        
          t1$ed= as.numeric(t1$Freq) * as.numeric(t1$nm2)
        
          t2 <- aggregate(ed ~ Hap, FUN = sum, data=t1)
          t2$ed = t2$ed/f_freq 
          t1 <- data.frame(gene_name=g, t(t2[,2]) )
          colnames(t1) <- c("gene_name", t( t2[,1] ) )
        
          all_genes_ed <- rbind(all_genes_ed, t1)
      }
      
      if(i_num ==1){
            tmp <- t1$Hap
            t1 <- data.frame(gene_name=g, t( t1$Freq/f_freq) )
            colnames(t1) <- c("gene_name", tmp )
            all_genes_ed <- rbind(all_genes_ed, t1)
      }       
    }
  }
  all_genes_ed<- all_genes_ed[order(as.character(all_genes_ed$gene_name)),]
  

  
  # only considering best editing distance ...
  
  bamNames=unique(as.character(read_to_gene$Hap))
  n_bams = length(bamNames)
 
  
  message(paste0( "Weighting reads for ", length(all_genes) , " genes: "),  Sys.time() )
  
  hap_combn <-  t(combn(bamNames,2))
  
  
  haps.weighed_reads_by_genes <- list()
  haps.weighed_reads_by_genes <- foreach(g=1:length(all_genes) ) %dopar% {
    
      gene <- as.character(all_genes[g])
      print(gene)
      #reads_in_a_gene <- subset(tmp.combine_all_reads, gene_name==gene, c(1,5:(n_bams+4))  ) 
    
    
      tmp <- subset(read_to_gene, as.character(gene_name)==gene)
    
    if(nrow(tmp) >=minDP_per_gene) {
      
        reads_in_a_gene <- data.frame(FragmentID=unique(as.character(tmp$FragmentID) ) )
      
        for(b in bamNames){
              tmp2 <- subset(tmp, as.character(Hap)==as.character(b) )
              #print(b)
              #print(table(tmp2$NH) )
              reads_in_a_gene <- merge(reads_in_a_gene, tmp2[,c(1:2)], by="FragmentID", all=TRUE)
              colnames(reads_in_a_gene)[ncol(reads_in_a_gene)] <- as.character(b)
        }
      
      
      get_the_same_hap <- c()
      for(i in bamNames){
            h1 <- reads_in_a_gene[,grep(i, colnames(reads_in_a_gene))]
            h1[is.na(h1)] <- 100
            
            get_the_same_hap_ <- c(i)
            for(i2 in bamNames){
                  h2 <- reads_in_a_gene[,grep(i2, colnames(reads_in_a_gene))]
                  h2[is.na(h2)] <- 100
                  if(all(h1==h2) & i!=i2){
                          get_the_same_hap_ <- c(get_the_same_hap_, i2)
                  }
            }
            if(length(get_the_same_hap_)>1){
                  get_the_same_hap_ <- get_the_same_hap_[order(as.character(get_the_same_hap_))]
                  get_the_same_hap_ <- paste(get_the_same_hap_, collapse=";")
                  get_the_same_hap <- rbind(get_the_same_hap, get_the_same_hap_)
            }
      }
      get_the_same_hap <- unique(get_the_same_hap)
      
      
      # to filter genes with the same edit distance across all haplotypes      
      
      tmp2 <- as.data.frame( t( reads_in_a_gene[,2:(n_bams+1)] ) )
      uniquelength <- sapply(tmp2, function(x) length(unique(x)))
      #tmp_sameED <- tmp[uniquelength==1,]
      reads_in_a_gene <- reads_in_a_gene[uniquelength!=1,]
      
      # to exclude the same haploytpe
      reads_in_a_gene_noDup <- reads_in_a_gene[!duplicated(as.list(reads_in_a_gene))]
      
      
      # to filter genes with the same edit distance across all haplotypes
      #a <- ifelse( reads_in_a_gene[,2]==reads_in_a_gene[,3:(n_bams+1)], 0,1 )
      #reads_in_a_gene <- reads_in_a_gene[rowSums(a) >0,]
      
     
      if(nrow(reads_in_a_gene_noDup) >0) {
        
            reads_in_a_gene_noDup <- reads_in_a_gene_noDup[!duplicated(reads_in_a_gene_noDup$FragmentID), ]
        
            #dat <- reads_in_a_gene[,2:(n_bams+1)]
            dat <- reads_in_a_gene_noDup[,2:(ncol(reads_in_a_gene_noDup))]
        
            weighed_reads <- apply( dat, 1, function(x) {
                ed <- as.numeric(x)
                ed[is.na(ed)] <- 100
                na_val <- which(ed==100)
                i_num <- length(unique(ed))
          
                for(i in i_num:1){
                    ord <- weightED(ed, i)
                    if(i==i_num){
                        weight <- 1/length(ord)
                        x[ord] <- weight
                    }
                    else{
                        x[ord] <- 0
                    }
                }
                x[na_val] <- 0
              x
            })
            weighed_reads <- as.data.frame(t(weighed_reads))
      }
      

      if(nrow(reads_in_a_gene_noDup) >0 & !is.null(get_the_same_hap)) {
            for(s in 1:nrow(get_the_same_hap)) {
                  identical_haps <- unlist(strsplit(get_the_same_hap[s,1], ";" ) )
                  to_copy <- colnames(weighed_reads)[colnames(weighed_reads) %in% identical_haps]
                  to_be_copied <- identical_haps[!(identical_haps %in% to_copy)] 
                  
                  #weighed_reads[,c(to_copy)]
                  weighed_reads[,c(to_be_copied)] <- weighed_reads[,c(to_copy)]
            }
      }
        
      if(nrow(reads_in_a_gene_noDup) ==0) {
            reads_in_a_gene_noDup <- data.frame(FragmentID=unique(as.character(tmp$FragmentID) ) )
        
            weighed_reads <- rep( 1/n_bams, n_bams )
            names(weighed_reads)  <- bamNames
            weighed_reads <- as.data.frame(weighed_reads)
            weighed_reads <- rep(weighed_reads, nrow(reads_in_a_gene_noDup))
            weighed_reads <- do.call(rbind.data.frame, weighed_reads)
            rownames(weighed_reads) <- c(1:nrow(reads_in_a_gene_noDup))
            names(weighed_reads)  <- bamNames
      }
      
        weighed_reads <- weighed_reads[,order(colnames(weighed_reads))]
      
        
          for(i in 1:nrow(hap_combn) ){
                tmp_ <- cbind( weighed_reads[,c(hap_combn[i,1])] , weighed_reads[,c(hap_combn[i,2])] )
                #tmp_ <- apply(tmp_, 1, max)
                tmp_ <- apply(tmp_, 1, sum)
          
                weighed_reads <- cbind( weighed_reads, tmp_)
                colnames(weighed_reads)[ncol(weighed_reads)] <- paste0( hap_combn[i,1] , "_", hap_combn[i,2])  
          }
        
          haps.weighed_reads_by_genes[[g]]  <- data.frame(FragmentID=reads_in_a_gene_noDup$FragmentID, weighed_reads)
    }
  }
  names(haps.weighed_reads_by_genes) <- all_genes
  message(paste0("Done:", Sys.time() ) )
  
  
  
  
  message("Get rates of weighted reads by genes...")
  haps.rates_of_weighed_reads_by_genes <- c()
  for(i in 1:length(haps.weighed_reads_by_genes)){
      gene <- names(haps.weighed_reads_by_genes)[i]
      #print(gene)
      tmp <- haps.weighed_reads_by_genes[[i]]
      if(!is.null(tmp)){
            tmp <- round( (colSums(tmp[,2:ncol(tmp)])/nrow(tmp) ), 5)  
            haps.rates_of_weighed_reads_by_genes <- rbind(haps.rates_of_weighed_reads_by_genes,
                                                          data.frame(gene_name=gene, t(tmp) ) )
      }
  }
  
  
  
  rst <- list( ed_summary=all_genes_ed,
               mapping_stats_to_genes=mapping_stats_to_genes,
               rates_of_weighed_reads_by_genes=haps.rates_of_weighed_reads_by_genes)
  return(rst)
  
}




#' To get mapping rates in all pairs of haplotypes and produce summerized figure 
#' 
#' @param ed output of 'predictRecombHapByGenes'
#' @param editing_distance number of mismatch
#' @param minDP_per_minimum read depth in genes
#' @param anno gtf version (ex, hg38)
#' @param sample_name output file name
#' @param read_length read length in bam files
#' @export
#' @import gplots
#' @import ggplot2
#' @import dplyr
#' @examples
#' output <- getMappingRatesFromPairs(ed, editing_distance=10, minDP_per_gene=30, anno="hg38", sample_name, read_length=50)
#' 
#' # Output data
#' # *mhc_heatmap.pdf : heatmap of selected haplotypes, a summerized figure of haplotype-specific mapping
#' # mapping rates and gene counts in all pairs of haplotypes
#' max_comb <- output$max_comb
#' # Mapping ranks in all pairs of haplotypes
#' mapping_ranks <- output$selected_genes


getMappingRatesFromPairs <- function(ed, editing_distance=10, minDP_per_gene=30, anno="hg38", sample_name, read_length=50){
     
        genes <- unique(as.character(ed$gene_name))
        genes <- genes[!is.na(genes)]
        genes <- data.frame(genes=genes)
        
        haps <- unique(as.character(ed$Hap) )
        hap_combn <-  as.data.frame( t(combn(haps,2)) )
        colnames(hap_combn) <- c("Hap1", "Hap2")
        hap_combn <- rbind(hap_combn, data.frame(Hap1=haps, Hap2=haps))
    
	      #genes <-  data.frame( genes=genes[grep("HLA", genes$genes), ] )
        #genes <- data.frame( genes=genes[1:10,] )
        
 
        # Get gene order
        hap_names <- c("apd", "cox", "dbb", "mann", "mcf", "pgf", "qbl", "ssto")
        hap_chrs<- c("GL000250.2", "GL000251.2", "GL000252.2", "GL000253.2", "GL000254.2", "chr6", "GL000255.2", "GL000256.2")
        haps <- data.frame(hap_names, hap_chrs)
        
        if(anno=="hg38"){          
              annoGenes <- read.table("/well/jknight/Wan/pipl/haplotype-mapping/summarizingAnno/annoGenes_ordered.txt", sep="\t", header=T)
              #annoGenes <- read.table("/Users/wanlee/Downloads/add_to_HapTag/annoGenes_ordered.txt", sep="\t", header=T)
       }
        

        z <- annoGenes[ (gsub( "(.*)_.*", "\\1", as.character(annoGenes$gene_from) ) %in% as.character(genes$genes) ) | (gsub( "(.*)_.*", "\\1", as.character(annoGenes$gene_to) ) %in% as.character(genes$genes) ),]  
                
        z$gene_from <- gsub("(.*)_.*", "\\1", z$gene_from)
        z$gene_to <- gsub("(.*)_.*", "\\1", z$gene_to)
        


          #########################################################
          ###
          ### get ordered genes in the MHC region
          ###
          #########################################################

        ordered_genes <- c()
        for(i in unique(as.character(z$hap_from)) ){
          hn <- subset(z, as.character(hap_from)==i & as.character(hap_to)==i)
          if(nrow(hn) >0 ){
            
            hn <- rbind( data.frame(gene_name=as.character(hn$gene_from), ord_bu=as.numeric(hn$ord_from) ),
                         data.frame(gene_name=as.character(hn$gene_to), ord_bu=as.numeric(hn$ord_to))  )
            hn <- hn[!is.na(hn$gene_name),]
            hn <- unique(hn)
            hn <- hn[order(as.numeric(hn$ord_bu)),]
            hn$ord <- c(1:nrow(hn))
            
            a0 <- as.character(subset(hn, ord==min(ord) )$gene_name )
            ordered_genes <- rbind(ordered_genes, 
                                   data.frame(gene_from="Start", gene_to=a0, hap_from="Start", hap_to=i, ord_from=0, ord_to=1) )
            
            for(i2 in 1: (length(hn$ord)-1) ){
              g1=hn$ord[i2]
              g2=hn$ord[i2+1]
              
              a1 <- as.character(subset(hn, ord==g1 )$gene_name )
              a2 <- as.character(subset(hn, ord==g2 )$gene_name )
              tmp <- data.frame(gene_from=a1, gene_to=a2, hap_from=i, hap_to=i, ord_from=g1, ord_to=g2)
              tmp <- unique(tmp)
              tmp <- subset(tmp, as.character(gene_from)!=as.character(gene_to) )
              #print(tmp)
              ordered_genes <- rbind(ordered_genes, tmp)
            }
            
            a3 <- as.character(subset(hn, ord==max(hn$ord) )$gene_name )
            ordered_genes <- rbind(ordered_genes, 
                                   data.frame(gene_from=a3, gene_to="End", hap_from=i, hap_to="End", 
                                              ord_from=max(hn$ord), ord_to=(max(hn$ord) +1)) )
          }
        }
        #ordered_genes <- subset(ordered_genes, gene_to!="End")
        


        num_genes_in_haps <- aggregate(ord_from ~ hap_from, data = ordered_genes, max)
        h1 <- as.character(num_genes_in_haps[which(num_genes_in_haps$ord_from==max(num_genes_in_haps$ord_from) ),]$hap_from)
        h2 <- as.character (haps$hap_names [!(haps$hap_names %in% c("apd", h1) )] )

        gene_order <- subset(ordered_genes, hap_from==h1)
        gene_order <- data.frame(gene_from=gene_order$gene_from, ord_from=gene_order$ord_from) 
        for(i in h2){
  
            tmp <- subset(ordered_genes, hap_from==i)
            tmp <- tmp[ !( tmp$gene_from %in% gene_order$gene_from), ]
  
            if(nrow(tmp) >0){
                  #print(tmp)
                  for( missing_g in 1:nrow(tmp) ){
                            print(missing_g)
                            print( as.character( tmp$gene_from[missing_g] ) )
                            
                            ord= gene_order[grep(as.character(tmp$gene_to[missing_g]), gene_order$gene_from),]$ord_from
                            if(length(ord) ==0){
                                    print( as.character( tmp$gene_from[missing_g] ) )
                            }
                            if(length(ord) ==1){
                                  gene_order[which(gene_order$ord_from>=ord),]$ord_from <- as.numeric( gene_order[which(gene_order$ord_from>=ord),]$ord_from +1 )
                                  gene_order <- rbind(gene_order, 
                                                      data.frame(gene_from=as.character( tmp$gene_from[missing_g] ),
                                                                ord_from=ord) )
                                  gene_order <- gene_order[order(gene_order$ord_from),]
                            }
      
                  }
              }   
      }
      gene_order=aggregate(ord_from ~ gene_from, data = gene_order, max)
      gene_order <- gene_order[order(gene_order$ord_from),]
      gene_order$ord=c(1:nrow(gene_order))



      ### filtering reads mapped to different genes in different haplotypes
      difGInDifHap <- data.frame(table(ed[, c("FragmentID", "gene_name")] ) )
      difGInDifHap <- subset(difGInDifHap , Freq >0)
      difGInDifHap <- data.frame(table(difGInDifHap[,"FragmentID"]))
      difGInDifHap <- subset(difGInDifHap, Freq >1)
      difGInDifHap <- merge(difGInDifHap, ed[, c("FragmentID", "gene_name")], by.x="Var1", by.y="FragmentID")
      difGInDifHap <- difGInDifHap[,c("Var1", "gene_name")]
      difGInDifHap <- unique(difGInDifHap)
      ed2 <- ed[!(as.character(ed$FragmentID) %in% as.character(difGInDifHap$Var1)), ]
      

      
      
        ############################################
        ###
        ### Get mapping rates by all possible pairs of haplotypes
        ###
        ############################################
           
      
      
      max_comb <- list()
      max_comb <- apply(genes, 1, function(g){
        
        message(g)
        s <- subset(ed2, as.character(gene_name)==as.character(g) ) 
        s <- subset(s, as.numeric(NM)<=as.numeric(editing_distance) )
        s$NM_FID <- paste0(s$NM, "_", s$FragmentID)
        
        # get reads with only minimum editing distance for each read
        s <- group_by(s, as.character(FragmentID) )
        #s <- group_by(s, as.character(NM_FID) )
        s <- filter(s, NM == min(NM))
        
        total_reads <- length( unique(as.character(s$FragmentID) ) )
        #total_reads <- length( unique(as.character(s$NM_FID) ) )
        
        
        num_reads <- apply(hap_combn, 1, function(x){
          #message(x)
          g1 <- as.character( subset(s, as.character(Hap)==as.character(x[1]) & as.numeric(NM)<=as.numeric(editing_distance) )$FragmentID )
          g2 <- as.character( subset(s, as.character(Hap)==as.character(x[2]) & as.numeric(NM)<=as.numeric(editing_distance) )$FragmentID ) 
          
          #g1 <- as.character( subset(s, as.character(Hap)=="pgf" & as.numeric(NM)<=as.numeric(editing_distance))$FragmentID )
          #g2 <- as.character( subset(s, as.character(Hap)=="qbl" & as.numeric(NM)<=as.numeric(editing_distance))$FragmentID )
          if(length(g1) >0 | length(g2) >0){
              input <- list(g1, g2)
              v <- venn(input, show.plot=FALSE)
              #attr(v, "dimnames")[[1]]
              head(v, n=4)
              s1=0
              s2=0
              s3=0
              s1_nm=0
              s2_nm=0
              s3_nm=0
              attr_hap1 <- attr(v, "intersections")$A
              attr_hap2 <- attr(v, "intersections")$B
              attr_intersect <- attr(v, "intersections")$`A:B`

	      attr_hap1 <- attr_hap1[!(attr_hap1 %in% attr_intersect)]
              attr_hap2 <- attr_hap2[!(attr_hap2 %in% attr_intersect)]
         

              if(length(attr_intersect) >0){
                  s3 <- s[as.character(s$FragmentID) %in% attr_intersect,]
                  s3_1 <- subset(s3, as.character(Hap)==as.character(x[1]) )
                  s3_2 <- subset(s3, as.character(Hap)==as.character(x[2]) )
                  tmp_s3 <- merge(s3_1, s3_2, by="FragmentID", all=TRUE)
                  s3 <- subset(tmp_s3, NM.x==NM.y)
                  s3_nm <- sum(as.numeric(s3$NM.x) )
                  s3=( nrow(s3) * read_length - sum(as.numeric(s3$NM.x) ))/(nrow(s3) * read_length )
                  s3_1 <- subset(tmp_s3, NM.x < NM.y)
                  s3_2 <- subset(tmp_s3, NM.x > NM.y)
                  if(nrow(s3_1) > 0){
                      attr_hap1 <- c(attr_hap1, s3_1$FragmentID)
                  }
                  if(nrow(s3_2) > 0){
                      attr_hap2 <- c(attr_hap2, s3_2$FragmentID)
                  }                                              
              }
            
              if(length(attr_hap1) >0){
                  s1 <- s[as.character(s$FragmentID) %in% attr_hap1,]
                  s1 <- subset(s1, as.character(Hap)==as.character(x[1]) )
                  s1_nm <- sum(as.numeric(s1$NM))
                  s1 <- ( nrow(s1) * read_length - sum(as.numeric(s1$NM)) )/(nrow(s1) * read_length )    
              }
              if(length(attr_hap2) >0) {
                  s2 <- s[as.character(s$FragmentID) %in% attr_hap2,]
                  s2 <- subset(s2, as.character(Hap)==as.character(x[2]) )
                  s2_nm <- sum(as.numeric(s2$NM))
                  s2 <- ( nrow(s2) * read_length - sum(as.numeric(s2$NM) ))/(nrow(s2) * read_length )
             } 
            
            if(as.character(x[1]) == as.character(x[2])){
                  s1=0
                  s2=0
                  s1_nm=0
                  s2_nm=0
                  attr_hap1 <- NULL
                  attr_hap2 <- NULL
            }
            
            all_rate  <- ( (length(attr_hap1) + length(attr_hap2) + length(attr_intersect)) * read_length -( s1_nm + s2_nm + s3_nm) ) / ( (length(attr_hap1) + length(attr_hap2) + length(attr_intersect)) * read_length )
            
            tmp=data.frame(gene_name=as.character(g), Hap1=x[1], Hap2=x[2], 
				only_h1= head(v)[3], only_h2=head(v)[2], intersect=head(v)[4],
                           	h1_mismatch=s1_nm, h2_mismatch=s2_nm, intersect_mismatch=s3_nm,
				only_h1_rate=s1, only_h2_rate=s2, intersect_rate=s3, mapping_rate=all_rate)
            tmp
          }
          else if(length(g1) ==0 & length(g2) ==0){
            tmp=data.frame(gene_name=as.character(g), Hap1=x[1], Hap2=x[2], 
				only_h1= 0, only_h2=0, intersect=0,
				h1_mismatch=0, h2_mismatch=0, intersect_mismatch=0,
                           	only_h1_rate=0, only_h2_rate=0, intersect_rate=0, mapping_rate=0)
            tmp 
          }
        })
        

        if( length(num_reads) >0 ) {
          num_reads <- do.call(rbind.data.frame, num_reads)
          num_reads$gene_count <- rowSums(num_reads[,4:6])
          
          if( max( num_reads$gene_count )>= minDP_per_gene ) {
            
              num_reads$rates_of_gene_counts <- as.numeric(num_reads$gene_count)/total_reads
              
              # num_reads <- num_reads[order(num_reads$gene_count, decreasing=T),]	
              # # changing name: best_pair -> putative_haps
              # num_reads$putative_haps <- as.numeric(factor(num_reads$gene_count,levels=unique(num_reads$gene_count)))
              # num_reads$putative_haps[num_reads$rates_of_gene_counts > 0.99] <- 1
              # #(1-144/145) < 0.1
              
               
              s1 <- subset(num_reads, rates_of_gene_counts >= (max(rates_of_gene_counts)-0.05) )
              s1$type="A"
              s2 <- subset(s1, mapping_rate!=max(mapping_rate) )
              
              s1 <- subset(s1, mapping_rate==max(mapping_rate) )
              s1$putative_haps <- 1
              putative_num <- s1
              
              s3 <- subset(num_reads, rates_of_gene_counts < (max(rates_of_gene_counts)-0.05) )
              
              if(nrow(s2)>0){
                s2 <- s2[order(s2$mapping_rate, s2$rates_of_gene_counts, decreasing=T),]
                s2$putative_haps <- c(2: (nrow(s2) +1) )
                putative_num <- rbind(putative_num, s2)
              }
              if(nrow(s3)>0){
                s3$type="B"
                s3 <- s3[order(s3$rates_of_gene_counts, s3$mapping_rate, decreasing=T),]
                s3$putative_haps <- c( (nrow(putative_num)+1): nrow(num_reads) )
                
                putative_num <- rbind(putative_num, s3)
              }
                   
              rownames(putative_num) <- NULL
                        
              max_comb[[g]] <-  data.frame(gene=as.character(g), putative_num ) 
            
          }
          
        }
      } )
      names(max_comb) <- genes$genes
      max_comb <- do.call(rbind.data.frame, max_comb)
      rownames(max_comb) <- NULL
      max_comb <- max_comb[,2:18]
      
        #max_comb <- max_comb[order(as.character(max_comb$gene)), ]
      
      	#write.table(max_comb, "1.0_1.0_max_comb.txt", sep="\t", quote=F, row.names=F)
      
      
        
        #######################################
        ###
        ### draw a heamap
        ###
        #######################################
        
        input_genes <- unique(as.character(max_comb$gene_name))
        
        #top1rates <- max_comb
        
        # to select haplotypes with only maximum gene count
        top1rates <- c()
        for(g in input_genes){
            s <- max_comb[grep(paste0("_", g, "_"), paste0("_", as.character(max_comb$gene_name), "_") ),]
          
            if( nrow(s) > 0){
                
                s1 <- subset(s, putative_haps==1)

                if(nrow(s1) >=1){
                      s1_1 <- subset(s1, as.numeric(s1$gene_count)==as.numeric(max(s1$gene_count) ) )
                      top1rates <- rbind(top1rates, s1_1)
                      s1_2 <- subset(s1, as.numeric(s1$gene_count)!=max(as.numeric(s1$gene_count) ) )
                      if(nrow(s1_2)>0){
                          s1_2$putative_haps <- 0
                          top1rates <- rbind(top1rates, s1_2)
                      }
                }
                s2 <- subset(s, putative_haps==0)
                top1rates <- rbind(top1rates, s2)
            }    
        }
        
        
        tmp <- merge(gene_order[,c("gene_from", "ord")], top1rates, by.x="gene_from", by.y="gene_name")

        tmp <- unique( tmp[, c("gene_from", "ord")] )
        tmp <- tmp[order(tmp$ord),]
        tmp$ord <- c(1:nrow(tmp))
        tmp <- merge(top1rates, tmp[,c("gene_from", "ord")], by.x="gene_name", by.y="gene_from") 
        
        
        tmp2 <- rbind( data.frame(gene_name=tmp$gene_name,
                                  hap=tmp$Hap1,
                                  ord=tmp$ord,
                                  Mapping_Rates=((tmp$only_h1 + tmp$intersect)/(tmp$gene_count/tmp$rates_of_gene_counts) ),
                                  putative_haps=tmp$putative_haps),
                       data.frame(gene_name=tmp$gene_name,
                                  hap=tmp$Hap2,
                                  ord=tmp$ord,
                                  Mapping_Rates=((tmp$only_h2 + tmp$intersect)/ (tmp$gene_count/tmp$rates_of_gene_counts) ),
                                  putative_haps=tmp$putative_haps) )

        tmp2 <- unique(tmp2)
        tmp2$gene_hap <- paste0(tmp2$gene_name, "_", tmp2$hap)
      
      

        top1rates <- c()
        for(g in unique(as.character(tmp2$gene_hap))){
                s <- subset(tmp2, gene_hap==g)
                if(nrow(s)==1){
                        top1rates <- rbind(top1rates, s)
                }
                if(nrow(s)>1){
                        s <- subset(s, putative_haps==1)
                        top1rates <- rbind(top1rates, s)
                }
        }


        top1rates <- top1rates[order(top1rates$ord),] 
        top1rates$ord <- top1rates$ord * -1

        gene_labels <- unique(top1rates$gene_name)
        gene_labels <- rev(gene_labels)


      
      top1rates_ <- subset(top1rates, putative_haps==1 & Mapping_Rates!=0)
      max_hap=length(unique(top1rates_$hap))
      a <- data.frame(table(top1rates_$gene_name),  stringsAsFactors=FALSE)

      
	    selected_genes <- c()
      # all haplotypes are identical
      top1rates_any <- subset(a, Freq==max_hap)
	    if(nrow(top1rates_any) >0){
      		    top1rates_any <- data.frame( top1rates_[as.character(top1rates_$gene_name) %in% as.character(top1rates_any$Var1),], type="any")
		          selected_genes <- data.frame(top1rates_any, n_hap=max_hap)
      }
      

      # more than 2 haplotypes are selected
      tmp <- subset(a, Freq>2 & Freq < max_hap )
      top1rates_3 <- top1rates_[as.character(top1rates_$gene_name) %in% as.character(tmp$Var1),]
      top1rates_3_1 <- c()
      top1rates_3_2 <- c()
      for(g in unique(as.character(top1rates_3$gene_name))){
            message(g)
            s1 <- subset(top1rates_3, as.character(gene_name)==g)
            s2 <- subset(max_comb, as.character(gene_name)==g & putative_haps==1)
            s2 <- s2[( s2$Hap1 %in% s1$hap ) & (s2$Hap2 %in% s1$hap),]
            if(sum(s2$only_h1)==0 & sum(s2$only_h2)==0){
                    top1rates_3_1 <- rbind(top1rates_3_1, data.frame( s1[1,], type="hm"))
                    top1rates_3_2 <- rbind(top1rates_3_2, data.frame( s1[2:nrow(s2),], type="hm" ) )
            }
            if(sum(s2$only_h1)!=0 & sum(s2$only_h2)!=0){
                      b <- data.frame( table( c( as.character(s2$Hap1), as.character(s2$Hap2) )) )
                      selected_hap=as.character( subset(b, Freq==max(b$Freq) )$Var1 )
                      #print(selected_hap)
                      if(length(selected_hap)==1){
                              selected_hap_=subset(s2, Hap1==selected_hap | Hap2==selected_hap)[1,] 
                      }
                      if(length(selected_hap)>1) {
                              selected_hap_=subset(s2, Hap1==selected_hap[1] )[1,]
                              if(nrow(subset(s2, Hap1==selected_hap[1] ))==0 ){
                                        selected_hap_=subset(s2, Hap2==selected_hap[1] )[1,]
                              }
                      }
                      selected_hap <- selected_hap_
                      top1rates_3_1 <- rbind(top1rates_3_1, 
                                             data.frame( s1[s1$hap %in% c(as.character(selected_hap$Hap1), as.character(selected_hap$Hap2)), ], 
                                                         type="ht") )
                      top1rates_3_2 <- rbind(top1rates_3_2, 
                                             data.frame( s1[! (s1$hap %in% c(as.character(selected_hap$Hap1), as.character(selected_hap$Hap2)) ), ],
                                                         type="ht") )                  
            }
      }
      top1rates_3_2 <- top1rates_3_2[!is.na(top1rates_3_2$gene_name),]
				# 3_selected: first one, 3_not: from 2nd ... basically random
				selected_genes <- rbind(selected_genes,
				                        data.frame(top1rates_3_1, n_hap="3_selected"),
				                        data.frame(top1rates_3_2, n_hap="3_not") )
	
      
      # 2 haplotypes are selected
      tmp <- subset(a, Freq==2)				## need to fix it. not all 2 is hetero... , 
      if(nrow(tmp)>0){
              top1rates_2 <- data.frame( top1rates_[as.character(top1rates_$gene_name) %in% as.character(tmp$Var1),], type="ht")
              selected_genes <- rbind(selected_genes,
                                      data.frame(top1rates_2, n_hap=2) )
              
      }
      # 1  haplotypes is selected
      tmp <- subset(a, Freq==1)
      if(nrow(tmp)>0){
              top1rates_1 <- data.frame( top1rates_[as.character(top1rates_$gene_name) %in% as.character(tmp$Var1),], type="hm" )
              selected_genes <- rbind(selected_genes,
                                      data.frame(top1rates_1, n_hap=1) )
              top1rates_2 <- rbind(top1rates_2, top1rates_1)
      }
      



	  top1rates_3_1_hm <- subset(top1rates_3_1, as.character(type)=="hm")
	  top1rates_3_1_ht <- subset(top1rates_3_1, as.character(type)=="ht")
	
    
	if(nrow(top1rates_any) >0 ){
      		mhc_heatmap <- ggplot(top1rates, aes(x = hap , y = ord, fill = (Mapping_Rates*100) )) + 
                            geom_tile() +
                            scale_fill_gradient() +
                            scale_y_continuous(expand=c(0,0),breaks=min(top1rates$ord):max(top1rates$ord),labels=gene_labels) +
                            xlab("Haplotypes") + ylab("Genes by the chromosomal order") +
                            geom_text(aes(label=round( (Mapping_Rates*100),2) ), size=3, color="white" ) +
                            geom_tile(data=top1rates_any, aes(x=hap, y=as.numeric(ord), color=putative_haps), size=1, fill=NA, color="grey") +
                            geom_tile(data=subset(top1rates_3_1, type=="hm"), aes(x=hap, y=ord, color=putative_haps), size=2,fill=NA, color="blue") +
                            geom_tile(data=subset(top1rates_3_1, type=="ht"), aes(x=hap, y=ord, color=putative_haps), size=1,fill=NA, color="blue") +
                            geom_tile(data=top1rates_2, aes(x=hap, y=ord, color=putative_haps),size=1,fill=NA, color="red") +
                            #geom_tile(data=top1rates_1, aes(x=hap, y=ord, color=putative_haps),size=2,fill=NA, color="red") +
                            ggtitle(sample_name )

	}
  
	if(nrow(top1rates_any) ==0 ){
		    mhc_heatmap <- ggplot(top1rates, aes(x = hap , y = ord, fill = (Mapping_Rates*100) )) +
                            geom_tile() +
                            scale_fill_gradient() +
                            scale_y_continuous(expand=c(0,0),breaks=min(top1rates$ord):max(top1rates$ord),labels=gene_labels) +
                            xlab("Haplotypes") + ylab("Genes by the chromosomal order") +
                            geom_text(aes(label=round( (Mapping_Rates*100),2) ), size=3, color="white" ) +
                            #geom_tile(data=top1rates_any, aes(x=hap, y=as.numeric(ord), color=putative_haps), size=1, fill=NA, color="grey") +
                            geom_tile(data=top1rates_3_1_hm, aes(x=hap, y=ord, color=putative_haps), size=2,fill=NA, color="blue") +
                            geom_tile(data=top1rates_3_1_ht, aes(x=hap, y=ord, color=putative_haps), size=1,fill=NA, color="blue") +
                            geom_tile(data=top1rates_2, aes(x=hap, y=ord, color=putative_haps),size=1,fill=NA, color="red") +
                            #geom_tile(data=top1rates_1, aes(x=hap, y=ord, color=putative_haps),size=2,fill=NA, color="red") +
                            ggtitle(sample_name )

	}



        ggsave(filename=paste0(sample_name, ".mhc_heatmap.pdf" ), mhc_heatmap, width=10, height=24)
      
        output2 <- list(selected_genes=selected_genes,
                       max_comb=max_comb)

	      return(output2)
}

    


#' To predict shortest path
#' 
#' @param paired_mapping_rates output 'predictRecombHapByGenes'  (*max_comb.txt)
#' @param anno gtf version (ex, "hg38")
#' @param penalty 1: penalty applied, 0: penalty not applied
#' @param sample_name output file name
#' @export
#' @import igraph
#' @import dplyr
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import plyr
#' @importFrom data.table data.table
#' @examples
#' output <- getHapsByShortestPaths(paired_mapping_rates, anno="hg38", penalty=1, sample_name)
#' 
#' # Output data
#' # *.shortestPath_heatmap.pdf : heatmap of selected haplotypes
#' # gene and haplotype to draw heatmap
#' heatmap_genes <- output$heatmap_genes
#' # predicted best pairs of haploytpes
#' summary_of_best_pairs <- output$summary_of_best_pairs

getHapsByShortestPaths <- function(paired_mapping_rates, anno="hg38", penalty=1, sample_name){
  
  unique_best_hits <- data.frame(table(paired_mapping_rates[,c("gene_name", "type")]) )
  
  paired_mapping_rates_noFT <-  paired_mapping_rates
  paired_mapping_rates <- subset(paired_mapping_rates, gene_count >= 50)   
  
  
  gene_list <- unique(as.character(paired_mapping_rates$gene_name)) 
  head(paired_mapping_rates)
  
  comb_haps <- unique(paste0( paired_mapping_rates[,c("Hap1")], "_", paired_mapping_rates[,c("Hap2")]) )
  
  v1 <- comb_haps
  v2 <- comb_haps
  tmp <- sort(apply(expand.grid(v1, v2), 1, paste, collapse = ",", sep = "")) 
  hap_combn <- data.frame(haps_from=gsub("(.*),(.*)", "\\1", tmp),
                          haps_to=gsub("(.*),(.*)", "\\2", tmp) )
  
  
  
  ### summary of best_pairs
  
  
  # additional filtering
  a_ <- subset(paired_mapping_rates, (only_h1>5 & only_h2>5 & Hap1!=Hap2)| (only_h1==0 & only_h2==0 & Hap1==Hap2) )       #wo.dp5
  

  a <- rbind(data.frame(gene_name=a_$gene_name, gene_hap=paste0(a_$gene_name, "_", a_$Hap1)),
             data.frame(gene_name=a_$gene_name, gene_hap=paste0(a_$gene_name, "_", a_$Hap2)))
  a <- unique(a)
  a2 <- data.frame(table(a$gene_name))
  
  
  n_hap1 <- a_[as.character(a_$gene_name) %in% as.character(a2[a2$Freq==1,]$Var1),]     # homozygous
  n_hap2 <- a_[as.character(a_$gene_name) %in% as.character(a2[a2$Freq==2,]$Var1),]     # heterozygous
  n_hap_any <- a_[as.character(a_$gene_name) %in% as.character(a2[a2$Freq>2,]$Var1),]   # more than two pairs
  n_hap_any_hm <- subset(n_hap_any, Hap1==Hap2)
  #n_hap_any_ht <- n_hap_any[!(n_hap_any$gene_name %in% n_hap_any_hm$gene_name), ]
  n_hap_any_ht <- subset(n_hap_any, Hap1!=Hap2)
  
  
  selected_genes <- c()
  if(nrow(n_hap1) >0){
    selected_genes <- rbind(selected_genes, data.frame(n_hap1, type="hm") )
  }
  if(nrow(n_hap2) >0){
    selected_genes <- rbind(selected_genes, data.frame(n_hap2, type="ht") )
  }
  if(nrow(n_hap_any_hm) >0){
    selected_genes <- rbind(selected_genes, data.frame(n_hap_any_hm, type="hm_multi") )
  }
  if(nrow(n_hap_any_ht) >0){
    selected_genes <- rbind(selected_genes, data.frame(n_hap_any_ht, type="ht_multi") )
  }
  
  
  
  input_structure <- unique(paste0( paired_mapping_rates[,c("Hap1")], "_", paired_mapping_rates[,c("Hap2")]) )
  input_structure <- data.frame(hap_pairs= input_structure[order(input_structure)])
  
  for(i in gene_list){
    s <- subset(paired_mapping_rates, as.character(gene_name)==i)
    s <- unique(s)
    s$tmp <- paste0(s$Hap1, "_", s$Hap2)
    s$mapping_rate[s$type=="B"] <- 0
    
    #input_structure <- merge(input_structure, s[,c("tmp", "rates_of_gene_counts")], by.x="hap_pairs", by.y="tmp", all=TRUE)
    input_structure <- merge(input_structure, s[,c("tmp", "mapping_rate")], by.x="hap_pairs", by.y="tmp", all=TRUE)
    
    colnames(input_structure)[ncol(input_structure)] <- i
    #print(i)
    #print(dim(input_structure))
  }
  input_structure <- unique(input_structure) 
  input_structure <- t(input_structure)
  colnames(input_structure) <- input_structure[1,]
  input_structure <- input_structure[2:nrow(input_structure),]
  
  
  
  
  #########################################################
  ###
  ### function: Combining two p-values
  ###
  #########################################################
  
  #combineTwoPvalues <- function(x,y) return(pchisq(-2*log(x/2)-2*log(y/2),4, low=F))
  
  
  
  #########################################################
  ###
  ### get annotation
  ###
  #########################################################
  
  hap_names <- c("apd", "cox", "dbb", "mann", "mcf", "pgf", "qbl", "ssto")
  
  if(anno=="hg38"){     
    hap_chrs<- c("GL000250.2", "GL000251.2", "GL000252.2", "GL000253.2", "GL000254.2", "chr6", "GL000255.2", "GL000256.2")
    haps <- data.frame(hap_names, hap_chrs)
    
    annoGenes <- read.table("/well/jknight/Wan/pipl/haplotype-mapping/summarizingAnno/annoGenes_ordered.txt", sep="\t", header=T)
    annoGenes$gene_from <- gsub("VARSL", "VARS2", as.character(annoGenes$gene_from) )
    annoGenes$gene_to <- gsub("VARSL", "VARS2", as.character(annoGenes$gene_to) )
    annoGenes$gene_to <- gsub("C6orf205", "MUC21", as.character(annoGenes$gene_to) )
    annoGenes$gene_to <- gsub("C6orf205", "MUC21", as.character(annoGenes$gene_to) )      
    annoGenes <- subset(annoGenes, as.character(gene_from)!="TRIM39-RPP21_pgf" & as.character(gene_to)!="TRIM39-RPP21_pgf" ) 	
    
    }
  
  
  
  #z <- annoGenes[ (gsub( "(.*)_.*", "\\1", as.character(annoGenes$gene_from) ) %in% as.character(rownames(input_structure)) ) | (gsub( "(.*)_.*", "\\1", as.character(annoGenes$gene_to) ) %in% as.character(rownames(input_structure)) ),]  
  z <- annoGenes
  z$gene_from <- gsub("(.*)_.*", "\\1", z$gene_from)
  z$gene_to <- gsub("(.*)_.*", "\\1", z$gene_to)
  # write.table(z, "gene_path_structure.txt", sep="\t", quote=F, row.names=F)  
  
  getMissingGenes <- unique(z[,c("gene_from", "hap_from")] ) 
  getMissingGenes <- data.frame(table(getMissingGenes$gene_from))
  getMissingGenes <- subset(getMissingGenes, Freq>0 & Freq<8)
  
  
  
  
  #########################################################
  ###
  ### get ordered genes in the MHC region
  ###
  #########################################################
  
  
  
  ordered_genes <- c()
  for(i in unique(as.character(z$hap_from)) ){
    hn <- subset(z, as.character(hap_from)==i & as.character(hap_to)==i)
    if(nrow(hn) >0 ){
      
      hn <- rbind( data.frame(gene_name=as.character(hn$gene_from), ord_bu=as.numeric(hn$ord_from) ),
                   data.frame(gene_name=as.character(hn$gene_to), ord_bu=as.numeric(hn$ord_to))  )
      hn <- hn[!is.na(hn$gene_name),]
      hn <- unique(hn)
      hn <- hn[order(as.numeric(hn$ord_bu)),]
      hn$ord <- c(1:nrow(hn))
      
      a0 <- as.character(subset(hn, ord==min(ord) )$gene_name )
      ordered_genes <- rbind(ordered_genes, 
                             data.frame(gene_from="Start", gene_to=a0, hap_from="Start", hap_to=i, ord_from=0, ord_to=1) )
      
      for(i2 in 1: (length(hn$ord)-1) ){
        g1=hn$ord[i2]
        g2=hn$ord[i2+1]
        
        a1 <- as.character(subset(hn, ord==g1 )$gene_name )
        a2 <- as.character(subset(hn, ord==g2 )$gene_name )
        tmp <- data.frame(gene_from=a1, gene_to=a2, hap_from=i, hap_to=i, ord_from=g1, ord_to=g2)
        tmp <- unique(tmp)
        tmp <- subset(tmp, as.character(gene_from)!=as.character(gene_to) )
        #print(tmp)
        ordered_genes <- rbind(ordered_genes, tmp)
      }
      
      a3 <- as.character(subset(hn, ord==max(hn$ord) )$gene_name )
      ordered_genes <- rbind(ordered_genes, 
                             data.frame(gene_from=a3, gene_to="End", hap_from=i, hap_to="End", 
                                        ord_from=max(hn$ord), ord_to=(max(hn$ord) +1)) )
    }
  }
  #ordered_genes <- subset(ordered_genes, gene_to!="End")
  
  
  
  
  num_genes_in_haps <- aggregate(ord_from ~ hap_from, data = ordered_genes, max)
  h1 <- as.character(num_genes_in_haps[which(num_genes_in_haps$ord_from==max(num_genes_in_haps$ord_from) ),]$hap_from)
  h2 <- as.character (haps$hap_names [!(haps$hap_names %in% c("apd", h1) )] )
  
  
  gene_order <- subset(ordered_genes, hap_from==h1)
  gene_order <- data.frame(gene_from=gene_order$gene_from, ord_from=gene_order$ord_from) 
  for(i in h2){
    
    tmp <- subset(ordered_genes, hap_from==i)
    tmp <- tmp[ !( tmp$gene_from %in% gene_order$gene_from), ]
    
    if(nrow(tmp) >0){
      
      for( missing_g in 1:nrow(tmp) ){
        #print(missing_g)
        #print( as.character( tmp$gene_from[missing_g] ) )
        
        ord= gene_order[grep(as.character(tmp$gene_to[missing_g]), gene_order$gene_from),]$ord_from
        if(length(ord) ==0){
          print( as.character( tmp$gene_from[missing_g] ) )
        }
        if(length(ord) ==1){
          gene_order[which(gene_order$ord_from>=ord),]$ord_from <- as.numeric( gene_order[which(gene_order$ord_from>=ord),]$ord_from +1 )
          gene_order <- rbind(gene_order, 
                              data.frame(gene_from=as.character( tmp$gene_from[missing_g] ),
                                         ord_from=ord) )
          gene_order <- gene_order[order(gene_order$ord_from),]
        }
        
      }
    }   
  }
  
  gene_order <- gene_order[gene_order$gene_from %in% rownames(input_structure), ]
  colnames(gene_order) <- c("gene_from", "ord_tmp")
  gene_order$ord_from=c(1:nrow(gene_order))
  
  gene_order=aggregate(ord_from ~ gene_from, data = gene_order, max)
  gene_order <- gene_order[order(gene_order$ord_from),]
  gene_order$ord=c(1:nrow(gene_order))
  
  
  #################################################################
  ###
  ### get ordered genes with pvlaues from predictRecombHapByGenes
  ###
  #################################################################
  
  
  input_structure <- data.frame(input_structure)
  input_structure$gene_name <- rownames(input_structure)
  mr <- merge(input_structure, gene_order, by.x="gene_name", by.y="gene_from")
  mr <- mr[order(as.numeric(mr$ord)),]
  
  
  t1 <- data.frame(c(rep(0, length(comb_haps)), 0,0) )
  rownames(t1) <- c(comb_haps, "ord_from", "ord")
  t1 <- data.frame(gene_name="Start", t(t1) )
  
  t2 <- data.frame( c(rep(0, length(comb_haps)), rep((max(mr$ord)+1), 2)) )
  rownames(t2) <- c(comb_haps, "ord_from", "ord")
  t2 <- data.frame(gene_name="End", t(t2) )
  
  mr <- rbind(t1, mr, t2)
  
 
  
  
  ### 
  tmp_haps <- data.frame(h1=gsub("(.*)_(.*)", "\\1", colnames(mr)[2:ncol(mr)] ),
                         h2=gsub("(.*)_(.*)", "\\2", colnames(mr)[2:ncol(mr)] ) )
  
  mr2 <- data.frame(haps=colnames(mr)[2:ncol(mr)])
  for(i in 1:nrow(mr)){
    tmp <- data.frame( tmp_haps, t(mr[i,2:ncol(mr)] ) ) 
    tmp$haps <- rownames(tmp)
    
    
    ### bu: to change mapping rate to 0 if both haplotypes are the same.
    a <- subset(tmp, as.character(h1)==as.character(h2) & as.numeric(as.character(tmp[,3]))==1 & h1!="ord")
    if(nrow(a) >0){
      i2 <-  as.character(a$h1)
      b2 <- tmp[as.character(tmp$h1) %in% i2 & as.character(tmp$h2) %in% i2 &  as.character(tmp$h1)!=as.character(tmp$h2), ] 
      b1 <- tmp[ !( as.character(tmp$h1) %in% i2 & as.character(tmp$h2) %in% i2 &  as.character(tmp$h1)!=as.character(tmp$h2) ), ] 
      
      if(nrow(b2) >0){
        b2[,3] <- 0
        b1[,3] <- as.numeric(as.character(b1[,3]))
        tmp <- rbind(b1, b2)
      }
      if(nrow(b2) ==0){
        tmp <- b1
      }

    }
    mr2 <- merge(mr2, tmp[,3:4], by="haps")
  }
  mr2 <- t(mr2)
  colnames(mr2) <- mr2[1,]
  mr2 <- mr2[2:nrow(mr2), ]
  mr2 <- data.frame(mr2)
  mr2=data.frame(gene_name=mr$gene_name, mr2)
  
  mr <- cbind( mr2[, -(grep("ord", colnames(mr2)) ) ], mr2[, grep("ord", colnames(mr2)) ] )
  
  mr_col <- colnames(mr)[2:(ncol(mr)-2)]
  rownames(mr) <- NULL
  
  
  
  
  
  path_frame <- c()
  for(ord in 1:(nrow(mr)-1) ){
    message(paste0(ord, ": ", mr[ord, c("gene_name")]) )
    
    g1 <- as.character(mr[ord, c("gene_name")])
    g2 <- as.character(mr[(ord+1), c("gene_name")])
    o1 <- as.numeric(as.character(mr[ord, c("ord")]))
    o2 <- as.numeric(as.character(mr[(ord+1), c("ord")]))
    
    gfrom <- data.frame(haps=mr_col, mr_from= as.numeric(t( mr[ord,2:(ncol(mr)-2)] )[,1] ) ) 
    gto  <-  data.frame(haps=mr_col,  mr_to=as.numeric(t( mr[(ord+1),2:(ncol(mr)-2)] )[,1] ) ) 
    tmp <- merge(hap_combn, gfrom, by.x="haps_from", by.y="haps")
    tmp <- merge(tmp, gto, by.x="haps_to", by.y="haps")
    tmp <- data.frame(gene_from=g1, gene_to=g2, ord_from=o1, ord_to=o2,  tmp)
    
    path_frame <- rbind(path_frame, tmp)             
  }
  
  path_frame <- path_frame[!is.na(path_frame$mr_from), ]
  path_frame <- path_frame[!is.na(path_frame$mr_to), ]
  
  path_frame$haps_to_1=gsub("(.*)_(.*)", "\\1", path_frame$haps_to)
  path_frame$haps_to_2=gsub("(.*)_(.*)", "\\2", path_frame$haps_to)
  path_frame$haps_from_1=gsub("(.*)_(.*)", "\\1", path_frame$haps_from)
  path_frame$haps_from_2=gsub("(.*)_(.*)", "\\2", path_frame$haps_from)
  
  path_frame$haps_from <- as.character( path_frame$haps_from)
  path_frame$haps_from[path_frame$gene_from=="Start"] <- "Start"
  path_frame$haps_from_1 <- as.character( path_frame$haps_from_1)
  path_frame$haps_from_2 <- as.character( path_frame$haps_from_2)
  path_frame$haps_from_1[path_frame$gene_from=="Start"] <- "Start"
  path_frame$haps_from_2[path_frame$gene_from=="Start"] <- "Start"
  
  path_frame$haps_to <- as.character( path_frame$haps_to)
  path_frame$haps_to[path_frame$gene_to=="End"] <- "End"
  path_frame$haps_to_1 <- as.character( path_frame$haps_to_1)
  path_frame$haps_to_2 <- as.character( path_frame$haps_to_2)
  path_frame$haps_to_1[path_frame$gene_to=="End"] <- "End"
  path_frame$haps_to_2[path_frame$gene_to=="End"] <- "End"
  
  path_frame <- unique(path_frame)
  
  
 
  
  ##############################################
  ### applying penalty
  ##############################################
  
  
  if(penalty==1){
    
    tmp <- path_frame
    unique_best_hits_A1 <- subset(unique_best_hits, type=="A" & Freq==1)
    
    sem <- 0.002783863

    
    tmp1 <- subset(tmp, as.character(haps_to)==as.character(haps_from) )
    tmp1 <- rbind(tmp1, subset(tmp, gene_from=="Start" | gene_to=="End"))
    
    tmp2 <- subset(tmp, as.character(haps_to)!=as.character(haps_from) )
    tmp2 <- subset(tmp2, gene_from!="Start" & gene_to!="End")

    tmp2_one <- subset(tmp2, as.character(haps_to_1)==as.character(haps_from_1) | as.character(haps_to_1)==as.character(haps_from_2) |  as.character(haps_to_2)==as.character(haps_from_1) | as.character(haps_to_2)==as.character(haps_from_2))
    tmp2_both <- subset(tmp2, ! (as.character(haps_to_1)==as.character(haps_from_1) | as.character(haps_to_1)==as.character(haps_from_2) |  as.character(haps_to_2)==as.character(haps_from_1) | as.character(haps_to_2)==as.character(haps_from_2)) )
    
    
    
    # no idea why done like below
    tmp2_with_penalty <- c()
    
    for(i in unique(as.character(tmp2_one$gene_from))){
      #print(i)
      s <- subset(tmp2_one, as.character(gene_from)==i)
      
      b1 <- subset(s,  as.numeric(mr_from)==1 & as.numeric(mr_to)==1 )
      b2 <- subset(s,  !( as.numeric(mr_from)==1 & as.numeric(mr_to)==1 ) )
      
      if(nrow(b1) >0){
        tmp2_with_penalty <- rbind(tmp2_with_penalty, b1)
      }
      if(nrow(b2) >0) {
        #                       # sperating: doesn't work well
        #                       b2_1 <- subset(b2, as.numeric(mr_to)==1 | haps_to_1==haps_to_2)
        #                       b2_2 <- subset(b2, as.numeric(mr_to)!=1 & haps_to_1!=haps_to_2  )
        #                       f(nrow(b2_1) >0){
        #                                tmp2_with_penalty <- rbind(tmp2_with_penalty, b2_1)
        #                       }
        #                       f(nrow(b2_2) >0){
        #                             b2_2$mr_to = b2_2$mr_to-sem
        #                             tmp2_with_penalty <- rbind(tmp2_with_penalty, b2_2)
        #                       }
        
        # together
        #b2$mr_from = b2$mr_from-sem
        b2$mr_to = b2$mr_to-sem
        tmp2_with_penalty <- rbind(tmp2_with_penalty, b2)
      }
      
    }
    
    for(i in unique(as.character(tmp2_both$gene_from))){
      #                 #print(i)
      
      #                   s <- subset(tmp2_both, as.character(gene_from)==i)
      #                   s$mr_from = s$mr_from-sem
      #                   s$mr_to = s$mr_to-sem
      #                   tmp2_with_penalty <- rbind(tmp2_with_penalty, s)
      
      
      s <- subset(tmp2_both, as.character(gene_from)==i)
      b1 <- subset(s,  as.numeric(mr_from)==1 & as.numeric(mr_to)==1 )
      if(nrow(b1) >0){
        b1$mr_from = b1$mr_from-sem
        tmp2_with_penalty <- rbind(tmp2_with_penalty, b1)
      }
      b2 <- subset(s, !( as.numeric(mr_from)==1 & as.numeric(mr_to)==1 ) )
      if(nrow(b2) >0) {
        b2$mr_from = b2$mr_from-sem  
        b2$mr_to = b2$mr_to-sem
        tmp2_with_penalty <- rbind(tmp2_with_penalty, b2)
      }
      
    }
    
    tmp2_with_penalty$mr_from[as.numeric(tmp2_with_penalty$mr_from) <0] <- 0
    tmp2_with_penalty$mr_to[as.numeric(tmp2_with_penalty$mr_to) <0] <- 0
    
    
    #path_frame <- rbind(tmp1, tmp2_with_penalty, tmp_missingGenes)
    path_frame <- rbind(tmp1, tmp2_with_penalty)
  }
  
  
  
  #path_frame$sum <- path_frame$mr_from + path_frame$mr_to
  path_frame$ave <- ( path_frame$mr_from + path_frame$mr_to ) /2
  #path_frame$mult <- path_frame$mr_from * path_frame$mr_to
  
  path_frame$ave[path_frame$ave==0] <- 0.0000000000000001
  path_frame$ave[path_frame$ave==1] <- 0.9999999999999999
  
  path_frame$weight <- 1/-(log( as.numeric( (1-path_frame$ave) ) )  )
  
  
  
  #########################################################
  ###
  ### Get the first shortest path
  ### 
  #########################################################
  
  path_frame2 <- path_frame
  path_frame2$gene_from2 <- paste0(path_frame2$gene_from, "_", path_frame2$haps_from)
  path_frame2$gene_to2 <- paste0(path_frame2$gene_to, "_", path_frame2$haps_to)
  path_frame2$gene_from2 <- gsub("(Start)_.*", "Start", path_frame2$gene_from2)
  path_frame2$gene_to2 <- gsub("(End)_.*", "End", path_frame2$gene_to2)
  
  path_frame2$haps_from2 <- as.character(path_frame2$haps_from)
  path_frame2$haps_from2[as.character(path_frame2$gene_from2)=="Start"] <- "Start"
  path_frame2$haps_to2 <- as.character(path_frame2$haps_to)
  path_frame2$haps_to2[as.character(path_frame2$gene_to2)=="End"] <- "End"
  
  path_frame2 <- path_frame2[,c(3:4,7:18)]
  #path_frame2 <- path_frame2[!is.na(path_frame2$sum),]
  path_frame2 <- unique(path_frame2)
  
  
  
  node_info <- rbind( data.frame(gene_hap=as.character(path_frame2$gene_from2),
                                 hap=as.character(path_frame2$haps_from2), stringsAsFactors=FALSE ), 
                      data.frame(gene_hap=as.character(path_frame2$gene_to2), 
                                 hap=as.character(path_frame2$haps_to2) , stringsAsFactors=FALSE) )
  node_info <- unique(node_info)
  
  
  # 1) from start to end
  g <- graph.data.frame(cbind(as.character(path_frame2$gene_from2), as.character(path_frame2$gene_to2) ), vertices=node_info, directed=TRUE)
  
  E(g)$weight <- as.numeric(path_frame2$weight)
  #distances(g)
  
  #ShortestPath=shortest_paths(g, "Start", "End", mode="out", output="both")
  ShortestPath=shortest_paths(g, "Start", "End", mode="all", output="both")   # all shows better results than out
  
  #ShortestPath=shortest_paths(g, "HLA-DQB1_cox_pgf", "End", mode="out", output="both")  
  top1 <- E(g)[E(g) %in% ShortestPath[[2]][[1]] ]
  selected_haps <- subgraph.edges(g, top1)
  top1_edges <- as.data.frame(get.edgelist(selected_haps) )
  
  
  
  
  # #2) from end to start
  # g2 <- graph.data.frame(cbind(as.character(path_frame2$gene_to2), as.character(path_frame2$gene_from2) ), vertices=node_info, directed=TRUE)
  # E(g2)$weight <- as.numeric(path_frame2$weight)
  # ShortestPath=shortest_paths(g2, "End", "Start", mode="all", output="both")
  # top1 <- E(g2)[E(g2) %in% ShortestPath[[2]][[1]] ]
  # selected_haps <- subgraph.edges(g2, top1)
  # top1_edges <- as.data.frame(get.edgelist(selected_haps) )
  
  
  
  #######################################
  ###
  ### draw a heamap
  ###
  #######################################
  
  #	all_mapping_rates <- rbind( 
  #				data.frame( gene_name=paired_mapping_rates$gene_name,
  #					gene_hap=paste0(paired_mapping_rates$gene_name, "_", paired_mapping_rates$Hap1),
  #					Mapping_Rates=(paired_mapping_rates$only_h1 + paired_mapping_rates$intersect) ...)
  
  
  
  
  dat <- rbind( data.frame(gene_name=gsub("(.*)_(.*)_(.*)", "\\1", top1_edges$V1),
                           hap=gsub("(.*)_(.*)_(.*)", "\\2", top1_edges$V1),
                           best_pair=1),
                data.frame(gene_name=gsub("(.*)_(.*)_(.*)", "\\1", top1_edges$V1),
                           hap=gsub("(.*)_(.*)_(.*)", "\\3", top1_edges$V1),
                           best_pair=1) )
  
  dat <- subset(dat, gene_name!="Start")
  dat$gene_hap <- paste0(dat$gene_name, "_", dat$hap)
  
  z2 <- z[,c("gene_from", "hap_from")]
  z2 <- unique(z2)
  z2 <- z2[z2$gene_from %in% dat$gene_name,]
  z2$gene_hap <- paste0(z2$gene_from, "_", z2$hap_from)
  
  dat <- merge(dat, z2, by="gene_hap", all=TRUE)
  dat$best_pair[is.na(dat$best_pair)] <- 0
  dat <- dat[,c("gene_hap", "gene_from", "hap_from", "best_pair")]
  
  dat <- merge(dat, gene_order, by="gene_from", all=TRUE)
  colnames(dat) <- c("gene_name", "gene_hap", "hap", "best_pair", "ord_from", "ord")
  dat <- dat[!is.na(dat$ord), ]
  dat <- dat[order(dat$ord),]
  
  
  mapping_rates <- paired_mapping_rates
  tmp1 <- aggregate(rates_of_gene_counts ~ gene_name, data = mapping_rates, max)
  tmp2 <- aggregate(gene_count ~ gene_name, data = mapping_rates, max)
  tmp <- merge(tmp1, tmp2, by="gene_name")
  tmp$total_reads <- tmp$gene_count / tmp$rates_of_gene_counts
  tmp <- tmp[,c("gene_name", "total_reads")]
  
  mapping_rates <- merge(mapping_rates, tmp, by="gene_name")
  mapping_rates$Mapping_Rates1 <- ( mapping_rates$only_h1 + mapping_rates$intersect ) / mapping_rates$total_reads
  mapping_rates$Mapping_Rates2 <- ( mapping_rates$only_h2 + mapping_rates$intersect ) / mapping_rates$total_reads
  mapping_rates2 <- rbind(data.frame(gene_hap=paste0(mapping_rates$gene_name, "_", mapping_rates$Hap1),
                                     Mapping_Rates=mapping_rates$Mapping_Rates1),
                          data.frame(gene_hap=paste0(mapping_rates$gene_name, "_", mapping_rates$Hap2),
                                     Mapping_Rates=mapping_rates$Mapping_Rates2) )
  
  mapping_rates2 <- unique(mapping_rates2)
  
  
  dat <- merge(dat, mapping_rates2, by.x="gene_hap")
  dat <- dat[order(dat$ord),]
  dat$ord <- dat$ord * -1
  dat <- unique(dat)
  
  
  gene_labels <- unique(dat$gene_name)
  gene_labels <- rev(gene_labels)
  #dat_1 <- subset(dat, best_pair==1)
  
  
  # head( selected_genes )
  
  ##########################################
  ### matching output based on top1path
  #######q###################################
  
  zyg <- subset(dat, best_pair==1)
  zyg_count <- data.frame(table(zyg$gene_name))
  
  
  #hm
  dat_1 <- dat[dat$gene_name %in% as.character(zyg_count$Var1[zyg_count$Freq==1]), ]
  dat_1 <- dat_1[dat_1$gene_hap %in% zyg$gene_hap,]
  
  #ht
  dat_2 <- dat[dat$gene_name %in% as.character(zyg_count$Var1[zyg_count$Freq==2]), ]
  dat_2 <- dat_2[dat_2$gene_hap %in% zyg$gene_hap,]
  
  
  # multiple hm
  dat_3 <- dat_1[dat_1$gene_name %in% n_hap_any_hm$gene_name, ]
  
  # multiple ht
  dat_4 <- dat_2[dat_2$gene_name %in% n_hap_any_ht$gene_name, ]
  
  
  dat_1 <- dat_1[!(dat_1$gene_name %in% dat_3$gene_name), ]
  dat_2 <- dat_2[!(dat_2$gene_name %in% dat_4$gene_name), ]
  
  
  
  
  paired_mapping_rates$g_h1_h2 <- paste0(paired_mapping_rates$gene_name, "_", paired_mapping_rates$Hap1, "_", paired_mapping_rates$Hap2) 
  
  selected_genes_by_top1path <- paired_mapping_rates[(paired_mapping_rates$g_h1_h2 %in% top1_edges$V1), ]
  type <- ifelse(as.character(selected_genes_by_top1path$Hap1)==as.character(selected_genes_by_top1path$Hap2), "homozygous", "heterozygous")     
  selected_genes_by_top1path$type= type
  
  

  
  selected_genes_readCount <- merge(unique(selected_genes_by_top1path[,c("gene_name", "only_h1", "only_h2", "intersect", "type")]), unique(dat[,c("gene_name", "ord")]), by="gene_name")       
  selected_genes_readCount <- selected_genes_readCount[order(selected_genes_readCount$ord, decreasing = TRUE), ]
  colnames(selected_genes_readCount) <- c("gene_name", "Hap1only", "Hap2only", "Hap1andHap2", "Zygosity", "ord")
  
  selected_genes_readCount <- melt(selected_genes_readCount, id=c("gene_name", "ord", "Zygosity") )
  colnames(selected_genes_readCount) <- c("gene_name", "ord",  "Zygosity", "Haplotype", "value")
  
  
  ### refine output
  
  dat_final <- c()
  selected_genes_readCount_final <- c()
  for(i in unique(as.character(selected_genes_readCount$gene_name)) ){
    
    a <- subset(dat, as.character(gene_name)==i)
    
    tmp <- subset(selected_genes_readCount, as.character(gene_name)==i)
    tmp <- subset(tmp, value!=0 )
    if(nrow(tmp)==1){
          tmp$Zygosity <- "homozygous"
          tmp$Haplotype <- "Hap1andHap2"
          #a$best_pair[which(a$best_pair==1)][2] <- 0
          if(length( a$best_pair[which(a$best_pair==1)] ) >1){
                print(a)
                a$best_pair[which(a$best_pair==1)] <- c (1,0)
                #print(a)
          }
    }
    
    if(nrow(tmp)==2){
      
      if( tmp$Haplotype[1]=="Hap1only" & tmp$Haplotype[2]=="Hap2only" ){
      }
      else if( any(grepl("Hap1only", tmp$Haplotype)) ){
        tmp <- data.frame(tmp[1,1:2], 
                          Zygosity="homozygous", Haplotype="Hap1andHap2",
                          value=sum(tmp$value))
        #a$best_pair[which(a$best_pair==1)][2] <- 0
        a$best_pair[which(a$best_pair==1)] <- c (1,0)
      }
      else if( any(grepl("Hap2only", tmp$Haplotype) )  ){
        tmp <- data.frame(tmp[1,1:2], 
                          Zygosity="homozygous", Haplotype="Hap1andHap2",
                          value=sum(tmp$value))
        #a$best_pair[which(a$best_pair==1)][1] <- 0
        a$best_pair[which(a$best_pair==1)] <- c (0,1)
      }
      
    }
    
    selected_genes_readCount_final <- rbind(selected_genes_readCount_final, tmp)
    
    dat_final <- rbind(dat_final, a)
  }
  
  
  paired_mapping_rates$Hap1 <- toupper(paired_mapping_rates$Hap1)
  paired_mapping_rates$Hap2 <- toupper(paired_mapping_rates$Hap2)
  
  dat_final$hap <- toupper(dat_final$hap)
  dat_0 <- subset(dat_final, best_pair==1)
  #tmp <- c( paste0(selected_genes$gene_name, "_", selected_genes$Hap1), paste0(selected_genes$gene_name, "_", selected_genes$Hap2) )
  #tmp <- unique(tmp)
  #dat_0 <- dat_0[dat_0$gene_hap %in% tmp,]
  
  
  final_set <- c()
  for(i in unique(as.character(dat_0$gene_name))){
    a <- subset(dat_0, gene_name==i)
    if(nrow(a)==2){
      b <- subset( paired_mapping_rates, as.character(gene_name)==i & as.character(Hap1)==as.character(a$hap)[1] & as.character(Hap2)==as.character(a$hap)[2] )
    }
    if(nrow(a)==1){
      b <- subset( paired_mapping_rates, as.character(gene_name)==i & as.character(Hap1)==as.character(a$hap)[1] & as.character(Hap2)==as.character(a$hap)[1] )
    }
    final_set <- rbind(final_set, b)
  }
  
  
  final_set <- merge(final_set, unique(selected_genes_readCount_final[,c("gene_name", "ord", "Zygosity")]), by="gene_name")
  
  
  summary_set <- final_set[, c("gene_name", "ord", "mapping_rate", "rates_of_gene_counts", "Zygosity")]
  summary_set$mapping_rate <- 1- summary_set$mapping_rate
  summary_set$rates_of_gene_counts <- summary_set$rates_of_gene_counts *100
  colnames(summary_set) <- c("gene_name", "ord", "Mismatching_rates", "Rates_of_gene_counts", "Zygosity")
  summary_set <- melt(summary_set, id.vars = c("gene_name", "ord", "Zygosity"))
  summary_set$variable <- factor(summary_set$variable, levels=c("Rates_of_gene_counts", "Mismatching_rates"))

  
  summary_heatmap <- ggplot(summary_set, aes(x = variable, y = ord , fill=Zygosity)) + 
    geom_tile() + scale_fill_brewer(palette="Pastel1") +
    geom_tile(data=summary_set, aes(x=variable, y=ord, color=Zygosity),size=0,fill=NA, color="grey") +
    geom_text(aes(label=round( value,2) ), size=5, color="black" ) +
    scale_y_continuous(expand=c(0,0),breaks=min(summary_set$ord):max(summary_set$ord),labels=NULL) +
    #theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
    ggtitle("Combined mapping rates") + 
    theme(legend.position="none") + ylab(NULL) + xlab("Mapping rates") +
    scale_x_discrete(breaks=c("Rates_of_gene_counts","Mismatching_rates"),
                     labels=c("Gene counts", "Mismatches"))
  
  
  # gene_readCount <-  ggplot(selected_genes_readCount_final, aes(x=ord, y=value, fill=Haplotype, order=desc(Haplotype ))) +
  gene_readCount <-  ggplot(selected_genes_readCount_final, aes(x=ord, y=value, fill=Haplotype, order=Haplotype )) + 
    geom_bar(stat="identity") + 
    scale_x_continuous(expand=c(0,0),breaks=min(selected_genes_readCount_final$ord):max(selected_genes_readCount_final$ord),labels=gene_labels) +
    #theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + 
    theme(axis.text.y=element_text(size = 12, face = "bold") ) +
    ylab("Read counts") + xlab(NULL) + coord_flip() +
    ggtitle(paste0( "Gene counts" ) )
  
  
  dat_final$Mapping_Rates <- dat_final$Mapping_Rates *100
  
  mhc_heatmap <- ggplot(dat_final, aes(x = hap , y = ord , fill=Mapping_Rates)) + 
    geom_tile() +
    #scale_fill_gradient() +
    #scale_fill_gradient2(low = "darkblue", mid = "blue", high = "lightblue") +
    #scale_fill_gradient2(high = "darkblue", mid = "lightskyblue1", low = "white")+
    scale_fill_gradientn(breaks=c(0,50,75,90,100), colours=c("lightblue1", "lightskyblue1", "steelblue1", "deepskyblue3", "darkblue"))+
    scale_y_continuous(expand=c(0,0),breaks=min(dat$ord):max(dat$ord),labels=gene_labels) +
    xlab("Haplotypes") + ylab("Genes by the chromosomal order") +
    geom_tile(data=dat_final, aes(x=hap, y=ord, color=best_pair),size=0,fill=NA, color="white") +
    geom_tile(data=dat_0, aes(x=hap, y=ord, color=best_pair),size=1,fill=NA, color="red") +    # was blue
    #geom_tile(data=dat_4, aes(x=hap, y=ord, color=best_pair),size=1,fill=NA, color="red") +    # was blue
    #geom_tile(data=dat_3, aes(x=hap, y=ord, color=best_pair),size=2,fill=NA, color="red") +    # was blue
    #geom_tile(data=dat_2, aes(x=hap, y=ord, color=best_pair),size=1,fill=NA, color="red") +
    #geom_tile(data=dat_1, aes(x=hap, y=ord, color=best_pair),size=2,fill=NA, color="red") +
    #geom_text(aes(label=round( (Mapping_Rates*100),2) ), size=5, color="white" ) +
    geom_text(aes(label=round( Mapping_Rates,2) ), size=5, color="white" ) +
    ggtitle(paste0("Predicted haplotypes and mapping rates") ) + 
    theme(legend.position="left", axis.text.y=element_text(size = 12, face = "bold"), axis.title.y = element_text(face="bold", size=18))
  
  gp1<- ggplot_gtable(ggplot_build( gene_readCount ))
  gp2<- ggplot_gtable(ggplot_build( mhc_heatmap ))
  gp3<- ggplot_gtable(ggplot_build( summary_heatmap )) 
  
  #maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
  #gp1$widths[2:3] <- maxWidth
  #gp2$widths[2:3] <- maxWidth
  
  #g_p12 <- grid.arrange(gp2, gp1, ncol=2, widths = c(2/3,1/3) )
  g_p123 <- grid.arrange(gp2, gp3, gp1, ncol=3, widths = c(8/15, 2/15, 5/15) )
  

  
  # testing
  test <- data.frame(table(paste0(selected_genes_by_top1path$Hap1, "_", selected_genes_by_top1path$Hap2)) )
  print(test)
  #write.table(test, paste0(sample_name, ".stats.tmp.txt"), quote=F)
  
 
  
  ggsave(filename=paste0(sample_name, ".shortestPath_heatmap.pdf" ), g_p123, width=18, height=22)
  
  
  output <- list(heatmap_genes=dat,
                 #summary_of_best_pairs=selected_genes_by_top1path,
                 summary_of_best_pairs=final_set)
  
}



#' To predict shortest path using results by permutation test
#' 
#' @param paired_mapping_rates output of 'predictRecombHapByGenes'  (*max_comb.txt)
#' @param permutation_result output of 'permutation_test'  (*permutation_test.txt)
#' @param anno gtf version (ex, "hg38")
#' @param penalty 1: penalty applied, 0: penalty not applied
#' @param sample_name output file name
#' @export
#' @import igraph
#' @import dplyr
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import plyr
#' @importFrom data.table data.table
#' @examples
#' output <- getHapsByShortestPaths_withPermutation_test(paired_mapping_rates, anno="hg38", penalty=1, sample_name)
#' 
#' # Output data
#' # *.shortestPath_heatmap.pdf : heatmap of selected haplotypes
#' # predicted best pairs of haploytpes
#' summary_of_best_pairs <- output$summary_of_best_pairs

getHapsByShortestPaths_withPermutation_test <- function(paired_mapping_rates, permutation_result, anno="hg38", penalty=1, sample_name){
  
  unique_best_hits <- data.frame(table(paired_mapping_rates[,c("gene_name", "type")]) )
  
  paired_mapping_rates_noFT <-  paired_mapping_rates
  paired_mapping_rates <- subset(paired_mapping_rates, gene_count >= 50)   
  

  gene_list <- unique(as.character(paired_mapping_rates$gene_name)) 
  head(paired_mapping_rates)
  
  comb_haps <- unique(paste0( paired_mapping_rates[,c("Hap1")], "_", paired_mapping_rates[,c("Hap2")]) )
  
  v1 <- comb_haps
  v2 <- comb_haps
  tmp <- sort(apply(expand.grid(v1, v2), 1, paste, collapse = ",", sep = "")) 
  hap_combn <- data.frame(haps_from=gsub("(.*),(.*)", "\\1", tmp),
                          haps_to=gsub("(.*),(.*)", "\\2", tmp) )
  
  
  
  ### summary of best_pairs


  # additional filtering
    a_ <- subset(paired_mapping_rates, (only_h1>5 & only_h2>5 & Hap1!=Hap2)| (only_h1==0 & only_h2==0 & Hap1==Hap2) )       #wo.dp5


  # Kat's and others for paper
  #a_ <- subset(paired_mapping_rates, putative_haps2==1)
  #a_ <- subset(a_, (only_h1>0 & intersect >0 & Hap1!=Hap2)| (only_h2>0 & intersect >0 & Hap1!=Hap2)| (only_h1==0 & only_h2==0 & Hap1==Hap2) )       #wo.dp5
  

  
  
  a <- rbind(data.frame(gene_name=a_$gene_name, gene_hap=paste0(a_$gene_name, "_", a_$Hap1)),
             data.frame(gene_name=a_$gene_name, gene_hap=paste0(a_$gene_name, "_", a_$Hap2)))
  a <- unique(a)
  a2 <- data.frame(table(a$gene_name))
  
  
  n_hap1 <- a_[as.character(a_$gene_name) %in% as.character(a2[a2$Freq==1,]$Var1),]     # homozygous
  n_hap2 <- a_[as.character(a_$gene_name) %in% as.character(a2[a2$Freq==2,]$Var1),]     # heterozygous
  n_hap_any <- a_[as.character(a_$gene_name) %in% as.character(a2[a2$Freq>2,]$Var1),]   # more than two pairs
  n_hap_any_hm <- subset(n_hap_any, Hap1==Hap2)
  #n_hap_any_ht <- n_hap_any[!(n_hap_any$gene_name %in% n_hap_any_hm$gene_name), ]
  n_hap_any_ht <- subset(n_hap_any, Hap1!=Hap2)
  
  
  selected_genes <- c()
  if(nrow(n_hap1) >0){
    selected_genes <- rbind(selected_genes, data.frame(n_hap1, type="hm") )
  }
  if(nrow(n_hap2) >0){
    selected_genes <- rbind(selected_genes, data.frame(n_hap2, type="ht") )
  }
  if(nrow(n_hap_any_hm) >0){
    selected_genes <- rbind(selected_genes, data.frame(n_hap_any_hm, type="hm_multi") )
  }
  if(nrow(n_hap_any_ht) >0){
    selected_genes <- rbind(selected_genes, data.frame(n_hap_any_ht, type="ht_multi") )
  }
  
  
  # the same results after changing mapping rates
  # selected_genes <- selected_genes[!( paste0(selected_genes$gene_name, "_", selected_genes$Hap1, "_", selected_genes$Hap2) %in% paste0(selected_genes_ft$gene_name, "_", selected_genes_ft$Hap1, "_", selected_genes_ft$Hap2) ), ]
  
  

  input_structure <- unique(paste0( paired_mapping_rates[,c("Hap1")], "_", paired_mapping_rates[,c("Hap2")]) )
  input_structure <- data.frame(hap_pairs= input_structure[order(input_structure)])
  
  for(i in gene_list){
        s <- subset(paired_mapping_rates, as.character(gene_name)==i)
        s <- unique(s)
        s$tmp <- paste0(s$Hap1, "_", s$Hap2)
        s$mapping_rate[s$type=="B"] <- 0
        
        #input_structure <- merge(input_structure, s[,c("tmp", "rates_of_gene_counts")], by.x="hap_pairs", by.y="tmp", all=TRUE)
        input_structure <- merge(input_structure, s[,c("tmp", "mapping_rate")], by.x="hap_pairs", by.y="tmp", all=TRUE)
        
        colnames(input_structure)[ncol(input_structure)] <- i
        #print(i)
        #print(dim(input_structure))
  }
  input_structure <- unique(input_structure) 
  input_structure <- t(input_structure)
  colnames(input_structure) <- input_structure[1,]
  input_structure <- input_structure[2:nrow(input_structure),]
  


  
  #########################################################
  ###
  ### function: Combining two p-values
  ###
  #########################################################
  
  #combineTwoPvalues <- function(x,y) return(pchisq(-2*log(x/2)-2*log(y/2),4, low=F))
  
  
  
  #########################################################
  ###
  ### get annotation
  ###
  #########################################################

  hap_names <- c("apd", "cox", "dbb", "mann", "mcf", "pgf", "qbl", "ssto")
  
  if(anno=="hg38"){     
        hap_chrs<- c("GL000250.2", "GL000251.2", "GL000252.2", "GL000253.2", "GL000254.2", "chr6", "GL000255.2", "GL000256.2")
        haps <- data.frame(hap_names, hap_chrs)
        
 	      annoGenes <- read.table("/well/jknight/Wan/pipl/haplotype-mapping/summarizingAnno/annoGenes_ordered.txt", sep="\t", header=T)
	      annoGenes$gene_from <- gsub("VARSL", "VARS2", as.character(annoGenes$gene_from) )
	      annoGenes$gene_to <- gsub("VARSL", "VARS2", as.character(annoGenes$gene_to) )
	      annoGenes$gene_to <- gsub("C6orf205", "MUC21", as.character(annoGenes$gene_to) )
        annoGenes$gene_to <- gsub("C6orf205", "MUC21", as.character(annoGenes$gene_to) )      
	      annoGenes <- subset(annoGenes, as.character(gene_from)!="TRIM39-RPP21_pgf" & as.character(gene_to)!="TRIM39-RPP21_pgf" ) 	
        
        identicalGeneList <- read.table("/well/jknight/Wan/git/haplotype-mapping/add_to_HapTag/identicalGenes2.txt", sep="\t", header=T)
  }

  
  
  #z <- annoGenes[ (gsub( "(.*)_.*", "\\1", as.character(annoGenes$gene_from) ) %in% as.character(rownames(input_structure)) ) | (gsub( "(.*)_.*", "\\1", as.character(annoGenes$gene_to) ) %in% as.character(rownames(input_structure)) ),]  
  z <- annoGenes
  z$gene_from <- gsub("(.*)_.*", "\\1", z$gene_from)
  z$gene_to <- gsub("(.*)_.*", "\\1", z$gene_to)
  # write.table(z, "gene_path_structure.txt", sep="\t", quote=F, row.names=F)  

	getMissingGenes <- unique(z[,c("gene_from", "hap_from")] ) 
	getMissingGenes <- data.frame(table(getMissingGenes$gene_from))
	getMissingGenes <- subset(getMissingGenes, Freq>0 & Freq<8)




  #########################################################
  ###
  ### get ordered genes in the MHC region
  ###
  #########################################################
  
  
  
  ordered_genes <- c()
  for(i in unique(as.character(z$hap_from)) ){
    hn <- subset(z, as.character(hap_from)==i & as.character(hap_to)==i)
    if(nrow(hn) >0 ){
      
      hn <- rbind( data.frame(gene_name=as.character(hn$gene_from), ord_bu=as.numeric(hn$ord_from) ),
                   data.frame(gene_name=as.character(hn$gene_to), ord_bu=as.numeric(hn$ord_to))  )
      hn <- hn[!is.na(hn$gene_name),]
      hn <- unique(hn)
      hn <- hn[order(as.numeric(hn$ord_bu)),]
      hn$ord <- c(1:nrow(hn))
      
      a0 <- as.character(subset(hn, ord==min(ord) )$gene_name )
      ordered_genes <- rbind(ordered_genes, 
                             data.frame(gene_from="Start", gene_to=a0, hap_from="Start", hap_to=i, ord_from=0, ord_to=1) )
      
      for(i2 in 1: (length(hn$ord)-1) ){
        g1=hn$ord[i2]
        g2=hn$ord[i2+1]
        
        a1 <- as.character(subset(hn, ord==g1 )$gene_name )
        a2 <- as.character(subset(hn, ord==g2 )$gene_name )
        tmp <- data.frame(gene_from=a1, gene_to=a2, hap_from=i, hap_to=i, ord_from=g1, ord_to=g2)
        tmp <- unique(tmp)
        tmp <- subset(tmp, as.character(gene_from)!=as.character(gene_to) )
        #print(tmp)
        ordered_genes <- rbind(ordered_genes, tmp)
      }
      
      a3 <- as.character(subset(hn, ord==max(hn$ord) )$gene_name )
      ordered_genes <- rbind(ordered_genes, 
                             data.frame(gene_from=a3, gene_to="End", hap_from=i, hap_to="End", 
                                        ord_from=max(hn$ord), ord_to=(max(hn$ord) +1)) )
    }
  }
  #ordered_genes <- subset(ordered_genes, gene_to!="End")
  
  
  
  
  num_genes_in_haps <- aggregate(ord_from ~ hap_from, data = ordered_genes, max)
  h1 <- as.character(num_genes_in_haps[which(num_genes_in_haps$ord_from==max(num_genes_in_haps$ord_from) ),]$hap_from)
  h2 <- as.character (haps$hap_names [!(haps$hap_names %in% c("apd", h1) )] )
  
  
  gene_order <- subset(ordered_genes, hap_from==h1)
  gene_order <- data.frame(gene_from=gene_order$gene_from, ord_from=gene_order$ord_from) 
  for(i in h2){
    
    tmp <- subset(ordered_genes, hap_from==i)
    tmp <- tmp[ !( tmp$gene_from %in% gene_order$gene_from), ]
    
    if(nrow(tmp) >0){

      for( missing_g in 1:nrow(tmp) ){
        #print(missing_g)
        #print( as.character( tmp$gene_from[missing_g] ) )
        
        ord= gene_order[grep(as.character(tmp$gene_to[missing_g]), gene_order$gene_from),]$ord_from
        if(length(ord) ==0){
          print( as.character( tmp$gene_from[missing_g] ) )
        }
        if(length(ord) ==1){
          gene_order[which(gene_order$ord_from>=ord),]$ord_from <- as.numeric( gene_order[which(gene_order$ord_from>=ord),]$ord_from +1 )
          gene_order <- rbind(gene_order, 
                              data.frame(gene_from=as.character( tmp$gene_from[missing_g] ),
                                         ord_from=ord) )
          gene_order <- gene_order[order(gene_order$ord_from),]
        }
        
      }
    }   
  }
  
  gene_order <- gene_order[gene_order$gene_from %in% rownames(input_structure), ]
  colnames(gene_order) <- c("gene_from", "ord_tmp")
  gene_order$ord_from=c(1:nrow(gene_order))
  
  gene_order=aggregate(ord_from ~ gene_from, data = gene_order, max)
  gene_order <- gene_order[order(gene_order$ord_from),]
  gene_order$ord=c(1:nrow(gene_order))
  
  
  #mhc_genes <- z
  
  
  
  #################################################################
  ###
  ### get ordered genes with pvlaues from predictRecombHapByGenes
  ###
  #################################################################
  

  input_structure <- data.frame(input_structure)
  input_structure$gene_name <- rownames(input_structure)
  mr <- merge(input_structure, gene_order, by.x="gene_name", by.y="gene_from")
  mr <- mr[order(as.numeric(mr$ord)),]
  
  
  t1 <- data.frame(c(rep(0, length(comb_haps)), 0,0) )
  rownames(t1) <- c(comb_haps, "ord_from", "ord")
  t1 <- data.frame(gene_name="Start", t(t1) )
  
  t2 <- data.frame( c(rep(0, length(comb_haps)), rep((max(mr$ord)+1), 2)) )
  rownames(t2) <- c(comb_haps, "ord_from", "ord")
  t2 <- data.frame(gene_name="End", t(t2) )
  
  mr <- rbind(t1, mr, t2)
  
  
#   
#   mr <- rbind( c("Start", rep(0, length(comb_haps)), 0,0), 
#                mr,
#                c("End", rep(0, length(comb_haps)), rep((max(mr$ord)+1), 2)) )
              

  
  ### 
  tmp_haps <- data.frame(h1=gsub("(.*)_(.*)", "\\1", colnames(mr)[2:ncol(mr)] ),
                         h2=gsub("(.*)_(.*)", "\\2", colnames(mr)[2:ncol(mr)] ) )
  
  mr2 <- data.frame(haps=colnames(mr)[2:ncol(mr)])
  for(i in 1:nrow(mr)){
        tmp <- data.frame( tmp_haps, t(mr[i,2:ncol(mr)] ) ) 
        tmp$haps <- rownames(tmp)
        
        
        ### bu: to change mapping rate to 0 if both haplotypes are the same.
        a <- subset(tmp, as.character(h1)==as.character(h2) & as.numeric(as.character(tmp[,3]))==1 & h1!="ord")
        if(nrow(a) >0){
              i2 <-  as.character(a$h1)
              b2 <- tmp[as.character(tmp$h1) %in% i2 & as.character(tmp$h2) %in% i2 &  as.character(tmp$h1)!=as.character(tmp$h2), ] 
              b1 <- tmp[ !( as.character(tmp$h1) %in% i2 & as.character(tmp$h2) %in% i2 &  as.character(tmp$h1)!=as.character(tmp$h2) ), ] 
              
              if(nrow(b2) >0){
                    b2[,3] <- 0
                    b1[,3] <- as.numeric(as.character(b1[,3]))
                    tmp <- rbind(b1, b2)
              }
              if(nrow(b2) ==0){
                tmp <- b1
              }
              
              # b2 <- tmp[as.character(tmp$h1) %in% i2 | as.character(tmp$h2) %in% i2, ]
              # b2[,3] <- 0
              # b3 <- tmp[!( as.character(tmp$h1) %in% i2) &  !(as.character(tmp$h2) %in% i2), ]
              # tmp <- rbind(b2, b3)
              # 
              # b1 <- tmp[as.character(tmp$h1) %in% i2 & as.character(tmp$h2) %in% i2 & as.character(tmp$h1)==as.character(tmp$h2), ]
              # b1[,3] <- 1
              # b3 <- tmp[!( as.character(tmp$h1) %in% i2 & as.character(tmp$h2) %in% i2 & as.character(tmp$h1)==as.character(tmp$h2) ), ]
              # tmp <- rbind(b1, b3)
        }
              mr2 <- merge(mr2, tmp[,3:4], by="haps")
  }
  mr2 <- t(mr2)
  colnames(mr2) <- mr2[1,]
  mr2 <- mr2[2:nrow(mr2), ]
  mr2 <- data.frame(mr2)
  mr2=data.frame(gene_name=mr$gene_name, mr2)

  mr <- cbind( mr2[, -(grep("ord", colnames(mr2)) ) ], mr2[, grep("ord", colnames(mr2)) ] )
  
  mr_col <- colnames(mr)[2:(ncol(mr)-2)]
  rownames(mr) <- NULL
  
  
  


  path_frame <- c()
  for(ord in 1:(nrow(mr)-1) ){
        message(paste0(ord, ": ", mr[ord, c("gene_name")]) )

	      g1 <- as.character(mr[ord, c("gene_name")])
	      g2 <- as.character(mr[(ord+1), c("gene_name")])
	      o1 <- as.numeric(as.character(mr[ord, c("ord")]))
	      o2 <- as.numeric(as.character(mr[(ord+1), c("ord")]))

	      gfrom <- data.frame(haps=mr_col, mr_from= as.numeric(t( mr[ord,2:(ncol(mr)-2)] )[,1] ) ) 
	      gto  <-  data.frame(haps=mr_col,  mr_to=as.numeric(t( mr[(ord+1),2:(ncol(mr)-2)] )[,1] ) ) 
	      tmp <- merge(hap_combn, gfrom, by.x="haps_from", by.y="haps")
 	      tmp <- merge(tmp, gto, by.x="haps_to", by.y="haps")
	      tmp <- data.frame(gene_from=g1, gene_to=g2, ord_from=o1, ord_to=o2,  tmp)
 
        path_frame <- rbind(path_frame, tmp)             
  }

  path_frame <- path_frame[!is.na(path_frame$mr_from), ]
  path_frame <- path_frame[!is.na(path_frame$mr_to), ]
  
  path_frame$haps_to_1=gsub("(.*)_(.*)", "\\1", path_frame$haps_to)
  path_frame$haps_to_2=gsub("(.*)_(.*)", "\\2", path_frame$haps_to)
  path_frame$haps_from_1=gsub("(.*)_(.*)", "\\1", path_frame$haps_from)
  path_frame$haps_from_2=gsub("(.*)_(.*)", "\\2", path_frame$haps_from)
 
  path_frame$haps_from <- as.character( path_frame$haps_from)
  path_frame$haps_from[path_frame$gene_from=="Start"] <- "Start"
  path_frame$haps_from_1 <- as.character( path_frame$haps_from_1)
  path_frame$haps_from_2 <- as.character( path_frame$haps_from_2)
  path_frame$haps_from_1[path_frame$gene_from=="Start"] <- "Start"
  path_frame$haps_from_2[path_frame$gene_from=="Start"] <- "Start"

  path_frame$haps_to <- as.character( path_frame$haps_to)
  path_frame$haps_to[path_frame$gene_to=="End"] <- "End"
  path_frame$haps_to_1 <- as.character( path_frame$haps_to_1)
  path_frame$haps_to_2 <- as.character( path_frame$haps_to_2)
  path_frame$haps_to_1[path_frame$gene_to=="End"] <- "End"
  path_frame$haps_to_2[path_frame$gene_to=="End"] <- "End"
  
  path_frame <- unique(path_frame)


 
    ##############################################
    ### applying penalty
    ##############################################


  if(penalty==1){
    
        tmp <- path_frame
        unique_best_hits_A1 <- subset(unique_best_hits, type=="A" & Freq==1)
        
        #tmp_missingGenes <- tmp[as.character(tmp$gene_from) %in% as.character(getMissingGenes$Var1) | as.character(tmp$gene_to) %in% as.character(getMissingGenes$Var1), ]
	      #tmp <- tmp[ !( as.character(tmp$gene_from) %in% as.character(getMissingGenes$Var1) | as.character(tmp$gene_to) %in% as.character(getMissingGenes$Var1) ), ]    

	     #tmp_missingGenes <- tmp[ as.character(tmp$gene_to) %in% as.character(getMissingGenes$Var1), ]
        #tmp <- tmp[ !( as.character(tmp$gene_to) %in% as.character(getMissingGenes$Var1) ), ]


       
 
        #data.summary<- read.table(paste0("/well/jknight/Wan/pipl/haplotype-mapping/sampling2/checkErrorRates/data.summary.txt"), sep="\t",  header=T )
        #sem <- subset(data.summary, gene=="all")$sem
        # 5 is a bit better than *10, 3 and 4 is similar. 3 and 5 is better than 1
        # after changing mapping rates, *1 is good
        sem <- 0.002783863
        #sem <- 0.03      # works well
        #sem <- 1- 0.002783863

        
          #tmp1 <- subset(tmp, as.character(haps_to)==as.character(haps_from) | (mr_from==1 & mr_to==1)  )
          #tmp2 <- subset(tmp, as.character(haps_to)!=as.character(haps_from) & (mr_from<1 | mr_to<1)  )
        
          tmp1 <- subset(tmp, as.character(haps_to)==as.character(haps_from) )
          tmp1 <- rbind(tmp1, subset(tmp, gene_from=="Start" | gene_to=="End"))

          tmp2 <- subset(tmp, as.character(haps_to)!=as.character(haps_from) )
          tmp2 <- subset(tmp2, gene_from!="Start" & gene_to!="End")
          #tmp2_1 <- subset(tmp2, mr_from==1 | mr_to==1)
          #tmp1 <- rbind(tmp1, tmp2_1)
          #tmp2 <- subset(tmp2, mr_from!=1 & mr_to!=1)
          
          
        
              tmp2_one <- subset(tmp2, as.character(haps_to_1)==as.character(haps_from_1) | as.character(haps_to_1)==as.character(haps_from_2) |  as.character(haps_to_2)==as.character(haps_from_1) | as.character(haps_to_2)==as.character(haps_from_2))
              tmp2_both <- subset(tmp2, ! (as.character(haps_to_1)==as.character(haps_from_1) | as.character(haps_to_1)==as.character(haps_from_2) |  as.character(haps_to_2)==as.character(haps_from_1) | as.character(haps_to_2)==as.character(haps_from_2)) )
            
    
              
              # no idea why done like below
            tmp2_with_penalty <- c()
        
            for(i in unique(as.character(tmp2_one$gene_from))){
                #print(i)
                s <- subset(tmp2_one, as.character(gene_from)==i)

                b1 <- subset(s,  as.numeric(mr_from)==1 & as.numeric(mr_to)==1 )
                b2 <- subset(s,  !( as.numeric(mr_from)==1 & as.numeric(mr_to)==1 ) )
          
                if(nrow(b1) >0){
                  tmp2_with_penalty <- rbind(tmp2_with_penalty, b1)
                }
                if(nrow(b2) >0) {

                      b2$mr_to = b2$mr_to-sem
                      tmp2_with_penalty <- rbind(tmp2_with_penalty, b2)
                }
          
              }
    
              for(i in unique(as.character(tmp2_both$gene_from))){

                  
                  s <- subset(tmp2_both, as.character(gene_from)==i)
                  b1 <- subset(s,  as.numeric(mr_from)==1 & as.numeric(mr_to)==1 )
                  if(nrow(b1) >0){
                        b1$mr_from = b1$mr_from-sem
                        tmp2_with_penalty <- rbind(tmp2_with_penalty, b1)
                  }
                  b2 <- subset(s, !( as.numeric(mr_from)==1 & as.numeric(mr_to)==1 ) )
                  if(nrow(b2) >0) {
     
                      b2$mr_from = b2$mr_from-sem    # *3 worse, *2 is the same with *1
                      b2$mr_to = b2$mr_to-sem
                      tmp2_with_penalty <- rbind(tmp2_with_penalty, b2)
                  }

              }
                       
              tmp2_with_penalty$mr_from[as.numeric(tmp2_with_penalty$mr_from) <0] <- 0
              tmp2_with_penalty$mr_to[as.numeric(tmp2_with_penalty$mr_to) <0] <- 0
             

             #path_frame <- rbind(tmp1, tmp2_with_penalty, tmp_missingGenes)
		        path_frame <- rbind(tmp1, tmp2_with_penalty)
  }
  


  #path_frame$sum <- path_frame$mr_from + path_frame$mr_to
  path_frame$ave <- ( path_frame$mr_from + path_frame$mr_to ) /2
  #path_frame$mult <- path_frame$mr_from * path_frame$mr_to
  
  path_frame$ave[path_frame$ave==0] <- 0.0000000000000001
  path_frame$ave[path_frame$ave==1] <- 0.9999999999999999
    
  path_frame$weight <- 1/-(log( as.numeric( (1-path_frame$ave) ) )  )

  

    #########################################################
    ###
    ### Get the first shortest path
    ### 
    #########################################################

    path_frame2 <- path_frame
    path_frame2$gene_from2 <- paste0(path_frame2$gene_from, "_", path_frame2$haps_from)
    path_frame2$gene_to2 <- paste0(path_frame2$gene_to, "_", path_frame2$haps_to)
    path_frame2$gene_from2 <- gsub("(Start)_.*", "Start", path_frame2$gene_from2)
    path_frame2$gene_to2 <- gsub("(End)_.*", "End", path_frame2$gene_to2)

    path_frame2$haps_from2 <- as.character(path_frame2$haps_from)
    path_frame2$haps_from2[as.character(path_frame2$gene_from2)=="Start"] <- "Start"
    path_frame2$haps_to2 <- as.character(path_frame2$haps_to)
    path_frame2$haps_to2[as.character(path_frame2$gene_to2)=="End"] <- "End"

    path_frame2 <- path_frame2[,c(3:4,7:18)]
    #path_frame2 <- path_frame2[!is.na(path_frame2$sum),]
    path_frame2 <- unique(path_frame2)



    node_info <- rbind( data.frame(gene_hap=as.character(path_frame2$gene_from2),
                                 hap=as.character(path_frame2$haps_from2), stringsAsFactors=FALSE ), 
                        data.frame(gene_hap=as.character(path_frame2$gene_to2), 
                                 hap=as.character(path_frame2$haps_to2) , stringsAsFactors=FALSE) )
    node_info <- unique(node_info)


  # 1) from start to end
    g <- graph.data.frame(cbind(as.character(path_frame2$gene_from2), as.character(path_frame2$gene_to2) ), vertices=node_info, directed=TRUE)
  
    E(g)$weight <- as.numeric(path_frame2$weight)
    #distances(g)

    #ShortestPath=shortest_paths(g, "Start", "End", mode="out", output="both")
    ShortestPath=shortest_paths(g, "Start", "End", mode="all", output="both")   # all shows better results than out

    #ShortestPath=shortest_paths(g, "HLA-DQB1_cox_pgf", "End", mode="out", output="both")  
    top1 <- E(g)[E(g) %in% ShortestPath[[2]][[1]] ]
    selected_haps <- subgraph.edges(g, top1)
    top1_edges <- as.data.frame(get.edgelist(selected_haps) )

    
    

    # #2) from end to start
    # g2 <- graph.data.frame(cbind(as.character(path_frame2$gene_to2), as.character(path_frame2$gene_from2) ), vertices=node_info, directed=TRUE)
    # E(g2)$weight <- as.numeric(path_frame2$weight)
    # ShortestPath=shortest_paths(g2, "End", "Start", mode="all", output="both")
    # top1 <- E(g2)[E(g2) %in% ShortestPath[[2]][[1]] ]
    # selected_haps <- subgraph.edges(g2, top1)
    # top1_edges <- as.data.frame(get.edgelist(selected_haps) )



    #######################################
    ###
    ### draw a heamap
    ###
    #######################################

#	all_mapping_rates <- rbind( 
#				data.frame( gene_name=paired_mapping_rates$gene_name,
#					gene_hap=paste0(paired_mapping_rates$gene_name, "_", paired_mapping_rates$Hap1),
#					Mapping_Rates=(paired_mapping_rates$only_h1 + paired_mapping_rates$intersect) ...)

 

 
    dat <- rbind( data.frame(gene_name=gsub("(.*)_(.*)_(.*)", "\\1", top1_edges$V1),
                           hap=gsub("(.*)_(.*)_(.*)", "\\2", top1_edges$V1),
                           best_pair=1),
                  data.frame(gene_name=gsub("(.*)_(.*)_(.*)", "\\1", top1_edges$V1),
                           hap=gsub("(.*)_(.*)_(.*)", "\\3", top1_edges$V1),
                           best_pair=1) )
                          
    dat <- subset(dat, gene_name!="Start")
    dat$gene_hap <- paste0(dat$gene_name, "_", dat$hap)
  
    z2 <- z[,c("gene_from", "hap_from")]
    z2 <- unique(z2)
    z2 <- z2[z2$gene_from %in% dat$gene_name,]
    z2$gene_hap <- paste0(z2$gene_from, "_", z2$hap_from)

    dat <- merge(dat, z2, by="gene_hap", all=TRUE)
    dat$best_pair[is.na(dat$best_pair)] <- 0
    dat <- dat[,c("gene_hap", "gene_from", "hap_from", "best_pair")]

    dat <- merge(dat, gene_order, by="gene_from", all=TRUE)
    colnames(dat) <- c("gene_name", "gene_hap", "hap", "best_pair", "ord_from", "ord")
    dat <- dat[!is.na(dat$ord), ]
    dat <- dat[order(dat$ord),]


    mapping_rates <- paired_mapping_rates
    tmp1 <- aggregate(rates_of_gene_counts ~ gene_name, data = mapping_rates, max)
    tmp2 <- aggregate(gene_count ~ gene_name, data = mapping_rates, max)
    tmp <- merge(tmp1, tmp2, by="gene_name")
    tmp$total_reads <- tmp$gene_count / tmp$rates_of_gene_counts
    tmp <- tmp[,c("gene_name", "total_reads")]

    mapping_rates <- merge(mapping_rates, tmp, by="gene_name")
    mapping_rates$Mapping_Rates1 <- ( mapping_rates$only_h1 + mapping_rates$intersect ) / mapping_rates$total_reads
    mapping_rates$Mapping_Rates2 <- ( mapping_rates$only_h2 + mapping_rates$intersect ) / mapping_rates$total_reads
    mapping_rates2 <- rbind(data.frame(gene_hap=paste0(mapping_rates$gene_name, "_", mapping_rates$Hap1),
                                       Mapping_Rates=mapping_rates$Mapping_Rates1),
                            data.frame(gene_hap=paste0(mapping_rates$gene_name, "_", mapping_rates$Hap2),
                                       Mapping_Rates=mapping_rates$Mapping_Rates2) )

    mapping_rates2 <- unique(mapping_rates2)


    dat <- merge(dat, mapping_rates2, by.x="gene_hap")
    dat <- dat[order(dat$ord),]
    dat$ord <- dat$ord * -1
    dat <- unique(dat)


    gene_labels <- unique(dat$gene_name)
    gene_labels <- rev(gene_labels)
    #dat_1 <- subset(dat, best_pair==1)


    ##########################################
    ### matching output with max_comb results 
    ##########################################

#     #hm
#     dat_0 <- subset(dat, best_pair==1)
#     dat_1 <- rbind( dat[dat$gene_hap %in% paste0(n_hap1$gene, "_", n_hap1$Hap1), ],
#                     dat[dat$gene_hap %in% paste0(n_hap1$gene, "_", n_hap1$Hap2), ] )
#     dat_1 <- unique(dat_1)
# 
#     #ht
#     dat_2 <- rbind( dat[dat$gene_hap %in% paste0(n_hap2$gene, "_", n_hap2$Hap1), ],
#                     dat[dat$gene_hap %in% paste0(n_hap2$gene, "_", n_hap2$Hap2), ] )
#     dat_2 <- unique(dat_2)
# 
#     # multiple hm
#     dat_3 <-  rbind( dat[dat$gene_hap %in% paste0(n_hap_any_hm$gene, "_", n_hap_any_hm$Hap1) & dat$best_pair==1, ],
#                      dat[dat$gene_hap %in% paste0(n_hap_any_hm$gene, "_", n_hap_any_hm$Hap2) & dat$best_pair==1, ] )
#     dat_3 <- unique(dat_3)
# 
#     # multiple ht
#     dat_4 <-  rbind( dat[dat$gene_hap %in% paste0(n_hap_any_ht$gene, "_", n_hap_any_ht$Hap1) & dat$best_pair==1, ],
#                      dat[dat$gene_hap %in% paste0(n_hap_any_ht$gene, "_", n_hap_any_ht$Hap2) & dat$best_pair==1, ] )
#     dat_4 <- unique(dat_4)
# 
# 
# 
#     selected_genes_readCount <- merge(unique(selected_genes[,c("gene_name", "only_h1", "only_h2", "intersect", "type")]), unique(dat[,c("gene_name", "ord")]), by="gene_name") 
#     selected_genes_readCount <- selected_genes_readCount[order(selected_genes_readCount$ord, decreasing = TRUE), ]
#     selected_genes_readCount$type <- gsub("(.*)_multi", "\\1", selected_genes_readCount$type)
#     selected_genes_readCount$type <- gsub("ht", "heterozygous", selected_genes_readCount$type)
#     selected_genes_readCount$type <- gsub("hm", "homozygous", selected_genes_readCount$type)
#     colnames(selected_genes_readCount) <- c("gene_name", "Hap1only", "Hap2only", "Hap1andHap2", "Zygosity", "ord")
#     selected_genes_readCount <- melt(selected_genes_readCount, id=c("gene_name", "ord", "Zygosity") )
#     colnames(selected_genes_readCount) <- c("gene_name", "ord",  "Zygosity", "Haplotype", "value")



    ##########################################
    ### matching output based on top1path
    #######q###################################
    
    zyg <- subset(dat, best_pair==1)
    zyg_count <- data.frame(table(zyg$gene_name))
    

    #hm
    dat_1 <- dat[dat$gene_name %in% as.character(zyg_count$Var1[zyg_count$Freq==1]), ]
    dat_1 <- dat_1[dat_1$gene_hap %in% zyg$gene_hap,]

    #ht
    dat_2 <- dat[dat$gene_name %in% as.character(zyg_count$Var1[zyg_count$Freq==2]), ]
    dat_2 <- dat_2[dat_2$gene_hap %in% zyg$gene_hap,]


    # multiple hm
    dat_3 <- dat_1[dat_1$gene_name %in% n_hap_any_hm$gene_name, ]

    # multiple ht
    dat_4 <- dat_2[dat_2$gene_name %in% n_hap_any_ht$gene_name, ]


    dat_1 <- dat_1[!(dat_1$gene_name %in% dat_3$gene_name), ]
    dat_2 <- dat_2[!(dat_2$gene_name %in% dat_4$gene_name), ]


 

        paired_mapping_rates$g_h1_h2 <- paste0(paired_mapping_rates$gene_name, "_", paired_mapping_rates$Hap1, "_", paired_mapping_rates$Hap2) 

        selected_genes_by_top1path <- paired_mapping_rates[(paired_mapping_rates$g_h1_h2 %in% top1_edges$V1), ]
        type <- ifelse(as.character(selected_genes_by_top1path$Hap1)==as.character(selected_genes_by_top1path$Hap2), "homozygous", "heterozygous")     
        selected_genes_by_top1path$type= type



        # to add all identical genes in the graph
#         identicalGenes_1 <- c()
#         identicalGenes_2 <- c()
#         for(i in 1:nrow(identicalGeneList)){
#               g1 <- as.character(identicalGeneList$gene_name[i])
#               tmp <- subset( dat, gene_name==as.character(identicalGeneList$gene_name[i]) )
#               tmp <- tmp[as.character(tmp$hap) %in% strsplit(as.character(identicalGeneList$haps[i]), ";")[[1]],]
#               if(g1==g2){
#                     identicalGenes_2 <- rbind(identicalGenes_2, tmp)
#               }             
#               if(g1!=g2){
#                     identicalGenes_1 <- rbind(identicalGenes_1, tmp)
#               }
#               g2 <- as.character(identicalGeneList$gene_name[i])
#         }

    


          # identicalGenes_1 <- c()
          # identicalGenes_2 <- c()
          # for(i in 1:nrow(identicalGeneList)){
          #       g1 <- as.character(identicalGeneList$gene_name)[i]
          #       
          #       tmp <- subset( dat, gene_name==g1)
          #       tmp <- tmp[as.character(tmp$hap) %in% strsplit(as.character(identicalGeneList$haps[i]), ";")[[1]],]
          #       
          #       if( nrow(subset(tmp, best_pair==1) ) >0 ){
          #             if(g1==g2){
          #                   identicalGenes_2 <- rbind(identicalGenes_2, tmp)
          #             }
          #             if(g1!=g2){
          #                   identicalGenes_1 <- rbind(identicalGenes_1, tmp)
          #             }
          #       }
          #       g2=g1
          # }



        

        selected_genes_readCount <- merge(unique(selected_genes_by_top1path[,c("gene_name", "only_h1", "only_h2", "intersect", "type")]), unique(dat[,c("gene_name", "ord")]), by="gene_name")       
        selected_genes_readCount <- selected_genes_readCount[order(selected_genes_readCount$ord, decreasing = TRUE), ]
        colnames(selected_genes_readCount) <- c("gene_name", "Hap1only", "Hap2only", "Hap1andHap2", "Zygosity", "ord")

        selected_genes_readCount <- melt(selected_genes_readCount, id=c("gene_name", "ord", "Zygosity") )
        colnames(selected_genes_readCount) <- c("gene_name", "ord",  "Zygosity", "Haplotype", "value")

        
        ### refine output
        dat_final <- c()
        selected_genes_readCount_final <- c()
        for(i in unique(as.character(selected_genes_readCount$gene_name)) ){
          
              a <- subset(dat, as.character(gene_name)==i)
          
              tmp <- subset(selected_genes_readCount, as.character(gene_name)==i)
              tmp <- subset(tmp, value!=0 )
              if(nrow(tmp)==1){
                    tmp$Zygosity <- "homozygous"
                    tmp$Haplotype <- "Hap1andHap2"
                    a$best_pair[which(a$best_pair==1)][2] <- 0
                    
              }
              if(nrow(tmp)==2){
                
                      if( tmp$Haplotype[1]=="Hap1only" & tmp$Haplotype[2]=="Hap2only" ){
                      }
                      else if( any(grepl("Hap1only", tmp$Haplotype)) ){
                            tmp <- data.frame(tmp[1,1:2], 
                                              Zygosity="homozygous", Haplotype="Hap1andHap2",
                                              value=sum(tmp$value))
                            a$best_pair[which(a$best_pair==1)][2] <- 0
                      }
                      else if( any(grepl("Hap2only", tmp$Haplotype) )  ){
                            tmp <- data.frame(tmp[1,1:2], 
                                              Zygosity="homozygous", Haplotype="Hap1andHap2",
                                              value=sum(tmp$value))
                            a$best_pair[which(a$best_pair==1)][1] <- 0
                      }
              }
              
              selected_genes_readCount_final <- rbind(selected_genes_readCount_final, tmp)
              dat_final <- rbind(dat_final, a)
        }

        
        permutation_result$gene_hap <- paste0(permutation_result$gene_name, "_", permutation_result$hap)
        dat_final <- merge(dat_final, permutation_result[,c("gene_hap", "pvalues")], by="gene_hap", all=TRUE)
        dat_final <- dat_final[!is.na(dat_final$best_pair), ]
        

        dat_0 <- subset(dat_final, best_pair==1)
        #tmp <- c( paste0(selected_genes$gene_name, "_", selected_genes$Hap1), paste0(selected_genes$gene_name, "_", selected_genes$Hap2) )
        #tmp <- unique(tmp)
        #dat_0 <- dat_0[dat_0$gene_hap %in% tmp,]
        

        #head(paired_mapping_rates)
        #head(dat_0)
        
        final_set <- c()
        for(i in unique(as.character(dat_0$gene_name))){
                a <- subset(dat_0, gene_name==i)
                if(nrow(a)==2){
                      b <- subset( paired_mapping_rates, as.character(gene_name)==i & as.character(Hap1)==as.character(a$hap)[1] & as.character(Hap2)==as.character(a$hap)[2] )
                }
                if(nrow(a)==1){
                  b <- subset( paired_mapping_rates, as.character(gene_name)==i & as.character(Hap1)==as.character(a$hap)[1] & as.character(Hap2)==as.character(a$hap)[1] )
                }
                final_set <- rbind(final_set, b)
        }
        
   
      final_set <- merge(final_set, unique(selected_genes_readCount_final[,c("gene_name", "ord", "Zygosity")]), by="gene_name")
      
      
       summary_set <- final_set[, c("gene_name", "ord", "mapping_rate", "rates_of_gene_counts", "Zygosity")]
      # #colnames(summary_set) <- c("gene_name", "ord", "mapping_rate", "rates_of_gene_counts", "Zygosity")
       summary_set <- melt(summary_set, id.vars = c("gene_name", "ord", "Zygosity"))
      
      
      summary_heatmap <- ggplot(summary_set, aes(x = variable, y = ord , fill=Zygosity)) + 
                                geom_tile() + scale_fill_brewer(palette="Pastel1") +
                                geom_tile(data=summary_set, aes(x=variable, y=ord, color=Zygosity),size=0,fill=NA, color="grey") +
                                geom_text(aes(label=round( (value*100),2) ), size=5, color="black" ) +
                                scale_y_continuous(expand=c(0,0),breaks=min(summary_set$ord):max(summary_set$ord),labels=NULL) +
                                #theme(axis.title.y=element_blank(), axis.text.y=element_blank()) +
                                ggtitle("Combined mapping rates") + 
                                theme(legend.position="none") + ylab(NULL) + xlab("Mapping rates") 
                          
      
      
   # gene_readCount <-  ggplot(selected_genes_readCount_final, aes(x=ord, y=value, fill=Haplotype, order=desc(Haplotype ))) +
    gene_readCount <-  ggplot(selected_genes_readCount_final, aes(x=ord, y=value, fill=Haplotype, order=Haplotype )) + 
                        geom_bar(stat="identity") + 
                        scale_x_continuous(expand=c(0,0),breaks=min(selected_genes_readCount_final$ord):max(selected_genes_readCount_final$ord),labels=gene_labels) +
                        #theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + 
                        theme(axis.text.y=element_text(size = 12, face = "bold") ) +
                        ylab("Read counts") + xlab(NULL) + coord_flip() +
                        ggtitle(paste0( "Gene counts" ) ) 
                        

    
    mhc_heatmap <- ggplot(dat_final, aes(x = hap , y = ord , fill=pvalues)) + 
                      geom_tile() +
                      #scale_fill_gradient() +
                      scale_fill_gradient2(low = "darkblue", mid = "blue", high = "lightblue") +
                      scale_y_continuous(expand=c(0,0),breaks=min(dat$ord):max(dat$ord),labels=gene_labels) +
                      xlab("Haplotypes") + ylab("Genes by the chromosomal order") +
                      geom_tile(data=dat_final, aes(x=hap, y=ord, color=best_pair),size=0,fill=NA, color="white") +
                      #geom_tile(data=identicalGenes_1, aes(x=hap, y=ord),size=2,fill=NA, color="lightgrey") + 
                      #geom_tile(data=identicalGenes_2, aes(x=hap, y=ord),size=2,fill=NA, color="darkgrey") + 
                      geom_tile(data=dat_0, aes(x=hap, y=ord, color=best_pair),size=1,fill=NA, color="red") +    # was blue
                      #geom_tile(data=dat_4, aes(x=hap, y=ord, color=best_pair),size=1,fill=NA, color="red") +    # was blue
                      #geom_tile(data=dat_3, aes(x=hap, y=ord, color=best_pair),size=2,fill=NA, color="red") +    # was blue
                      #geom_tile(data=dat_2, aes(x=hap, y=ord, color=best_pair),size=1,fill=NA, color="red") +
                      #geom_tile(data=dat_1, aes(x=hap, y=ord, color=best_pair),size=2,fill=NA, color="red") +
                      geom_text(aes(label=round( (Mapping_Rates*100),2) ), size=5, color="white" ) +
                      ggtitle(paste0("Predicted haplotypes and mapping rates") ) + 
                      theme(legend.position="left", axis.text.y=element_text(size = 12, face = "bold"), axis.title.y = element_text(face="bold", size=18)) 

  
    
        gp1<- ggplot_gtable(ggplot_build( gene_readCount ))
        gp2<- ggplot_gtable(ggplot_build( mhc_heatmap ))
        gp3<- ggplot_gtable(ggplot_build( summary_heatmap )) 
        
        #maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
        #gp1$widths[2:3] <- maxWidth
        #gp2$widths[2:3] <- maxWidth
        
        #g_p12 <- grid.arrange(gp2, gp1, ncol=2, widths = c(2/3,1/3) )
        g_p123 <- grid.arrange(gp2, gp3, gp1, ncol=3, widths = c(8/15, 2/15, 5/15) )
        
 
               
        
        # testing
        #test <- data.frame(table(paste0(selected_genes_by_top1path$Hap1, "_", selected_genes_by_top1path$Hap2)) )
        #print(test)
	      #write.table(test, paste0(sample_name, ".stats.tmp.txt"), quote=F)


        #subset(selected_genes_by_top1path, gene_name=="RPP21")
        #subset(selected_genes_by_top1path, gene_name=="C6orf136")
        #subset(selected_genes_by_top1path, gene_name=="AGER")
        
        
        
    ggsave(filename=paste0(sample_name, ".mhc_heatmap.pdf" ), g_p123, width=18, height=22)
    #ggsave(filename=paste0(sample_name, ".mhc_heatmap.pdf" ), g_p12)


    output <- list(heatmap_genes=dat,
                   #summary_of_best_pairs=selected_genes_by_top1path,
                   summary_of_best_pairs=final_set)

}

