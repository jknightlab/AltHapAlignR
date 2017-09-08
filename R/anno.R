#' get exon data
#'
#' It produces min.info <- c("gene_id", "transcript_id", "gene_type", "gene_status", "gene_name", "transcript_type", "transcript_status",  "transcript_name", "exon_number", "exon_id", "level")
#' @param gtf_granges GRanges object
#' @export


getExons <- function(gtf_granges){
  
  # get subset with exons
  requiredMeta <- c("type", "group")
  presentMeta <- requiredMeta %in% names(mcols(gtf_granges)) 
  if(!all(presentMeta)){
	  idx <- which(!presentMeta)
	  stop("Required columns ", paste(presentMeta[idx], collapse=", "), " are missing.")
  }
  min.info <- c("gene_id", "transcript_id", "gene_type", "gene_status", "gene_name", "transcript_type", "transcript_status",  "transcript_name", "exon_number", "exon_id", "level")
  gtf.exons <- subset(gtf_granges, gtf_granges$type=="exon") 
  
  if(length(gtf.exons) >0){
    # get exon attributes
    message("Parsing gene/transcript/exon ids.")
    exon.info <- gsub("\"|;", "", gtf.exons$group)   
    exon.info <- strsplit(gsub("  ", " ", exon.info), split=" ")
	if(any(sapply(exon.info, length) < length(min.info))) 
		stop("GTF file is missing annotations")
    exon.info <- lapply(exon.info, function(x) {
      data <- x[seq(2, length(x), 2)]
      names(data) <- x[seq(1, length(x), 2)]
      data })
    
    attribs <- names(exon.info[[1]])
    if (!all(min.info %in% attribs)) stop("Not all required annotations are in this GTF file.")
    exon.info <- do.call(data.frame, c(lapply(attribs, function(x) sapply(exon.info, "[", x)), stringsAsFactors = FALSE))
    colnames(exon.info) <- attribs
    
    gtf.exons <- cbind(ranges(gtf.exons), strand=strand(gtf.exons), exon.info)
  } else {
    warning("no exon found ...")
	gtf.exons <- data.frame(matrix(nrow=0, ncol=4+length(min.info)))
	names(gtf.exons) <- c("start", "end", "width", "strand", min.info)
  }
  return(gtf.exons)
}



#' Convert a GTF file to an exon table
#'
#' This will produce a table of parsed exon data from gtf file. It uses a function 'getExons'. 
#' @param gtf Transcript annotations. Can either be provided as a GRanges object or as the name of
#' a GTF file containing the annotations.
#' @param refver human genome reference version
#' @param grange Character vector of the form \code{c(chromosome, start, end)}, 
#' indicating the genomic region from which exons should be considered. Set \code{start}
#' or \code{end} to 0 to indicate that the region should extend to the respective end 
#' of the chromosome. Setting this to \code{NULL} will disable subsetting altogether. 
#'
#' @export
#' @import rtracklayer
#' @import GenomicRanges
#' @import GenomicFeatures
#' 
#' @examples
#' gtf=system.file("extdata", "gencode.v21.hla.annotation.gtf", package = "HapTag")
#' gtf2ExonTable(gtf, refver="hg38", grange=c("chr6", "28528216", "33389373") )
#' 
#' # to use HapTagAnno containing human reference genomes and parsed data.
#' library(HapTagAnno.hg38)
#' gtf2ExonTable(gtf, refver="hg38", grange=c("chr6", "28528216", "33389373") )

### hg19 cox : SCAND3 ~ DAXX
### hg19 pgf : from SCAND3 ( right after GPX5) to DAXX (right before KIFC1 ) : 28502730 ~ 33359312

### hg38 cox : GPX5 ~ LYPLA2P1
### hg38 pgf : from GPX5 ( right after GPX6 ) to LYPLA2P1 (right before RPL35AP4) : 28524925 ~ 33389373

gtf2ExonTable <- function(gtf, refver="hg38", grange=c("chr6", "28524925", "33389373") ){
  chr <- ""
  start <- 0
  end <- 0
  if(!is.null(grange)){
  	chr=as.character(grange[1])
  	start=as.numeric(grange[2])
  	end=as.numeric(grange[3])
  }
  
  message("The reference genome version ", refver, "(", chr, ":", start, "-", end, ")") 

  if(is(gtf, "connection") || is.character(gtf)){
      gtf <- import.gff(gtf, asRangedData = FALSE)
  }
  genome(gtf)[is.na(genome(gtf))] <- refver
  if(any(genome(gtf) != refver)){
	  warning("Reference version mismatch. Expected ", refver, " found ", 
			  genome(gtf)[which(genome(gtf) != refver)[1]], ".")
  }

  if(length(gtf) >0){
	if(length(chr) > 0 && chr != ""){
		gtf <- gtf[seqnames(gtf) == chr]
	}
	if(start > 0){
		gtf <- gtf[start(gtf) >= start]
	}
	if(end > 0){
		gtf <- gtf[end(gtf) <= end]
	}  	
  }else{
	warning("no annotation provided ...")
  }
  exonTable <- getExons(gtf)

  return(exonTable)
}



#' Filter and merge a exon table
#'
#' This filers exon data by annotation level and merges exons across multiple isoforms. 
#' @param exon_table output from 'gtf2ExonTable'
#' @param seqname chromosome name, it uses only single chromosome
#' @param levelBelow, evidence level of annotation. If levelBelow=3, it uses only level 1 and 2 data. levels: 1(verified loci), 2(manually annotated loci), 3(automatically annotated loci)
#' @param transcript_type filter exons by specific biotypes. transcript_type="all" uses all data. default: 'protein_coding', Further details about transcript biotypes are here. http://www.gencodegenes.org/gencode_biotypes.html
#' @export
#' @import GenomicRanges
#' @examples
#' # to use your own reference genome
#' gtf=system.file("extdata", "gencode.v21.hla.annotation.gtf", package = "HapTag")
#' exTab <- gtf2ExonTable(gtf, refver="hg38", grange=c("chr6", "28528216", "33389373"))
#' filtered <- filterAndMergeExons(exTab)
#' 
#' # to use HapTagAnno containing human reference genomes and parsed data.
#' library(HapTagAnno.hg38)
#' gtf <- system.file("extdata/gencode.v21.hla.annotation.gtf", package="HapTag")
#' exTab <- gtf2ExonTable(gtf, refver="hg38", grange=c("chr6", "28528216", "33389373"))
#' filtered <- filterAndMergeExons(exTab)
#' 

filterAndMergeExons <- function( exon_table, seqname="chr6" , levelBelow=3, transcript_type="protein_coding"){
  # exon_table: output from 'gtf2ExonTable', column names should be the same with output from 'gtf2ExonTable'
  # levelBelow=3: fileter exons with level 1 or 2
  # transcript_type: single type or all option
  
  num_exons_before_filter <- nrow(exon_table)
  num_genes_before_filter <- length(unique(exon_table$gene_id))
  
  if(transcript_type=="all"){
    exon_table <- subset(exon_table, exon_table$level < levelBelow)
  }else {
    exon_table <- subset(exon_table, exon_table$level < levelBelow & transcript_type==transcript_type )
  }
  
  num_exons_after_filter <- nrow(exon_table)
  num_genes_after_filter <- length(unique(exon_table$gene_id))
  geneids <- unique(exon_table$gene_id)
  geneids <- as.character(geneids)
  message("Number of unique genes: ", length(geneids))
  
  stats <- cbind(c(num_exons_before_filter, num_exons_after_filter),
                 c(num_genes_before_filter, num_genes_after_filter) )
  
  all_exonsMerged <- list()
  for(g in geneids){
    sub <- unique(cbind(seqnames=seqname, subset(exon_table, exon_table$gene_id==g, c("start", "end", "strand", "gene_name", "exon_id" )) ) )
    
    gene_granges <- makeGRangesFromDataFrame(sub, start.field="start", end.field="end", strand.field="strand", keep.extra.columns=TRUE)
    g_name <- as.character(gene_granges$gene_name[1])
    g_strand <- as.character(gene_granges@strand@values)
    message(g_name)
    
    if(length(gene_granges) ==1){
      tmp <- data.frame( gene_granges@ranges)
      exonsMerged <- cbind(1, tmp$start, tmp$end, tmp$width, g_strand, g_name, 1, as.character(gene_granges$exon_id))
    }
    else{
      by_exons <- unlist(split(gene_granges, gene_granges@elementMetadata$exon_id))
      overlappedExons <- as.matrix(findOverlaps(by_exons[1:length(by_exons)], by_exons[1:length(by_exons)], type="any", select="all" ) )
      overlappedExons <- data.frame(overlappedExons)
      
      exonsMerged <- sapply( names(table(overlappedExons$subjectHits)), function(x) {			
        tmp <- subset(overlappedExons, overlappedExons$queryHits==x, c("subjectHits") )
        hits_group <- unique(overlappedExons[overlappedExons$queryHits %in% tmp$subjectHits, c("subjectHits") ])
        hits_group <- hits_group[order(hits_group)]
        
        hits_group <- unique(overlappedExons[overlappedExons$queryHits %in% hits_group, c("subjectHits") ])
        hits_group <- hits_group[order(hits_group)]
        
        hits_group <- unique(overlappedExons[overlappedExons$queryHits %in% hits_group, c("subjectHits") ])
        hits_group <- hits_group[order(hits_group)]
        
        tmp <- by_exons[hits_group]
        tmp <- data.frame(tmp@ranges, tmp$exon_id)
        tmp <- c( min(tmp$start), max(tmp$end), (max(tmp$end) - min(tmp$start) +1), g_strand, g_name, nrow(tmp),
                  do.call(paste, c(as.list(tmp$tmp.exon_id) , sep=";") ) )
      })
      exonsMerged <- unique(data.frame(t(exonsMerged)) )
      exonsMerged <- exonsMerged[order(exonsMerged[,1]),]
      
      exonsMerged <- data.frame(exn=as.character(c(1:nrow(exonsMerged))),
                                start=as.character(exonsMerged[,1]),
                                end=as.character(exonsMerged[,2]),
                                width=as.character(exonsMerged[,3]),
                                gene_strand=as.character(exonsMerged[,4]),
                                gene_name=as.character(exonsMerged[,5]),
                                n_of_merged_exons=as.character(exonsMerged[,6]),
                                merged_exon_ids=as.character(exonsMerged[,7]), stringsAsFactors=FALSE )
      
    }
    colnames(exonsMerged) <- c("exn", "start", "end", "width", "gene_strand", "gene_name", "n_of_merged_exons", "merged_exon_ids")
    gid_name <- paste(exonsMerged[,6][1], g, sep="_")
    message(gid_name)
    all_exonsMerged[[gid_name]] <- exonsMerged
  }
  return(all_exonsMerged)
} 


#' fasta2DNAStringSet
#'
#' converting fasta format to DNAStringSet
#' @param fasta fastq format data
#' @export
#' @import Biostrings

fasta2DNAStringSet <- function(fasta){
  x23 <- readDNAStringSet(fasta, format="fasta", use.names=TRUE)
  return(x23)
}

