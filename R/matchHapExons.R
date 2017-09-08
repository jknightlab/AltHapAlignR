
#' get variants between two aligned sequences
#'
#' @param alignment alignment
#' @param snp snp
#' @export
#' @import Biostrings

getVariantsByAlignment__ <- function(alignment, snp){ 
  seq1aln <- pattern(alignment)
  seq2aln <- subject(alignment)
  
  output <- c()
  alnlen  <- nchar(seq1aln)
  for (i in 1:alnlen) {
    n1 <- substring(seq1aln, i, i)
    n2 <- substring(seq2aln, i, i)
    
    if(n1 != n2){
      if(n1 == "-"){
        output <- rbind(output, data.frame(alnpos=i, hap1_nt=n1, hap2_nt=n2, hap1="del", hap2="in"))
      }
      else if(n2 == "-"){
        output <- rbind(output, data.frame(alnpos=i, hap1_nt=n1, hap2_nt=n2, hap1="in", hap2="del"))
      }
      else{
        output <- rbind(output, data.frame(alnpos=i, hap1_nt=n1, hap2_nt=n2, hap1="snp", hap2="snp"))
      }
    }
  }
  return(output)
}


#' load snp data to find variants by aligned sequences (need to change to load dbSNP for only two haplotypes, not entire dbSNP)
#'
#' @param refver a version of reference sequence, hg38 and hg19 versions are available from HapTagAnno.hg38 and HapTagAnno.hg19, respectively.
#' @param seq1 a chromosome name of sequence 1
#' @param seq2 a chromosome name of sequence 2
#' @export
#' @examples
#' # for reference hg38
#' require(HapTagAnno.hg38)
#' # for PGF and COX haplotypes
#' loadDataForVariantsByAlignment(refver="hg38", seq1="chr6", seq2="GL000251.2")
#'
#' # for reference hg19
#' \dontrun{
#' require(HapTagAnno.hg19)
#' # for PGF and COX haplotypes
#' loadDataForVariantsByAlignment(refver="hg19", seq1="chr6", seq2="GL000251.1")
#' }

loadDataForVariantsByAlignment <- function(refver="hg38", seq1="chr6", seq2="GL000251.2"){
  if(refver=="hg38"){
    if(!exists("refseq")) {
      #require(HapTagAnno_hg38)
      message("loading reference sequences...")
      data(ref_hg38)
    }
    if(!exists("dbSNP")) {
      #require(HapTagAnno_hg38)     
      message("loading dbSNP141 data...")
      data(dbSNP141_hg38)
    }
  }
  message("Done...")
}


#' get variants between two sequences after pairwiseAlignment
#'
#' @param refver a version of human genome reference (it may not need here...)
#' @param refseq reference sequences. If no reference sequence is supplied the \code{ref_MHC} object
#' from the \code{HapTagAnno} package that was last loaded will be used.
#' @param hap1_dbSNP dbSNP data of haplotype 1, should be loaded before running the function
#' @param hap2_dbSNP dbSNP data of haplotype 2, should be loaded before running the function
#' @param seq1 position information of first sequence, c(chromsome name, start position, end position, strand)
#' @param seq2 position information of second sequence, c(chromsome name, start position, end position, strand)
#' @param hap_names names of two haplotypes
#' @export
#' @import Biostrings
#' @import sqldf
#' @importFrom dplyr left_join
#' @examples
#' # load reference sequences and  
#' library(HapTagAnno.hg38)
#' # load snp data
#' hap1_dbSNP <- subset(dbSNP_MHC, dbSNP_MHC$chr=="chr6")
#' hap2_dbSNP <- subset(dbSNP_MHC, dbSNP_MHC$chr=="chr6_GL000251v2_alt")
#' # load reference sequences and dbSNP by a function 'loadDataForVariantsByAlignment'
#' loadDataForVariantsByAlignment(refver="hg38", seq1="chr6", seq2="GL000251.2")
#'

getVariantsByAlignment <- function(refver="hg38", refseq, hap1_dbSNP, hap2_dbSNP, seq1, seq2, hap_names=c("hap1", "hap2")){

  message("starting...")
  
  if(missing(refseq)) refseq <- ref_MHC
  sameSeq <- data.frame()
  hapSpe <- data.frame()
  hap_mismatch <- data.frame()
  hap_indel <- data.frame()
  
  col_names <- c("h1_start", "h1_end", "h1_seq", "h2_start", "h2_end", "h2_seq", "h1_type", "h2_type",
                   "h1_chr", "h1_rs_start", "h1_rs_end", "h1_rsID", "h1_strand", "h1_refNCBI", "h1_refUCSC", "h1_observed", "h1_class",
                   "h2_chr", "h2_rs_start", "h2_rs_end", "h2_rsID", "h2_strand", "h2_refNCBI", "h2_refUCSC", "h2_observed", "h2_class")

  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -2, baseOnly = TRUE)
  
  dnastring1 <- subseq(refseq[[seq1[1]]], start=as.numeric(seq1[2]), end=as.numeric(seq1[3]) )
  dnastring2 <- subseq(refseq[[seq2[1]]], start=as.numeric(seq2[2]), end=as.numeric(seq2[3]) )
  
  # 1) the same merged exons        
  if(dnastring1 == dnastring2){
    message("Two exons are identical.")
    sameSeq <- data.frame(h1_start=seq1[2], h1_end=seq1[3], h1_seq="",
                          h2_start=seq2[2], h2_end=seq2[3], h2_seq="",
                          h1_type="s", h2_type="s" )
  }
  
  if(dnastring1 != dnastring2){
    message("Two exons are not identical.")
    
    if( length(dnastring1)==length(dnastring2) ){
      align <- pairwiseAlignment(dnastring1, dnastring2,  type = "global", substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
    }
    if( length(dnastring1)!=length(dnastring2) ){
      align <- pairwiseAlignment(dnastring1, dnastring2,  type = "local", substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
    }
    
    seq1aln <- pattern(align)
    seq2aln <- subject(align)
    
    tmp <- data.frame()
    alnlen  <- nchar(seq1aln)
    tmp_h1_pos=0
    tmp_h2_pos=0
    for (i in 1:alnlen) {
      n1 <- substring(seq1aln, i, i)
      n2 <- substring(seq2aln, i, i)
      if(n1 != n2){
        if(n1 == "-"){
          tmp_h2_pos=tmp_h2_pos+1
          h1_pos=as.numeric(as.numeric(seq1[2]) + tmp_h1_pos-1 + seq1aln@range@start -1)
          h2_pos=as.numeric(as.numeric(seq2[2]) + tmp_h2_pos -1 + seq2aln@range@start -1)
          
          tmp <- rbind(tmp, data.frame(alnpos=i, seq1aln_pos=as.numeric(tmp_h1_pos), seq2aln_pos=as.numeric(tmp_h2_pos),
                                       h1_start=h1_pos, h1_end=h1_pos, h1_seq=n1, 
                                       h2_start=h2_pos, h2_end=h2_pos, h2_seq=n2,
                                       h1_type=paste("del",as.numeric(tmp_h1_pos), sep="_"), h2_type="in"))
        }
        else if(n2 == "-"){
          tmp_h1_pos=tmp_h1_pos +1
          h1_pos=as.numeric(as.numeric(seq1[2]) + tmp_h1_pos-1 + seq1aln@range@start -1)
          h2_pos=as.numeric(as.numeric(seq2[2]) + tmp_h2_pos-1 + seq2aln@range@start -1)
          
          tmp <- rbind(tmp, data.frame(alnpos=i, seq1aln_pos=as.numeric(tmp_h1_pos), seq2aln_pos=as.numeric(tmp_h2_pos),
                                       h1_start=h1_pos, h1_end=h1_pos, h1_seq=n1, 
                                       h2_start=h2_pos, h2_end=h2_pos, h2_seq=n2,
                                       h1_type="in", h2_type=paste("del", as.numeric(tmp_h2_pos), sep="_")))
        }
        else{
          tmp_h1_pos=tmp_h1_pos+1
          tmp_h2_pos=tmp_h2_pos+1
          h1_pos=as.numeric(as.numeric(seq1[2]) + tmp_h1_pos-1 + seq1aln@range@start -1)
          h2_pos=as.numeric(as.numeric(seq2[2]) + tmp_h2_pos-1 + seq2aln@range@start -1)
          
          tmp <- rbind(tmp, data.frame(alnpos=i, seq1aln_pos=as.numeric(tmp_h1_pos), seq2aln_pos=as.numeric(tmp_h2_pos),
                                       h1_start=h1_pos, h1_end=h1_pos, h1_seq=n1,
                                       h2_start=h2_pos, h2_end=h2_pos, h2_seq=n2,
                                       h1_type="snp", h2_type="snp"))
        }
      }
      else{
        tmp_h1_pos=tmp_h1_pos+1
        tmp_h2_pos=tmp_h2_pos+1
      }
    }
    
    if(nrow(tmp) >0){
      # 2) mismatches
    
      snps <- subset(tmp, tmp$h1_type=="snp")
      if(nrow(snps) >0){
        message("single variants...")
        h <- subset(hap1_dbSNP, hap1_dbSNP$class=="single")
        a1 <- left_join(snps, h, by=c("h1_end" = "end"))
        a1 <- data.frame(a1[,4:13], end=a1$h1_end, a1[,14:19])
        h <- subset(hap2_dbSNP, hap2_dbSNP$class=="single")
        a2 <- left_join(a1, h, by=c("name"="name"))
        hap_mismatch <- data.frame(a2[,1:20], name=a2$name, a2[,21:25])
        names(hap_mismatch) <- col_names
     }
    
      # 3)indels
    
      ids <- unique ( c(as.character(tmp$h1_type), as.character(tmp$h2_type) ) )  
      ids <- ids[grep("del", ids)]
    
      if(length(ids) >0 ){
        message("indels...")
        a <- data.frame()
        for( i in ids){
          x <- subset(tmp, tmp$h1_type==i | tmp$h2_type==i)
          a <- rbind( a, 
                  data.frame(  h1_start=min(x$h1_start), h1_end=max(x$h1_end), h1_seq=paste(x$h1_seq, collapse=""),
                               h2_start=min(x$h2_start), h2_end=max(x$h2_end), h2_seq=paste(x$h2_seq, collapse=""),
                               h1_type=x$h1_type[1], h2_type=x$h2_type[1] ) )	
        }
    
        h <- subset(hap1_dbSNP, hap1_dbSNP$class!="single") 
        a1 <- " select a.*, h.* from a LEFT JOIN h ON ( h.start=(a.h1_start-1) or h.start=a.h1_start ) and h.end=a.h1_end and h.class!='single' "
        a2 <- sqldf(a1)
        h <- subset(hap2_dbSNP, hap2_dbSNP$class!="single") 
        a3 <- " select a2.*, h.* from a2 LEFT JOIN h ON ( h.start=(a2.h2_start-1) or h.start=a2.h2_start ) and h.end=a2.h2_end and h.class!='single' "
        hap_indel <- sqldf(a3)
        colnames(hap_indel) <- col_names
      }
    
    }
    
    # 4) hapSpecificSequence
    message("haplotype-specific sequences...")
    spe <- rbind(
		data.frame( 	exonStart=1, exonEnd=(seq1aln@range@start-1), 
				start=as.numeric(seq1[2]), end=(seq1aln@range@start-1 + as.numeric(seq1[2])-1), type="h1S"), 
                data.frame(	exonStart=(seq1aln@range@start+seq1aln@range@width), exonEnd=length(dnastring1),
				start=(seq1aln@range@start+seq1aln@range@width + as.numeric(seq1[2])-1), end=as.numeric(seq1[3]), type="h1S"),
                data.frame( 	exonStart=1, exonEnd=(seq2aln@range@start-1),
				start=seq2[2], end=(seq2aln@range@start-1 + as.numeric(seq2[2])-1), type="h2S"),
                data.frame( 	exonStart=(seq2aln@range@start+seq2aln@range@width), exonEnd=length(dnastring2),
				start=(seq2aln@range@start+seq2aln@range@width + as.numeric(seq2[2])-1), end=as.numeric(seq2[3]), type="h2S" ) )
    hapSpe <- subset(spe, spe$start < spe$end)
    
  }
  
  if(nrow(sameSeq) >0){
      sameSeq <- data.frame(h1_name=hap_names[1], h2_name=hap_names[2], sameSeq)
  }
  if(nrow(hapSpe) >0){
      hapSpe <- data.frame(h1_name=hap_names[1], h2_name=hap_names[2], hapSpe)
  }
  if(nrow(hap_mismatch) >0){
      hap_mismatch <- data.frame(h1_name=hap_names[1], h2_name=hap_names[2], hap_mismatch)
  }
  if(nrow(hap_indel) >0){
      hap_indel <- data.frame(h1_name=hap_names[1], h2_name=hap_names[2], hap_indel)
  }
  
  message("done.")  
  output <- list(sameSeq=sameSeq, hapSpe=hapSpe, hap_mismatch=hap_mismatch, hap_indel=hap_indel)
}





#' pairwise sequence alignment
#'
#' It produce best hits of exon sequences between two haplotypes by pairwiseAlignment
#' @param hap1_exon an exon data.frame, an output from a function 'filterAndMergeExons', which column names are c("exn", "start", "end", "width", "gene_strand", "gene_name", "n_of_merged_exons", "merged_exon_ids")
#' @param hap2_exon another exon data.frame to compaire with hap1_exon table
#' @param chrs chromosome names of two haplotypes to compare
#' @param refseq 'DNAStringSet' of reference sequences, 
#' @import Biostrings
#' @export
#' @examples
#' # to use fasta file from 'HapTagAnno.hg38'
#' library(HapTagAnno.hg38)
#' fasta=ref_MHC   # human reference version hg38
#' 
#' # to use your own reference genome
#' fasta="/path/to/custom/reference/genome.fa"
#' 
#' gtf <- system.file("extdata", "gencode.v21.hla.annotation.gtf", package = "HapTag")
#' gtf2ExonTable(gtf, refver="hg38", grange=c("chr6", "28528216", "33389373"))


seqCompare <- function(hap1_exon, hap2_exon, chrs=c("chr6", "GL000251.2"), refseq){
	# hap1_exon, hap2_exon : output from a function 'filterAndMergeExons', columns include "exn","start","end","width","gene_strand","gene_name","n_of_merged_exons","merged_exon_ids"


        mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -2, baseOnly = TRUE)

	hap1_exon$start <- as.numeric(hap1_exon$start)
	hap2_exon$start <- as.numeric(hap2_exon$start)
        hap1_exon$end <- as.numeric(hap1_exon$end)
	hap2_exon$end <- as.numeric(hap2_exon$end)

        output <- c()
        if(nrow(hap1_exon) >0 & nrow(hap2_exon) >0 ){
                message("starting pairwiseAlignment")

                for( i in c(1:nrow(hap1_exon)) ){
                        rst <- c()
                        seq1 <- subseq(refseq[[chrs[1]]], start=hap1_exon$start[i], end=hap1_exon$end[i])
                        if(hap1_exon$gene_strand[i] == "-"){
                                seq1 <- reverseComplement(seq1)
                        }
                        for(j in c(1:nrow(hap2_exon)) ) {
                                seq2 <- subseq(refseq[[chrs[2]]], start=hap2_exon$start[j], end=hap2_exon$end[j])
                                if(hap2_exon$gene_strand[j] == "-"){
                                        seq2 <- reverseComplement(seq2)
                                }

				if( seq1==seq2) {
				     type="s"
			        }
				if( seq1!=seq2) {
                                        type="d"
                                }

				if( hap1_exon$width[i]==hap2_exon$width[j] ){
					align <- pairwiseAlignment(seq1, seq2,  type = "global", substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)
				 }
				 if( hap1_exon$width[i]!=hap2_exon$width[j] ){
					align <- pairwiseAlignment(seq1, seq2,  type = "local", substitutionMatrix = mat, gapOpening = -5, gapExtension = -2)	# not for HLA-A exon 4 
				 }
        
                                hap1_len <- align@pattern@unaligned@ranges@width
                                hap2_len <- align@subject@unaligned@ranges@width
                                hap1_mis <- length(align@pattern@mismatch[[1]])
                                hap2_mis <- length(align@subject@mismatch[[1]])
                                hap1_indel <- length(align@pattern@indel@unlistData@width)
                                hap1_sum <- sum(align@pattern@indel@unlistData@width)
                                hap2_indel <- length(align@subject@indel@unlistData@width)
                                hap2_sum <- sum(align@subject@indel@unlistData@width)

	                        hap1_range <- data.frame(align@pattern@range)
	                        hap2_range <- data.frame(align@subject@range)

                                rst_ <- data.frame(seq1_name=as.character(rownames(hap1_exon)[i]), seq2_name=as.character(rownames(hap2_exon)[j]), 
	                        			alignScore=as.numeric(align@score),
                                                    	seq1_exonNumber=as.numeric(hap1_exon$exn[i]), seq2_exonNumber=as.numeric(hap2_exon$exn[j]),
                                                    	seq1_len=as.numeric(hap1_len), seq2_len=as.numeric(hap2_len),
                                                    	seq1_misMatch=as.numeric(hap1_mis), seq2_misMatch=as.numeric(hap2_mis),
                                                    	seq1_indel=as.numeric(hap1_indel), seq1_sum=as.numeric(hap1_sum),
                                                    	seq2_indel=as.numeric(hap2_indel), seq2_sum=as.numeric(hap2_sum),
	                        			seq1_strand=as.character(hap1_exon$gene_strand[i]),
				                        seq2_strand=as.character(hap2_exon$gene_strand[j]),
							
			                        	seq1_aligned_start=as.numeric(hap1_range[,1]),
		                        		seq1_aligned_end=as.numeric(hap1_range[,2]), 
						        seq1_aligned_len=as.numeric(hap1_range[,3]),
			                        	seq2_aligned_start=as.numeric(hap2_range[,1]),
                                                    	seq2_aligned_end=as.numeric(hap2_range[,2]),  
                                                    	seq2_aligned_len=as.numeric(hap2_range[,3]),
							seq1_unaligned_start=as.numeric(hap1_exon$start[i]),
						        seq1_unaligned_end=as.numeric(hap1_exon$end[i]),
		              		                seq2_unaligned_start=as.numeric(hap2_exon$start[j]),
	              					seq2_unaligned_end=as.numeric(hap2_exon$end[j]),
	              					type=as.character(type) )

                                                    rst <- rbind(rst, rst_)
                        }
                        rst <- subset(rst, rst[,3]==max(rst[,3]) )
                        output <- rbind(output, rst)
                }
	    }
  
      return(output)
}

