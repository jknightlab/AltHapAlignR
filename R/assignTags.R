#' Replacing multiple patterns 
#'
#' Allow the use of multiple pattern/replacement pairs in a single function call.
#' @param pattern Character vector containing regular expressions.
#' @param replacement Character vector containing one replacement string for each 
#' regular expression in \code{pattern}.
#' @param x Character vector where matches are sought.
#' @param ... Further arguments passed to gsub.
#' @export

mgsub <- function(pattern, replacement, x, ...) {
	if (length(pattern)!=length(replacement)) {
		stop("pattern and replacement do not have the same length.")
	}
	result <- x
	for (i in 1:length(pattern)) {
		result <- gsub(pattern[i], replacement[i], result, ...)
	}
	result
}


#' To assign haplotype-specific tags
#' 
#' @param hap1 a bam file of haplotype 1
#' @param hap2 a bam file of haplotype 2
#' @param mapper mapping program used, only tophat or bwa format available 
#' @export
#' @import Rsamtools
#' @import GenomicAlignments
#' @importFrom data.table data.table
#' @examples
#' # Input data: bam files from mapping to two haplotype references  
#' hap1 <- system.file("extdata", "MHC_PGF.bam", package = "HapTag")
#' hap2 <- system.file("extdata", "MHC_COX.bam", package = "HapTag")
#' output <- tagReads(hap1, hap2, mapper="tophat")
#' 
#' # Output data
#' # unique number of all reads and paired-end reads
#' stat_uniqNumber <- output$stats$uniqueNumber
#' # Frequency of reads for each group 
#' stat_groups <- output$stats$groups
#' # Frequency of reads for each tag
#' stat_tagXG <- output$stats$tagXG
#' # Fragment ids with mapping tags and assigned tags.
#' fragmentIDs <- output$tagXG

tagReads <- function(hap1, hap2, mapper="tophat"){
	# both bam files are all filtered out with specific region, ex, mhc region
	# hap1: bam file of hap1
	# hap2: bam file of hap2
  message("Processing BAM files ", hap1, " and ", hap2, " (mapped with ", mapper, ")")
 
  # set parameters
  what= c("qname", "flag", "cigar")
  
  if(mapper=="tophat"){
  	tags <- c("NM", "NH", "HI")
  }
  if(mapper=="bwa"){
	tags <- c("NM", "X0", "X1")
  }
  flags <- scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE)
  param_hap1=ScanBamParam(what=what, tag=tags, flag=flags)
  param_hap2=ScanBamParam(what=what, tag=tags, flag=flags)


  # 1. read bam: all by scanBam
  
  message("Load all reads in the MHC region from haplotype_1 sample ")
  hap1_mhc  <- scanBam(hap1, param=param_hap1)
  message("Load all reads in the MHC region from haplotype_2 sample ")
  hap2_mhc <- scanBam(hap2, param=param_hap2)
  message("all fragments of haplotype_1 and haplotype_2 data loaded...")

  message(" Stats: unfiltered hits ...")
  temp_ <- data.frame(table(hap1_mhc[[1]]$qname))
  all_hap1_qnames <- data.frame(table(temp_$Freq))
  all_hap1_qnames <- data.table(hitsOfQname=all_hap1_qnames[,1], Freq.hap1=all_hap1_qnames[,2], key="hitsOfQname")
  temp_ <- data.frame(table(hap2_mhc[[1]]$qname))
  all_hap2_qnames <- data.frame(table(temp_$Freq))
  all_hap2_qnames <- data.table(hitsOfQname=all_hap2_qnames[,1], Freq.hap2=all_hap2_qnames[,2], key="hitsOfQname")
 
  
  # 2. read bam : paired reads by GenomicAlignments 

  message("Load paired reads from haplotype_1 sample...")
  hap1_paired <- readGAlignmentPairs(hap1, use.names=FALSE, param=param_hap1)
  hap1_dumped <- getDumpedAlignments()						# alignments with ambiguous are dropped because the seqname or strand of the alignments in the pair were not concordant
  message("Load paired reads from haplotype_2 sample...")
  hap2_paired <- readGAlignmentPairs(hap2, use.names=FALSE, param=param_hap2)
  hap2_dumped <- getDumpedAlignments()                             

  if(mapper=="bwa"){
	hap1_paired@first@elementMetadata$NH <- hap1_paired@first@elementMetadata$X0
	hap1_paired@last@elementMetadata$NH <- hap1_paired@last@elementMetadata$X0
	hap2_paired@first@elementMetadata$NH <- hap2_paired@first@elementMetadata$X0
	hap2_paired@last@elementMetadata$NH <- hap2_paired@last@elementMetadata$X0

	hap1_paired@first@elementMetadata$HI <- hap1_paired@first@elementMetadata$X1
	hap1_paired@last@elementMetadata$HI <- hap1_paired@last@elementMetadata$X1	
	hap2_paired@first@elementMetadata$HI <- hap2_paired@first@elementMetadata$X1
	hap2_paired@last@elementMetadata$HI <- hap2_paired@last@elementMetadata$X1
  }

  if(all(hap1_paired@first@elementMetadata$NH == hap1_paired@last@elementMetadata$NH) == FALSE){
    warning("Paried reads of haplotype_1 sample have different NH data.")
  }
  if(all(hap2_paired@first@elementMetadata$NH == hap2_paired@last@elementMetadata$NH) == FALSE){
    warning("Paried reads of haplotype_2 sample have different NH data.")
  }


  # get data.table of all paired set...
  message("Get data table of all paired reads...")
  dt.hap1.paired <- data.table(FragmentID=hap1_paired@first@elementMetadata$qname, 
                       hap1_cigar_1=hap1_paired@first@elementMetadata$cigar, 
                       hap1_nm_1=hap1_paired@first@elementMetadata$NM, 
                       hap1_nh_1=hap1_paired@first@elementMetadata$NH, 
                       hap1_hi_1=hap1_paired@first@elementMetadata$HI, 
                       hap1_cigar_2=hap1_paired@last@elementMetadata$cigar, 
                       hap1_nm_2=hap1_paired@last@elementMetadata$NM, 
                       hap1_nh_2=hap1_paired@last@elementMetadata$NH, 
                       hap1_hi_2=hap1_paired@last@elementMetadata$HI, 
                       key="FragmentID" )

  dt.hap2.paired <- data.table(FragmentID=hap2_paired@first@elementMetadata$qname, 
                       hap2_cigar_1=hap2_paired@first@elementMetadata$cigar, 
                       hap2_nm_1=hap2_paired@first@elementMetadata$NM, 
                       hap2_nh_1=hap2_paired@first@elementMetadata$NH, 
                       hap2_hi_1=hap2_paired@first@elementMetadata$HI, 
                       hap2_cigar_2=hap2_paired@last@elementMetadata$cigar, 
                       hap2_nm_2=hap2_paired@last@elementMetadata$NM, 
                       hap2_nh_2=hap2_paired@last@elementMetadata$NH, 
                       hap2_hi_2=hap2_paired@last@elementMetadata$HI, 
                       key="FragmentID" )


  temp_paired <- data.table(FragmentID=unique(hap1_paired@first@elementMetadata$qname), hap1.paired="1", key="FragmentID")
  temp_all_mhc <- data.table(FragmentID=unique(hap1_mhc[[1]]$qname), hap1.allmhc="1", key="FragmentID")
  hap1_all <- merge(temp_paired, temp_all_mhc, all=TRUE)

  temp_paired <- data.table(FragmentID=unique(hap2_paired@first@elementMetadata$qname), hap2.paired="1", key="FragmentID")
  temp_all_mhc <- data.table(FragmentID=unique(hap2_mhc[[1]]$qname), hap2.allmhc="1", key="FragmentID")
  hap2_all <- merge(temp_paired, temp_all_mhc, all=TRUE)
  
  # stat1
  message("Get unique number of reads from haplotype_1 and haplotype_2 sample...")
  uniqueNumber <- data.frame(haplotype_1=c(length(unique(hap1_mhc[[1]]$qname)), 
                                   length(unique(hap1_paired@first@elementMetadata$qname)) ),
                             haplotype_2=c(length(unique(hap2_mhc[[1]]$qname)),
                                   length(unique(hap2_paired@first@elementMetadata$qname)) ))
  rownames(uniqueNumber) <- c("all", "paired")
  message(uniqueNumber)
  

  # stat 2
  message("Compare reads between haplotype_1 and haplotype_2 sample...")
  # Summerize data ...
  all_fragments <- merge(hap1_all, hap2_all, all=TRUE)  # "ID", "hap1.paired", "hap1.allmhc", "hap2.paired", "hap2.allmhc"  
  all_fragments[is.na(all_fragments)] <- 0
  id.matrix <- data.frame(all_fragments, Num=c(1:dim(all_fragments)[1]), 
                          Gr=paste(all_fragments$hap1.paired, all_fragments$hap1.allmhc, all_fragments$hap2.paired, all_fragments$hap2.allmhc, sep="_"))
  
  Gr.types <- data.frame(table(id.matrix$Gr))
  gr <- c("1_1_1_1", "1_1_0_0", "1_1_0_1", "0_0_1_1", "0_1_1_1", "0_0_0_1", "0_1_0_0", "0_1_0_1")
  t1 <- Gr.types$Var1
  t1 <- mgsub(gr, c("haplotype_1.paired_haplotype_2.paired", "haplotype_1.paired_haplotype_2.NA....", 
                    "haplotype_1.paired_haplotype_2.improp", "haplotype_1.NA...._haplotype_2.paired", 
                    "haplotype_1.improp_haplotype_2.paired", "haplotype_1.NA...._haplotype_2.improp", 
                    "haplotype_1.improp_haplotype_2.NA....", "haplotype_1.improp_haplotype_2.improp") ,t1)
  t2 <- Gr.types$Var1
  t2 <- mgsub(gr, c("Both or by groups", 
                    "haplotype_1_specific", "haplotype_1_specific", 
                    "haplotype_2_specific", "haplotype_2_specific", 
                    "Dump2", "Dump3", "Dump1"), t2)
  t3 <- Gr.types$Var1
  t3 <- mgsub(gr, c("A0,A0X,A1hap1,A1hap2,A2,A2X,B1hap1,B1hap2,BZ,C1hap1,C1hap2,CZ", 
                    "haplotype_1", "haplotype_1", "haplotype_2", "haplotype_2", "D2", "D3", "D1"), t3)


  # stat 3
  Gr.types <- data.frame(Group=Gr.types$Var1, Freq=Gr.types$Freq, Desc=t1, Anno=t2, tagXG=t3)
  #message(Gr.types)

  
  # 1.  Group 1_1_1_1  : haplotype_1.paired, haplotype_2.paired
  # 1.1 tag NH1
  both_paired.ids <- subset(id.matrix, id.matrix$Gr=="1_1_1_1")
  both_paired.hap1 <- merge(dt.hap1.paired, both_paired.ids)
  both_paired.hap2 <- merge(dt.hap2.paired, both_paired.ids)

  hap1_nh1.ids <- subset(both_paired.hap1, both_paired.hap1$hap1_nh_1==1 & both_paired.hap1$hap1_nh_2==1, c(1,4,8)) 
  hap2_nh1.ids <- subset(both_paired.hap2, both_paired.hap2$hap2_nh_1==1 & both_paired.hap2$hap2_nh_2==1, c(1,4,8)) 
  
  all_nh1.ids <- merge(hap1_nh1.ids, hap2_nh1.ids, all=TRUE)
  both_nh1.ids <- merge(hap1_nh1.ids, hap2_nh1.ids)
  both_nh1.hap1 <- merge(both_nh1.ids, dt.hap1.paired)
  both_nh1.hap2 <- merge(both_nh1.ids, dt.hap2.paired)
  

  fr_nh1 <- merge(both_nh1.hap1, both_nh1.hap2, all=TRUE, by="FragmentID")
  #fr_nh1 <- fr_nh1[order(fr_nh1$hap1_pos_1),]
  
  # A0 and A0X
  sCigar.sNM <- subset(fr_nh1, 
                          as.character(fr_nh1$hap1_cigar_1)==as.character(fr_nh1$hap2_cigar_1) &
                          as.character(fr_nh1$hap1_cigar_2)==as.character(fr_nh1$hap2_cigar_2) &
                          fr_nh1$hap1_nm_1==fr_nh1$hap2_nm_1 & fr_nh1$hap1_nm_2==fr_nh1$hap2_nm_2)
                         
  # A0
  sCigar.sNM_0 <- subset(sCigar.sNM, 
                          sCigar.sNM$hap1_nm_1==0 & sCigar.sNM$hap1_nm_2==0 & sCigar.sNM$hap2_nm_1==0 & sCigar.sNM$hap2_nm_2==0)
  
  # A0X
  sCigar.sNM_0X <- subset(sCigar.sNM, 
                          !(sCigar.sNM$hap1_nm_1==0 & sCigar.sNM$hap1_nm_2==0 & sCigar.sNM$hap2_nm_1==0 & sCigar.sNM$hap2_nm_2==0))

  # A2 (potentially structrual variation like indels between two hoplotypes)
  dCigar.sNM <- subset(fr_nh1, 
                          (as.character(fr_nh1$hap1_cigar_1)!=as.character(fr_nh1$hap2_cigar_1) | as.character(fr_nh1$hap1_cigar_2)!=as.character(fr_nh1$hap2_cigar_2)) &
                          fr_nh1$hap1_nm_1==fr_nh1$hap2_nm_1 & fr_nh1$hap1_nm_2==fr_nh1$hap2_nm_2)
  
  # A1*
  sCigar.dNM <- subset(fr_nh1, 
                          as.character(fr_nh1$hap1_cigar_1)==as.character(fr_nh1$hap2_cigar_1) & 
                          as.character(fr_nh1$hap1_cigar_2)==as.character(fr_nh1$hap2_cigar_2) &
                          ( fr_nh1$hap1_nm_1!=fr_nh1$hap2_nm_1 | fr_nh1$hap1_nm_2!=fr_nh1$hap2_nm_2) )
  
  
  # A1hap2
  nh1.sd.hap2 <- subset(sCigar.dNM, 
                          (sCigar.dNM$hap1_nm_1 + sCigar.dNM$hap1_nm_2) > (sCigar.dNM$hap2_nm_1 + sCigar.dNM$hap2_nm_2) )
  message(dim(nh1.sd.hap2))
  # A1hap1
  nh1.sd.hap1 <- subset(sCigar.dNM, 
                          (sCigar.dNM$hap1_nm_1 + sCigar.dNM$hap1_nm_2) < (sCigar.dNM$hap2_nm_1 + sCigar.dNM$hap2_nm_2) )
  message(dim(nh1.sd.hap1))

  # B1*
  dCigar.dNM <- subset(fr_nh1, 
                       (as.character(fr_nh1$hap1_cigar_1)!=as.character(fr_nh1$hap2_cigar_1) | as.character(fr_nh1$hap1_cigar_2)!=as.character(fr_nh1$hap2_cigar_2))  & 
                         ( fr_nh1$hap1_nm_1!=fr_nh1$hap2_nm_1 | fr_nh1$hap1_nm_2!=fr_nh1$hap2_nm_2) )
  
  # B1hap2
  nh1.dd.hap2 <- subset(dCigar.dNM, (dCigar.dNM$hap1_nm_1 + dCigar.dNM$hap1_nm_2) > (dCigar.dNM$hap2_nm_1 + dCigar.dNM$hap2_nm_2) )
  # B1hap1
  nh1.dd.hap1 <- subset(dCigar.dNM, (dCigar.dNM$hap1_nm_1 + dCigar.dNM$hap1_nm_2) < (dCigar.dNM$hap2_nm_1 + dCigar.dNM$hap2_nm_2) )
  # BZ
  amb <- subset(dCigar.dNM, (dCigar.dNM$hap1_nm_1 + dCigar.dNM$hap1_nm_2) == (dCigar.dNM$hap2_nm_1 + dCigar.dNM$hap2_nm_2) )
  
  # C1hap1
  hap1.nh1_hap2.na <- subset(all_nh1.ids, is.na(all_nh1.ids$hap2_nh_1) )
  # C1hap2
  hap1.na_hap2.nh1 <- subset(all_nh1.ids, is.na(all_nh1.ids$hap1_nh_1) )

  if ( dim(merge(hap1.nh1_hap2.na, both_paired.ids))[1] != dim(hap1.nh1_hap2.na)[1] )
    {message("Fragments mapped on haplotype_1 are not in a list of both paried... may be out of 1_1_1_1 group...")}
  if ( dim(merge(hap1.na_hap2.nh1, both_paired.ids))[1] != dim(hap1.na_hap2.nh1)[1] )
    {message("Fragments mapped on haplotype_2 are not in a list of both paired... may be out of 1_1_1_1 group...") }


  message("get reads ids with assigned group ids...")
  # 1: NH1, 2: paired, NH >1, 
	ids.NH1 <- data.frame(FragmentID=sCigar.sNM_0$FragmentID, haplotype_1.NH="1", haplotype_2.NH="1", CIGAR="S", NM="0", XG="A0")		# both, no mismatch
	if( nrow(sCigar.sNM_0X) >0 ){
		ids.NH1 <- rbind( ids.NH1,
				data.frame(FragmentID=sCigar.sNM_0X$FragmentID, haplotype_1.NH="1", haplotype_2.NH="1", CIGAR="S", NM="S", XG="A0X") )	# both, with mismatch
	}
	if( nrow(nh1.sd.hap1) >0){
		ids.NH1 <- rbind( ids.NH1,
				data.frame(FragmentID=nh1.sd.hap1$FragmentID, haplotype_1.NH="1", haplotype_2.NH="1", CIGAR="S", NM="D:hap1<hap2", XG="A1hap1") )	# both or haplotype_1, hap1 had small NM
	}
	if( nrow(nh1.sd.hap2) >0){
		ids.NH1 <- rbind( ids.NH1,
				data.frame(FragmentID=nh1.sd.hap2$FragmentID, haplotype_1.NH="1", haplotype_2.NH="1", CIGAR="S", NM="D:hap1>hap2", XG="A1hap2") )	# both or haplotype_2, hap2 has small NM
	}
	if( nrow(dCigar.sNM) >0){
                ids.NH1 <- rbind( ids.NH1,
				data.frame(FragmentID=dCigar.sNM$FragmentID, haplotype_1.NH="1", haplotype_2.NH="1", CIGAR="D", NM="S", XG="A2") )	# both (structural variation)
        }
	if( nrow(nh1.dd.hap1) >0){
                ids.NH1 <- rbind( ids.NH1,
				data.frame(FragmentID=nh1.dd.hap1$FragmentID, haplotype_1.NH="1", haplotype_2.NH="1", CIGAR="D", NM="D:hap1<hap2", XG="B1hap1") )	# different cigar, hap1 has small NM
        }
	if( nrow(nh1.dd.hap2) >0){
                ids.NH1 <- rbind( ids.NH1,
				data.frame(FragmentID=nh1.dd.hap2$FragmentID, haplotype_1.NH="1", haplotype_2.NH="1", CIGAR="D", NM="D:hap1>hap2", XG="B1hap2") )	# different cigar, hap2 has small NM
        }
	if( nrow(amb) >0){
                ids.NH1 <- rbind( ids.NH1,
				data.frame(FragmentID=amb$FragmentID, haplotype_1.NH="1", haplotype_2.NH="1", CIGAR="D", NM="D", XG="BZ") )		# both or either or none
        }
	if( nrow(hap1.nh1_hap2.na) >0){
                ids.NH1 <- rbind( ids.NH1,
				data.frame(FragmentID=hap1.nh1_hap2.na$FragmentID, haplotype_1.NH="1", haplotype_2.NH=">=2", CIGAR=".", NM=".", XG="C1hap1") )	# haplotype_1, haplotype_1 has NH1
        }
	if( nrow(hap1.na_hap2.nh1) >0){
                ids.NH1 <- rbind( ids.NH1,
				data.frame(FragmentID=hap1.na_hap2.nh1$FragmentID, haplotype_1.NH=">=2", haplotype_2.NH="1", CIGAR=".", NM=".", XG="C1hap2") )	# haplotype_2, haplotype_2 has NH1
        }

  temp1 <- both_paired.hap1[!(both_paired.hap1$FragmentID %in% ids.NH1$FragmentID),]
  temp2 <- both_paired.hap2[!(both_paired.hap2$FragmentID %in% ids.NH1$FragmentID),]
  
  ids.NH.over2 <- rbind( data.frame(FragmentID=unique(temp1$FragmentID), haplotype_1.NH=">=2", haplotype_2.NH=">=2", CIGAR=".", NM=".", XG="CZ" ) )	# both have NH>=2

  spe_hap1_fid=subset(id.matrix, id.matrix$Gr=="1_1_0_0" | id.matrix$Gr=="1_1_0_1", c(1) )
  spe_hap2_fid=subset(id.matrix, id.matrix$Gr=="0_0_1_1" | id.matrix$Gr=="0_1_1_1", c(1) )

  spe_hap1_fid <- dt.hap1.paired[dt.hap1.paired$FragmentID %in% spe_hap1_fid$FragmentID, ]
  ids.haplotype_1_spe <- data.frame(
			FragmentID=spe_hap1_fid$FragmentID,
			haplotype_1.NH=as.character(spe_hap1_fid$hap1_nh_1),
			haplotype_2.NH="NA", 
			CIGAR=".", NM=".", XG="haplotype_1" )
  spe_hap2_fid <- dt.hap2.paired[dt.hap2.paired$FragmentID %in% spe_hap2_fid$FragmentID, ]
  ids.haplotype_2_spe <- data.frame(
                        FragmentID=spe_hap2_fid$FragmentID,
                        haplotype_1.NH="NA",
			haplotype_2.NH=as.character(spe_hap2_fid$hap2_nh_1),
                        CIGAR=".", NM=".", XG="haplotype_2" )
  

  tagXG <- rbind(ids.NH1, ids.NH.over2, ids.haplotype_1_spe, ids.haplotype_2_spe)

  if( nrow(subset(id.matrix, id.matrix$Gr=="0_1_0_1")) >0 ){
	tagXG <- rbind( tagXG, data.frame(FragmentID=subset(id.matrix, id.matrix$Gr=="0_1_0_1", c(1) ), haplotype_1.NH="0", haplotype_2.NH="0", CIGAR=".", NM=".", XG="D1" ) )
  }
  if( nrow(subset(id.matrix, id.matrix$Gr=="0_0_0_1")) >0 ){
	tagXG <- rbind( tagXG, data.frame(FragmentID=subset(id.matrix, id.matrix$Gr=="0_0_0_1", c(1) ), haplotype_1.NH="NA", haplotype_2.NH="0", CIGAR=".", NM=".", XG="D2" ) )
  }
  if( nrow(subset(id.matrix, id.matrix$Gr=="0_1_0_0")) >0 ){
	tagXG <- rbind( tagXG, data.frame(FragmentID=subset(id.matrix, id.matrix$Gr=="0_1_0_0", c(1) ), haplotype_1.NH="0", haplotype_2.NH="NA", CIGAR=".", NM=".", XG="D3" )  )
  }

  temp <- data.frame(table(tagXG[,2:6]))
  tagXG.table <- subset(temp, temp$Freq >0)
  tagXG.desc <- tagXG.table$XG
  tagXG.desc <- mgsub(c("A0X", "A0", "A1hap1", "A1hap2", "A2", "B1hap1", "B1hap2", "BZ", 
                        "C1hap1", "C1hap2", "CZ", "haplotype_2", "haplotype_1", "D1", "D2", "D3"),
                      c("sCigar.sNM>0", "sCigar.sNM0", "sCigar.dNM_haplotype_1", "sCigar.dNM_haplotype_2", "dCigar.sNM", 
                        "dCigar.dNM_haplotype_1", "dCigar.dNM_haplotype_2", "Both_ambiguous", 
                        "NH1inP.NH>1inC", "NH>1inP.NH1inC", "NH>1inBoth", 
                        "haplotype_2", "haplotype_1", "Dump1", "Dump2", "Dump3"), tagXG.desc)
  tagXG.gr <- tagXG.table$XG
  tagXG.gr <- mgsub(c("A0X", "A0", "A1hap1", "A1hap2", "A2", "B1hap1", "B1hap2", "BZ", 
                      "C1hap1", "C1hap2", "CZ", "haplotype_2", "haplotype_1", "D1", "D2", "D3"),
                    c("1_1_1_1", "1_1_1_1","1_1_1_1", "1_1_1_1", "1_1_1_1", "1_1_1_1", "1_1_1_1", "1_1_1_1", 
                      "1_1_1_1", "1_1_1_1", "1_1_1_1", 
                      "0_0_1_1,0_1_1_1", "1_1_0_0,1_1_0_1", "0_1_0_1", "0_0_0_1", "0_1_0_0" ) ,tagXG.gr)

  
  getRates.p <- function(dt, totalNumber) {
    result <- c()
    for (i in 1:length(dt)) {
      if(dt[i] ==0 | dt[i] == "NA"){ result <- c(result, 0)  }
      else {  result <- c(result, round( ((tagXG.table$Freq[i]/totalNumber)*100),2 ))  }
    }
    result
  }
  getRates.a <- function(dt, totalNumber) {
    result <- c()
    for (i in 1:length(dt)) {
      if(dt[i] == "NA" ){ result <- c(result, 0)  }
      else {  result <- c(result, round( ((tagXG.table$Freq[i]/totalNumber)*100),2 ))  }
    }
    result
  }

  tagXG.rates.hap1.paired <- getRates.p(tagXG.table$haplotype_1.NH, uniqueNumber$haplotype_1[2])
  tagXG.rates.hap2.paired <- getRates.p(tagXG.table$haplotype_2.NH, uniqueNumber$haplotype_2[2])
  tagXG.rates.hap1.all <- getRates.a(tagXG.table$haplotype_1.NH, uniqueNumber$haplotype_1[1])
  tagXG.rates.hap2.all <- getRates.a(tagXG.table$haplotype_2.NH, uniqueNumber$haplotype_2[1])
  tagXG.table <- cbind(tagXG.gr, tagXG.desc, tagXG.table, 
                       tagXG.rates.hap1.paired, tagXG.rates.hap2.paired, 
                       tagXG.rates.hap1.all, tagXG.rates.hap2.all)
 
  stats <- list(uniqueNumber=uniqueNumber, groups=Gr.types, tagXG=tagXG.table)

  output <- list(stats=stats, tagXG=tagXG)
  message("Return output...")
  return(output)
}






