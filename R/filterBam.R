

#' Flitering a bam file by read names (fragment ids)
#'
#' This filters a bam by a table with read names and tags(NH, HI, etc.). 
#' An input table can be an output from a function, 'bamByFeatures'. 
#' @param bfl A bam file
#' @param FragmentID A table of read names with additional colnums for tag information depending on user filtering option. The first column should be fragment ids. currently only NH and HI tags are availbe. 
#' @param outputName An output file name
#' @param readFiltering "unique": it will produce only uniquely mapped reads with 'NH==1' tag from the input 'FragmentID'. "all": all reads from the input 'FragmentID'.  "FALSE": currently not available.
#' @param tagFiltering When an 'readFiltering' is 'FLASE', it will use a tagFilering. Currently, only "NH" and "HI" tags are usable. To use this option, the 'FragmentID' table should have columns with the same names of tags.
#' @param mapper a name of mapping tool. currently, tophat and bwa are available
#' @export
#' @import GenomicAlignments
#' @import rtracklayer
#' @examples
#' # Input data: a bam file and fragment ids
#' bfl <- system.file("extdata", "MHC_PGF.bam", package = "HapTag")
#' # FragmentID: Input table, It's an output from a function 'tagReads'. 
#' # The first column should be read names. 
#' FragmentID <- system.file("extdata", "assignedTags.txt", package = "HapTag")
#' # To filter a bam with uniquely mapped reads.
#' bamFilterByFragmentIDs(bfl, FragmentID, outputName="output_bam", 
#' 		readFiltering="unique", tagFiltering=c("NH", "HI"), mapper="tophat" )


bamFilterByFragmentIDs<- function(bfl, FragmentID, outputName="output_bam", readFiltering="all", tagFiltering=c("NH", "HI"), mapper=c("tophat", "bwa")){
	# output: bam or grange
	# first column sould be fragment ids
	
	mapper <- match.arg(mapper)
	FragmentID <- data.frame(FragmentID)
	colnames(FragmentID)[1] <- "ID"
 
	message('Reading the bam file...')
	if(mapper=="tophat"){
		tags <- c("NM", "NH", "HI") 
 	}
	if(mapper=="bwa"){
		tags <- c("NM", "X0", "X1")
	}

	aln <- readGAlignments(bfl, param=ScanBamParam(what=scanBamWhat(), tag=tags))
 	names(aln) <- mcols(aln)$qname 

	# to avoid wrong filtering...
	ftrd_bam <- aln[paste("", mcols(aln)$qname, "", sep="_") %in% paste("", FragmentID$ID, "", sep="_")]

	if(mapper=="bwa"){
		ftrd_bam@elementMetadata$NH <- ftrd_bam@elementMetadata$X0
		ftrd_bam@elementMetadata$HI <- ftrd_bam@elementMetadata$X1
	}

	if(readFiltering == "all"){
		message('No filtering, all hits with the fragment ids in the input file...')
	}

	if(readFiltering == "unique"){
		message('Filtering multiple hits and get only uniquely mapped reads...')
                tmp <- as.data.frame(table(ftrd_bam@elementMetadata$NH==1))
		ftrd_bam <- ftrd_bam[!is.na(ftrd_bam@elementMetadata$NH)]
		ftrd_bam <- ftrd_bam[ftrd_bam@elementMetadata$NH==1]
	}

	#if(uniqueRead == "FALSE"){
        #        message('Filtering reads by a tagFiltering option ...')
	#	tmp <- paste(mcols(ftrd_bam)$qname, ftrd_bam@elementMetadata$NH, ftrd_bam@elementMetadata$HI, sep="_")
	#	ft <- paste(FragmentID$ID, FragmentID$NH, FragmentID$HI, sep="_")
	#	ftrd_bam <- ftrd_bam[paste("", tmp, "", sep="_") %in% paste("", ft, "", sep="_")]
        #}
	message("Reads after filtering: ", length(ftrd_bam))

	if(mapper=="bwa"){
                ftrd_bam@elementMetadata$X0 <- ftrd_bam@elementMetadata$NH
                ftrd_bam@elementMetadata$X1 <- ftrd_bam@elementMetadata$HI
		ftrd_bam@elementMetadata$NM <- as.integer(ftrd_bam@elementMetadata$NM)
		ftrd_bam@elementMetadata$X0 <- as.integer(ftrd_bam@elementMetadata$X0)
		ftrd_bam@elementMetadata$X1 <- as.integer(ftrd_bam@elementMetadata$X1)
        }

	if(length(ftrd_bam) >0){
			message('Writing the filtered bam file...')
			export(ftrd_bam, BamFile(paste(outputName, ".bam", sep="")))
	}
	if(length(ftrd_bam)==0){
		message('No read in the bam file...')
	}
}




#' filtering a bam file by feature information
#'
#' It deals with only singe chromosome so a bam file should be filtered by one chromosome (with a certain location). 
#' Using a filtered gtf with single chromosome or a certain region of the chromosome will make it run faster. 
#' @param bam a bam file
#' @param readLen a size of a single read
#' @param gtf gtf file
#' @param chr chromosome
#' @param feature feature type, currently filtering with either 'gene' or 'transcript' feature is available. 
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import sqldf

# readsToGenes : previous name of a function

bamByFeatures <- function(bam, readLen=51, gtf, chr="chr6", feature="gene"){

	# gtf file required
	# only for one chromosome, so a bam file shoud be first filtered by one chromosome (with specific locations) 
	# geneList(locus) is IRange format filtered by specific genomic feature (here, gene feature), use gene name with location due to the same gene names
        getOverlappingFeature  <- function(read_names, readLen, pos, ops, cigar, geneList){
                ranges_on_ref <- cigarRangesAlongReferenceSpace(cigar=cigar, pos=pos, ops=ops)
                names(ranges_on_ref) <- read_names
                unlisted_ranges_on_ref <- unlist(ranges_on_ref)

                # get read length with M cigar before findOverlaps with genes

                tmp1 <- unlisted_ranges_on_ref[as.numeric(unlisted_ranges_on_ref@width) < readLen, ]
                tmp1 <- data.frame(tmp1)[,3:4]
                ops_M_before=aggregate(width ~ names, data=tmp1, FUN=sum)

                # find overlaps reads with gene regions

                hits <- findOverlaps(unlisted_ranges_on_ref, geneList)
                read2gene <- cbind(hits@queryHits, as.data.frame(unlisted_ranges_on_ref[hits@queryHits]),
                                hits@subjectHits, as.data.frame(geneList[hits@subjectHits]) )

                # get read length with M cigar after findOverlaps with genes

                tmp2 <- subset(read2gene, as.numeric(read2gene$width) < readLen)
                tmp2 <- cbind(names2=paste(tmp2[,5], tmp2[,10], sep="_tmp2_"), tmp2)
                ops_M_after=aggregate(width ~ names2, data=tmp2, FUN=sum)
                ops_M_after$names <- gsub("(.*)_tmp2_.*", "\\1", ops_M_after$names2, perl=TRUE)

                # compare ops_M_before with ops_M_after

                tmp <- " select * from ops_M_after LEFT JOIN ops_M_before ON (ops_M_after.names=ops_M_before.names and ops_M_after.width=ops_M_before.width)"
                tmp <- sqldf(tmp)
                tmp <- as.character(tmp[!is.na(tmp[,4]),c(1)])

                a <- read2gene[paste(read2gene[,5],read2gene[,10], sep="_tmp2_") %in% tmp,]
                read2gene_ft <- rbind(subset(read2gene, as.numeric(read2gene$width)==readLen), a)

                return(read2gene_ft)
        }



        gtf_feature <- subset(gtf, gtf[,3]==as.character(feature) & gtf[,1]==chr)
        message(nrow(gtf_feature), " ", feature, "s found")

	if( feature=="transcript"){
		gtf_names <- paste(gsub(".*transcript_name (.*); level .*", "\\1", gtf_feature[,9], perl=TRUE), gtf_feature[,4], gtf_feature[,5], sep="_")
		length(gtf_names)
	}

	if( feature=="gene"){
		gtf_names <- paste(gsub(".*gene_name (.*); transcript_type .*", "\\1", gtf_feature[,9], perl=TRUE), gtf_feature[,4], gtf_feature[,5], sep="_")
		length(gtf_names)
	}
	
	geneList <- IRanges(gtf_feature[,4], width=(gtf_feature[,5]-gtf_feature[,4]+1 ), names=gtf_names)


	ops <- c("M")
	what=scanBamWhat()
	tags <- c("NM", "NH", "HI")
	flags <- scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE)
	param=ScanBamParam(what=what, flag=flags, tag=tags)
	paired <- readGAlignmentPairs(bam, use.names=FALSE, param=param)      # use.names=TRUE takes longer
	dumped <- getDumpedAlignments()

	read_names <-  paste(paired@first@elementMetadata$qname, paste(paired@first@elementMetadata$NH, paired@first@elementMetadata$HI, sep="."), sep="_tmp_" )

	first_seq <- paired@first@elementMetadata$seq
	first_pos <-  paired@first@elementMetadata$pos
	first_cigar <- paired@first@elementMetadata$cigar
	first_strand <- paired@first@elementMetadata$strand

	last_seq <- paired@last@elementMetadata$seq
	last_pos <-  paired@last@elementMetadata$pos
	last_cigar <- paired@last@elementMetadata$cigar
	last_strand <- paired@last@elementMetadata$strand

	##### paired readds
	# 1) link reads to genes

	first <- getOverlappingFeature(read_names=read_names, readLen=readLen, pos=first_pos, ops=ops, cigar=first_cigar, geneList=geneList)
	colnames(first) <- c("queryHits", "queryHits_start", "queryHits_end", "queryHits_width", "queryHits_names", "subjectHits", "subjectHits_start", "subjectHits_end", "subjectHits_width", "subjectHits_names")

	last <- getOverlappingFeature(read_names=read_names, readLen=readLen, pos=last_pos, ops=ops, cigar=last_cigar, geneList=geneList)
	colnames(last) <- c("queryHits", "queryHits_start", "queryHits_end", "queryHits_width", "queryHits_names", "subjectHits", "subjectHits_start", "subjectHits_end", "subjectHits_width", "subjectHits_names")

	# 2) get not splitted reads and get paried reads mapped to the same gene

	first_noSplittedRead <- subset(first, first[,4]==readLen)
	last_noSplittedRead <- subset(last, last[,4]==readLen)

	sql_paired_reads_mapped_to_the_same_gene_1 <- " select * from first_noSplittedRead LEFT JOIN last_noSplittedRead ON (first_noSplittedRead.queryHits_names=last_noSplittedRead.queryHits_names)"
	paired_reads_mapped_to_genes <- sqldf(sql_paired_reads_mapped_to_the_same_gene_1)
	paired_reads_mapped_to_the_same_gene_1 <- unique( subset(paired_reads_mapped_to_genes, paired_reads_mapped_to_genes[,10]==paired_reads_mapped_to_genes[,20], c(5,10)) ) 

	groups_read2genes <- data.frame( read_NH_HI=paired_reads_mapped_to_the_same_gene_1$queryHits_names, readSplitted=0, mapped2TheSameGene=1, linkedGene=as.character(paired_reads_mapped_to_the_same_gene_1$subjectHits_names) )	

	# 3) get splitted reads and get paried reads mapped to the same gene

	first_splittedRead <- subset(first, first[,4] < readLen)
	last_splittedRead <- subset(last, last[,4] < readLen)

	sql_paired_reads_mapped_to_the_same_gene_2_1 <- " select * from first_splittedRead LEFT JOIN last ON (first_splittedRead.queryHits_names=last.queryHits_names and first_splittedRead.subjectHits_names=last.subjectHits_names)"  
	paired_reads_mapped_to_genes_2_1 <- sqldf(sql_paired_reads_mapped_to_the_same_gene_2_1)
	paired_reads_mapped_to_genes_2_1 <-  paired_reads_mapped_to_genes_2_1[!is.na(paired_reads_mapped_to_genes_2_1[,11]), ]

	sql_paired_reads_mapped_to_the_same_gene_2_2 <- " select * from first LEFT JOIN last_splittedRead ON (first.queryHits_names=last_splittedRead.queryHits_names and first.subjectHits_names=last_splittedRead.subjectHits_names)"
        paired_reads_mapped_to_genes_2_2 <- sqldf(sql_paired_reads_mapped_to_the_same_gene_2_2)
	paired_reads_mapped_to_genes_2_2 <- paired_reads_mapped_to_genes_2_2[!is.na(paired_reads_mapped_to_genes_2_2[,11]), ]

	paired_reads_mapped_to_the_same_gene_2 <- unique(rbind(paired_reads_mapped_to_genes_2_1[,c(5,10)], paired_reads_mapped_to_genes_2_2[,c(5,10)]))
	groups_read2genes <- rbind( groups_read2genes,
			data.frame( read_NH_HI=paired_reads_mapped_to_the_same_gene_2$queryHits_names, readSplitted=1, mapped2TheSameGene=1, linkedGene=paired_reads_mapped_to_the_same_gene_2$subjectHits_names) )


	# 4) others: one of both has a gene overlapped, both have no gene overlapped

	others <- unique( read_names[ !(read_names %in% groups_read2genes$read_NH_HI)] )
	groups_read2genes <- rbind( groups_read2genes, data.frame( read_NH_HI=others, readSplitted="0or1", mapped2TheSameGene=0, linkedGene=as.character("-") ) )


	##### dumped reads

	if(length(dumped) >0) {	
		d_read_names <- paste( dumped@elementMetadata$qname, paste(dumped@elementMetadata$NH, dumped@elementMetadata$HI, sep="."), sep="_tmp_")
		groups_read2genes <- rbind( groups_read2genes,
                        data.frame( read_NH_HI=d_read_names, readSplitted="d_0or1", mapped2TheSameGene=0, linkedGene=as.character("-")) )
	}

	groups_read2genes <- data.frame(read_name=gsub("(.*)\\_tmp.*", "\\1", groups_read2genes$read_NH_HI, perl=TRUE),
					read_NH=gsub(".*\\_tmp\\_(\\d+)\\..*", "\\1", groups_read2genes$read_NH_HI, perl=TRUE),
					read_HI=gsub(".*\\_tmp\\_\\d+\\.(.*)", "\\1", groups_read2genes$read_NH_HI, perl=TRUE),
					groups_read2genes[,2:4])
	return(groups_read2genes)
}
