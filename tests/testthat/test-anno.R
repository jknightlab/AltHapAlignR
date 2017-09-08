context("Importing exons")

gtfFile <- system.file("extdata", "gencode.v21.hla.annotation.gtf", package="HapTag")
gr <- rtracklayer::import.gff(gtfFile, asRangedData = FALSE)

test_that("annotated GRanges object is parsed properly", {
			exons <- getExons(gr)
			
			expect_is(exons, "data.frame")
			expect_is(getExons(subset(gr, seqnames(gr) == "foo")), "data.frame")
			expect_equal(nrow(exons), length(subset(gr, type == "exon")))
			expect_true(all(c("start", "end", "width", "strand", "gene_id", 
									"transcript_id", "gene_type", "gene_status", 
							"gene_name", "transcript_type", "transcript_status",  
							"transcript_name", "exon_number", "exon_id", "level") %in% names(exons)))
			
			gr2 <- gr
			gr2$group <- ""
			expect_error(getExons(granges(gr)), "missing")
			expect_error(getExons(gr2), "missing annotations")
			expect_warning(getExons(subset(gr, type != "exon")), "no exon")
		}
)

test_that("GTF data can be imported",{
			expect_identical(gtf2ExonTable(gtfFile), gtf2ExonTable(gr))
			expect_identical(gtf2ExonTable(file(gtfFile)), gtf2ExonTable(gr))
		}
)

test_that("subsetting of annotation works",{
			gr6 <- subset(gr, seqnames(gr) == "chr6")
			genome(gr6) <- "hg38"
			allExons <- gtf2ExonTable(gr6, grange=NULL)
			expect_identical(gtf2ExonTable(gr6, grange=c("", "0", "0")), allExons)
			expect_identical(max(gtf2ExonTable(gr6, grange=c("chr6", "28000000", "0"))$end), 
					max(allExons$end))
			expect_identical(min(gtf2ExonTable(gr6, grange=c("chr6", "0", "34000000"))$start), 
					min(allExons$start))
			expect_equal(nrow(gtf2ExonTable(gr, 
									grange=c("chr6", allExons$start[1], allExons$end[1]))), 1)
			
			expect_warning(gtf2ExonTable(subset(gr, seqnames(gr) == "foo")), "no annotation")
			expect_warning(gtf2ExonTable(gr6, refver="hg19"), "Reference version")
		}
)