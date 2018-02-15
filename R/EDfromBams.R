#' To get gene feature from a GTF file
#' 
#' @param gtf GTF file
#' @param type selecting gene feature
#' @export
#' @import rtracklayer
#' @importFrom rtracklayer blocks
#' @importFrom rtracklayer path
#' @examples
#' gtf=system.file("extdata", "gencode.v21.chr_patch_hapl_HLA.annotation.gtf", 
#'                  package = "AltHapAlignR")
#' output <- getGeneList(gtf, type="protein_coding")
#' 
#' # Output data
#' # summary of gene information, format: data.frame


getGeneList <- function(gtf, type="protein_coding"){
        g <- import(gtf, "gtf")
        g <- subset(g, type=="gene")
        g <- data.frame(chr=as.character(g@seqnames),
                        start=as.numeric(g@ranges@start),
                        end=as.numeric(g@ranges@start + g@ranges@width),
                        gene_name=as.character(g$gene_name),
                        gene_type=as.character(g$gene_type),
                        gene_status=as.character(g$gene_status),
                        transcript_type=as.character(g$transcript_type) )
        g <-subset(g, gene_type==type) 
}


#' To get editing distance from multiple bam files
#' 
#' @param gtf gtf file <path/to/GTF_file> 
#' @param bamFiles bam files <path/to/BAM_file_1> <path/to/BAM_file_2>
#' @param output_name output name
#' @param r combine redundnat gene names, example, VARSL-VARS2,C6orf205-MUC21, optional
#' @param virtualenv path of virtualenv, optional
#' @export
#' @examples
#' bamFiles=system.file("extdata", "mapping2*bam", package = "AltHapAlignR")
#' gtf=system.file("extdata", "gencode.v21.chr_patch_hapl_HLA.annotation.gtf", 
#'                  package = "AltHapAlignR")
#' EDfromBams(bamFiles, gtf, output_name="output.txt", virtualenv=0, r="VARSL-VARS2,C6orf205-MUC21")

          
EDfromBams <- function(bamFiles, gtf, output_name="ed_table.txt", virtualenv=0, r="VARSL-VARS2,C6orf205-MUC21"){
    
    py_path=system.file("scripts", "make_a_table.py", package = "AltHapAlignR")

    #py_path="./inst/scripts/make_a_table.py"
    bamFiles <- paste(bamFiles, collapse = " ")
    
    command <- paste("python", py_path,  gtf, bamFiles, ">", output_name )

    
    if(!is.na(r)){
        r <- gsub("-", " ", r)
        r <- gsub(",", " -r ", r)
        r <- paste0("-r ", r)
        command <- paste("python ", py_path, r,  gtf, bamFiles, ">", output_name )
    }
    
    
    if( virtualenv!=0 ){
          command <- paste0("./inst/scripts/wrapper.sh ", virtualenv, " ", command )
          #print(command)
    }
    if( virtualenv==0 ){
          #print(command)
    }
    
    output= system(command, intern=T)
    return(output)
}




