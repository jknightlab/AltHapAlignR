#' To get gene feature from a GTF file
#' 
#' @param gtf GTF file
#' @param type selecting gene feature
#' @export
#' @import rtracklayer
#' @examples
#' output <- getGeneList(gtf, type="protein_coding")
#' 
#' # Output data
#' # summary of gene information, format: data.frame

gtf="inst/extdata/gencode.v21.chr_patch_hapl_HLA.annotation.gtf"

getGeneList <- function(gtf, type="protein_coding"){
        g <- import(gtf, "gtf")
        g <- subset(g, type=="gene")
        g <- data.frame(chr=seqnames(g),
                        start=g@ranges@start,
                        end=g@ranges@start + g@ranges@width,
                        gene_name=g$gene_name,
                        gene_type=g$gene_type)
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
#' EDfromBams(bamFiles="./inst/extdata/mapping2*bam", gtf="./inst/extdata/gencode.v25.chr_patch_hapl_HLA.annotation.gtf", output_name="output.txt", virtualenv=NA, r="VARSL-VARS2,C6orf205-MUC21")


          
EDfromBams <- function(bamFiles, gtf, output_name="ed_table.txt", virtualenv=NA, r="VARSL-VARS2,C6orf205-MUC21"){

    py_path=paste0(system.file(package="AltHapAlignR"), "./inst/scripts/make_a_table.py")
    bamFiles <- paste(bamFiles, collapse = " ")
    
    command <- paste("python", py_path,  gtf, bamFiles, ">", output_name )

    
    if(!is.na(r)){
        r <- gsub("-", " ", r)
        r <- gsub(",", " -r ", r)
        r <- paste0("-r ", r)
        command <- paste("python ", py_path, r,  gtf, bamFiles, ">", output_name )
    }
    
    
    if(!is.na(virtualenv) ){
          command <- paste0("./inst/scripts/wrapper.sh ", virtualenv, " ", command )
          print(command)
    }
    if(is.na(virtualenv) ){
      print(command)
    }
    
    output= system(command, intern=T)
    return(output)
}


# ./make_a_table.py -r VARSL VARS2 -r C6orf205 MUC21 gencode.v21.only_MHC.annotation.gtf */*picard*bam > output_file
# example: 
gtf="./inst/extdata/gencode.v25.chr_patch_hapl_HLA.annotation.gtf"
#bamFiles="mapping2apd.bam,mapping2apd.bam,mapping2cox.bam,mapping2mann.bam,mapping2mcf.bam,mapping2qbl.bam,mapping2ssto.bam,mapping2pgf.bam"
bamFiles="./inst/extdata/mapping2*"
output_name ="test.txt"
#export PYTHONPATH=/well/jknight/Irina/Programs/HTSeq-0.6.1/build/lib.linux-x86_64-2.7:/well/jknight/software/rescomp/lib/python2.7/site-packages/${PYTHONPATH:+:$PYTHONPATH}
#export PYTHONPATH=/Users/wanlee/Documents/well/test_althap:/Library/Python/2.7/site-packages/${PYTHONPATH:+:$PYTHONPATH}


#EDfromBams(bamFiles, gtf, output_name="with_virtual.txt", virtualenv="/Users/wanlee/Documents/well/test_althap/althapalign_virtualenv", r="VARSL-VARS2,C6orf205-MUC21")


#EDfromBams(bamFiles, gtf, output_name="without_virtual.txt", r="VARSL-VARS2,C6orf205-MUC21")


#EDfromBams(bamFiles, gtf, output_name="test_with_virtual.txt", virtualenv=1, r="VARSL-VARS2,C6orf205-MUC21")

EDfromBams(bamFiles, gtf, output_name="without_virtual_21.txt", r="VARSL-VARS2,C6orf205-MUC21")
EDfromBams(bamFiles, gtf, output_name="without_virtual_25_r2_1.txt", r="VARSL-VARS2,C6orf205-MUC21")


gtf="./inst/extdata/gencode.v21.chr_patch_hapl_HLA.annotation.gtf"
bamFiles="./inst/extdata/mapping2*"
EDfromBams(bamFiles, gtf, output_name="without_virtual_v21_r1_1.txt", r="VARSL-VARS2,C6orf205-MUC21")
