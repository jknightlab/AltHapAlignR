


#' To get mapping rates in all pairs of haplotypes and produce summerized figure 
#' 
#' @param ed_table a file with a table of editing distances generated from a function 'EDframBams'
#' @param hap_names haplotype names matching to column names of the input file 'ed_table' (eg. c("apd", "cox", "dbb", "mann", "mcf", "pgf", "qbl", "ssto") )
#' @param read_length read length in bam files
#' @export
#' @import reshape2
#' @import dplyr
#' @importFrom dplyr first slice filter
#' @importFrom utils combn read.table
#' @examples
#' ed_table=system.file("extdata", "example_ed_table.txt", package = "AltHapAlignR")
#' ed_table <- read.table(ed_table, sep="\t", header=TRUE)
#' output <- getMappingRatesFromPairs(ed_table, 
#'                hap_names=c("apd", "cox", "dbb", "mann", "mcf", "pgf", "qbl", "ssto"), 
#'                read_length=50)


getMappingRatesFromPairs <- function(ed_table, hap_names=c("apd", "cox", "dbb", "mann", "mcf", "pgf", "qbl", "ssto"),  read_length=50){

        message("Loading a file...")
        f <- ed_table
        colnames(f)[3:ncol(f)] <- hap_names
        f[is.na(f)] <- 100


        hap_combn <-  t(combn(hap_names,2))
        hap_combn <- rbind(hap_combn, cbind(hap_names, hap_names))
        
        genes <- data.frame(genes=unique(as.character(f$gene_name) ) )
          
        message("Calculating paired mapping rates...")
        paired_mr <- list()
        paired_mr <- apply(genes, 1, function(g) {
                message(g)
                reads_in_a_gene <- subset(f, gene_name==g)
                total_number_of_reads_in_a_gene <- nrow(reads_in_a_gene)
                total_number_of_reads_in_a_gene
        
                a <- melt(reads_in_a_gene, id.vars = c("read_name", "gene_name") )
                a$value[is.na(a$value)] <- 100
                a <- group_by(a, as.character(read_name) )
                a <- as.data.frame(a)
                
                a <- subset(a, value==min(a$value))
                a <- a[,1:4]
        
                reads_in_a_gene <- dcast(a, read_name ~ variable, value.var = "value")
                reads_in_a_gene[is.na(reads_in_a_gene)] <- read_length
        

                #num_reads <- list()
                num_reads <- apply(hap_combn, 1, function(x){
                        #message(as.character(x[1]))
                        #message(as.character(x[2]))
                        
                        only_h1=0
                        only_h2=0
                        intersect=0
                        rate_of_gene_count=0
                    
                        col_h1 <- grep ( as.character(x[1]), colnames(reads_in_a_gene) )
                        col_h2 <- grep ( as.character(x[2]), colnames(reads_in_a_gene) )
                  
                        if(length(col_h1)==1 & length(col_h2)==1 ){
                              if( as.character(x[1])==as.character(x[2])) {
                                paired_haps <- reads_in_a_gene[, c(1, col_h1) ]
                                paired_haps <- subset(paired_haps, paired_haps[, 2]!=read_length) 
                                total_number_of_reads_in_a_pair <- nrow(paired_haps)
                                
                                intersect_mismatch = subset(paired_haps, paired_haps[,2]!=read_length)
                                intersect_mismatch = sum(intersect_mismatch[,2])
                                
                                tmp <- data.frame(gene_name=g, 
                                                  Hap1=as.character(x[1]), Hap2=as.character(x[2]),
                                                  gene_count=total_number_of_reads_in_a_pair,
                                                  rate_of_gene_count=total_number_of_reads_in_a_pair/total_number_of_reads_in_a_gene,
                                                  only_h1=0, 
                                                  only_h2=0, 
                                                  intersect=nrow(paired_haps) ,
                                                  h1_mismatch=0, 
                                                  h2_mismatch=0, 
                                                  intersect_mismatch=intersect_mismatch,
                                                  h1_mismatching_rate= 0,
                                                  h2_mismatching_rate= 0,
                                                  intersect_mismatching_rate= intersect_mismatch/(nrow(paired_haps)*read_length),
                                                  mismatching_rate= (intersect_mismatch)/(total_number_of_reads_in_a_pair*read_length)  )
                              }

                        
                            if( as.character(x[1])!=as.character(x[2]) ) {
                                    paired_haps <- reads_in_a_gene[, c(1, col_h1, col_h2) ]

                                    paired_haps <- subset(paired_haps, paired_haps[, 2]!=read_length | paired_haps[, 3]!=read_length ) 
                  
                                    total_number_of_reads_in_a_pair <- nrow(paired_haps)
                                    #total_number_of_reads_in_a_pair

                                    only_h1 <- subset(paired_haps, paired_haps[, 2] < paired_haps[, 3] )

                                    only_h2 <- subset(paired_haps, paired_haps[, 2] > paired_haps[, 3]  )

                                    both_haps <- subset(paired_haps, paired_haps[,2] == paired_haps[,3] )

                  
                                  if(!is.null(only_h1)){
                                        h1_mismatch = subset(only_h1, only_h1[,2]!=read_length)
                                        h1_mismatch = sum(h1_mismatch[,2])
                                  }

                                  if(!is.null(only_h2)){
                                        h2_mismatch = subset(only_h2, only_h2[,3]!=read_length)
                                        h2_mismatch = sum(h2_mismatch[,3])
                                  }

                                  if(!is.null(both_haps)){
                                        intersect_mismatch = subset(both_haps, both_haps[,2]!=both_haps[,3])
                                        intersect_mismatch = sum(intersect_mismatch[,2])
                                  }

                                  tmp <- data.frame(gene_name=g, 
                                                Hap1=as.character(x[1]), Hap2=as.character(x[2]),
                                                gene_count=total_number_of_reads_in_a_pair,
                                                rate_of_gene_count=total_number_of_reads_in_a_pair/total_number_of_reads_in_a_gene,
                                                only_h1=nrow(only_h1), 
                                                only_h2=nrow(only_h2), 
                                                intersect=nrow(both_haps) ,
                                                h1_mismatch=h1_mismatch, 
                                                h2_mismatch=h2_mismatch, 
                                                intersect_mismatch=intersect_mismatch,
                                                h1_mismatching_rate= h1_mismatch/(nrow(only_h1)*read_length),
                                                h2_mismatching_rate= h2_mismatch/(nrow(only_h2)*read_length),
                                                intersect_mismatching_rate= intersect_mismatch/(nrow(both_haps)*read_length),
                                                mismatching_rate= (h1_mismatch+h2_mismatch+intersect_mismatch)/(total_number_of_reads_in_a_pair*read_length)  )
                                #tmp
                 
                            }
                              tmp
                        }
                } )
        
               if( length(num_reads) >0 ) {
                          num_reads <- do.call(rbind.data.frame, num_reads)
                          num_reads <- unique(num_reads)
                          num_reads <- num_reads[ order(-num_reads$rate_of_gene_count, num_reads$mismatching_rate), ]
                
                          r1 <- num_reads$rate_of_gene_count[1]
                          r2 <- num_reads$mismatching_rate[1]
                          putative_haps_ <- c(1)
                          putative_haps=1
                
                          if(nrow(num_reads)>1){
                                    for(i in 2:nrow(num_reads)){
                                            if( num_reads$rate_of_gene_count[i]==r1 & num_reads$mismatching_rate[i]==r2){
                                                    putative_haps = putative_haps
                                            }
                                             else{
                                                    putative_haps = putative_haps +1
                                            }
                                            r1=num_reads$rate_of_gene_count[i]
                                            r2=num_reads$mismatching_rate[i]
                                            putative_haps_ <- c(putative_haps_, putative_haps)
                                        }
                            }
                            group_a <- subset(num_reads, rate_of_gene_count > (max(num_reads$rate_of_gene_count)-0.03) ) 
                            group_b <- subset(num_reads, rate_of_gene_count <= (max(num_reads$rate_of_gene_count)-0.03) ) 
                
                            if(nrow(group_b)>0){
                                      num_reads <- rbind(data.frame(group_a, type="A"),
                                                         data.frame(group_b, type="B") )
                            }
                            if(nrow(group_b)==0){
                                      num_reads <- rbind(data.frame(group_a, type="A") )
                            }
                            num_reads$putative_haps <- putative_haps_

                    }
        
                paired_mr[[g]] <- num_reads
        } )

        paired_mr <- do.call(rbind.data.frame, paired_mr)
        return(paired_mr)
}





#' To predict shortest path
#' 
#' @param paired_mapping_rates output from 'predictRecombHapByGenes'  
#' @param gtf a gtf file name containing gene feature of all haplotypes
#' @param hap_names haplotype names matching to column names (eg. c("apd", "cox", "dbb", "mann", "mcf", "pgf", "qbl", "ssto") )
#' @param penalty 1: penalty applied, 0: penalty not applied
#' @param sample_name output file name
#' @import igraph
#' @importFrom igraph union simplify path as_data_frame groups
#' @importFrom dplyr first slice
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @importFrom gridExtra combine
#' @importFrom plyr desc
#' @importFrom data.table data.table
#' @importFrom stats aggregate filter
#' @export
#' @examples
#' paired_mapping_rates=system.file("extdata", "example_paired_mapping_rates.txt", 
#'                  package = "AltHapAlignR")
#' paired_mapping_rates <- read.table(paired_mapping_rates, sep="\t", header=TRUE)
#' gtf=system.file("extdata", "gencode.v21.chr_patch_hapl_HLA.annotation.gtf", 
#'                  package = "AltHapAlignR")
#' output <- heatmapByShortestPaths(paired_mapping_rates, gtf, 
#'                  hap_names= c("apd", "cox", "dbb", "mann", "mcf", "pgf", "qbl", "ssto"), 
#'                  penalty=1, sample_name="test)
#' # Input data: max_comb from predictRecombHapByGenes
#' # Output data
#' # *.shortestPath_heatmap.pdf : heatmap of selected haplotypes
#' # gene and haplotype to draw heatmap
#' heatmap_genes <- output$heatmap_genes
#' # predicted best pairs of haploytpes
#' summary_of_best_pairs <- output$summary_of_best_pairs


heatmapByShortestPaths <- function(paired_mapping_rates, gtf, hap_names= c("apd", "cox", "dbb", "mann", "mcf", "pgf", "qbl", "ssto"), penalty=1, sample_name="sample_name"){
  
        ########################################################################################
        # Load annotation from a gtf file 
        ########################################################################################
      genes_from_gtf <- getGeneList(gtf, type="protein_coding")
      
      if(nrow(genes_from_gtf)==0){
            message("no gtf file loaded.")
            
            genes_from_gtf <- read.table( system.file("extdata", "example_gencode.v21.pc_genes.txt", package = "AltHapAlignR"), sep="\t", header=TRUE)
            message("example gtf file loaded.")
      }
      
      
      unique_best_hits <- data.frame(table(paired_mapping_rates[,c("gene_name", "type")]) )

      haps <- data.frame(hap_chrs= c("GL000250", "GL000251", "GL000252", "GL000253", "GL000254", "chr6", "GL000255", "GL000256"),
                         hap_names=c("apd", "cox", "dbb", "mann", "mcf", "pgf", "qbl", "ssto") )
      
      #paired_mapping_rates <- subset(paired_mapping_rates, gene_count >= 50)   
  
  
      gene_list <- unique(as.character(paired_mapping_rates$gene_name)) 
      
      
      ########################################################################################
      # get genes in the input data
      ########################################################################################
      
      listed_genes <- genes_from_gtf[which(as.character(genes_from_gtf$gene_name) %in% gene_list),]
      listed_genes$chr <- gsub("\\.\\d+", "", listed_genes$chr)
      listed_genes <- listed_genes[which(listed_genes$chr %in% haps$hap_chrs), ]
      listed_genes <- merge(listed_genes, haps, by.x="chr", by.y="hap_chrs")
      colnames(listed_genes)[8] <- "hap_n"
      # hap_n <- c()
      # for(a in as.character(listed_genes$chr) ){
      #   a1 <- as.character(haps$hap_names)[ grep(a, haps$hap_chrs) ]
      #   hap_n <- c(hap_n, a1)
      # }
      # listed_genes$hap_n <- hap_n
      
      
  
      comb_haps <- unique(paste0( paired_mapping_rates[,c("Hap1")], "_", paired_mapping_rates[,c("Hap2")]) )
  
      v1 <- comb_haps
      v2 <- comb_haps
      tmp <- sort(apply(expand.grid(v1, v2), 1, paste, collapse = ",", sep = "")) 
      hap_combn <- data.frame(haps_from=gsub("(.*),(.*)", "\\1", tmp),
                          haps_to=gsub("(.*),(.*)", "\\2", tmp) )
  
  
  
      ### summary of best_pairs
  
  
      # additional filtering
      a_ <- subset(paired_mapping_rates, (only_h1>5 & only_h2>5 & as.character(Hap1)!=as.character(Hap2) ) | (only_h1==0 & only_h2==0 & as.character(Hap1)==as.character(Hap2)) )       #wo.dp5
      #a__ <- subset(paired_mapping_rates, (only_h1>5 & Hap1!=Hap2) | (only_h2>5 & Hap1!=Hap2) |(only_h1==0 & only_h2==0 & Hap1==Hap2 & intersect >5) )       #wo.dp5
  
  
      a <- rbind(data.frame(gene_name=a_$gene_name, gene_hap=paste0(a_$gene_name, "_", a_$Hap1)),
                 data.frame(gene_name=a_$gene_name, gene_hap=paste0(a_$gene_name, "_", a_$Hap2)))
      a <- unique(a)
      a2 <- data.frame(table(a$gene_name))
  
  
      n_hap1 <- a_[as.character(a_$gene_name) %in% as.character(a2[a2$Freq==1,]$Var1),]     # homozygous
      n_hap2 <- a_[as.character(a_$gene_name) %in% as.character(a2[a2$Freq==2,]$Var1),]     # heterozygous
      n_hap_any <- a_[as.character(a_$gene_name) %in% as.character(a2[a2$Freq>2,]$Var1),]   # more than two pairs
      n_hap_any_hm <- subset(n_hap_any, as.character(Hap1)==as.character(Hap2) )
      #n_hap_any_ht <- n_hap_any[!(n_hap_any$gene_name %in% n_hap_any_hm$gene_name), ]
      n_hap_any_ht <- subset(n_hap_any, as.character(Hap1)!=as.character(Hap2) )
  
  
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
  
  
  
        paired_mapping_rates$mapping_rate <- 1- paired_mapping_rates$mismatching_rate
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
 
      tmp <- c()
      for(i in as.character(haps$hap_names)){

              a1 <- subset(listed_genes, hap_n==i)
              redu <- data.frame( table(a1[,c("gene_name")] ) )
              redu <- subset(redu, Freq>1)
              
              if(nrow(redu) >0){
                      b1 <- a1[ !(a1$gene_name %in% redu$Var1 ), ]
                      b2 <- c()
                      for(i2 in unique(as.character(redu$Var1)) ) {
                                #print(i2)
                                a2 <- subset(redu, Var1==i2)
                                a3 <- a1[a1$gene_name %in% a2$Var1, ]
                                a3 <- a3[1,]
                                b2 <- rbind(b2, a3)
                      }
                      a1 <- rbind(b1, b2)
              }

              a1 <- a1[order(as.numeric(a1$start) ), ]
              a1$ord=c(1:nrow(a1))
              tmp <- rbind(tmp, a1)
      }
      listed_genes <- tmp
      
      listed_genes$gene_name <- gsub("VARSL", "VARS2", listed_genes$gene_name)
      listed_genes$gene_name <- gsub("C6orf205", "MUC21", listed_genes$gene_name)
      ng_in_haps <- aggregate(ord ~ hap_n, listed_genes, max)
      
      tmp <- subset(listed_genes, hap_n=="pgf")
      tmp2 <- listed_genes[!(listed_genes$gene_name %in% tmp$gene_name), ]
      
      gene_order <- c()
      if(nrow(tmp2)>0){
              missing_genes <- unique(as.character(tmp2$gene_name))
              for(i2 in missing_genes){
                      a1 <- subset(tmp2, gene_name==i2)
                      #print(a1)
                      a1 <- merge(a1, ng_in_haps, by="hap_n")
                      a1 <- subset(a1, ord.y==max(a1$ord.y))
                      g_bf <- subset(listed_genes, ord==(a1$ord.x[1] - 1) & hap_n==a1$hap_n[1], c("gene_name"))
                      #g_af <- subset(listed_genes, ord==(a1$ord.x[1] + 1) & hap_n==a1$hap_n[1], ("gene_name") )
                      
                      tmp_ord=subset(listed_genes, hap_n=="pgf" & gene_name==as.character(g_bf), c('ord') )
                      tmp_ord <- as.numeric(tmp_ord) +0.5
                      #subset(listed_genes, hap_n=="pgf" & gene_name==as.character(g_af) )
                      gene_order <- rbind(gene_order,
                                          data.frame(gene_name=as.character(a1$gene_name[1]),
                                                     ord=tmp_ord) )
                      
              }
      }
      
      gene_order <-  rbind(gene_order,
                           unique(tmp[,c("gene_name", "ord")]) )
      
      gene_order <- gene_order[order(gene_order$ord), ]
      gene_order$ord2 <- c(1:nrow(gene_order))
      gene_order <- gene_order[, c("gene_name", "ord2")]
      colnames(gene_order)[2] <- "ord"
      
 
  
  
  #################################################################
  ###
  ### get ordered genes with pvlaues from predictRecombHapByGenes
  ###
  #################################################################
  
  
  input_structure <- data.frame(input_structure)
  input_structure$gene_name <- rownames(input_structure)
  mr <- merge(input_structure, gene_order, by="gene_name")
  mr <- mr[order(as.numeric(mr$ord)),]
  
  
  t1 <- data.frame(c(rep(0, length(comb_haps)), 0) )
  rownames(t1) <- c(comb_haps, "ord")
  t1 <- data.frame(gene_name="Start", t(t1) )
  
  t2 <- data.frame( c(rep(0, length(comb_haps)), rep((max(mr$ord)+1), 1)) )
  rownames(t2) <- c(comb_haps, "ord")
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
  colnames(mr)[ncol(mr)] <- "ord"
  mr_col <- colnames(mr)[2:(ncol(mr)-1)]
  rownames(mr) <- NULL
  
  
  
  
  message("Making path structure...")
  path_frame <- c()
  for(ord in 1:(nrow(mr)-1) ){
    #message(paste0(ord, ": ", mr[ord, c("gene_name")]) )
    
    g1 <- as.character(mr[ord, c("gene_name")])
    g2 <- as.character(mr[(ord+1), c("gene_name")])
    o1 <- as.numeric(as.character(mr[ord, c("ord")]))
    o2 <- as.numeric(as.character(mr[(ord+1), c("ord")]))
    
    gfrom <- data.frame(haps=mr_col, mr_from= as.numeric(t( mr[ord,2:(ncol(mr)-1)] )[,1] ) ) 
    gto  <-  data.frame(haps=mr_col,  mr_to=as.numeric(t( mr[(ord+1),2:(ncol(mr)-1)] )[,1] ) ) 
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
  message("Constructing the shortest path...")
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
  

  
  z2 <- listed_genes[, c("gene_name", "hap_n")]
  z2$gene_hap <- paste0(z2$gene_name, "_", z2$hap_n)
  
  

  #dat <- merge(dat, z2, by="gene_hap", all=TRUE)
  dat <- merge(dat, z2, by="gene_hap", all=TRUE)
  dat$best_pair[is.na(dat$best_pair)] <- 0
  dat$gene_name.y <- gsub("(.*)_(.*)", "\\1", dat$gene_hap)
  dat$hap_n <- gsub("(.*)_(.*)", "\\2", dat$gene_hap)
  
  dat <- dat[,c("gene_hap", "gene_name.y", "hap_n", "best_pair")]
  colnames(dat) <- c("gene_hap", "gene_name", "hap_from", "best_pair")
  dat$hap_from <- toupper(dat$hap_from)
  dat$gene_hap <- toupper(dat$gene_hap)
  
  
  dat <- merge(dat, gene_order, by="gene_name", all=TRUE)
  colnames(dat) <- c("gene_name", "gene_hap", "hap", "best_pair", "ord")
  
  dat <- dat[!is.na(dat$ord), ]
  dat <- dat[order(dat$ord), ]
  
  mapping_rates <- paired_mapping_rates
  tmp1 <- aggregate(rate_of_gene_count ~ gene_name, data = mapping_rates, max)
  tmp2 <- aggregate(gene_count ~ gene_name, data = mapping_rates, max)
  tmp <- merge(tmp1, tmp2, by="gene_name")
  tmp$total_reads <- tmp$gene_count / tmp$rate_of_gene_count
  tmp <- tmp[,c("gene_name", "total_reads")]
  
  mapping_rates <- merge(mapping_rates, tmp, by="gene_name")
  mapping_rates$Mapping_Rates1 <- ( mapping_rates$only_h1 + mapping_rates$intersect ) / mapping_rates$total_reads
  mapping_rates$Mapping_Rates2 <- ( mapping_rates$only_h2 + mapping_rates$intersect ) / mapping_rates$total_reads
  mapping_rates2 <- rbind(data.frame(gene_hap=paste0(mapping_rates$gene_name, "_", mapping_rates$Hap1),
                                     Mapping_Rates=mapping_rates$Mapping_Rates1),
                          data.frame(gene_hap=paste0(mapping_rates$gene_name, "_", mapping_rates$Hap2),
                                     Mapping_Rates=mapping_rates$Mapping_Rates2) )
  
  mapping_rates2 <- unique(mapping_rates2)
  mapping_rates2$gene_hap <- toupper(mapping_rates2$gene_hap)
  
  dat <- merge(dat, mapping_rates2, by="gene_hap")
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
  
  
  # #hm
  # dat_1 <- dat[dat$gene_name %in% as.character(zyg_count$Var1[zyg_count$Freq==1]), ]
  # dat_1 <- dat_1[dat_1$gene_hap %in% zyg$gene_hap,]
  # 
  # #ht
  # dat_2 <- dat[dat$gene_name %in% as.character(zyg_count$Var1[zyg_count$Freq==2]), ]
  # dat_2 <- dat_2[dat_2$gene_hap %in% zyg$gene_hap,]
  # 
  # 
  # # multiple hm
  # dat_3 <- dat_1[dat_1$gene_name %in% n_hap_any_hm$gene_name, ]
  # 
  # # multiple ht
  # dat_4 <- dat_2[dat_2$gene_name %in% n_hap_any_ht$gene_name, ]
  # 
  # 
  # dat_1 <- dat_1[!(dat_1$gene_name %in% dat_3$gene_name), ]
  # dat_2 <- dat_2[!(dat_2$gene_name %in% dat_4$gene_name), ]
  
  
  
  
  paired_mapping_rates$g_h1_h2 <- paste0(paired_mapping_rates$gene_name, "_", paired_mapping_rates$Hap1, "_", paired_mapping_rates$Hap2) 
  
  selected_genes_by_top1path <- paired_mapping_rates[(paired_mapping_rates$g_h1_h2 %in% top1_edges$V1 ), ]
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
    #message(i)
    a <- subset(dat, as.character(gene_name)==i)
    
    tmp <- subset(selected_genes_readCount, as.character(gene_name)==i)
    tmp <- subset(tmp, value!=0 )
    if(nrow(tmp)==1){
      tmp$Zygosity <- "homozygous"
      tmp$Haplotype <- "Hap1andHap2"
      #a$best_pair[which(a$best_pair==1)][2] <- 0
      if(length( a$best_pair[which(a$best_pair==1)] ) >1){
        #message("step1")
        a$best_pair[which(a$best_pair==1)] <- c (1,0)
      }
    }
    
    if(nrow(tmp)==2){
      
      if( tmp$Haplotype[1]=="Hap1only" & tmp$Haplotype[2]=="Hap2only" ){
      }
      else if( tmp$Haplotype[1]=="Hap1only" & tmp$Haplotype[2]=="Hap1andHap2" ){
      }
      else if( tmp$Haplotype[1]=="Hap2only" & tmp$Haplotype[2]=="Hap1andHap2" ){
      }
      else if( any(grepl("Hap1only", tmp$Haplotype)) ){
        tmp <- data.frame(tmp[1,1:2], 
                          Zygosity="homozygous", Haplotype="Hap1andHap2",
                          value=sum(tmp$value))
        #a$best_pair[which(a$best_pair==1)][2] <- 0
        #message("step2")
        a$best_pair[which(a$best_pair==1)] <- c (1,0)
      }
      else if( any(grepl("Hap2only", tmp$Haplotype) )  ){
        tmp <- data.frame(tmp[1,1:2], 
                          Zygosity="homozygous", Haplotype="Hap1andHap2",
                          value=sum(tmp$value))
        #a$best_pair[which(a$best_pair==1)][1] <- 0
        #message("step3")
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
  
  
  summary_set <- final_set[, c("gene_name", "ord", "mapping_rate", "rate_of_gene_count", "Zygosity")]
  summary_set$mapping_rate <- 1- summary_set$mapping_rate
  summary_set$rate_of_gene_count <- summary_set$rate_of_gene_count *100
  colnames(summary_set) <- c("gene_name", "ord", "Mismatching_rates", "Rates_of_gene_counts", "Zygosity")
  summary_set <- melt(summary_set, id.vars = c("gene_name", "ord", "Zygosity"))
  summary_set$variable <- factor(summary_set$variable, levels=c("Rates_of_gene_counts", "Mismatching_rates"))
  
  
  summary_heatmap <- ggplot(summary_set, aes(x = variable, y = ord , fill=Zygosity)) + 
    geom_tile() + scale_fill_brewer(palette="Pastel1") +
    geom_tile(data=summary_set, aes(x=variable, y=ord, color=Zygosity),size=0,fill=NA, color="grey") +
    geom_text(aes(label=round( value,3) ), size=5, color="black" ) +
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
  #test <- data.frame(table(paste0(selected_genes_by_top1path$Hap1, "_", selected_genes_by_top1path$Hap2)) )

  if(length(unique(dat_final$gene_name) )  > 15){
            height_length= length(unique(dat_final$gene_name) ) * 0.28
  }
  if(length(unique(dat_final$gene_name) )  < 15){
            height_length= length(unique(dat_final$gene_name) ) * 0.32
  }
  
  ggsave(filename=paste0(sample_name, ".shortestPath_heatmap.pdf" ), g_p123, width=18, height=height_length)
  
  
  output <- list(heatmap_genes=dat,
                 #summary_of_best_pairs=selected_genes_by_top1path,
                 summary_of_best_pairs=final_set)
  
}

