#'@title gene position counts
#'@description Inputs needed are: recoded genetic data formatted in PLINK format, SNP name with BP (position) and gene name with START and END position. The first six columns of the input genetic data follow standard PLINK .raw format. Column names as FID, IID, PAT, MAT, SEX and PHENOTYPE followed by SNP information as recoded by the PLINK software. The function returns allelic counts per gene per sample (where each row represents a gene and each column represents an individual starting with the second column where first column contains gene information). 

gene_pos_counts<-function(dt_gen,dt_snp,dt_gene){ 

##07 20 2020

#' @export
#' @import data.table
#' @importFrom dplyr %>%
#' @importFrom data.table rowid
#' @importFrom data.table .SD
#' @importFrom data.table :=
#' @param dt_gen recoded genetic data from PLINK 
#' @param dt_gene with CHR START END GENE as column names. Where CHR should be integer 1-22. START and END column should be integer. GENE column contains gene names
#' @param dt_snp with SNP BP as column names
#'
#' @examples
#' #Package provides sample data that are loaded with package loading. 
#' #not RUN
#' data(recodedgen) #PLINK raw formatted data of 10 individiduals with 10 SNPs
#'
#' data(genecoord) #gene coordinates with START, END, CHR and GENE names. 
#' #Five genes with start and end genomic coordinates
#'
#' data(snppos) #SNP and BP column names with SNP names and SNP genomic location in BP. 
#' #10 SNPs with genomic location
#'
#' gene_pos_counts(recodedgen, snppos, genecoord) #run the function
#' #end not RUN
#'
#' @return Returns an object of data.table class as an output with allelic gene counts within each sample where each row corresponds to gene and column to individual IDs from column second. The first column contains gene names.
#'
#' @author Sanjeev Sariya
#'

    dt_gen<-data.table::as.data.table(dt_gen) # convert into data.table
    dt_gene<-data.table::as.data.table(dt_gene)
    dt_snp<-data.table::as.data.table(dt_snp)
    START<-END<-GENE<-BP<-NULL ## binding the variable locally to the function

    if(all(garcom_check_column_names(dt_gene, c("START","END","GENE")))){
        ## all good with gene data
    }else{
        stop("column names don't match for gene data")
    }

    if(all(garcom_check_column_names(dt_snp, c("SNP","BP")))){
        ## all good with SNP data
    }else{
        stop("column names don't match for snp data")
    }
    ##Check ends 

    colnames(dt_gen) <- gsub("_.*","",colnames(dt_gen)) ##Remove _ from recode format

    ## https://gist.github.com/nacnudus/ef3b22b79164bbf9c0ebafbf558f22a0

    snp_withingenes<-dt_snp[dt_gene, c("SNP","GENE","START","END"), on=list(BP>=START , BP<=END), nomatch=0] # inner join ##https://stackoverflow.com/questions/63290994/foverlaps-data-table-error-ys-key-must-be-identical-to-the-columns-specified

    if(nrow(snp_withingenes) == 0){
        stop("No snps within any gene boundaries provided")	
    }
    ##if gene sum is zero. Stop and return


    if(isFALSE( (any( colnames(dt_gen) %in% unique(snp_withingenes$SNP)) )) ){
        stop("No SNPs overlapping between genetic data and SNP annotation with Gene boundaries")
    }
    ##if nothing matches then Stop and error out


    ##dt_gen_subset<- dt_gen[,.SD,.SDcols=intersect(colnames(dt_gen),unique(snp_withingenes$SNP ))] %>%  .[, rowid := dt_gen$IID ] %>% data.table::transpose(keep.names = "SNP", make.names="rowid")

    dt_gen_subset<- dt_gen[,.SD,.SDcols=intersect(colnames(dt_gen),unique(snp_withingenes$SNP ))] ## get columns that intersect 
    dt_gen_subset[, rowid := dt_gen$IID ] ## assign row names as IIDs
    dt_gen_subset<-data.table::transpose(dt_gen_subset,keep.names = "SNP", make.names="rowid") ## tranpose data and have some fancy settings

    ## we perform inner join. SNPs that are found in input boundaries-annotation as well as .raw data
#### subsetsnps_genes_lefted_join <- snp_withingenes[dt_gen_subset, on="SNP", nomatch=0] %>% .[, c("START","END","SNP"):=NULL] %>% data.table::setcolorder(.,c("GENE")) ## remove START, END, BP and SNP column, and in the put GENE column and then order 

    subsetsnps_genes_lefted_join <- snp_withingenes[dt_gen_subset, on="SNP", nomatch=0] 
    subsetsnps_genes_lefted_join[, c("START","END","SNP"):=NULL]  
    
    ##https://stackoverflow.com/a/32277135/2740831
    ##  matrix_withallelecount_withinGene <-subsetsnps_genes_lefted_join[,lapply(.SD,sum,na.rm=TRUE),by=GENE] %>% .[ rowSums(.[,-c("GENE")]) > 0,]

    matrix_withallelecount_withinGene <-subsetsnps_genes_lefted_join[,lapply(.SD,sum,na.rm=TRUE),by=GENE] 
    matrix_withallelecount_withinGene<-matrix_withallelecount_withinGene[ rowSums(matrix_withallelecount_withinGene[,-c("GENE")]) > 0,]
    
    if(nrow(matrix_withallelecount_withinGene)>0){
        return(matrix_withallelecount_withinGene)
    }
    else{
        return(NULL) ## if nothing is left after gene >0 
    }

} ## function ends
