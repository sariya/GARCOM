#' @title gene annotation counts using VCF data
#' @description Function returns a matrix with allelic (reference) counts per gene per individual for SNP-gene annotation
#' @details Inputs needed are a vcf data and a data frame of SNP-gene annotation. The function returns a matrix of allelic counts (reference) per gene per sample (where each row represents a gene and each column represents an individual starting with the second column where first column contains gene information).
#' 

vcf_counts_annot<-function(vcf_data,df_snpgene,keep_indiv=NULL,extract_SNP=NULL,filter_gene=NULL){
    ## added on 08 28 2020
    #' @export
    #' @import vcfR
    #' @import data.table
    #'
    #' @param vcf_data an object of vcfR class
    #'
    #' @param df_snpgene a data frame that contains SNP and annotated gene with SNP and GENE as column name
    #'
    #' @param keep_indiv an option to specify individuals to retain. Mutation counts will be provided for individuals included in the list only. Default is all individuals. Provide list of individuals in a vector.
    #'
    #' @param extract_SNP an option to specify SNPs for which mutation counts are needed. Mutation counts will be provided for SNPs included in the list only. Default is all SNPs.
    #'
    #' @param filter_gene an option to filter in a list of Genes. Mutation counts will be provided for genes specifed in the list only. Default is all genes. Provide list of genes in a vector.
    #'@examples 
    #'\dontrun{
    #' vcf_counts_annot(vcf,df_snpgene_test)
    #' }
    #'
    #' @return Returns an matrix of data.table class as an output with allelic (reference) gene counts within each sample where each row corresponds to gene and column to individual IDs from column second. The first column contains gene names.
    #'
    #' @author Sanjeev Sariya
    #'
    SNP<-GENE<-NULL ## binding the variable locally to the function    
    if( class(vcf_data)[1] !="vcfR" ){
        print("VCF annot: vcfR class not found")
    }
    
    df_snpgene<-data.table::as.data.table(df_snpgene)
    if(all(garcom_check_column_names(df_snpgene,c("SNP","GENE")))){
        ## all good with SNP data
    }else{
        stop("VCF annot: column names don't match for snp-gene data")
    }
    ## check ends
    if(FALSE == isTRUE(garcom_check_unique(df_snpgene) )){
        stop("VCF annot: Duplicate SNP-Gene annotation values")
    }
    
    if(is.null(keep_indiv)==TRUE){
        
        genotyped_extracted<-vcfR::extract.gt(vcf_data,element = "GT",as.numeric=TRUE,convertNA=TRUE) 
    } else{
        genotyped_extracted <- tryCatch({
            vcfR::extract.gt(vcf_data,element = "GT",as.numeric=TRUE,convertNA=TRUE)[,keep_indiv]
            
        }, warning = function(w) {
            message(paste("warning vcf counts annot in subsetting individuals", w))
            
        }, error =function(e) {
            message(paste(" vcf counts annot: error subsetting individuals ", e))
            stop("Exiting vcf counts annot ")
            
        }) 
    } ## else ends for subsetting individuals 
    
    df_genotyped_extracted<-data.table::data.table(genotyped_extracted,keep.rownames="SNP") ## assign row names as SNP to perform a join
    
    if(is.null(extract_SNP) == FALSE){
        
        extract_SNP<-as.character(extract_SNP)
        
        df_genotyped_extracted<-tryCatch({
            df_genotyped_extracted[SNP %in% extract_SNP,]
            
        }, warning = function(w) {
            message(paste("warning vcf counts annot in subsetting SNPs", w))
            
        }, error =function(e) {
            message(paste(" vcf counts annot: error subseting SNPs ", e))
            stop("Exiting vcf counts annot ")
            
        } ) 
    } ## else ends for subseting SNPs 
    ## subsetting based on SNPs ends

    if(is.null(filter_gene) == FALSE){
        ##Start process to filter genes. The list provided by user is what we'd like to keep
        filter_gene<-as.character(filter_gene) ## turn into character
        df_snpgene<-garcom_filter_gene(df_snpgene,filter_gene) ##filter SNP-gene annotation based on Gene list
    }
    ##sub-setting ends for gene filter
    
    jointed_gene_VCFGT<-df_snpgene[df_genotyped_extracted,on=c(SNP="SNP"),nomatch=0L]
    jointed_gene_VCFGT<-jointed_gene_VCFGT[,SNP:=NULL] ### remove SNP cols
    
    if(nrow(jointed_gene_VCFGT)==0){
        message("VCF annot: No SNPs match with the annotation")
        return(NULL)
    }

    jointed_gene_VCFGT<-jointed_gene_VCFGT[, lapply(.SD,as.numeric),by="GENE"] ## convert into numeric
    jointed_genesSNP_filtered<-jointed_gene_VCFGT[,lapply(.SD,sum,na.rm=TRUE),by=GENE] # get sum within a gene
    jointed_genesSNP_filtered<-jointed_genesSNP_filtered[rowSums(jointed_genesSNP_filtered[,-c("GENE")]) > 0,]
    
    if(nrow(jointed_genesSNP_filtered) == 0){
        stop("VCF annot: No genes with any reference genotyped in the VCF with supplied SNP gene annotation")
    }
    return(jointed_genesSNP_filtered)
    
}
## function ends
