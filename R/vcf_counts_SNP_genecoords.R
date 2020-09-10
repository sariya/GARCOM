## added on 09 08 2020

#'@title VCF gene position counts
#'@description Function returns a matrix with allelic counts per gene per individual for SNP and gene coordinates as inputs
#'

vcf_counts_SNP_genecoords<-function(vcf_data,df_snppos,df_genecoords,keep_indiv=NULL,extract_SNP=NULL,filter_gene=NULL){
    
    #' @export
    #' @import data.table
    #' @import vcfR
    
    #' @param vcf_data an object of vcfR class
    #' @param keep_indiv an option to specify individuals to retain. Mutation counts will be provided for individuals included in the list only. Default is all individuals. Provide list of individuals in a vector.
    
    #' @param df_snppos a dataframe for SNP information with SNP BP as column names.
    #' @param df_genecoords a dataframe for gene boundaries with CHR START END GENE as column names. Where CHR should be integer 1-22. START and END column should be integer. GENE column contains gene names
    #'
    #' @param extract_SNP an option to specify SNPs for which mutation counts are needed. Mutation counts will be provided for SNPs included in the list only. Default is all SNPs.
    #'@examples 
    #'\dontrun{
    #' vcf_counts_SNP_genecoords(vcf_data_test,df_snppos_test,df_genecoords_test)
    #' }
    #'
#' @param filter_gene an option to filter in Genes. Mutation counts will be provided for genes included in the list only. Default is all genes.

    #' @return Returns an matrix of data.table class as an output with allelic (reference) gene counts within each sample where each row corresponds to gene and column to individual IDs from column second. The first column contains gene names.
    #'
    #' @author Sanjeev Sariya
    #'
    
    SNP<-START<-END<-GENE<-BP<-NULL ## bind variable locally to the function
    df_genecoords<-data.table::as.data.table(df_genecoords)
    df_snppos<-data.table::as.data.table(df_snppos)
    
    if(is.null(keep_indiv)==TRUE){
        genotyped_extracted<-vcfR::extract.gt(vcf_data,element = "GT",as.numeric=TRUE,convertNA=TRUE) 
    }else{
        
        genotyped_extracted <- tryCatch({
            vcfR::extract.gt(vcf_data,element = "GT",as.numeric=TRUE,convertNA=TRUE)[,keep_indiv]
            
        }, warning = function(w) {
            message(paste("warning vcf_counts_SNP_genecoords ", w))
            
        }, error =function(e) {
            message(paste(" vcf counts SNP genecoords: error subsetting individuals ", e))
            stop("Exiting vcf counts SNP genecoords ")
            
        }) 
        
    } ## else ends for checking individual sub-setting 
    
    df_genotyped_extracted<-data.table::data.table(genotyped_extracted,keep.rownames="SNP") ## assign row names as SNP to perform a join
    if(is.null(extract_SNP) == FALSE){
        
        extract_SNP<-as.character(extract_SNP)
        
        df_genotyped_extracted<-tryCatch({
            df_genotyped_extracted[SNP %in% extract_SNP,]
            
        }, warning = function(w) {
            message(paste("warning vcfcounts SNPgenecoords in subsetting SNPs", w))
            
        }, error =function(e) {
            message(paste("vcf countsSNP genecoords: error subseting SNPs ", e))
            stop("Exiting vcf counts annot ")
            
        } ) 
    } ## else ends for subseting SNPs 
    ## subsetting based on SNPs ends

    if(is.null(filter_gene) == FALSE){
        filter_gene<-as.character(filter_gene) ## turn into character
        df_genecoords<-garcom_filter_gene(df_genecoords,filter_gene) ##filter genes based on Gene list vector
    }
    ##subsetting ends for gene filtering
    
    snp_withingenes<-df_snppos[df_genecoords, c("SNP","GENE","START","END"), on=list(BP>=START, BP<=END), nomatch=0]
    
    if(nrow(snp_withingenes) == 0){
        stop("VCF counts SNP pos: No snps within any gene boundaries provided")	
    }
    
    if(isFALSE( (any( (df_genotyped_extracted$SNP) %in% unique(snp_withingenes$SNP))))){
        stop("VCF counts SNP pos: No SNPs overlapping between genetic data and SNP annotation with Gene boundaries")
    }
    ## initial checks end
    
    jointed_gene_VCFGT<-snp_withingenes[df_genotyped_extracted, on=c(SNP="SNP"),nomatch=0L]
    jointed_gene_VCFGT[,c("START","END","SNP"):=NULL]  
    
    ##https://stackoverflow.com/a/32277135/2740831
    matrix_withallelecount_withinGene <-jointed_gene_VCFGT[,lapply(.SD,sum,na.rm=TRUE),by=GENE] 
    matrix_withallelecount_withinGene<-matrix_withallelecount_withinGene[ rowSums(matrix_withallelecount_withinGene[,-c("GENE")]) > 0,]
    
    if(nrow(matrix_withallelecount_withinGene)>0){
        return(matrix_withallelecount_withinGene)
    }
    else{
        return(NULL) ## if nothing is left after gene >0 
    }
    
} 
## function ends
