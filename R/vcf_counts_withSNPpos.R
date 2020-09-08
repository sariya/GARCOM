## added on 09 08 2020

#'@title VCF gene position counts
#'@description Function returns a matrix with allelic counts per gene per individual for SNP and gene coordinates as inputs
#'

vcf_counts_SNP_genecoords<-function(vcf_data,df_snppos,df_genecoords){
    
    #' @export
    #' @import data.table
    #' @import vcfR
    
    #' @param vcf_data an object of vcfR class
    
    #' @param df_snppos a dataframe for SNP information with SNP BP as column names.
    #' @param df_gene_coords a dataframe for gene boundaries with CHR START END GENE as column names. Where CHR should be integer 1-22. START and END column should be integer. GENE column contains gene names
    #'
    #'@examples 
    #'\dontrun{
    #' vcf_counts_SNP_genecoords(vcf_data_test,df_snppos_test,df_genecoords_test)
    #' }
    #'
    #' @return Returns an matrix of data.table class as an output with allelic (reference) gene counts within each sample where each row corresponds to gene and column to individual IDs from column second. The first column contains gene names.
    #'
    #' @author Sanjeev Sariya
    #'
    
    IID<-START<-END<-GENE<-BP<-NULL ## bind variable locally to the function
    
    genotyped_extracted<-vcfR::extract.gt(vcf_data,element = "GT",as.numeric=TRUE,convertNA=TRUE) 
    
    df_genotyped_extracted<-data.table::data.table(genotyped_extracted,keep.rownames=TRUE) #3 use rn later while merging
    
    snp_withingenes<-df_snppos[df_genecoords, c("SNP","GENE","START","END"), on=list(BP>=START, BP<=END), nomatch=0] 
    if(nrow(snp_withingenes) == 0){
        stop("VCF counts SNP pos: No snps within any gene boundaries provided")	
    }
    
    if(isFALSE( (any( (df_genotyped_extracted$rn) %in% unique(snp_withingenes$SNP))))){
        stop("VCF counts SNP pos: No SNPs overlapping between genetic data and SNP annotation with Gene boundaries")
    }
    ## initial checks end
    
    jointed_gene_VCFGT<-snp_withingenes[df_genotyped_extracted, on=c(SNP="rn"), nomatch=0L]
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

