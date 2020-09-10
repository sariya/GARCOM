#'@title gene position counts
#'@description Function returns matrix with allelic counts per gene per individual for SNP and gene coordinates as inputs
#'
#'@details Inputs needed are: recoded genetic data formatted in PLINK format, SNP name with BP (position) and gene name with START and END position. The first six columns of the input genetic data follow standard PLINK .raw format. Column names as FID, IID, PAT, MAT, SEX and PHENOTYPE followed by SNP information as recoded by the PLINK software. The function returns allelic counts per gene per sample (where each row represents a gene and each column represents an individual starting with the second column where first column contains gene information). 

#' @usage gene_pos_counts(dt_gen,dt_snp,dt_gene,keep_indiv=NULL,
#' extract_SNP=NULL,filter_gene=NULL,
#' impute_missing=FALSE,impute_method="mean")
#'

gene_pos_counts<-function(dt_gen,dt_snp,dt_gene,keep_indiv=NULL,extract_SNP=NULL,filter_gene=NULL,impute_missing=FALSE,impute_method="mean"){ 

    ##07 20 2020
    
    #' @export
    #' @import data.table
    #' @importFrom data.table rowid
    #' @importFrom data.table .SD
    #' @importFrom data.table :=
    #' @param dt_gen a dataframe for genetic data that follows PLINK format (.raw)
    #'
    #' @param dt_gene a dataframe for gene boundaries with CHR START END GENE as column names. Where CHR should be integer 1-22. START and END column should be integer. GENE column contains gene names
    #' @param dt_snp a dataframe for SNP information with SNP BP as column names.
    #'
    #' @param keep_indiv an option to specify individuals to retain. Mutation counts will be provided for individuals provided in the list only. Default is all individuals. 
    #' @param extract_SNP an option to specify SNPs for which mutation counts are needed. Mutation counts will be provided for SNPs included in the list only. Default is all SNPs.
    #' @param filter_gene an option to filter in Genes. Mutation counts will be provided for genes included in the list only. Default is all genes.
    #'
    #' @param impute_missing an option to impute missing genotypes. Default is FALSE. 
    #'
    #' @param impute_method an option to specify method to specify imptuation method. Default method is impute to the mean. Alternatively imputation can be carried out by median. Function accepts method in quotes: "mean" or "median". Data are rounded to the second decimal places (e.g. 0.1234 will become 0.12.).
    #'
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
    #'
    #' #subset individuals
    #' gene_pos_counts(recodedgen, snppos, genecoord,keep_indiv=c("IID_sample2","IID_sample4"))
    #'
    #' #subset genes
    #' gene_pos_counts(recodedgen,snppos,genecoord,filter_gene=c("GENE1","GENE2")) 
    #'
    #' #subset genes and individual iids
    #' gene_pos_counts(recodedgen,snppos,genecoord,filter_gene=c("GENE1","GENE2"),
    #' keep_indiv=c("IID_sample10","IID_sample4")) 
    #'
    #' ##impute by mean
    #' gene_pos_counts(recodedgen,snppos,genecoord,impute_missing=TRUE,impute_method="mean")
    #'
    #' #end not RUN
    #'
    #' @return Returns an object of data.table class as an output with allelic gene counts within each sample where each row corresponds to gene and column to individual IDs from column second. The first column contains gene names.
    #'
    #' @author Sanjeev Sariya
    #'
    
    IID<-START<-END<-GENE<-BP<-NULL ## bind variable locally to the function
    dt_gen<-data.table::as.data.table(dt_gen) # convert into data.table
    dt_gen[, IID:=as.character(IID)] ## convert into character in case IIDs are integer values. 
    dt_gene<-data.table::as.data.table(dt_gene)
    dt_snp<-data.table::as.data.table(dt_snp)

    if(all(garcom_check_column_names(dt_gene, c("START","END","GENE")))){
        ## all good with gene data
    }else{
        stop("gene pos: column names don't match for gene data")
    }

    if(all(garcom_check_column_names(dt_snp, c("SNP","BP")))){
        ## all good with SNP data
    }else{
        stop("gene pos: column names don't match for snp data")
    }
    ##Check ends 

    if(FALSE == isTRUE(garcom_check_duplicates(dt_snp,"SNP"))){

        stop("gene pos: duplicate SNP names")
    }
    ##check ends for SNP data.table
    if(FALSE == isTRUE(garcom_check_duplicates(dt_gene,"GENE"))){

        stop("gene pos: duplicate GENE names")
    }
    ##check ends for GENE data.table

    if(is.null(keep_indiv) == FALSE ){
        keep_indiv<-as.character(keep_indiv) ## convert them into character
        dt_gen<-garcom_subsetIIDs(dt_gen,keep_indiv) ## it returned a sub-setted data with iids of interest

    }
    ## sub-setting complete for individuals interested

    if(is.null(extract_SNP) == FALSE){
        extract_SNP<-as.character(extract_SNP)
        dt_snp<-garcom_subsetSNPs(dt_snp,extract_SNP) ## returns data with overlapping SNPs
    }
    ##check ends for sub-setting SNPs

    if(is.null(filter_gene) == FALSE){
        filter_gene<-as.character(filter_gene) ## turn into character
        dt_gene<-garcom_filter_gene(dt_gene,filter_gene) ##filter SNP-gene annotation based on Gene list
    }
    ##check ends for filtering Genes

    if(isTRUE(impute_missing)){

        ##we pass impute method and genetic data frame
        dt_gen<-garcom_impute(dt_gen,impute_method)
    }

    ##impute check ends
    colnames(dt_gen) <- gsub("_.*","",colnames(dt_gen)) ##Remove _ from recode format

    ## https://gist.github.com/nacnudus/ef3b22b79164bbf9c0ebafbf558f22a0

    snp_withingenes<-dt_snp[dt_gene, c("SNP","GENE","START","END"), on=list(BP>=START, BP<=END), nomatch=0] # inner join ##https://stackoverflow.com/questions/63290994/foverlaps-data-table-error-ys-key-must-be-identical-to-the-columns-specified

    if(nrow(snp_withingenes) == 0){
        stop("gene pos: No snps within any gene boundaries provided")	
    }
    ##if gene sum is zero. Stop and return

    if(isFALSE( (any( colnames(dt_gen) %in% unique(snp_withingenes$SNP))))){
        stop("gene pos: No SNPs overlapping between genetic data and SNP annotation with Gene boundaries")
    }
    ##if nothing matches then Stop and error out

    dt_gen_subset<- dt_gen[,.SD,.SDcols=intersect(colnames(dt_gen),unique(snp_withingenes$SNP))] ## get columns that intersect 
    dt_gen_subset[, rowid := dt_gen$IID ] ## assign row names as IIDs
    dt_gen_subset<-data.table::transpose(dt_gen_subset,keep.names = "SNP", make.names="rowid") ## tranpose data and have some fancy settings

    ## we perform inner join. SNPs that are found in input boundaries-annotation as well as .raw data

    subsetsnps_genes_lefted_join <- snp_withingenes[dt_gen_subset,on="SNP",nomatch=0] 
    subsetsnps_genes_lefted_join[,c("START","END","SNP"):=NULL]  
    
    ##https://stackoverflow.com/a/32277135/2740831
    matrix_withallelecount_withinGene <-subsetsnps_genes_lefted_join[,lapply(.SD,sum,na.rm=TRUE),by=GENE] 

    matrix_withallelecount_withinGene<-matrix_withallelecount_withinGene[ rowSums(matrix_withallelecount_withinGene[,-c("GENE")]) > 0,]
    
    if(nrow(matrix_withallelecount_withinGene)>0){
        return(matrix_withallelecount_withinGene)
    }
    else{
        return(NULL) ## if nothing is left after gene >0 
    }

} ## function ends
