#'@title gene annotation counts
#'@description The function returns a matrix with allelic counts per gene per individual for annotated SNPs
#'
#' @usage gene_annot_counts(dt_gen,dt_snpgene,keep_indiv=NULL,
#' extract_SNP=NULL,filter_gene=NULL,
#' impute_missing=FALSE,impute_method="mean")

gene_annot_counts<-function(dt_gen,dt_snpgene,keep_indiv=NULL,extract_SNP=NULL,filter_gene=NULL,impute_missing=FALSE,impute_method="mean"){
    ##07 10 2020
    
    #' @export
    #' @import data.table
    #' @importFrom data.table :=
    #' @param dt_gen a dataframe for genetic data that follows PLINK format (.raw) 
    #' @param dt_snpgene a dataframe that contains SNP and annotated gene with SNP and GENE as column name
    #'
    #' @param keep_indiv an option to specify individuals to retain. Mutation counts will be provided for individuals included in the list only. Default is all individuals. Provide list of individuals in a vector.
    #'
    #' @param extract_SNP an option to specify SNPs for which mutation counts are needed. Mutation counts will be provided for SNPs provided in the list only. Default all SNPs are used. Provide list of SNPs in a vector.
    #' @param filter_gene an option to filter in a list of Genes. Mutation counts will be provided for genes specifed in the list only. Default is all genes. Provide list of genes in a vector.
    #'
    #' @param impute_missing an option to impute missing genotypes. Default is FALSE. 
    #'
    #' @param impute_method an option to specify imptuation method. Default method is imputation to the mean. Alternatively imputation can be carried out by median. Function accepts method in quotes: "mean" or "median". Data are rounded to the second decimal places (e.g. 0.1234 will become 0.12).
    #'
    #'
    #' @details Inputs needed are recoded genetic data formatted in PLINK format (.raw) and SNP-gene annotation data. The first six columns of the input genetic data follow standard PLINK .raw format. Column names as FID, IID, PAT, MAT, SEX and PHENOTYPE followed by SNP information as recoded by the PLINK software. SNP-gene data has two columns: GENE and SNP names. The function returns allelic counts per gene per sample (where each row represents a gene and each column represents an individual starting with the second column where first column contains gene information). 
    #'
    #'
    #' @examples
    #'
    #' #Package provides sample data that are loaded with package loading. 
    #'
    #' data(recodedgen) #PLINK raw formatted data of 10 individiduals with 10 SNPs
    #'
    #' data(snpgene) #SNP and its respective GENE annotated. 
    #' #Here 10 SNPs are shown annotated in five genes. 
    #' #A SNP can be annotated in multiple genes. 
    #'
    #' gene_annot_counts(recodedgen,snpgene) #run the function
    #'
    #' #subset Genes
    #' gene_annot_counts(recodedgen,snpgene,filter_gene=c("GENE1","GENE2"))
    #'
    #' #Subset individuals
    #' gene_annot_counts(recodedgen, snpgene,keep_indiv=c("IID_sample1","IID_sample8"))
    #'
    #' #subset with genes and samples
    #' gene_annot_counts(recodedgen,snpgene,filter_gene=c("GENE1","GENE2"),
    #' keep_indiv=c("IID_sample1","IID_sample8"))
    #'
    #' #impute missing using default method. 
    #'
    #' gene_annot_counts(recodedgen,snpgene,impute_missing=TRUE)
    #'
    #' #Subset on individuals and impute for missing values. Default as mean
    #' gene_annot_counts(recodedgen,snpgene,impute_missing=TRUE,
    #' keep_indiv=c("IID_sample1","IID_sample2","IID_sample10"))
    #'
    #' #impute using median method
    #' gene_annot_counts(recodedgen,snpgene,impute_missing=TRUE,impute_method="median")
    #'
    #' #end not RUN
    #'
    #' @return Returns an object of data.table class as an output with allelic gene counts within each sample where each row corresponds to gene and column to individual IDs from column second. The first column contains gene names.
    #'
    #' @author Sanjeev Sariya
    #'
    
    IID<-SNP<-GENE<-NULL ## binding the variable locally to the function    
    dt_gen<-data.table::as.data.table(dt_gen) ## make data.table format for higher speed
    dt_gen[, IID:=as.character(IID)] ## convert into character in case IIDs are integer values. 

    dt_snpgene<-data.table::as.data.table(dt_snpgene)

    ##https://www.r-bloggers.com/no-visible-binding-for-global-variable/
    
    if(all(garcom_check_column_names(dt_snpgene,c("SNP","GENE")))){
        ## all good with SNP data
    }else{
        stop("gene annot: column names don't match for snp-gene data")
    }
    ## check ends

    if(FALSE == isTRUE(garcom_check_unique(dt_snpgene) )){
        stop("gene annot: Duplicate SNP-Gene annotation values")
    }

    if(is.null(keep_indiv) == FALSE ){
        keep_indiv<-as.character(keep_indiv) ## convert them into character
        dt_gen<-garcom_subsetIIDs(dt_gen,keep_indiv) ## it returned a sub-setted data with iids of interest

    }
    ##check ends for sub-setting IIDs

    if(is.null(extract_SNP) == FALSE){
        extract_SNP<-as.character(extract_SNP)
        dt_snpgene<-garcom_subsetSNPs(dt_snpgene,extract_SNP)
    }
    ##check ends for sub-setting SNPs

    if(is.null(filter_gene) == FALSE){
        ##
        ##Start process to filter genes. The list provided by user is what we'd like to keep
        filter_gene<-as.character(filter_gene) ## turn into character
        dt_snpgene<-garcom_filter_gene(dt_snpgene,filter_gene) ##filter SNP-gene annotation based on Gene list
    }
    ##check ends for sub-setting Genes

    if(isTRUE(impute_missing)){

        ##we pass impute method and genetic data frame
        dt_gen<-garcom_impute(dt_gen,impute_method)
    }
    ##check ends for imputing genetic data

    colnames(dt_gen) <- gsub("_.*","",colnames(dt_gen)) ## remove underscore generate from plink

    IID_samples<-data.frame("IID"=dt_gen[,"IID"])
    IID_samples$IID<-as.character(IID_samples$IID) ##make character we can have IIDs as numbers

    SNP_names<-colnames(dt_gen)[c(7:length(colnames(dt_gen)))] # use this to assign SNP column when piping

    dt_gen_transposed<-data.table(data.table::transpose(dt_gen) )
    
    dt_gen_transposed<-dt_gen_transposed[,.SD[-1:-6]] ## remove first six rows

    data.table::setnames(dt_gen_transposed,IID_samples$IID) ## set column names
    dt_gen_filtered<-dt_gen_transposed[,c("SNP") := SNP_names]

    ##https://gist.github.com/nacnudus/ef3b22b79164bbf9c0ebafbf558f22a0
  
    jointed_genesSNP<-dt_snpgene[dt_gen_filtered,on="SNP", nomatch=0] ## do a left join on the data.table RAW
    jointed_genesSNP<-jointed_genesSNP[,SNP:=NULL] ## remove SNP column

    if(nrow(jointed_genesSNP)==0){
        message("gene annot: No SNPs match with the annotation")
        return(NULL)
    }
    ##check if join gives any rows. If not return NULL

    jointed_genesSNP<-jointed_genesSNP[, lapply(.SD,as.numeric), by="GENE"] ## convert into numeric
    ##https://stackoverflow.com/a/62959318/2740831

    jointed_genesSNP_filtered<-jointed_genesSNP[,lapply(.SD,sum,na.rm=TRUE),by=GENE] # get sum within a gene
    jointed_genesSNP_filtered<-jointed_genesSNP_filtered[rowSums(jointed_genesSNP_filtered[,-c("GENE")]) > 0,] ##get count minus gene column and keep only genes with sum more than 0. test with a test case here

    if(nrow(jointed_genesSNP_filtered)==0){
        message("gene annot: All genes with zero count")
        return(NULL)
    }
    else{
        return(jointed_genesSNP_filtered)
    }

    ##
    ##https://stackoverflow.com/questions/50768717/failure-using-data-table-inside-package-function
    ## https://stackoverflow.com/questions/10527072/using-data-table-package-inside-my-own-package

} ## function ends 
###
