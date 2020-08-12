#'@title gene annotation counts
#'@description Inputs needed are recoded genetic data formatted in PLINK format and SNP-gene annotation data . The first six columns of the input genetic data follow standard PLINK .raw formt. Column names as FID, IID, PAT, MAT, SEX and PHENOTYPE followed by SNP information as recoded by the PLINK software. SNP-gene data has two columns: GENE and SNP names. The function returns allelic counts per gene per sample (where each row represents a gene and each column represents an individual starting with the second column where first column contains gene information). 

gene_annot_counts<-function(dt_gen,dt_snpgene, keep_indiv=NULL, extract_SNP=NULL,filter_gene=NULL,impute_missing=FALSE){
##07 10 2020

#' @export
#' @import data.table
#' @importFrom dplyr %>%
#' @importFrom data.table :=
#' @param dt_gen recoded genetic data from PLINK
#' @param dt_snpgene with SNP and GENE as column names
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
#' #end not RUN
#'
#' @return Returns an object of data.table class as an output with allelic gene counts within each sample where each row corresponds to gene and column to individual IDs from column second. The first column contains gene names.
#'
#' @author Sanjeev Sariya
#'
    
    dt_gen<-data.table::as.data.table(dt_gen) ## make data.table format for higher speed
dt_gen[, IID:=as.character(IID)] ## convert into character in case IIDs are integer values. 

    dt_snpgene<-data.table::as.data.table(dt_snpgene)

    SNP<-GENE<-NULL ## binding the variable locally to the function
    ##https://www.r-bloggers.com/no-visible-binding-for-global-variable/
    
    if(all(garcom_check_column_names(dt_snpgene, c("SNP","GENE")))){
        ## all good with SNP data
    }else{
        stop("column names don't match for snp-gene data")
    }
    ## check ends

if(FALSE == isTRUE(garcom_check_unique(dt_snpgene) )){
stop("Duplicate SNP-Gene annotation values")
}

if(is.null(keep_indiv) == FALSE ){
keep_indiv<-as.character(keep_indiv) ## convert them into character
dt_gen<-garcom_subsetIIDs(dt_gen,keep_indiv) ## it returned a sub-setted data with iids of interest

}
###check ends for sub-setting IIDs

if(is.null(extract_SNP) == FALSE){
extract_SNP<-as.character(extract_SNP)
dt_snpgene<-garcom_subsetSNPs(dt_snpgene,extract_SNP)

}
###check ends for sub-setting SNPs

if(is.null(filter_gene) == FALSE){
filter_gene<-as.character(filter_gene) ## turn into character
dt_snpgene<-garcom_filter_gene(dt_snpgene,filter_gene) ##filter SNP-gene annotation based on Gene list

}
###check ends for sub-setting Genes

    colnames(dt_gen) <- gsub("_.*","",colnames(dt_gen)) ## remove underscore generate from plink

    IID_samples<-data.frame("IID"=dt_gen[,"IID"])
    IID_samples$IID<-as.character(IID_samples$IID) ##make character we can have IIDs as numbers

    SNP_names<-colnames(dt_gen)[c(7:length(colnames(dt_gen)))] # use this to assign SNP column when piping

    dt_gen_transposed<- data.table(data.table::transpose(dt_gen) )
    
    dt_gen_transposed<-dt_gen_transposed[,.SD[-1:-6]] ## remove first six rows

    data.table::setnames(dt_gen_transposed,IID_samples$IID) ## set column names
    dt_gen_filtered<-dt_gen_transposed[, c("SNP") := SNP_names ]

    ##https://gist.github.com/nacnudus/ef3b22b79164bbf9c0ebafbf558f22a0
  
    jointed_genesSNP<-dt_snpgene[dt_gen_filtered , on="SNP", nomatch=0] ## do a left join on the data.table RAW
    jointed_genesSNP<- jointed_genesSNP[,SNP:=NULL]   ## remove SNP column

    if(nrow(jointed_genesSNP)==0){
        message("No SNPs match with the annotation")
        return(NULL)
    }
    ##check if join gives any rows. If not return NULL

    jointed_genesSNP<-jointed_genesSNP[, lapply(.SD, as.numeric), by="GENE"] ## convert into numeric
    ##https://stackoverflow.com/a/62959318/2740831

    jointed_genesSNP_filtered<-jointed_genesSNP[,lapply(.SD,sum,na.rm=TRUE),by=GENE] # get sum within a gene
    jointed_genesSNP_filtered<-jointed_genesSNP_filtered[ rowSums(jointed_genesSNP_filtered[,-c("GENE")]) > 0,]  # get count minus gene column and keep only genes with sum more than 0. test with a test case here

    if(nrow(jointed_genesSNP_filtered) ==0){
        message("All genes with zero count")
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
