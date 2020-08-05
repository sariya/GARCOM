#'@title gene annotation counts
#'@description Inputs as recoded genetic data formatted in PLINK format and SNP-gene annotation data as parameters. The first six columns of the input genetic data follow standard PLINK .raw formt. Column names as FID, IID, PAT, MAT, SEX and PHENOTYPE followed by SNP information as recoded by the PLINK software. The function returns allelic counts per gene per sample (where each row represents a gene and each column represents an individual starting with the second column where first column contains gene information). 

gene_annot_counts<-function(dt_gen,dt_snpgene){
##07 10 2020

#' @export
#' @param dt_gen recoded genetic data from PLINK
#' @param dt_snpgene with SNP and GENE as column names
#' @examples
#'
#' Package provides sample data that are loaded with package loading. 
#'
#' data(recodedgen) #PLINK raw formatted data of 10 individiduals with 10 SNPs
#'
#' data(snpgene) #SNP and its respective GENE annotated. 
#' #Here 10 SNPs are shown annotated in five genes. 
#' #A SNP can be annotated in multiple genes. 
#'
#' small_output<-gene_annot_counts(recodedgen,snpgene) #run the function
#'
#' #end not RUN
#'
#' @return Returns an object of data.table class as an output with allelic gene counts within each sample where each row corresponds to gene and column to individual IDs from column second. The first column contains gene names.
#'
#' @author Sanjeev Sariya
#'
    
    dt_gen<-data.table::as.data.table(dt_gen) ## make data.table format for higher speed
    dt_snpgene<-data.table::as.data.table(dt_snpgene)

    if(all(garcom_check_column_names(dt_snpgene, c("SNP","GENE")))){
        ## all good with SNP data
    }else{
        stop("column names don't match for snp-gene data")
    }
    ## check ends

    colnames(dt_gen) <- colnames(dt_gen) %>% gsub("_.*","",.) ## remove underscore genearte from plink

    IID_samples<-as.data.frame(dt_gen[,2]) %>% `colnames<-` (c("IID")) ## use this later

    SNP_names<-colnames(dt_gen)[c(7:length(colnames(dt_gen)))] # use this to assign SNP column when piping

    dt_gen_filtered<- data.table::transpose(dt_gen) %>% .[,.SD[-1:-6]] %>% data.table::setnames(.,IID_samples$IID) %>% .[, c("SNP") := SNP_names ]

    ##https://gist.github.com/nacnudus/ef3b22b79164bbf9c0ebafbf558f22a0
    jointed_genesSNP<-dt_snpgene[dt_gen_filtered , on="SNP", nomatch=0] %>% .[,SNP:=NULL]   ## do a left join on the data.table RAW and remove SNP column

    jointed_genesSNP<-jointed_genesSNP[, lapply(.SD, as.numeric), by="GENE"] ## convert into numeric
    ##https://stackoverflow.com/a/62959318/2740831

    jointed_genesSNP_filtered<-jointed_genesSNP[,lapply(.SD,sum,na.rm=TRUE),by=GENE] %>% .[ rowSums(.[,-c("GENE")]) > 0,] 

    if(nrow(jointed_genesSNP_filtered) ==0){
        print("All genes with zero count")
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
