#'@title Provide SNP index using genetic boundaries and PLINK bim file
#'
#'@description Function returns a vector of unique integer values for SNP position/index using genetic boundaries and plink bim file
#'
#'@details PLINK (.bim) file will be read
#'
#' @usage get_plinkSNP_index_geneboundaries(plinkbim_file,genecoord_dt)
#'

get_plinkSNP_index_geneboundaries<-function(plinkbim_file,genecoord_dt){
#'
#' @export
#'
#' @import data.table
#'
#' @param plink_file plink .bed file with path 
#'
#' @param genecoord_dt a dataframe for gene boundaries with CHR START END GENE as column names. Where CHR should be integer 1-22. START and END column should be integer. GENE column contains gene names
#'
#' @return a vector of integer that represents SNP index that are found within the gene boundaries and SNP bim file based on the interection. Returned values are unique if found, otherwise NULL is returned.
#'
#' @author Sanjeev Sariya
#'
    START<-END<-GENE<-BP<-NULL ## bind variable locally to the function
    genecoord_dt <- data.table::as.data.table(genecoord_dt)

    if(all(garcom_check_column_names(genecoord_dt,c("START","END","GENE")))){
        ## all good with gene data
    }else{
        stop("column names don't match for gene data\n")
    }
    if(FALSE == isTRUE(garcom_check_duplicates(genecoord_dt,"GENE"))){

        stop("Duplicate GENE names\n")
    }
    ##check ends for GENE data.table

    bimsnps_dt<-data.table::fread(plinkbim_file,header=FALSE)
    
    ##snp_withingenes<- bimsnps_dt[genecoord_dt,c("V2","GENE","START","END"),on=list(V4>=START,V4<=END),nomatch=0] # inner join ##https://stackoverflow.com/questions/63290994/foverlaps-data-table-error-ys-key-must-be-identical-to-the-columns-specified

    index_pergenelist<- with(genecoord_dt[bimsnps_dt,.(.I,GENE),on=.(START<=V4,END>= V4),nomatch=0],split(I,GENE))

    #unique_snpindex<-unique(match(snp_withingenes$V2,bimsnps_dt$V2)) # provide only unique # V2 in the bim file is SNP/rsid

    if(length(index_pergenelist)==0){
        cat("No intersection of SNPs found between Gene and SNP within the BIM file. No index\n")
        return (NULL)
    }
    else{
        return(index_pergenelist)
    }
     
}
## function ends added on Jan 10 2021


