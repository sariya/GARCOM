#07 21 2020

#' @importFrom stats median 

garcom_check_column_names<-function(temp_data,col_names){
## 07 21 2020

    ##function check column names for any dataframe and supplied vector of column names
    ##two params: dataframe and vector of columns we are looking for.
    ##return TRUE or FALSE
    
    return(col_names %in% colnames(temp_data))

}
##function ends

garcom_check_duplicates<-function(temp_data, column_name){

    ##08/11/2020
    ##two params: dataframe and second is column name for which we
    ##we'd like to ensure no duplicates are present

    if( length(unique(temp_data[,get(column_name)])) !=length(temp_data[,get(column_name)])){
        return(FALSE)
    } else{
        return(TRUE)
    }
    ##check ends

}
##function ends

garcom_check_unique<-function(temp_data){

    ##08/11/2020
    ##one param: dataframe 
    ##we'd like to ensure no duplicates are present
    ## This is conducted with all columns. to ensure no row duplicates are present
    if( nrow(unique(temp_data))!=nrow(temp_data)) {
        return(FALSE)
    } else{
        return(TRUE)
    }
##check ends

}
##function ends

garcom_subsetIIDs<-function(tempdata,iids_to_keep){
    ##08/11/2020

    ##two params: data and individual ids to keep
    ##we return genetic data with subset IIDs 
    ##We stop if IIDs in the vector aren't found
    ##we subset provided dataframe with vector of IIDs supplied  

    if(isTRUE(anyNA(iids_to_keep))){

        stop("There are individuals with IIDs as NA in the list provided to extract. Exiting...")
    }
    ##if check ends for any NA in the IIDs to extract

    if( isTRUE( any("" == iids_to_keep)) | isTRUE( any('' == iids_to_keep)) ){
        stop("There is an IID with length zero name. Cannot extract it. Exiting... ")
    }
    ##check ends for any ids with length as zero

    row_index_gendata_subset<-match(iids_to_keep,tempdata[,get("IID")]) ## store indices to return data.table with selected individuals
    intersected_IIDs<-intersect(iids_to_keep,tempdata[,get("IID")]) ## check length later with this

    if(length(intersected_IIDs) == 0){
        stop("no iids intersect. Exiting...")
    }
    ##if to check length of IIDs intersect is zero ends

    return(tempdata[row_index_gendata_subset,])
}
##function ends

garcom_subsetSNPs<-function(tempdata,snps_to_keep){
    ##08/11/2020

    ##two params: SNP data and second parameter is vector of SNPs that are to be subsetted
    ##we return subset data
    ##Stop if SNPs aren't subsetted
    ##we subset provided dataframe with vector of SNPs supplied
    
    if(isTRUE(anyNA(snps_to_keep))){

        stop("There are SNP names as NA in the list provided to extract. Exiting...")
    }
    ##if check ends for any NA in SNP

    if( isTRUE( any("" == snps_to_keep)) | isTRUE( any('' == snps_to_keep)) ){
        stop("There is SNP with length zero name. Cannot extract it. Exiting... ")
    }
    ##check ends for any SNPs with length as zero

    index_SNPs_subset<- which(tempdata[,get("SNP")]  %in% snps_to_keep )
    
    if(length(index_SNPs_subset) ==0){
        stop("SNPs not found to sub-set. Exiting...")
    }

    return(tempdata[index_SNPs_subset,]) ## return data that are extracted with list of interest
}
##function ends

garcom_filter_gene<-function(tempdata,filter_gene){

    ##08/11/2020

    ##Two parameters: gene data and vector of genes that are to be subsetted
    ##return if genes are subsetted
    ##Stop if no genes are found
    ##we subset provided dataframe with vector of genes supplied
    
    if(isTRUE(anyNA(filter_gene))){

        stop("There are Gene names as NA in the list provided to extract. Exiting...")
    }
    ##if check ends for any NA in GENE

    if( isTRUE( any("" == filter_gene)) | isTRUE( any('' == filter_gene)) ){
        stop("There is GENE with length zero name. Cannot extract it. Exiting... ")

    } ##check ends for any GENE with length as zero

    keep_gene_index<- which(tempdata[,get("GENE")]  %in% filter_gene )
    
    if(length(keep_gene_index) ==0){
        stop("Genes not found to sub-set. Exiting...")
    }

    return(tempdata[keep_gene_index,])
}
##function ends

garcom_impute<-function(temp_genetic,temp_impute_method){

    ##08 17 2020
    ##two params: genetic data and method to impute
    ##we exit if mean and median isn't provided.
    ##we impute data and then return imputed genetic data
    ## return rounded to two decimal places
    
    if(isTRUE(!(temp_impute_method %in% c("mean","median")))){
        stop("impute method doesn't match mean or median. Exiting")
    }

    if(temp_impute_method=="mean"){
        temp_genetic[]<-temp_genetic[, c(7:ncol(temp_genetic)):= lapply(.SD,function(x) ifelse(is.na(x),mean(x, na.rm = TRUE),x)),.SDcols = c(7:ncol(temp_genetic))]

    }
    ## end for mean imputation
    
    if(temp_impute_method=="median"){
        temp_genetic[]<-temp_genetic[, c(7:ncol(temp_genetic)):= lapply(.SD, function(x) ifelse(is.na(x),stats::median(x, na.rm = TRUE),x)),.SDcols = c(7:ncol(temp_genetic))]
    }
    
    ## end for median imputation
    
    ##round to two decimal places: 0.1234 will become 0.12 only. 0.00 will be 0.00
    ## https://stackoverflow.com/a/12135122/2740831
    temp_genetic[]<-temp_genetic[,c(7:ncol(temp_genetic)):=lapply(.SD,function(x) format(round(x,2),nsmall=2)), .SDcols=c(7:ncol(temp_genetic))]

    temp_genetic[,(c(seq(from=7,to=ncol(temp_genetic)))) := lapply(.SD,as.numeric),.SDcols = c(seq(from=7,to=ncol(temp_genetic)))] ## convert into numeric. For some weird reason everything is turned in character class when rounding up

    return(temp_genetic) ##return imputed data
}
##function ends
