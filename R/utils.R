#07 21 2020

garcom_check_column_names<-function(temp_data,col_names){
## 07 21 2020

##two params: dataframe and vector of columns we are looking for.
##return TRUE or FALSE
return(col_names %in% colnames(temp_data))

}
##function ends

garcom_check_duplicates<-function(temp_data, column_name){

##08/11/2020
##two params: dataframe and second is column name for which we
##we'd like to ensure no duplicates are present
##
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
if( nrow(unique(temp_data)) !=nrow(temp_data)) {
return(FALSE)
} else{
return(TRUE)
}
##check ends

}
##function ends

garcom_subsetIIDs<-function(tempdata,iids_to_keep){
##08/11/2020

if(isTRUE(anyNA(iids_to_keep))){

stop("There are individuals with IIDs as NA in the list provided to extract")
}
##if check ends for any NA

if( isTRUE( any("" == iids_to_keep)) | isTRUE( any('' == iids_to_keep)) ){
stop("There is an IID with length zero. Cannot extract it. Exiting... ")

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
print("we are in the subset of SNPs")

}
##function ends