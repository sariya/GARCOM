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

}
##function ends

garcom_subsetSNPs<-function(tempdata,snps_to_keep){

}
##function ends