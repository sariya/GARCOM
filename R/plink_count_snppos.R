#'@title plink gene counts using SNP coords and gene boundaries
#'
#'@description Function returns matrix with allelic counts per gene per individual for SNP and gene coordinates as inputs for PLINK (.bed) file format
#'
#'@details PLINK (.bed) file will be read
#'
#' @usage plink_count_snppos(plink_file,genecoord_dt,
#' snp_index,individuals_index)
#'

plink_count_snppos<-function(plink_file,genecoord_dt,snp_index=NULL,individual_index=NULL){

#' @export
#'
#' @import data.table
#'
#' @importFrom bigsnpr snp_attach
#' 
#' @importFrom bigsnpr snp_readBed2

#' @param plink_file plink .bed file with path 
#'
#' @param genecoord_dt a dataframe for gene boundaries with CHR START END GENE as column names. Where CHR should be integer 1-22. START and END column should be integer. GENE column contains gene names
#'
#' @param snp_index a vector of integer that specifies SNPs to read. Default all SNPs will be read.
#'
#' @param individual_index a vector of integer that specifies individuals to select. Default all individuals will be read.
#'
#' @examples 
#' \dontrun{
#' plink_count_snppos(path_plinkbed_file,data_genecoord,snp_index,individual_select)
#' }
#'
#' @author Sanjeev Sariya
#'
    START<-END<-GENE<-BP<-NULL ## bind variable locally to the function
    plink_rds <-NULL
    genecoord_dt <- data.table::as.data.table(genecoord_dt) #save CHR, BP, START END as data.table
    
    if(all(garcom_check_column_names(genecoord_dt,c("START","END","GENE")))){
        ## all good with gene data
    }else{
        stop("column names don't match for gene data")
    }
    
    ##check ends for SNP data.table
    if(FALSE == isTRUE(garcom_check_duplicates(genecoord_dt,"GENE"))){

        stop("duplicate GENE names")
    }
    ##check ends for GENE data.table

    if((is.null(individual_index)==FALSE ) & (is.null(snp_index)==FALSE) ){
        cat("User provided snp and individuals to select\n")
        plink_rds <- bigsnpr::snp_readBed2(plink_file,backingfile=tempfile(),ind.col=snp_index,ind.row=individual_index)
        
    }
    if( (is.null(individual_index)==TRUE) & (is.null(snp_index)==FALSE) ){
        cat("User provided snp  to select\n")
        plink_rds <- bigsnpr::snp_readBed2(plink_file,backingfile=tempfile(),ind.col=snp_index)
    }
    
    if( (is.null(individual_index)==FALSE) & (is.null(snp_index)==TRUE) ){
        cat("User provided individuals to select\n")
        plink_rds <- bigsnpr::snp_readBed2(plink_file,backingfile=tempfile(),ind.row=individual_index)
    }    
    
    if( (is.null(individual_index)==TRUE) & (is.null(snp_index)==TRUE) ){
        cat("No user no SNPs selected. Load complete data\n")
        plink_rds <- bigsnpr::snp_readBed2(plink_file,backingfile=tempfile())
    }    
    
    ## Loading the data from backing files
    data_plink <- bigsnpr::snp_attach(plink_rds)
    
    plink_snp_information_dt<-as.data.table(data_plink$map)
    plink_fam_information<-as.data.table(data_plink$fam)
    
    plink_genotype_dt<- as.data.table(data_plink$genotypes[])
    colnames(plink_genotype_dt) <-plink_snp_information_dt$marker.ID
    cat("Genotypes have been loaded from plink file\n")
    
    snp_withingenes<- plink_snp_information_dt[genecoord_dt,c("marker.ID","GENE","START","END"),on=list(physical.pos>=START,physical.pos<=END),nomatch=0] # inner join ##https://stackoverflow.com/questions/63290994/foverlaps-data-table-error-ys-key-must-be-identical-to-the-columns-specified
    if(nrow(snp_withingenes) == 0){
        stop("No snps within any gene boundaries provided")	
    }
    
    ## we check if the SNPs that are found in plink map data overlap with inner-joined data 
    if(isFALSE( (any( plink_snp_information_dt$marker.ID %in% unique(snp_withingenes$marker.ID))))){
        stop("No SNPs overlapping between genetic data and SNP annotation with Gene boundaries")
    }
    
    ##get index of SNPs and store in a vector to be used while subsetting plink genotype data.table
    snps_intersect_index<- base::match(intersect(plink_snp_information_dt$marker.ID,unique(snp_withingenes$marker.ID)),plink_snp_information_dt$marker.ID) 
    
    ## subset plink data
    plink_genotype_dt_subset<-plink_genotype_dt[,.SD,.SDcols=snps_intersect_index ] ## get columns that intersect 
    plink_genotype_dt_subset[,rowid := plink_fam_information$sample.ID ] ## assign row names as IIDs
    
    plink_genotype_dt_subset_transposed<-data.table::transpose(plink_genotype_dt_subset,keep.names="SNP",make.names="rowid") ## tranpose data and have some fancy settings
    
    subsetsnps_genes_lefted_join <-snp_withingenes[plink_genotype_dt_subset_transposed,on=c("marker.ID"="SNP"),nomatch=0] 
    subsetsnps_genes_lefted_join[,c("START","END","marker.ID"):=NULL]  
    matrix_withallelecount_withinGene <-subsetsnps_genes_lefted_join[,lapply(.SD,sum,na.rm=TRUE),by=GENE] 
    matrix_withallelecount_withinGene<-matrix_withallelecount_withinGene[rowSums(matrix_withallelecount_withinGene[,-c("GENE")]) > 0,]
    
    if(nrow(matrix_withallelecount_withinGene)>0){
        return(matrix_withallelecount_withinGene)
    }
    else{
        return(NULL) ## if nothing is left after gene >0 
    }
    
}
## function ends


