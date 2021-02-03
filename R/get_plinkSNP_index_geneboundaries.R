#'@title Provide SNP index using genetic boundaries and PLINK bim file
#'
#'@description Function returns data frame for SNP position/index using genetic boundaries and plink bim file
#'
#'@details PLINK (.bim) file will be read
#'
#' @usage get_plinkSNP_index_geneboundaries(plinkbim_file,genecoord_dt)
#'

get_plinkSNP_index_geneboundaries<-function(plinkbim_file,genecoord_dt,freq_file_plink=NULL,threshold_freq_min=NULL,threshold_freq_max=NULL){
#'
#' @export
#'
#' @import data.table
#'
#' @param plink_file plink .bim file with path 
#' 
#' @param freq_file_plink .frq file generated from the PLINK file. It should contain same number of SNPs as .bim
#' 
#' @param threshold_freq_max frequency threshold to keep or exclude SNPs. value between 0 and 0.50 as PLINK default. The number of SNPs must match provided in the plink BIM file. The provided threshold will include value provided. For example, if threshold is 0.01, SNPs filtered will include MAF of 0.01 and above. Default no SNPs are filtered.
#'
#' @param threshold_freq_min frequency threshold to keep or exclude SNPs. value between 0 and 0.50 as PLINK default. The number of SNPs must match provided in the plink BIM file. The provided threshold will include value provided. For example, if threshold is 0.01, SNPs filtered will include MAF of 0.01 and above. Default no SNPs are filtered.
#'
#' @param genecoord_dt a dataframe for gene boundaries with CHR START END GENE as column names. Where CHR should be integer 1-22. START and END column should be integer. GENE column contains gene names
#'
#' @return a data frame that contains SNP index found within the gene boundaries and plink bim file, with SNP, base pair and gene information. NULL is returned if no values are found intersecting.
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

    ##index_pergenelist<- with(genecoord_dt[bimsnps_dt,.(.I,GENE),on=.(START<=V4,END>= V4),nomatch=0],split(I,GENE))     ##unique_snpindex<-unique(match(snp_withingenes$V2,bimsnps_dt$V2)) # provide only unique # V2 in the bim file is SNP/rsid

    matched_snps_withingene<-genecoord_dt[bimsnps_dt,.(V2, V4,GENE),on=.(START<=V4,END>=V4),nomatch=0] 
    
    matched_snps_withingene$index<-chmatch(matched_snps_withingene[,V2],bimsnps_dt[,V2])
    colnames(matched_snps_withingene ) <-c("SNP","BP","GENE","index")
    
    if( (is.null(freq_file_plink)==FALSE) & (is.null(threshold_freq_max)==FALSE) ){
        
if(is.null(threshold_freq_max)){
threshold_freq_max=0.50
} else{
        threshold_freq_max<-as.numeric(threshold_freq_max)

}

if(is.null(threshold_freq_min)){
threshold_freq_min=0
} else{
        threshold_freq_min<-as.numeric(threshold_freq_min)

}

        if(threshold_freq_max >0.5 ){
            
            cat("Issue with provided frequency value maximum. It can only take until 0.5\n")
            stop("Exiting code")
        }
        ##if ends for frequency threshold
        
        freq_plink_dt<-data.table::fread(freq_file_plink,header=TRUE) ##read data 
        freq_plink_dt_filtered<-freq_plink_dt[(MAF>=threshold_freq_min & MAF<=threshold_freq_max ),] ## filter based on the threshold provided
        
        filtered_maf_matched<-freq_plink_dt_filtered[matched_snps_withingene,c("SNP","BP","GENE","index"), on=c("SNP"),nomatch=0] ## do a join on MAF and already index SNP, gene BP position

         ## Do check on row count
        if(nrow(filtered_maf_matched)==0){
            cat("No intersection of SNPs found between Gene and SNP within the BIM file after MAF filtering. No index\n")
            return (NULL)
        }
        else{
            return(filtered_maf_matched)
        }

        ##check ends for row count
    }

    ## If both are not null then return already jointed, filtered data
    
    else{
        if(nrow(matched_snps_withingene)==0){
            cat("No intersection of SNPs found between Gene and SNP within the BIM file. No index\n")
            return (NULL)
        }
        else{
            return(matched_snps_withingene)
        }
    }
    ##else ends for no NULL in freq and threshold
    ########
}
## function ends added on Jan 10 2021


