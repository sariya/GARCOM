
plink_count_snppos<-function(plink_file,genecoord_dt,snp_pos_dt,snp_index,individuals_index){
## added on 12 23 2020

.libPaths(c("/home/ss5505/libraries_R/R_LIB4.0",.libPaths()))
library(bigsnpr)
library(data.table)

plink_rds <- bigsnpr::snp_readBed2(plink_file,backingfile = tempfile(), ind.col=snp_index, ind.row=individuals_index)

     # Loading the data from backing files
     data_plink <- snp_attach(plink_rds)

plink_snp_information_dt<- as.data.table(data_plink$map)
plink_fam_information<-as.data.table(data_plink$fam)

plink_genotype_dt<-  as.data.table(data_plink$genotypes[] )
colnames(plink_genotype_dt) <-plink_snp_information_dt$marker.ID

#     str(data_plink)     # dim(G <- data_plink$genotypes)    #  G[1:8, 1:8]
cat("Genotypes have been loaded from plink file\n")

  snp_withingenes<- snp_pos_dt[genecoord_dt,c("SNP","GENE","START","END"),on=list(BP>=START,BP<=END),nomatch=0] # inner join ##https://stackoverflow.com/questions/63290994/foverlaps-data-table-error-ys-key-must-be-identical-to-the-columns-specified
    if(nrow(snp_withingenes) == 0){
        stop("No snps within any gene boundaries provided")	
    }

## we check if the SNPs that are found in plink map data overlap with inner-joined data 
    if(isFALSE( (any( plink_snp_information_dt$marker.ID %in% unique(snp_withingenes$SNP))))){
        stop("No SNPs overlapping between genetic data and SNP annotation with Gene boundaries")
    }

 ##get index of SNPs and store in a vector to be used while subsetting plink genotype data.table
snps_intersect_index<- base::match(intersect( plink_snp_information_dt$marker.ID ,unique(snp_withingenes$SNP)),plink_snp_information_dt$marker.ID) 

 ## subset plink data
plink_genotype_dt_subset<- plink_genotype_dt[,.SD,.SDcols=snps_intersect_index ] ## get columns that intersect 
    plink_genotype_dt_subset[, rowid := plink_fam_information$sample.ID ] ## assign row names as IIDs

    plink_genotype_dt_subset_transposed<-data.table::transpose(plink_genotype_dt_subset,keep.names = "SNP", make.names="rowid") ## tranpose data and have some fancy settings

   subsetsnps_genes_lefted_join <- snp_withingenes[plink_genotype_dt_subset_transposed,on="SNP",nomatch=0] 
    subsetsnps_genes_lefted_join[,c("START","END","SNP"):=NULL]  
    matrix_withallelecount_withinGene <-subsetsnps_genes_lefted_join[,lapply(.SD,sum,na.rm=TRUE),by=GENE] 
    matrix_withallelecount_withinGene<-matrix_withallelecount_withinGene[ rowSums(matrix_withallelecount_withinGene[,-c("GENE")]) > 0,]

    if(nrow(matrix_withallelecount_withinGene)>0){
        return(matrix_withallelecount_withinGene)
    }
    else{
        return(NULL) ## if nothing is left after gene >0 
    }


}
## function ends


snppos_dt<-data.table::fread("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/SOFTWARES/R_package_development/test_vcffiles/test_snppos.txt", header=TRUE)
genecoordpass_dt<-data.table::fread("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/SOFTWARES/R_package_development/test_vcffiles/sample_genes.txt", header=TRUE)
file_plink<-"/mnt/mfs/statgen/UKBiobank/data/exome_files/ukb23155_c22_b0_v1.bed"


output_test<-plink_count_snppos(plink_file=file_plink,genecoord_dt=genecoordpass_dt,snp_pos=snppos_dt,snp_index=c(1:4000))

output_testfiltered<-plink_count_snppos(plink_file=file_plink,genecoord_dt=genecoordpass_dt,snp_pos=snppos_dt,snp_index=c(1:4000),individuals_index=c(1:3000))

output_testfiltered<-plink_count_snppos(plink_file=file_plink,genecoord_dt=genecoordpass_dt,snp_pos=snppos_dt,snp_index=c(1:4000),individuals_index=c(1:3000))

output_testfiltered_complete<-plink_count_snppos(plink_file=file_plink,genecoord_dt=genecoordpass_dt,snp_pos=snppos_dt,snp_index=c(1:414980),individuals_index=c(1:3000))


