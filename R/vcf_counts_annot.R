## added on 08 28 2020

##
##SNP gene annotation based
##
##

vcf_counts_annot<-function(vcf_data,df_snpgene){

if( class(vcf_data)[1] !="vcfR" ){
print("vcfR class not found")
}

genotyped_extracted<-vcfR::extract.gt(vcf_data,element = "GT",as.numeric=TRUE,convertNA=TRUE) 

df_genotyped_extracted<-data.table::data.table(genotyped_extracted,keep.rownames=TRUE)
jointed_gene_VCFGT<-df_snpgene[df_genotyped_extracted, on=c(SNP="rn"), nomatch=0L]
jointed_gene_VCFGT<-jointed_gene_VCFGT[,SNP:=NULL] ### remove SNP cols



if(nrow(jointed_gene_VCFGT)==0){
        message("No SNPs match with the annotation")
        return(NULL)
}

jointed_gene_VCFGT<-jointed_gene_VCFGT[, lapply(.SD,as.numeric), by="GENE"] ## convert into numeric
jointed_genesSNP_filtered<-jointed_gene_VCFGT[,lapply(.SD,sum,na.rm=TRUE),by=GENE] # get sum within a gene
jointed_genesSNP_filtered<-jointed_genesSNP_filtered[rowSums(jointed_genesSNP_filtered[,-c("GENE")]) > 0,]

if(nrow(jointed_genesSNP_filtered) == 0){
stop("No genes with any reference genotyped in the VCF with supplied SNP gene annotation")
}
return(jointed_genesSNP_filtered)

}

vcf_outputannotated<-vcf_counts_annot(vcf,df_snpgene)

.libPaths(c( "/home/ss5505/libraries_R/R_LIB4.0",.libPaths()))

#Loading required namespace: adegenet

library(vcfR)
library(data.table)


##df_snpgene<- data.table(read.table("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/SOFTWARES/R_package_development/test_vcffiles/CHR22_snp.gene",header=TRUE))

##pubvcf <- vcfR::read.vcfR("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/SOFTWARES/R_package_development/test_vcffiles/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",verbose =TRUE)

##vcf <- vcfR::read.vcfR("/mnt/mfs/hgrcgrid/shared/GT_ADMIX/SOFTWARES/R_package_development/test_vcffiles/CHR22_ADC.vcf.gz", verbose=TRUE)




