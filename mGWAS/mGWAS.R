library(data.table)
library(dplyr)

j <- Sys.getenv("SGE_TASK_ID")
j <- as.numeric(j)

psam = fread("GCAT_transformed_metabolites.psam")

platforms = fread("platforms.csv")

platform_trait = NULL

trait = colnames(psam)[j]

trait = gsub("_log10","",trait)
trait = gsub("_logn","",trait)
trait = gsub("_rank_inv_norm","",trait)

if(as.character(platforms[platforms$Fenotype %in% trait,1]) == "NMRr"){platform_trait="NMRr"}
if(as.character(platforms[platforms$Fenotype %in% trait,1]) == "LCms"){platform_trait="LCms"}
if(as.character(platforms[platforms$Fenotype %in% trait,1]) == "GCms"){platform_trait="GCms"}

## mGWAS for log10 transformed traits

if(length(grep("_log10",colnames(psam)[j]))==1){

system(paste0("plink2 ",
              " --psam GCAT_transformed_metabolites.psam ",
              " --pvar GCAT_imputed.pvar ",
              " --pgen GCAT_imputed.pgen ",
              " --glm cols=chrom,pos,ref,alt,a1count,nobs,a1freq,test,orbeta,se,ci,tz,p hide-covar ",
              " --chr 1-22 ",
              " --threads 4 ",
	            " --memory 15000 ",
              " --ci 0.95 ",
              " --pheno-name ",colnames(psam)[j],
              " --covar-name GENDER,AGE,PC1,PC2,PC3,PC4,Chylomicrons,",platform_trait,
              " --covar-variance-standardize ",
              " --out outputs/log10/",colnames(psam)[j]))


system(paste0("plink2 ",
              " --psam GCAT_transformed_metabolites.psam ",
	            " --pvar GCAT_imputed.pvar ",
              " --pgen GCAT_imputed.pgen ",
              " --glm cols=chrom,pos,ref,alt,a1count,nobs,a1freq,test,orbeta,se,ci,tz,p hide-covar ",
              " --chr 23 ",
              " --threads 4 ",
              " --memory 15000 ",
              " --ci 0.95 ",
	            " --xchr-model 1 ",
              " --pheno-name ",colnames(psam)[j],
              " --covar-name AGE,PC1,PC2,PC3,PC4,Chylomicrons,",platform_trait,
              " --covar-variance-standardize ",
              " --out outputs/log10/",colnames(psam)[j],"_chrX"))

system(paste0("gzip outputs/log10/",colnames(psam)[j],".",colnames(psam)[j],".glm.linear"))
system(paste0("gzip outputs/log10/",colnames(psam)[j],"_chrX.",colnames(psam)[j],".glm.linear"))


}


## mGWAS for log transformed traits


if(length(grep("_logn",colnames(psam)[j]))==1){

system(paste0(" mGWAS_PLINK/plink2 ",
              " --psam  outputs/GCAT.psam ",
              " --pvar GCAT_imputed.pvar ",
              " --pgen GCAT_imputed.pgen ",
              " --glm cols=chrom,pos,ref,alt,a1count,nobs,a1freq,test,orbeta,se,ci,tz,p hide-covar ",
              " --chr 1-22 ",
              " --threads 4 ",
              " --ci 0.95 ",
              " --pheno-name ",colnames(psam)[j],
              " --covar-name GENDER,AGE,PC1,PC2,PC3,PC4,Chylomicrons,",platform_trait,
              " --covar-variance-standardize ",
              " --out outputs/logn/",colnames(psam)[j]))


system(paste0(" mGWAS_PLINK/plink2 ",
              " --psam  outputs/GCAT.psam ",
              " --pvar  GCAT_imputed.pvar ",
              " --pgen  GCAT_imputed.pgen ",
              " --glm cols=chrom,pos,ref,alt,a1count,nobs,a1freq,test,orbeta,se,ci,tz,p hide-covar ",
              " --chr 23 ",
              " --threads 4 ",
              " --ci 0.95 ",
              " --xchr-model 1 ",
              " --pheno-name ",colnames(psam)[j],
              " --covar-name AGE,PC1,PC2,PC3,PC4,Chylomicrons,",platform_trait,
              " --covar-variance-standardize ",
              " --out outputs/logn/",colnames(psam)[j],"_chrX"))

system(paste0("gzip outputs/logn/",colnames(psam)[j],".",colnames(psam)[j],".glm.linear"))
system(paste0("gzip outputs/logn/",colnames(psam)[j],"_chrX.",colnames(psam)[j],".glm.linear"))


}


## mGWAS for rank inverse normal transformed traits


if(length(grep("rank_inv_norm",colnames(psam)[j]))==1){


system(paste0(" mGWAS_PLINK/plink2 ",
              " --psam outputs/GCAT.psam ",
              " --pvar GCAT_imputed.pvar ",
              " --pgen GCAT_imputed.pgen ",
              " --glm cols=chrom,pos,ref,alt,a1count,nobs,a1freq,test,orbeta,se,ci,tz,p hide-covar ",
              " --chr 1-22 ",
              " --threads 4 ",
              " --ci 0.95 ",
              " --pheno-name ",colnames(psam)[j],
              " --out outputs/rank_inv_norm/",colnames(psam)[j]))


system(paste0(" mGWAS_PLINK/plink2 ",
              " --psam outputs/GCAT.psam ",
              " --pvar GCAT_imputed.pvar ",
              " --pgen GCAT_imputed.pgen ",
              " --glm cols=chrom,pos,ref,alt,a1count,nobs,a1freq,test,orbeta,se,ci,tz,p hide-covar ",
              " --chr 23 ",
              " --ci 0.95 ",
              " --threads 4 ",
              " --xchr-model 1 ",
              " --pheno-name ",colnames(psam)[j],
              " --out  outputs/rank_inv_norm/",colnames(psam)[j],"_chrX"))

system(paste0("gzip outputs/rank_inv_norm/",colnames(psam)[j],".",colnames(psam)[j],".glm.linear"))
system(paste0("gzip outputs/rank_inv_norm/",colnames(psam)[j],"_chrX.",colnames(psam)[j],".glm.linear"))


}
