library(dplyr)
library(data.table)
library(RNOmni)

setwd("github")

# 3 metabolite transformations #########

# log10
# log
# rank based inverse normal transformation (INT)  https://cran.r-project.org/web/packages/RNOmni/vignettes/RNOmni.html

psam = fread("GCAT.psam") # reading .psam PLINK file

# the columns 14:204 are the metabolites
# column 3 SEX
# columns 4:7 PCAs
# column 9 AGE
# column 10 Chylomicrons
# column 11 NMRr
# column 12 LCms
# column 13 GCms


platform_trait = NULL

for(i in 14:204){
  
  aux = log10(psam[,i])
  aux2 = log(psam[,i])

  aux[aux==-Inf] = NA
  aux2[aux2==-Inf] = NA
  
  if(platforms[platforms$Fenotype %in% colnames(psam)[i],1] == "NMRr"){platform_trait=11}
  if(platforms[platforms$Fenotype %in% colnames(psam)[i],1] == "LCms"){platform_trait=12}
  if(platforms[platforms$Fenotype %in% colnames(psam)[i],1] == "GCms"){platform_trait=13}
  
  my_lm = lm(psam[,i] ~ factor(psam[,3]) + psam[,4] + psam[,5] + psam[,6] + psam[,7] +
                        psam[,9] + psam[,10] + factor(psam[,platform_trait]))
  
  my_residuals = rankNorm(residuals(my_lm))
  
  rank_inv_normal = rep(NA,4988) # use a vector to fill the residuals since there are NA in the residuals
  
  rank_inv_normal[as.numeric(names(my_residuals))] = my_residuals
  
  aux3 = rank_inv_normal
  
  psam = cbind(psam,aux,aux2,aux3)
  
  colnames(psam)[ncol(psam)] = paste0(colnames(psam)[i],"_rank_inv_norm")
  colnames(psam)[ncol(psam)-1] = paste0(colnames(psam)[i],"_logn")
  colnames(psam)[ncol(psam)-2] = paste0(colnames(psam)[i],"_log10")
  
  print(i)
  
}

psam %>% fwrite("GCAT_transformed_metabolites.psam",row.names = F,quote = F,sep = " ")
