
##################################################################################
#this function allows to add two columns to an ROH file (pos1, pos2)             #
#                                                                                #
#!!! ROH : an ROh file computed by plink                                         #                                       
#!!! submaps : the object submaps containing all the submaps and summarys        #
#                                                                                #
#*** return the new dataframe modified                                           #
##################################################################################

addCM <- function (ROH, submaps) {
  
  ROH$POS1_cM <- 0
  ROH$POS2_cM <- 0
  ROH$cM <- 0
  map <- submaps@bedmatrix@snps
  
  #add cM information
  SNP1 <- map$dist[ match(ROH$SNP1, map$id) ]
  SNP2 <- map$dist[ match(ROH$SNP2, map$id) ]
      
  ROH$POS1_cM <- SNP1
  ROH$POS2_cM <-  SNP2
  ROH$cM      <- SNP2 - SNP1
  
  ROH
  
}




