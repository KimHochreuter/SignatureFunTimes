################################################################################
##
##  WGS_PCAWG_96 (from the Pan-Cancer Analysis of Whole Genomes)
##
################################################################################


WGS_PCAWG_96 <- read_csv("Data/WGS_PCAWG.96.csv")

#Determine the number of patient within each cancer
MARTA <- WGS_PCAWG_96
MARTA <- rbind(colnames(MARTA), MARTA)
colnames(MARTA) <- c("Mutation type","Trinucleotide", 
                     gsub("\\-.*","", colnames(WGS_PCAWG_96)[3:length(WGS_PCAWG_96[1,])]))
table(colnames(MARTA))

#Extract the patient with liver cancer for further analysis
Liver <- cbind(WGS_PCAWG_96[,c(1,2)],
               WGS_PCAWG_96[,gsub("\\-.*","",colnames(WGS_PCAWG_96)) == "Liver"])


rownames(Liver) = paste(substr(Liver$Trinucleotide, start = 1, stop = 1),"[",Liver$`Mutation type`,"]", substr(Liver$Trinucleotide, start = 3, stop = 3), sep = "")
Liver = Liver[,3:ncol(Liver)]
save(Liver, file = "DATA/Liver326.RData")
