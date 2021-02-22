#R PROJECT FOR DR OGHBAIE

#Author: Laia Fernández Calvo

library(ggplot2)


#Read in file
proteinGroups <- read.delim("~/Documents/GitHub/Mehrnoosh_Rproject/proteinGroups.txt")

names(proteinGroups)


####
#3.Remove every record that seem to be a Potential Contaminant or Reverse sequence from the table.
####

levels(proteinGroups$Reverse)
levels(proteinGroups$Potential.contaminant)


cleandata <- subset(proteinGroups, Reverse != "+")
cleandata <- subset(proteinGroups, Potential.contaminant != "+")



####
#4. Separate different tissue LFQ intensity columns. Don’t forget to include gene.name and protein IDs.
#There are three different tissues [144/159/163], they could be either tumor or normal [T/N].
#The name format is like [144/159/163] [T/N]_[ORF1/IgG]_[1-10][a/b/c].
#Different antibodies were used for immuno-precipitation [ORF1/IgG]. 
#The number following by the antibody refer to specific condition (combination of tissue and antibody).
#In total we have 10 conditions here. The last lower alphabet is the number of replicates.
#A comparison could be either (159T_ORF1 against 159T_IgG) or (159T_ORF1 against 159N_ORF1).
####


t144 <- subset(cleandata, select = c(Protein.IDs, Gene.names, grep("LFQ.intensity.144", names(cleandata))))
t159 <- subset(cleandata, select = c(Protein.IDs, Gene.names, grep("LFQ.intensity.159", names(cleandata))))
t163 <- subset(cleandata, select = c(Protein.IDs, Gene.names, grep("LFQ.intensity.163", names(cleandata))))



####
#5. Log2 transform LFQ intensity.
####

#"Zero values are missing ( the intensity doesn’t reach a threshold ) or they just didn’t exist." 

#convert zero to NA 

t144[t144 == 0] <- NA
t159[t159 == 0] <- NA
t163[t163 == 0] <- NA

#"If only one replicate is non-zero or all replicates were zero, it didn’t exist."

#if within a condition, only one replicate is not NA - convert that to NA

changetoNA <- function(x){
  xt=t(x)
  for (i in 1:nrow(x)){
    NAcount=0
    for(j in 1:3){
      if(is.na(xt[j,i])) NAcount = NAcount+1
    }
    if(NAcount==2) {
      for(j in 1:3){
        xt[j,i]=NA
      }
    }
  }
  return(t(xt))
}


t144[3:5] <- changetoNA(t144[3:5])
t144[6:8] <- changetoNA(t144[6:8])
t144[9:11] <- changetoNA(t144[9:11])
t144[12:14] <- changetoNA(t144[12:14])


t159[3:5] <- changetoNA(t159[3:5])
t159[6:8] <- changetoNA(t159[6:8])
t159[9:11] <- changetoNA(t159[9:11])


t163[3:5] <- changetoNA(t163[3:5])
t163[6:8] <- changetoNA(t163[6:8])
t163[9:11] <- changetoNA(t163[9:11])


logt144 <- t144
logt144[, 3:14] <- log2(t144[3:14])

logt159 <- t159
logt159[, 3:11] <- log2(t159[3:11])

logt163 <- t163
logt163[, 3:11] <- log2(t163[3:11])



####
#6. Calculate the average log2 LFQ intensities in each condition.
####

#for tissue 144
logt144$mean10 <- rowMeans(logt144[,c(3,4,5)], na.rm = TRUE)
logt144$mean2 <- rowMeans(logt144[,c(6,7,8)], na.rm = TRUE)
logt144$mean1 <- rowMeans(logt144[,c(9,10,11)], na.rm = TRUE)
logt144$mean9 <- rowMeans(logt144[,c(12,13,14)], na.rm = TRUE)

#for tissue 159
logt159$mean5 <- rowMeans(logt159[,c(3,4,5)], na.rm = TRUE)
logt159$mean4 <- rowMeans(logt159[,c(6,7,8)], na.rm = TRUE)
logt159$mean3 <- rowMeans(logt159[,c(9,10,11)], na.rm = TRUE)


#for tissue 163
logt163$mean8 <- rowMeans(logt163[,c(3,4,5)], na.rm = TRUE)
logt163$mean7 <- rowMeans(logt163[,c(6,7,8)], na.rm = TRUE)
logt163$mean6 <- rowMeans(logt163[,c(9,10,11)], na.rm = TRUE)



####
#7. Calculate log2foldchange by subtracting the avg log2 LFQ intensity of control from case.
####

#To compare:
# 163T_ORF1 vs. 163N_ORF1
# 163T_ORF1 vs. 163T_IgG
# 
# 159T_ORF1 vs. 159N_ORF1
# 159T_ORF1vs. 159T_IgG
# 
# And two version of 
# 144T_ORF1 vs. 144T_IgG


#for tissue 144
logt144$cond_9_10<- logt144$mean9 - logt144$mean10
logt144$cond_1_2<- logt144$mean1 - logt144$mean2

#for tissue 159
logt159$TvsN<- logt159$mean3 - logt159$mean5
logt159$ORFvsIgG<- logt159$mean3 - logt159$mean4

#for tissue 163
logt163$TvsN<- logt163$mean6 - logt163$mean8
logt163$ORFvsIgG<- logt163$mean6 - logt163$mean7


####
#8. Perform t.test() between every records in case and control and extract the p.value from it.
#You can use try() to continue your calculation whenever not enough replicates exist.
####

ttest <- function(x,y){
  p_values=vector()
  for(i in 1:nrow(x)){
    if (all(is.na(t(x[i,])) | all(is.na(t(y[i,])))))
      p_values[i] <- NA
    else
      try(
        p_values[i] <- t.test(t(x[i,]),t(y[i,]))$p.value
      )
  }
  return(p_values)

}


#for tissue 144
# 144T_ORF1 vs. 144T_IgG
logt144$p_val9_10 <- ttest(logt144[,3:5], logt144[,12:14])
logt144$p_val1_2 <-ttest(logt144[,6:8], logt144[,9:11])


#for tissue 159
logt159$p_valT_N <- ttest(logt159[,3:5], logt159[,12:14])
logt159$p_valORF_IgG <- ttest(logt159[,3:5], logt159[,6:8])

#for tissue 163
logt163$p_valT_N <- ttest(logt163[,3:5], logt163[,12:14])
logt163$p_valORF_IgG <- ttest(logt163[,6:8], logt163[,12:14])



####
#9. Calculate adjusted p.value from the p.value distribution.
####

#for tissue 144
# 144T_ORF1 vs. 144T_IgG
logt144$p_adjust_9_10 <- p.adjust(logt144$p_val9_10 , method = "BH")
logt144$p_adjust_1_2 <- p.adjust(logt144$p_val1_2 , method = "BH")


#for tissue 159
logt159$p_adjust_T_N <- p.adjust(logt159$p_valT_N , method = "BH")
logt159$p_adjust_ORF_IgG <- p.adjust(logt159$p_valORF_IgG , method = "BH")


#for tissue 163
logt163$p_adjust_T_N <- p.adjust(logt163$p_valT_N , method = "BH")
logt163$p_adjust_ORF_IgG <- p.adjust(logt163$p_valORF_IgG , method = "BH")



####
#10. Draw the volcano plot using -log10(adjusted p.value) on the y-axis and log2fold change on the x-axis.
####

#volcano plot function

#Decided to choose significant values with p.value < 0.05 and to put threshold for the fold change at
#-1/1 log2fold change. These could be changed if desired.


volplot <- function(dataframe, x, y, title, label){
  fn=paste(getwd(),"/",title,".png",sep="")
  volcano_plot <- ggplot(data=dataframe, aes(x=x, y=-log10(y), col= label)) +
    theme_bw()+
    geom_point(size = 1) + 
    geom_vline(xintercept=c(-1, 1), col="blue") +
    geom_hline(yintercept=-log10(0.05), col="green")+
    ylim(-1,4)+
    xlim(-7,10)+
    xlab("log2fold") +
    ylab("-log10(adj. p.value)") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.4, face="bold.italic", size=14),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.x = element_text(face="bold",size=10),
          axis.text.y = element_text(face="bold", size=10), 
          legend.position = "right")+
    scale_color_manual(values=c("downregulated" = "blue", "non_signif." = "grey", "upregulated" = "red"))
  ggsave(fn, plot=volcano_plot)
}


#Adding labels

#tissue 159
# add column of for expression
logt159$diffexpres <- "non_signif."
# if log2Foldchange > 1 and pvalue < 0.05, set as "upregulated"
logt159$diffexpres[logt159$ORFvsIgG > 1 & logt159$p_adjust_ORF_IgG < 0.05] <- "upregulated"
# if log2Foldchange < -1 and pvalue < 0.05, set as "downregulated"
logt159$diffexpres[logt159$ORFvsIgG < -1 & logt159$p_adjust_ORF_IgG < 0.05] <- "downregulated"

#same thing with the rest of data

logt159$diffexpres2 <- "non_signif."
logt159$diffexpres2[logt159$TvsN > 1 & logt159$p_adjust_T_N < 0.05] <- "upregulated"
logt159$diffexpres2[logt159$TvsN < -1 & logt159$p_adjust_T_N < 0.05] <- "downregulated"


#tissue 144
logt144$diffexpres <- "non_signif."
logt144$diffexpres[logt144$cond_9_10 > 1 & logt144$p_adjust_9_10 < 0.05] <- "upregulated"
logt144$diffexpres[logt144$cond_9_10 < -1 & logt144$p_adjust_9_10 < 0.05] <- "downregulated"

logt144$diffexpres2 <- "non_signif."
logt144$diffexpres2[logt144$cond_1_2 > 1 & logt144$p_adjust_1_2 < 0.05] <- "upregulated"
logt144$diffexpres2[logt144$cond_1_2 < -1 & logt144$p_adjust_1_2 < 0.05] <- "downregulated"


#tissue 163
logt163$diffexpres <- "non_signif."
logt163$diffexpres[logt163$ORFvsIgG > 1 & logt163$p_adjust_ORF_IgG < 0.05] <- "upregulated"
logt163$diffexpres[logt163$ORFvsIgG < -1 & logt163$p_adjust_ORF_IgG < 0.05] <- "downregulated"

logt163$diffexpres2 <- "non_signif."
logt163$diffexpres2[logt163$TvsN > 1 & logt163$p_adjust_T_N < 0.05] <- "upregulated"
logt163$diffexpres2[logt163$TvsN < -1 & logt163$p_adjust_T_N < 0.05] <- "downregulated"


##
#Volcano plots
##
volplot(logt159, logt159$ORFvsIgG, logt159$p_adjust_ORF_IgG, "tissue_159_ORF_vs_IgG", logt159$diffexpres)
volplot(logt159, logt159$TvsN, logt159$p_adjust_T_N, "tissue_159_Tumor_vs_Normal_(ORF)", logt159$diffexpres2)

volplot(logt163, logt163$ORFvsIgG, logt163$p_adjust_ORF_IgG, "tissue_163_ORF_vs_IgG", logt163$diffexpres)
volplot(logt163, logt163$TvsN, logt163$p_adjust_T_N, "tissue_163_Tumor_vs_Normal_(ORF)", logt163$diffexpres2)

volplot(logt144, logt144$cond_9_10, logt144$p_adjust_9_10, "tissue_144_ORF_vs_IgG_(cond 9vs10)", logt144$diffexpres)
volplot(logt144, logt144$cond_1_2, logt144$p_adjust_1_2, "tissue_144_ORF_vs_IgG_(cond 1vs2)", logt144$diffexpres2)



