#Read in file
proteinGroups <- read.delim("~/Documents/GitHub/Mehrnoosh_Rproject/proteinGroups.txt")

names(proteinGroups)

#3.Remove every record that seem to be a Potential Contaminant or Reverse sequence from the table.

levels(proteinGroups$Reverse)
levels(proteinGroups$Potential.contaminant)


cleandata <- subset(proteinGroups, Reverse != "+")
cleandata <- subset(proteinGroups, Potential.contaminant != "+")


#4. Separate different tissue LFQ intensity columns. Don’t forget to include gene.name and protein IDs.
#There are three different tissues [144/159/163], they could be either tumor or normal [T/N].
#The name format is like [144/159/163] [T/N]_[ORF1/IgG]_[1-10][a/b/c].
#Different antibodies were used for immuno-precipitation [ORF1/IgG]. 
#The number following by the antibody refer to specific condition (combination of tissue and antibody).
#In total we have 10 conditions here. The last lower alphabet is the number of replicates.
#A comparison could be either (159T_ORF1 against 159T_IgG) or (159T_ORF1 against 159N_ORF1).

#Ask me questions if you have problem figuring the experiment design.


t144 <- subset(cleandata, select = c(Protein.IDs, Gene.names, grep("LFQ.intensity.144", names(cleandata))))
t159 <- subset(cleandata, select = c(Protein.IDs, Gene.names, grep("LFQ.intensity.159", names(cleandata))))
t163 <- subset(cleandata, select = c(Protein.IDs, Gene.names, grep("LFQ.intensity.163", names(cleandata))))


#5. Log2 transform LFQ intensity.

logt144 <- t144
logt144[, 3:14] <- log2(t144[3:14])

logt159 <- t159
logt159[, 3:14] <- log2(t159[3:11])

logt163 <- t163
logt163[, 3:14] <- log2(t163[3:11])


#6. Calculate the average log2 LFQ intensities in each condition.


View(proteinGroups$Intensity)
View(proteinGroups$Intensity.144T_IgG_10a)



#7. Calculate log2foldchange by subtracting the avg log2 LFQ intensity of control from case.





#8. Perform t.test() between every records in case and control and extract the p.value from it.
#You can use try() to continue your calculation whenever not enough replicates exist.






#9. Calculate adjusted p.value from the p.value distribution.








#10. Draw the volcano plot using -log10(adjusted p.value) on the y-axis and log2fold change on the x-axis.













