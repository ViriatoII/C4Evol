#install.packages("ape")
#install.packages("caper")
#install.packages("sommer")


library("sommer")
library("ape")
library("caper")

###### Phylo, ortho and SNP matrices  #####################
dir="~/mounts/hilbert/project/qggp/C34_PS/experiments/annotation/orthofinder/abinitio_ann/pre_blasts/OrthoFinder/Results_msa_final/" #/OrthoFinder/Results_msa_newGG/"

# Create Covariance table from Phylo tree
tree <- ape::read.tree(paste(dir,"Species_Tree/SpeciesTree_rooted.txt", sep = "")) 
phyl_covar <- ape::vcv.phylo(tree)
#write.csv(x = phyl_covar, file = "01_C4breed/phylo_covar.csv", sep = ",")

# Orthogroups (created by orthofinder, separating ancient paralogs) -> NEEDS PYTHON EDITING: /annotation/orthofinder/03_orthotable_transform.py
ogroups = read.csv("mounts/hilbert/project/qggp/C34_PS/experiments/annotation/orthofinder/orthogroups.csv")


###### SNP matrix #####

SNPs = data.frame(read.csv("mounts/hilbert/project/qggp/C34_PS/experiments/sebastian/GLDP1/GLDP1_SNP_table_new2.csv", stringsAsFactors=FALSE, sep = ",", quote = ""))
SNPs <- data.frame(SNPs) # stringsAsFactors=FALSE)

rownames(SNPs) <- SNPs$POS

SNPs <- SNPs[,c(4:13,21:25)]
SNPs <- t(SNPs)  # transpose
          # replaced b_repanda fy b_gravinae1 TEMPORARILY
rownames(SNPs) <- c("b_gravinae1", "b_tournefortii",          
  "d_viminea", "e_sativa", "b_rapa", "h_incana",
  "d_erucoides", "b_oleracea" ,"d_tenuisiliqua" ,
  "m_moricandioides", "m_nitens", "m_suffruticosa",
  "m_arvensis" , "d_tenuifolia", "g_gynandra")

SNPs <- SNPs[sort(rownames(SNPs)),]

###### CCP values #############

species = c("b_gravinae1","b_juncea", "b_napus",
            "b_nigra","b_oleracea", "b_rapa",
            "b_repanda", "b_tournefortii", "d_erucoides",
            "d_tenuifolia", "d_tenuisiliqua",
            "d_viminea","e_sativa", #"g_gynandra",
            "h_incana","m_arvensis","m_moricandioides",
            "m_nitens","m_suffruticosa","r_sativus")

ccp = c(37.217, 46.286, 46.155, 47.32, 50.516, 49.674, 55.603, 47.461,
        30.3, 12.411, 49.376, 51.083, 51.819, # 4.278,
        38.856, 24.425, 51.587, 21.043, 24.873, 45.187)

DT2 <- data.frame(species=species, ccp= ccp)

  # filter DT2 and phyl covar matrix to have exactly the same species
#species <- rownames(phyl_covar)[rownames(phyl_covar) %in% species ]
#DT2 <- DT2[DT2$species %in% species,] 
#rownames(DT2) <- NA

phyl_covar <- phyl_covar[species,species]

###### Ortho Model #############

# Table for modeling. Filter to only have desired species
model_table <- t(ogroups[,species])
colnames(model_table) <- ogroups[,"HOG"]    # REMOVE the 1:6

model_table <- as.data.frame(model_table)# , stringsAsFactors=FALSE)
#model_table <- model_table[1:15,] # fix stupid bug

# Eliminate single gene families
model_table= model_table[colSums(model_table)>1]


###### GWAS copy number #####################
result <- GWAS(ccp ~ 1, 
               random =~vs(species, Gu=phyl_covar),
               rcov=~units, gTerm = "u:species",
               M = model_table, data = DT2 )

# store the most relevant family names
#relevant <- names(which((10^(-result$scores[2,]) < 0.05))) # exp(2.5)

p_vals <- round(10^(-result$scores[2,((10^(-result$scores[2,]) < 0.05))]), 4)
relevant <- names(p_vals)
write.csv(data.frame("p_val"=p_vals, "HOG"=relevant), file = "interesting_copyNUM.csv", quote = FALSE, row.names = FALSE)

# Filter and save relevant families
min_sps <- c("d_tenuifolia", "m_arvensis", "m_nitens" ,  "m_suffruticosa")
med_sps <- c("d_erucoides", "b_gravinae1", "h_incana")
C3_sps <- species [species %in% c(min_sps, med_sps)==FALSE]

out_copyNUM <- data.frame("HOG" = 0, "p-vals" = 0 ,"C3-C4_avg_genes"=0, "C3-C4_sd_genes"=0,
                          "C3-C4*_avg_genes"=0, "C3-C4*_sd_genes"=0,
                          "C3_avg_genes"=0, "C3_sd_genes"=0, check.names = F)

for (rel in relevant){
  #get species present in OG. 
  present <- rownames(model_table[which(model_table[rel]>0),])
  
  # Ignore OGs related to high ccp values
  if (mean(DT2$ccp[DT2$species %in% present])< 40){
    
    ##  Obtain gene number per species
    DT3 <- DT2["ccp" ]
    rownames(DT3) <- DT2$species
    DT3$copy_number = 0
    DT3[which(model_table[rel]>0,), "copy_number"] <- model_table[which(model_table[rel]>0),rel]
    
    # Discriminate gene number per group
    min_genes <- DT3[min_sps,]$copy_number
    med_genes <- DT3[med_sps,]$copy_number
    C3_genes  <- DT3[C3_sps,]$copy_number
    
    if (sum(min_genes>0, med_genes>0)>1) {
      out_copyNUM[rel,c("C3-C4_avg_genes", "C3-C4_sd_genes")] <- c(round(mean(min_genes),2), round(sd(min_genes),2))
      out_copyNUM[rel,c("C3-C4*_avg_genes", "C3-C4*_sd_genes")] <- c(round(mean(med_genes),2), round(sd(med_genes),2))
      out_copyNUM[rel,c("C3_avg_genes", "C3_sd_genes")] <- c(round(mean(C3_genes),2), round(sd(C3_genes),2))
      
    }
  }
}

out_copyNUM["p-vals"] <- p_vals[row.names(out_copyNUM)]
out_copyNUM["HOG"] <- row.names(out_copyNUM)
out_copyNUM <- out_copyNUM[-1,]
write.csv(out_copyNUM, file = "interesting_copyNUM.csv", quote = FALSE, row.names = FALSE)




## Copy number Plots (ccp values) ############

plot(10^(-result$scores[2,][1:100]))
# it's the negative logarithm of the p-value. The higher the better (smaller the p-value is)
#10^-1.33  


library("ggplot2")
# Relevant GWAS scores
qplot(seq_along(result$scores[2,relevant]), result$scores[2,relevant]) +
  geom_point(colour ="navyblue", size= 1.2) +
        ggtitle("GWAS scores per gene family") + 
        theme(axis.text = element_text(size = rel(1.4), face = "bold"),
              plot.title = element_text(hjust = 0.5, face = "bold", size = rel(2)),
              axis.title = element_text(size=rel(1.8))) +
      xlab("Gene families") + ylab("-log10(p-val)")

# For copy number table
for (rel in relevant[1:10]){
  
  #get species present in OG. 
  present <- rownames(model_table[which(model_table[rel]>0),])
  
  # Ignore OGs related to high ccp values
  #print(mean(DT2$ccp[DT2$species %in% present]))
  if (mean(DT2$ccp[DT2$species %in% present])< 30){
    
    DT3 <- DT2 
    DT3$copy_number = 0
    
    DT3[which(model_table[rel]>0,), "copy_number"] = model_table[which(model_table[rel]>0),rel]
    #DT3$copy_number = factor(DT3$copy_number)
    
    #DT3 <- DT3[-14,] # remove Cleome
    
    print(ggplot(DT3, aes(reorder(species, ccp), ccp,colour=copy_number)) +
            geom_point(aes(fill=copy_number), 
                       colour="black",pch=21, size=rel(7)) +
            ylab("CCP") +
            theme(axis.title = element_text(size=rel(1.8)),
                  axis.text.x = element_text(size = rel(1.6), angle = 50, hjust = 1, face = "bold"),
                  axis.text.y = element_text(size = rel(1.4), hjust = 1, face = "bold"),
                  plot.title = element_text(size = rel(2), hjust = 0.5, face = "bold")) +
            labs(color = "Present in\ngene family", size = rel(1.5)) + 
            xlab("Species") + 
            scale_color_gradient(low = "gray10", high="steelblue1") +
            ggtitle(rel)
    )
  }
}

## Save plots of specific HOGs selected by Urte #####

selected <- c("N0.HOG0003237","N0.HOG0025962",
              "N0.HOG0028581","N0.HOG0034511",
              "N0.HOG0034527","N0.HOG0021276",
              "N0.HOG0026267")

for (rel in selected){
  
  #get species present in OG. 
  present <- rownames(model_table[which(model_table[rel]>0),])
  
  # Ignore OGs related to high ccp values
    
  DT3 <- DT2 
  DT3$copy_number = 0
    
  DT3[which(model_table[rel]>0,), "copy_number"] = model_table[which(model_table[rel]>0),rel]
  
  ploti <- ggplot(DT3, aes(reorder(species, ccp), ccp,colour=copy_number)) +
          geom_point(aes(fill=copy_number),
                     colour="black",pch=21, size=rel(7)) +
            ylab("CCP") +
          theme(axis.title = element_text(size=rel(1.8)),
                axis.text.x = element_text(size = rel(1), angle = 50, hjust = 1, face = "bold"),
                axis.text.y = element_text(size = rel(1.4), hjust = 1, face = "bold"),
                plot.title = element_text(size = rel(2), hjust = 0.5, face = "bold")) +
          labs(color = "Present in\ngene family", size = rel(1.5)) + 
          xlab("Species") + 
          scale_color_gradient(low = "gray10", high="steelblue1") +
          ggtitle(rel)
  
  ggsave(ploti, filename = rel, device = "png")
  
        
}


###### GWAS Presence/absence ###########

################ Presence-Absence table #
PA_table <- model_table
PA_table[model_table > 1] = 1 

result_PA <- GWAS(ccp ~ 1, 
               random =~vs(species, Gu=phyl_covar),
               rcov=~units, gTerm = "u:species",
               M = PA_table, data = DT2 )


# store the most relevant family names
#relevantPA <- names(which((10^(-result_PA$scores[2,]) < 0.05))) # exp(2.5)

p_vals_PA <- round(10^(-result_PA$scores[2,((10^(-result_PA$scores[2,]) < 0.05))]), 4)
relevantPA <- names(p_vals_PA)
write.csv(data.frame("p_val"=p_vals_PA, "HOG"=relevantPA), file = "interesting_PA.csv", quote = FALSE, row.names = FALSE)


## Filter and save most relevant families
out_PA <- data.frame("HOG" = 0, "p-vals" = 0 ,"C3-C4_avg_genes"=0, "C3-C4_sd_genes"=0,
                          "C3-C4*_avg_genes"=0, "C3-C4*_sd_genes"=0,
                          "C3_avg_genes"=0, "C3_sd_genes"=0, check.names = F)

for (rel in relevantPA){
  #get species present in OG. 
  present <- rownames(model_table[which(model_table[rel]>0),])
  
  # Ignore OGs related to high ccp values
  if (mean(DT2$ccp[DT2$species %in% present])< 40){
    
    ##  Obtain gene number per species
    DT3 <- DT2["ccp" ]
    rownames(DT3) <- DT2$species
    DT3$copy_number = 0
    DT3[which(model_table[rel]>0,), "copy_number"] <- model_table[which(model_table[rel]>0),rel]
    
    # Discriminate gene number per group
    min_genes <- DT3[min_sps,]$copy_number
    med_genes <- DT3[med_sps,]$copy_number
    C3_genes  <- DT3[C3_sps,]$copy_number
    
    if (sum(min_genes>0, med_genes>0)>1) {
      out_PA[rel,c("C3-C4_avg_genes", "C3-C4_sd_genes")] <- c(round(mean(min_genes),2), round(sd(min_genes),2))
      out_PA[rel,c("C3-C4*_avg_genes", "C3-C4*_sd_genes")] <- c(round(mean(med_genes),2), round(sd(med_genes),2))
      out_PA[rel,c("C3_avg_genes", "C3_sd_genes")] <- c(round(mean(C3_genes),2), round(sd(C3_genes),2))
      
    }
  }
}

out_PA["p-vals"] <- p_vals_PA[row.names(out_PA)]
out_PA <- out_PA[-1,]
out_PA["HOG"] <- row.names(out_PA)
write.csv(out_PA, file = "interesting_PA.csv", quote = F, row.names = F)


## PA plots (ccp values) #######

require(ggplot2)
# For P/Absense table
for (rel in relevantPA[1:10]){
  
  #get species present in OG
  present <- rownames(PA_table[rel]>0)
  
  #print(mean(DT2$ccp[DT2$species %in% present]))
  
  # Ignore OGs related to high ccps
  if (mean(DT2$ccp[DT2$species %in% present])< 42){
  
  DT3 <- DT2 
  DT3$pres = "NOT"
  
  DT3[which(model_table[rel]>0), "pres"] = "PRES"
  
  #DT3 <- DT3[-14,] # remove Cleome
  
  print(ggplot(DT3, aes(reorder(species, ccp), ccp,colour=pres)) +
    geom_point(size = rel(7)) + 
    ylab("CCP") +
    theme(axis.title = element_text(size=rel(1.8)),
          axis.text.x = element_text(size = rel(1.6), angle = 50, hjust = 1, face = "bold"),
          axis.text.y = element_text(size = rel(1.4), hjust = 1, face = "bold"),
          plot.title = element_text(size = rel(2), hjust = 0.5, face = "bold")) +
    labs(color = "Present in\ngene family", size = rel(1.5)) + 
    xlab("Species") +
    ggtitle(rel)
  )
  
  }
}


## Save plots of specific HOGs selected by Urte #####
selected <- c("N0.HOG0030986","N0.HOG0025956","N0.HOG0021770",
              "N0.HOG0022963","N0.HOG0024157","N0.HOG0024382")

for (rel in selected){
  
  #get species present in OG. 
  present <- rownames(PA_table[which(PA_table[rel]>0),])
  
  # Ignore OGs related to high ccp values
  
  DT3 <- DT2 
  DT3$pres = "NOT"
  
  DT3[which(model_table[rel]>0), "pres"] = "PRES"
  
  #DT3 <- DT3[-14,] # remove Cleome
  
  ploti <- ggplot(DT3, aes(reorder(species, ccp), ccp,colour=pres)) +
          geom_point(size = rel(7)) + 
          ylab("CCP") +
          theme(axis.title = element_text(size=rel(1.8)),
                axis.text.x = element_text(size = rel(1), angle = 50, hjust = 1, face = "bold"),
                axis.text.y = element_text(size = rel(1.4), hjust = 1, face = "bold"),
                plot.title = element_text(size = rel(2), hjust = 0.5, face = "bold")) +
          labs(color = "Present in\ngene family", size = rel(1.5)) + 
          xlab("Species") +
          ggtitle(rel,)
  ggsave(ploti, filename = rel, device = "png")
}



###### SNP model #############

model_table[1:300,1:5]
SNPs[1:4,1:5]


# remove positions with more than 5 NAs?
SNPs <- SNPs[,which(colSums(is.na(SNPs))<5 )]


result <- GWAS(ccp ~ 1, 
               random =~vs(species, Gu=phyl_covar),
               rcov=~units, gTerm = "u:species",
               M = SNPs, data = DT2 )

result$scores[4,]  #[4] #[3,result$scores[3,]>20]
result$scores[3,]>20


#### Others #################

# Only 98 families have single protein
sum(rowSums(ogroups[,4:length(colnames(ogroups))])<2)

# Plot this
#hist(rowSums(ogroups[,4:length(colnames(ogroups))]), 500, 
#  xlim = range(0,100), main = "Number of genes per Orthogroup",
#  xlab = "Number of Proteins in OG (Max is 1400)")

sp_pres <- ogroups[,4:length(colnames(ogroups))]
sp_pres[sp_pres>1]=1

