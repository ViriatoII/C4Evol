

### FOR BENJAMIN  ####

# We use a var-covar matrix. But:
   # Min and max of this are close to each other. 
   # The diagonal line is NOT one
   # Just looks different from a Kinship matrix 


#install.packages("ape")
#install.packages("caper")
#install.packages("sommer")
#install.packages("devtools")
#install.packages("phytools")
#install.packages("heatmap.plus")
#install.packages("RColorBrewer")


library(devtools)
#install_github("vqv/ggbiplot")
#devtools::install_github("slowkow/ggrepel")

library("sommer")
library("ape")
library("caper")
library("phytools")
library("ggplot2")
library("ggbiplot")
library("heatmap.plus")
library("RColorBrewer")


################ Phylo covar matrix and contrasts  #####################

#dir="~/mounts/hilbert/project/qggp/C34_PS/experiments/annotation/orthofinder/db_ann/pre_blasts0/OrthoFinder/Results_pre_blasts9_msa_final/" #/OrthoFinder/Results_msa_final/"
#dir="01_C4breed/orthofinder/db_ann/pre_blasts/OrthoFinder/Results_pre_blasts8_msa/"
#dir="01_C4breed/orthofinder/db_ann_filt/pre_blasts/OrthoFinder/Results_pre_blasts7_tetra_raxml/"

dir="C:/unix/01_C4breed/"

# Create Covariance table from Phylo tree
#tree <- ape::read.tree(paste(dir,"Species_Tree/SpeciesTree_rooted.txt", sep = "")) 

tmp = "C:/unix/01_C4breed/singleOGs.final.tree" #qggp/C34_PS/experiments/phylogeny/astral/single_OGs/singleOGs.final.tree"
tree <- ape::read.tree(tmp)
#tree$edge.length
tree <- compute.brlen(tree)
tree <- ape::root(tree, outgroup = "g_gynandra", resolve.root = TRUE)

phyl_covar  <- ape::vcv.phylo(tree)

heatmap(phyl_covar, labCol = colSums(ogroups[-1]), symm = T )

#write.csv(x = phyl_covar, file = "01_C4breed/phylo_covar.csv", sep = ",")

##################### Orthogroup matrix            ##################################
#(created by orthofinder, separating ancient paralogs) 
# NEEDS PYTHON EDITING: /annotation/orthofinder/03_orthotable_transform.py

ogroups = read.csv(paste(
            dir,  "ogroups.csv", sep ="")) #"Phylogenetic_Hierarchical_Orthogroups/ogroups.csv"

##################### Heatmap of gene clusters     #####

#install.packages("corrplot")
#library(corrplot)

sp_counts <- ogroups [, -1]
#rownames(sp_counts) <- ogroups$HOG
sp_counts[sp_counts>0] <- 1 

# Create correlation matrix and Heatmap  (filter sp_counts for species correlation, use ogroups for all genes)
min_sps <- 9
filtered <- sp_counts[rowSums(sp_counts)> min_sps,]  # not filtering now
filtered.cor = cor(as.matrix(filtered))

#png(file=paste("~/Desktop/meeting/images/C34/orthogroups/heatmap_HOGs_", min_sps, "Sps"), 
    #width = 1100, height = 714, res = 115)

heatmap(as.matrix(filtered.cor),  symm = T,
        labCol = colSums(filtered), keep.dendro = FALSE ,
        main = paste ("HOGs with >", min_sps, " Sps"))
#dev.off()


## Plot SP num distribution across HOGs
t <- data.frame(Sps=rowSums(sp_counts))
row.names(t) <- NULL

ggplot(t, aes(x=Sps))+ 
  geom_histogram(binwidth=1) + xlim(c(0,29)) +
  scale_x_continuous(breaks=1:max(t)) + 
  ggtitle("Distribution of \n species NUM in HOGs") +
  ylab("N° HOGs") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20)
  )

colSums(tmp)
for (i in 1:10){
  print(tmp[i,tmp[i,]>0])}

##################### Plot genes per HOG           ########################

# Get only HOGs with more than > species
sp_filt=9
HOG_filt <- rowSums(ogroups[,-c(1)]!=0)> sp_filt # filter for groups with more than 3 species

gens_per_sp <- sort(colSums(ogroups[HOG_filt, -c(1)]))  #1:4

#gens_per_sp <- sort(colSums(ogroups[-1]))  # ERASE LATER
sp_num <- length(gens_per_sp)

par(mar=c(9,5,3,1))

# Boxplot
boxplot(ogroups[HOG_filt,names(gens_per_sp)], 
        ylim=c(0,8), xlab =NA, names = NA, 
        ylab= "N° Genes per HOG", 
        main = paste("Genes per HOG per species \n (HOG>", sp_filt, "sps)")) 

axis(1,at=1:sp_num,labels=names(gens_per_sp), las =2)

# Normal plot
plot(gens_per_sp/sum(HOG_filt),xaxt="n", ylab = "Avg. genes per HOG", xlab = NA,  cex.lab = 1.4, cex.axis = 1.5)
axis(1,at=1:sp_num,labels=names(gens_per_sp), las =2, cex.axis = 1.1,)


# Stacked plot Single-Copy Genes + Copy_num genes (paralogs)
tmp <- ogroups[,-1] #[HOG_filt,-1]
t <- data.frame("unassigned" = c("443","365","387","623","1649","442","607","300","911","674","794","798","379","739","304","841","994","828","2165","757","648","619","575","805","551","3848","2887","4556","1291","586","646","2381"),
                "total" =colSums(tmp)-colSums(tmp==1),
                "single_copy" = colSums(tmp==1))


# Same thing but stacked plot is in %
t <- data.frame("total" =1-(colSums(tmp==1)/colSums(tmp)),
                "single_copy" = colSums(tmp==1)/colSums(tmp))

options(scipen=3)                

barplot(as.matrix(t(t)), las = 2)#,
        #main = 'Total Genes and Single Copy genes per species', las = 2)

##################### Correlation between geneNUM and BUSCO stats ###########

for_corr <- data.frame(row.names = c("a_alpina",	"a_arabicum", "a_thaliana",	"b_gravinae", "b_juncea", "b_napus", "b_nigra", "b_oleracea", "b_rapa", "b_repanda", "b_tournefortii","c_annua", "d_acris",	"d_erucoides", "d_harra", "d_muralis", "d_tenuifolia",  "d_tenuisiliqua",	"d_viminea","e_sativa",	"g_gynandra",	"h_incana",	"h_incana3", "m_arvensis",	"m_moricandioides",	"m_nitens",	"m_sinaica", "m_spinosa", "m_suffruticosa", "r_raphanistrum", "r_sativus", "s_alba"),
           dups = c(0.03	,0.03,	0.01, 0.09,	0.73, 0.75,	0.16,	0.17,	0.16,	0.15,	0.12, 0.12, 0.29,	0.20, 0.15, 0.88,	0.14, 0.14, 0.13, 0.20, 0.07, 0.15, 0.14, 0.27, 0.18, 0.14, 0.12, 0.15, 0.17, 0.13, 0.27, 0.18),
           L50 = c(4,311,3,7248,9,9,5,5,6,232,8,36,3714,38,13174,117,86,316,23,120,151,74,382,2402,1110,3576,1831,11820,121,6724,5,77),
           L90 = c(7,2884,5,28717,242,54,134,9,183,3059,21,132,36273,216,59884,1558,447,3557,93,630,564,1266,11058,8526,4268,13039,19129,80352,1239,29702,801,2170),
           compl = c(0.90,0.89,0.94,0.65,0.95,0.95,0.90,0.94,0.93,0.91,0.90,0.91,0.94,0.97,0.79,0.98,0.95,0.90,0.92,0.96,0.92,0.91,0.90,0.95,0.95,0.86,0.94,0.74,0.95,0.93,0.92,0.94))

# Add the avg. genes per species per HOG to for_corr table
for_corr$genes <- colMeans(ogroups[rownames(for_corr)])  #colSums
#for_corr <- for_corr[-c(1:3),]  # filter species out of cor analysis

View(for_corr)
plot(for_corr$genes, for_corr$dups)
plot(for_corr$genes, for_corr$L50)
plot(for_corr$genes, for_corr$compl)

cor.test(for_corr$genes, for_corr$dups)
cor.test(for_corr$genes, for_corr$L50)
cor.test(for_corr$genes, for_corr$compl)

for_corr$genes <- for_corr$genes / for_corr$dups
#for_corr$genes <- for_corr$genes / for_corr$compl

# Normalize gene counts in Ogroups!
for (sp in colnames(ogroups[-1])){
  ogroups[sp] <- ogroups[sp]/
    for_corr[sp,"dups"] / for_corr[sp,"L50"]  # EXPERIMENTAL
}

#View(for_corr)
##################### SNP matrix    #############################

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


##################### CCP values    ####################################

# Filter out these:
filt_out = c(2,3,11,15,16,24,25) # 1,19,31) -> always out

species = c(#"a_thaliana",
            "b_gravinae","b_juncea",  "b_napus",
            "b_nigra","b_oleracea", "b_rapa",
            "b_repanda", "b_tournefortii", "c_annua",
            "d_erucoides",   "d_muralis",
            "d_tenuifolia", "d_tenuisiliqua","d_viminea", 
            "d_acris", "d_harra",
            "e_sativa", #"g_gynandra",
            "h_incana", "h_incana3","m_arvensis","m_moricandioides",
            "m_nitens","m_suffruticosa", "m_sinaica", "m_spinosa", 
            "r_raphanistrum", "r_sativus",
            "s_alba")#, "s_arvensis" )

ccp = c(#60.0, 
        37.217, 46.286, 46.155,
        47.32, 50.516, 49.674, 
        55.603, 47.461, 55.33, # c_annua
        30.3,   35.04,  # d_muralis
        12.411, 49.376, # d_tenuisiliqua
        51.083,  55.57, 53.25, # d_acris d_harra
        51.819, # 4.278, # Gg
        50.50, 38.856,  24.425, # ma
        51.587, 21.043, # Mn
        24.873, 23.70, 17.75, # Mp
        54.16, 45.187,
        49.96) #, 54.60 )

species <- species[-filt_out]
ccp <- ccp[-filt_out]

DT2 <- data.frame(species=species, ccp= ccp)

  # filter DT2 and phyl covar matrix to have exactly the same species
#species <- rownames(phyl_covar)[rownames(phyl_covar) %in% species ]
#DT2 <- DT2[DT2$species %in% species,] 
#rownames(DT2) <- NA
phyl_covar <- phyl_covar[species,species]

## Erase this line
#for (i in species){ if (i %in% row.names(phyl_covar) == FALSE ){print(i)}}


##################### Model table (filtered)  ###################################

# Table for modeling. Filter to only have desired species
model_table   <-   t(ogroups[,species])
colnames(model_table) <- ogroups[,"HOG"]    # REMOVE HOG COLUMN

#model_table <- as.data.frame(model_table) #, stringsAsFactors=F)

# Eliminate HOGs with < 10 species
model_table = model_table[,colSums(model_table !=0)>9]

# Eliminate families with supermultiplicated genes in any species (Noticed in PCA analysis: B. olereaceae proteome is way off)
model_table = model_table[,-which(colSums(model_table>100)>0)]   # was 200 genes, now it's 100 (normalized value())

# Eliminate families with too many genes 
#model_table = model_table[,-which(colSums(model_table)>800)]   

##################### PCA (Principal Component Analysis) ###########

# Normalize with scale (-avg / sd) 
PCA1 <- prcomp(scale(t(model_table))) #species visual
PCA2 <- prcomp(model_table) #orthogroups variance

## Put main components in model table (Wrong?)
#model_table[c("PC1","PC2", "PC3"),] <- PCA$rotation[,c(1,2,3)]

# PCA plot   (ggbiplot doesn't work with so many)
#ggbiplot(PCA, center=T,scale.=T) #, labels=rownames(forPCA))

plot(PCA2$rotation[,c(1,2)], 
     xlab = paste("PC1 - Explained variance: ", 
                  summary(PCA2)$importance[,1][2]),
     ylab = paste("PC2 - Explained variance: ", 
                  summary(PCA2)$importance[,2][2]))

relevant <- names(sort(PCA2$rotation[,2], decreasing = T))

head(relevant)
View(model_table[,relevant[1:10]])

library("ggplot2")
library("ggrepel")
qplot(PCA1$rotation[,1],PCA1$rotation[,2], 
      label = names(PCA1$rotation[,2])) + 
  geom_point(colour ="navyblue", size= 1.2) +
  ggtitle("PCA of HOG table") + 
  theme(axis.text = element_text(size = rel(1.4), face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = rel(2)),
        axis.title = element_text(size=rel(1.8))) +
  xlab(paste("PC1 - Explained variance: ", 
             summary(PCA1)$importance[,1][2])) +
  ylab(paste("PC2 - Explained variance: ", 
             summary(PCA1)$importance[,2][2])) +
  geom_text_repel()



##################### MMER    CopyNUM MMER     #################

# GLDP1 is "N0.HOG0015932"

library(sommer)

p_vals <- c()

# APPLY MMER TO EVERY MARKER
for(i in 15001:ncol(model_table)){ 
  #browser() 
  # Prepare DT table
  DT3 <- cbind(DT2, model_table[,i], PCA1$rotation[,1], 
               PCA1$rotation[,2], PCA1$rotation[,3])
  colnames(DT3) <- c("species","ccp", "marker", "PC1", "PC2", "PC3")          
  DT3$dups <- for_corr[species, "dups"]
  
  # Apply Mixed Modelling 
  result <- mmer(ccp~1 + marker + PC1 + PC2,
                 random=~vs(species), #, Gu=phyl_covar),
                 rcov=~units, data=DT3)
  
  result2 <- mmer(ccp~1 + PC1 + PC2,
                  random=~vs(species), #, Gu=phyl_covar),
                  rcov=~units, data=DT3, verbose = F)
  
  u.hat <- result2$U$`u:species`$ccp
  u.hat <- u.hat[match(DT3$species,names(u.hat))]
  
  #plot(DT3$ccp, u.hat)
  
  # Anova 
  aov <- anova.mmer(result, result2) #type = 2)["marker",5]
  
  p_vals <- append(p_vals, 
                   as.numeric(strsplit(
                     paste(aov["mod2","PrChisq"]), " ")[[1]][1])
  )
}

head(p_vals)
which(p_vals<0.05)

plot(-log10(p_vals), xlab = "Orthogroup (HOG)",
     ylab = "-log10(p_val)") +
  abline(h=1.30103, col="red")


#colnames(model_table[,p_vals<0.05])

write(p_vals, file = "p_vals_copyNUM_tetras.txt.txt", ncolumns = 1)
relevant <- colnames(model_table[,p_vals<0.05])
write(relevant, file = "relevantHOGs_copyNUM_tetras.txt", ncolumns = 1)

############ Plots for RELEVANT HOGs ##############


# See HOGs that are flagged
View(model_table[,which(p_vals<0.05)])

sp <-  c("d_tenuifolia", "d_muralis", "d_viminea")
sp <-  c("h_incana3", "h_incana", "b_tournefortii")
sp <-  c("m_arvensis", "m_suffruticosa", "m_moricandioides")
sp <-  c("b_napus", "b_nigra", "d_acris")

barplot(model_table[sp,which(p_vals<0.05)], beside = T,
        col = c(3,7,2),cex.axis=2, xaxt="n")
legend(x = "topright",
       fill = c(3,7,2),
       legend = sp)


mycols = c()

for (i in ccp){
  if (i < 30){mycols <- c(mycols,"green")}
  else if (i <40){
    mycols <- c(mycols,"yellow")}
  else{mycols <- c(mycols,"grey")}}


# Boxplot
par(mar=c(8, 4, 2, 2) + 0.1)
boxplot(t(model_table[species,which(p_vals<0.05)]), 
        ylim=c(0,8), xlab =NA, names = NA, 
        ylab= "N° Genes per HOG", 
        col = mycols,
        main = "Gene distribution in flagged HOGs" ) 

axis(1,at=1:nrow(model_table),
     labels=row.names(model_table), 
     las =2)
par(op)


model_table[,which(p_vals<0.05)]


############ Copy number Plots (GWAS and ccp values) #######################

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
for (rel in relevant){
  
  #get species present in OG 
  present <- rownames(model_table[which(model_table[,rel]>0),])
  
  # Ignore OGs related to high ccp values
  #print(mean(DT2$ccp[DT2$species %in% present]))
  if (mean(DT2$ccp[DT2$species %in% present])< 35){
    
    print(rel)
    DT3 <- DT2 
    DT3$copy_number <-0
    
    DT3[which(model_table[,rel]>0,), "copy_number"] <-  model_table[which(model_table[,rel]>0),rel]

    
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

############ Save plots of specific HOGs selected by Urte #####

selected <- c("N0.HOG0014922","N0.HOG0015728","N0.HOG0018254",
              "N0.HOG0018651","N0.HOG0018903","N0.HOG0018997",
              "N0.HOG0019858")

for (rel in relevant){
  
  #get species present in OG. 
  present <- rownames(model_table[which(model_table[,rel]>0),])
  
  # Ignore OGs related to high ccp values # REMEMBER TO REMOVE IT IF URTE WANTS IT
  if (mean(DT2$ccp[DT2$species %in% present])< 35){
  DT3 <- DT2 
  DT3$copy_number = 0
    
  DT3[which(model_table[,rel]>0,), "copy_number"] <- model_table[
                                    which(model_table[,rel]>0),rel]
  
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
          scale_fill_gradientn(colours = c("black","steelblue1")) +
          ggtitle(rel)
  #ploti
  ggsave(ploti, filename = rel, device = "png")
  }
}


##################### MMER Presence/absence ###########

PA_table <- model_table
PA_table[model_table > 1] = 1 


p_vals_PA <- c()

# APPLY MMER TO EVERY MARKER
for(i in 5001:12000){#1:ncol(PA_table)){ 
  #browser() 
  # Prepare DT table
  DT3 <- cbind(DT2, PA_table[,i], PCA1$rotation[,1], 
               PCA1$rotation[,2], PCA1$rotation[,3])
  colnames(DT3) <- c("species","ccp", "marker", "PC1", "PC2", "PC3")          
  DT3$dups <- for_corr[species, "dups"]
  
  # Apply Mixed Modelling 
  result <- mmer(ccp~1 + marker + PC1 + PC2,
                 random=~vs(species), #, Gu=phyl_covar),
                 rcov=~units, data=DT3)
  
  result2 <- mmer(ccp~1 + PC1 + PC2,
                  random=~vs(species), #, Gu=phyl_covar),
                  rcov=~units, data=DT3, verbose = F)
  
  u.hat <- result2$U$`u:species`$ccp
  u.hat <- u.hat[match(DT3$species,names(u.hat))]
  
  #plot(DT3$ccp, u.hat)
  
  # Anova 
  aov <- anova.mmer(result, result2) #type = 2)["marker",5]
  
  p_vals_PA <- append(p_vals_PA, 
                   as.numeric(strsplit(
                     paste(aov["mod2","PrChisq"]), " ")[[1]][1])
  )
}

head(p_vals_PA)
which(p_vals_PA<0.05)

plot(-log10(p_vals_PA), xlab = "Gene Families (HOGs)") +
  abline(h=1.30103, col="red")


#colnames(model_table[,p_vals<0.05])

write(p_vals_PA, file = "p_vals_PA_no_tetras.txt", ncolumns = 1)
relevant <- colnames(model_table[,p_vals_PA<0.05])
write(relevant, file = "relevantHOGs_PA_no_tetras.txt", ncolumns = 1)















# store the most relevant family names
relevantPA <- names(which((10^(-result_PA$scores[2,]) < 0.05))) # exp(2.5)

p_vals_PA <- round(10^(-result_PA$scores[2,((10^(-result_PA$scores[2,]) < 0.05))]), 4)
relevantPA <- names(p_vals_PA)
write.csv(data.frame("p_val"=p_vals_PA, "HOG"=relevantPA), file = "interesting_PA.csv", quote = FALSE, row.names = FALSE)


## Filter and save most relevant families
out_PA <- data.frame("HOG" = 0, "p-vals" = 0 ,"C3-C4_avg_genes"=0, "C3-C4_sd_genes"=0,
                          "C3-C4*_avg_genes"=0, "C3-C4*_sd_genes"=0,
                          "C3_avg_genes"=0, "C3_sd_genes"=0, check.names = F)
library("ggplot2")
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


############ PA plots (scores and ccp values) #######

require(ggplot2)

##### Relevant GWAS scores
qplot(seq_along(result_PA$scores[2,relevantPA]), result_PA$scores[2,relevantPA]) +
  geom_point(colour ="navyblue", size= 1.2) +
  ggtitle("GWAS scores per gene family (PA)") + 
  theme(axis.text = element_text(size = rel(1.4), face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = rel(2)),
        axis.title = element_text(size=rel(1.8))) +
  xlab("Gene families") + ylab("-log10(p-val)")

####


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


############ Save plots of specific HOGs selected by Urte #####

selected <- c("N0.HOG0003237",  "N0.HOG0025962",  "N0.HOG0028581", 
              "N0.HOG0034511",  "N0.HOG0034527",  "N0.HOG0021276", 
              "N0.HOG0026267")


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



############ SNP model #############

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


## Others #################

# Only 98 families have single protein
sum(rowSums(ogroups[,4:length(colnames(ogroups))])<2)

# Plot this
#hist(rowSums(ogroups[,4:length(colnames(ogroups))]), 500, 
#  xlim = range(0,100), main = "Number of genes per Orthogroup",
#  xlab = "Number of Proteins in OG (Max is 1400)")

sp_pres <- ogroups[,4:length(colnames(ogroups))]
sp_pres[sp_pres>1]=1



############ Experimental segment ########
## Feldenstein contrasts exploratory analysis
# https://lukejharmon.github.io/ilhabela/instruction/2015/07/02/phylogenetic-independent-contrasts/
ttree <- drop.tip(tree,setdiff(tree$tip.label, species))
ape::pic(x = 1: 19, phy=ttree, var.contrasts = T)

x <- fastBM(tree)
y <- fastBM(tree)
plot(x,y)

phylomorphospace(tree,cbind(x,y),label="off",node.size=c(0.5,0.7))


############ GWAS copy number #####################

library("sommer")

model_table <- t(matrix(model_table))

result <- GWAS(ccp ~ 1, 
               random =~vs(species, Gu=phyl_covar),
               rcov=~units, gTerm = "u:species",
               M = model_table, data = DT2, n.PC = 3) #, n.core = 4 )

# store the most relevant family names
relevant <- names(which((10^(-result$scores[2,]) < 0.05))) # exp(2.5)

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


