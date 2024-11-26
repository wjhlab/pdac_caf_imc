################################
#    Neighborhood analysis     #
################################
# Fig2C, 2D, 3K, 3L, 5H, 6A, 6B 

library(ComplexHeatmap)
library(reshape2)
library(dittoSeq)
library(dplyr)
library(tibble)
library(cowplot)

### Load backup data ===========================================================

#CAF
cafsubset_data<- readRDS('./backup/CAFsubset_phenographoutput_01_095_fixed.RDS')
data_caf<- cafsubset_data[[1]]
#Tumor
tumorsubset_data<- readRDS('./backup/tumorsubset_phenographoutput_01_095_fixed.RDS')
data_tumor<-tumorsubset_data[[1]]
#Immune cells 
cd45subset_data<- readRDS('./backup/CD45subset_phenographoutput_01_095_fixed.RDS')
data_cd45<- cd45subset_data[[1]]



tumorcell_ind <- as.numeric(rownames(data_tumor))

data_t_basal <- data_tumor[data_tumor$cluster1m_refined == "Basal",]
data_t_classical <- data_tumor[data_tumor$cluster1m_refined == "Classical",]
data_t_mixed <- data_tumor[data_tumor$cluster1m_refined == "Mixed",]


basalcell_ind <- paste(data_t_basal$ImageId,data_t_basal$CellId,sep="_")
classicalcell_ind <- paste(data_t_classical$ImageId,data_t_classical$CellId,sep="_")
mixedcell_ind <- paste(data_t_mixed$ImageId,data_t_mixed$CellId,sep="_")

data_caf_clean_ext <- data_caf %>%
  mutate(Id_ext = paste(ImageId, CellId, sep = "_"),
         NN1_ext = paste(ImageId, NN1, sep = "_"),
         NN2_ext = paste(ImageId, NN2, sep = "_"),
         NN3_ext = paste(ImageId, NN3, sep = "_"))


caf_ind <- c(data_caf_clean_ext$NN1_ext,data_caf_clean_ext$NN2_ext,data_caf_clean_ext$NN3_ext)
caf_ind_basal <- caf_ind[caf_ind %in% basalcell_ind]
caf_ind_classical <- caf_ind[caf_ind %in% classicalcell_ind]
caf_ind_mixed <- caf_ind[caf_ind %in% mixedcell_ind]

cafs_near_basal <- c(data_caf_clean_ext[data_caf_clean_ext$NN1_ext %in% caf_ind_basal,]$cluster1m_refined,
                     data_caf_clean_ext[data_caf_clean_ext$NN2_ext %in% caf_ind_basal,]$cluster1m_refined,
                     data_caf_clean_ext[data_caf_clean_ext$NN3_ext %in% caf_ind_basal,]$cluster1m_refined)

cafs_near_classical <- c(data_caf_clean_ext[data_caf_clean_ext$NN1_ext %in% caf_ind_classical,]$cluster1m_refined,
                         data_caf_clean_ext[data_caf_clean_ext$NN2_ext %in% caf_ind_classical,]$cluster1m_refined,
                         data_caf_clean_ext[data_caf_clean_ext$NN3_ext %in% caf_ind_classical,]$cluster1m_refined)

cafs_near_mixed <- c(data_caf_clean_ext[data_caf_clean_ext$NN1_ext %in% caf_ind_mixed,]$cluster1m_refined,
                     data_caf_clean_ext[data_caf_clean_ext$NN2_ext %in% caf_ind_mixed,]$cluster1m_refined,
                     data_caf_clean_ext[data_caf_clean_ext$NN3_ext %in% caf_ind_mixed,]$cluster1m_refined)


caf_near_basal<-as.data.frame(table(cafs_near_basal))
caf_near_classical<-as.data.frame(table(cafs_near_classical))
caf_near_mixed<-as.data.frame(table(cafs_near_mixed))

caf_near_epi_sum <- data.frame(basal = caf_near_basal$Freq, classical = caf_near_classical$Freq, mixed = caf_near_mixed$Freq)
rownames(caf_near_epi_sum) <- caf_near_basal$cafs_near_basal

##Heatmap of CAF subtypes by tumor cell subtypes
### Fig 2C ####
pdf('./output/Fig2C.pdf')
Heatmap(t(scale(t(as.matrix(caf_near_epi_sum)))),
        name = "Scaled\nFreq",
        width = ncol(caf_near_epi_sum)*unit(7, "mm"), 
        height = nrow(caf_near_epi_sum)*unit(7, "mm"),
        heatmap_legend_param = list(
          labels_gp=gpar(fontsize=6)))
dev.off()

# Convert data to format for ggplot
plot_data_long<- melt(caf_near_epi_sum%>% rownames_to_column())
names(plot_data_long)<-c("TumorCellType", "Subtype", "Frequency")


colorassigned<- dittoColors()[1:length(unique(plot_data_long$TumorCellType))]
clusterlevels=c("iCAF",
                "myCAF",
                "CD105+ CAF",
                "ApCAF",
                "CXCL12+ CAF",
                "FAP+ CAF",
                "CAF undefined")
names(colorassigned)<- clusterlevels
plot_data_long$TumorCellType<- factor(plot_data_long$TumorCellType, levels=clusterlevels)

# Create a stacked bar chart with custom colors
### Fig 2D ####
pdf('./output/Fig2D.pdf', height=3.5, width=3)
ggplot(plot_data_long, aes(x = Subtype, y = Frequency, fill = TumorCellType)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90,vjust =0.5, hjust = 1, color="black"))+
  labs(x = "", y = "Frequency", fill = "Tumor Cell Type")+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = colorassigned) 
dev.off()




### Figure 3K, 3L =====
data_tumor_clean_ext2 <- data_tumor %>%
  mutate(Id_ext = paste(ImageId, CellId, sep = "_"),
         NN1_ext = paste(ImageId, NN1, sep = "_"),
         NN2_ext = paste(ImageId, NN2, sep = "_"),
         NN3_ext = paste(ImageId, NN3, sep = "_"),
         NN4_ext = paste(ImageId, NN4, sep = "_"),
         NN5_ext = paste(ImageId, NN5, sep = "_"),
         NN6_ext = paste(ImageId, NN6, sep = "_"),
         NN7_ext = paste(ImageId, NN7, sep = "_"),
         NN8_ext = paste(ImageId, NN8, sep = "_"),
         NN9_ext = paste(ImageId, NN9, sep = "_"),
         NN10_ext = paste(ImageId, NN10, sep = "_"),
         NN11_ext = paste(ImageId, NN11, sep = "_"),
         NN12_ext = paste(ImageId, NN12, sep = "_"))
keymarkerlist <- c("KI67","HLADR","VIM","NCAD","Ecad")

CAF_ind <- paste0(data_caf$ImageId, "_", data_caf$CellId)

##for CAF near/not near Basal
basalnearCAF <- c()

data_tumor_clean_ext2_subtype<- data_tumor_clean_ext2[data_tumor_clean_ext2$cluster1m_refined=="Basal", ]
#save the cell index for any basal cell that has a given all CAF as its nearest neighbor (any one of 12)
#if any of the 12 neighbors belong within the CAF indices and the number of TRUEs are >0, then return the cell index
for(i in 1:nrow(data_tumor_clean_ext2_subtype)){
  if(
    sum(data_tumor_clean_ext2_subtype[i,c('NN1_ext','NN2_ext','NN3_ext','NN4_ext','NN5_ext','NN6_ext',
                                  'NN7_ext','NN8_ext','NN9_ext','NN10_ext','NN11_ext','NN12_ext')] %in% CAF_ind) >0){
    basalnearCAF <- c(basalnearCAF,data_tumor_clean_ext2_subtype[i,]$Id_ext)
  }}

#create a ggdf based on the cell indices returned from the loop  
ggdfvio_caf<- rbind(data_tumor_clean_ext2[data_tumor_clean_ext2$Id_ext %in% basalnearCAF,c("Id_ext",keymarkerlist)],
                    data_tumor_clean_ext2[!data_tumor_clean_ext2$Id_ext %in% basalnearCAF,c("Id_ext",keymarkerlist)])
ggdfvio_caf$near <- 1
ggdfvio_caf[ggdfvio_caf$Id_ext %in% basalnearCAF,]$near<-"Near"
ggdfvio_caf[!ggdfvio_caf$Id_ext %in% basalnearCAF,]$near<-"Not_Near"
ggdfvio_caf[,keymarkerlist] <- asinh(ggdfvio_caf[,keymarkerlist] / 0.8)


pdf("./output/3K_basal_CAF_NN_Ecad.pdf",width=3,height=3)
ggplot(ggdfvio_caf, aes(x=near, y=Ecad, fill=near))+
  geom_violin(linewidth=0.25)+
  geom_boxplot(outlier.color = NA, width=0.25, alpha=0.2, linewidth=0.25)+
  theme(axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        axis.line = element_line(color="black",linewidth = 0.25),
        axis.text = element_text(color="black", size=7),
        panel.background = element_blank())+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  stat_compare_means(aes(label = paste0("p=", after_stat(p.format))))
dev.off()

pdf("./output/3L_basal_CAF_NN_Vim.pdf",width=3,height=3)
ggplot(ggdfvio_caf, aes(x=near, y=VIM, fill=near))+
  geom_violin(linewidth=0.25)+
  geom_boxplot(outlier.color = NA, width=0.25, alpha=0.2, linewidth=0.25)+
  theme(axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        axis.line = element_line(color="black",linewidth = 0.25),
        axis.text = element_text(color="black", size=7),
        panel.background = element_blank())+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  stat_compare_means(aes(label = paste0("p=", after_stat(p.format))))
dev.off()



### Figure 5 (IL8) =====

data_caf<- cafsubset_data[[1]]

data_caf1 <- data.matrix(data_caf[,Stromamarkers])
data_caf2 <- asinh(data_caf1 / 0.8)

#phenograph clustering of data
rng <- colQuantiles(data_caf2, probs = c(0.01, 0.95))
data_caf01 <- t((t(data_caf2) - rng[, 1]) / (rng[, 2] - rng[, 1]))
data_caf01[data_caf01 < 0] <- 0; data_caf01[data_caf01 > 1] <- 1
#set.seed(1234)
#phenographout<-Rphenograph(data_caf01)
phenographout<- readRDS('./backup/phenograph_01_095_fixed_caf.rds')
data_caf$cluster_str<-factor(membership(phenographout[[2]]))

#save phenographoutput for subset
#save as RDS file
#cafsubset_data <- list(data_caf, data, data01, csv_full)
#saveRDS(cafsubset_data, "CAFsubset_phenographoutput_01_095_fixed.RDS")

cluster_mean <- data.frame(data_caf01, cluster = data_caf$cluster_str, check.names = FALSE) %>%
  group_by(cluster) %>% summarize_all(list(mean))
cluster_mean_mat<-cluster_mean[,Stromamarkers] #here you are adding this subsetting
rownames(cluster_mean_mat)<-1:nrow(cluster_mean_mat)

## Load merge file
clusterMergeFile = paste0(work,"/Config/merge_CAF.xlsx") 
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c("FAP+ IL6+",
                "VIM+",
                "VIM+ CD105+",
                "HLADR+",
                "SMA+ CD105+",
                "SMA+ IL8+",
                "SMA+ CXCL12+",
                "FAP+",
                "SMA+",
                "CD105+",
                "SMA+ PDPN+",
                "FAP+ CD105+",
                "CAF undefined",
                "SMA+ CXCL12+ FAP+ PDPN+",
                "CD105+ HLADR+",
                "FAP+ CXCL12+ IL6+ IL8+")
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
names(colorassigned)<-clusterlevels
mm1 <- match(data_caf$cluster_str, cluster_merging$original_cluster)
data_caf$cluster2m_refined <- cluster_merging$new_cluster[mm1]


data_caf_clean_ext <- data_caf %>%
  mutate(Id_ext = paste(ImageId, CellId, sep = "_"),
         NN1_ext = paste(ImageId, NN1, sep = "_"),
         NN2_ext = paste(ImageId, NN2, sep = "_"),
         NN3_ext = paste(ImageId, NN3, sep = "_"))


caf_ind <- c(data_caf_clean_ext$NN1_ext,data_caf_clean_ext$NN2_ext,data_caf_clean_ext$NN3_ext)
caf_ind_basal <- caf_ind[caf_ind %in% basalcell_ind]
caf_ind_classical <- caf_ind[caf_ind %in% classicalcell_ind]
caf_ind_mixed <- caf_ind[caf_ind %in% mixedcell_ind]

cafs_near_basal <- c(data_caf_clean_ext[data_caf_clean_ext$NN1_ext %in% caf_ind_basal,]$cluster2m_refined,
                     data_caf_clean_ext[data_caf_clean_ext$NN2_ext %in% caf_ind_basal,]$cluster2m_refined,
                     data_caf_clean_ext[data_caf_clean_ext$NN3_ext %in% caf_ind_basal,]$cluster2m_refined)

cafs_near_classical <- c(data_caf_clean_ext[data_caf_clean_ext$NN1_ext %in% caf_ind_classical,]$cluster2m_refined,
                         data_caf_clean_ext[data_caf_clean_ext$NN2_ext %in% caf_ind_classical,]$cluster2m_refined,
                         data_caf_clean_ext[data_caf_clean_ext$NN3_ext %in% caf_ind_classical,]$cluster2m_refined)

cafs_near_mixed <- c(data_caf_clean_ext[data_caf_clean_ext$NN1_ext %in% caf_ind_mixed,]$cluster2m_refined,
                     data_caf_clean_ext[data_caf_clean_ext$NN2_ext %in% caf_ind_mixed,]$cluster2m_refined,
                     data_caf_clean_ext[data_caf_clean_ext$NN3_ext %in% caf_ind_mixed,]$cluster2m_refined)


caf_near_basal<-as.data.frame(table(cafs_near_basal))
caf_near_classical<-as.data.frame(table(cafs_near_classical))
caf_near_mixed<-as.data.frame(table(cafs_near_mixed))

caf_near_epi_sum <- data.frame(basal = caf_near_basal$Freq, classical = caf_near_classical$Freq, mixed = caf_near_mixed$Freq)
rownames(caf_near_epi_sum) <- caf_near_basal$cafs_near_basal

##Heatmap of CAF subtypes by tumor cell subtypes
pdf('./output/Fig5H.pdf')
Heatmap(t(scale(t(as.matrix(caf_near_epi_sum)))),
        name = "Scaled\nFreq",
        width = ncol(caf_near_epi_sum)*unit(7, "mm"), 
        height = nrow(caf_near_epi_sum)*unit(7, "mm"),
        heatmap_legend_param = list(
          labels_gp=gpar(fontsize=6)))
dev.off()


### Figure 6 (immune tumor nearest neighbor) =====
get_immune_indices <- function(cluster_name) {
  data_subset <- data_cd45[data_cd45$cluster1m_refined == cluster_name, ]
  paste(data_subset$ImageId, data_subset$CellId, sep = "_")
}

# Define the cluster names of interest
cluster_names <- c("T cell", 
                   "unassigned", 
                   "CD57+", 
                   "B cell", 
                   "Myeloid")
indices_list <-lapply(cluster_names, get_immune_indices)
names(indices_list) <- cluster_names

# Function to find neighboring cells for specific indices
get_tumor_near <- function(ind) {
  tumor_ind <- c(data_tumor_clean_ext$NN1_ext, data_tumor_clean_ext$NN2_ext, data_tumor_clean_ext$NN3_ext)
  tumor_ind <- tumor_ind[tumor_ind %in% ind]
  c(data_tumor_clean_ext[data_tumor_clean_ext$NN1_ext %in% tumor_ind, ]$cluster1m_refined,
    data_tumor_clean_ext[data_tumor_clean_ext$NN2_ext %in% tumor_ind, ]$cluster1m_refined,
    data_tumor_clean_ext[data_tumor_clean_ext$NN3_ext %in% tumor_ind, ]$cluster1m_refined)
}



# Extend data_tumor with necessary columns for neighbor identification
data_tumor_clean_ext <- data_tumor %>%
  mutate(Id_ext = paste(ImageId, CellId, sep = "_"),
         NN1_ext = paste(ImageId, NN1, sep = "_"),
         NN2_ext = paste(ImageId, NN2, sep = "_"),
         NN3_ext = paste(ImageId, NN3, sep = "_"))

# Function to get neighboring clusters for a specific index
get_tumor_near <- function(ind) {
  tumor_ind <- c(data_tumor_clean_ext$NN1_ext, data_tumor_clean_ext$NN2_ext, data_tumor_clean_ext$NN3_ext)
  tumor_ind <- tumor_ind[tumor_ind %in% ind]
  c(data_tumor_clean_ext[data_tumor_clean_ext$NN1_ext %in% tumor_ind, ]$cluster1m_refined,
    data_tumor_clean_ext[data_tumor_clean_ext$NN2_ext %in% tumor_ind, ]$cluster1m_refined,
    data_tumor_clean_ext[data_tumor_clean_ext$NN3_ext %in% tumor_ind, ]$cluster1m_refined)
}

tumor_near_results <- lapply(indices_list, get_tumor_near)

# Convert each result to a data frame with frequency counts
tumor_near_tables <- lapply(tumor_near_results, function(near) as.data.frame(table(near)))

# Combine the tables into a summary data frame
tumor_near_immune_sum <- do.call(cbind, lapply(tumor_near_tables, `[[`, "Freq"))
colnames(tumor_near_immune_sum) <- cluster_names
rownames(tumor_near_immune_sum) <- tumor_near_tables$`B cell`$near  

# Remove rows where row name is "Unassigned"
tumor_near_immune_sum <- tumor_near_immune_sum[rownames(tumor_near_immune_sum) != "Unassigned", ]


pdf('./output/Fig6A.pdf')
Heatmap(t(scale(t(as.matrix(tumor_near_immune_sum)))),
        name = "Scaled\nFreq",
        width = ncol(tumor_near_immune_sum)*unit(7, "mm"), 
        height = nrow(tumor_near_immune_sum)*unit(7, "mm"),
        heatmap_legend_param = list(
          labels_gp=gpar(fontsize=6)))
dev.off()


## Tcell subtypes 
tcellsubset_data<- readRDS('./backup/tcellsubset_phenographoutput_01_095_fixed.RDS')
data_tcell<- tcellsubset_data[[1]]

# Function to get indices based on cluster type
get_Tcell_indices <- function(cluster_name) {
  data_subset <- data_tcell[data_tcell$cluster1m_refined2== cluster_name, ]
  paste(data_subset$ImageId, data_subset$CellId, sep = "_")
}


# Retrieve indices for each T-cell cluster
clusterlevels=c("CD8+",
                "CD4",
                "activated CD8",
                "activated CD4",
                "memory CD4",
                "Naïve CD4",
                "naïve CD8",
                "dysfunctional CD8")
indices_list <- lapply(clusterlevels, get_Tcell_indices)
names(indices_list) <- clusterlevels

# Function to find neighboring cells for specific indices
get_tumor_near <- function(ind) {
  tumor_ind <- c(data_tumor_clean_ext$NN1_ext, data_tumor_clean_ext$NN2_ext, data_tumor_clean_ext$NN3_ext)
  tumor_ind <- tumor_ind[tumor_ind %in% ind]
  c(data_tumor_clean_ext[data_tumor_clean_ext$NN1_ext %in% tumor_ind, ]$cluster1m_refined,
    data_tumor_clean_ext[data_tumor_clean_ext$NN2_ext %in% tumor_ind, ]$cluster1m_refined,
    data_tumor_clean_ext[data_tumor_clean_ext$NN3_ext %in% tumor_ind, ]$cluster1m_refined)
}

# Apply the function for each T-cell type and store in a list
tumor_near_results <- lapply(indices_list, get_tumor_near)

# Convert to data frames and count occurrences
tumor_near_tables <- lapply(tumor_near_results, function(near) as.data.frame(table(near)))

tumor_near_Tcell_sum  <- data.frame(CD8 = tumor_near_tables$`CD8+`$Freq, 
                                    CD4 = tumor_near_tables$CD4$Freq, 
                                    activatedCD8 = tumor_near_tables$`activated CD8`$Freq, 
                                    activatedCD4 = tumor_near_tables$`activated CD4`$Freq, 
                                    memoryCD4 = tumor_near_tables$`memory CD4`$Freq,
                                    naiveCD8 = tumor_near_tables$`naïve CD8`$Freq,
                                    naiveCD4 = tumor_near_tables$`Naïve CD4`$Freq,
                                    dysfunctionalCD8 = tumor_near_tables$`dysfunctional CD8`$Freq)

rownames(tumor_near_Tcell_sum) <- tumor_near_tables$`CD8+`$near

#remove Unassigned tumor cell column
tumor_near_Tcell_sum <- tumor_near_Tcell_sum[rownames(tumor_near_Tcell_sum)[rownames(tumor_near_Tcell_sum) != "Unassigned"],]


# Convert data to format for ggplot
plot_data_long<- melt(tumor_near_Tcell_sum%>% rownames_to_column())
names(plot_data_long)<-c("TumorCellType", "ImmuneType", "Frequency")

clusterlevels=c("CD4",
                "CD8",
                "activatedCD8",
                "memoryCD4",
                "activatedCD4",
                "naiveCD4",
                "naiveCD8",
                "dysfunctionalCD8")
colorassigned<- c("#4D4D4D",
                  "#FFBE2D",
                  "#80C7EF",
                  "#00F6B3",
                  "#F4EB71",
                  "#06A5FF", 
                  "#FF8320",
                  "#D99BBD")
names(colorassigned)<- clusterlevels
plot_data_long$ImmuneType<- factor(plot_data_long$ImmuneType, levels=clusterlevels)

# Create a stacked bar chart with custom colors
pdf('./output/Fig6B.pdf', height=3.5, width=3)
ggplot(plot_data_long, aes(x = TumorCellType, y = Frequency, fill = ImmuneType)) +
  geom_bar(stat = "identity") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90,vjust =0.5, hjust = 1, color="black"))+
  labs(x = "", y = "Frequency", fill = "T cell subtypes")+
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = colorassigned) 
dev.off()















