#################### PURPOSE ####################################
# This code takes in features that have been processed (nans removed, z-scored..
# sign-adjusted, harmonized, and filtered on kurtosis/skewness/iqr) and applies
# consensus clustering and pca to make nice dendrograms of the cases clustered 
# into phenotypes along with the features group together through clustering
#
##################### Add libraries ##############################

# for R markdown use ```{r blah for .Rmd}

library(gplots)
library(ggplot2)
library("readxl")
library(heatmap3)
library(sigclust)
library(ConsensusClusterPlus)
library(colorRamps)
library(plotrix)
library(gmodels)
#library(vioplot)
library(scales)
library(plyr)
library(lattice)
library(cluster)
#library(fpc)
library(survC1)     # c-statistic
library("tidyverse")  # (Includes tibble, magrittr, tidyr, ggplot2, others)
library("survival")
#library("survminer")  # survival curves in ggplot
library("survC1")
library(survival)
library(survMisc)
#library("xlsx")
#library("ggpubr")
library("dendextend") #for dendrograms
library("factoextra") #for visualizing feature contributions to principle components
# for autoplot:
#install_github('sinhrks/ggfortify')
library(ggfortify)



################### DEFINE ALL FUNCTIONS: #############################

## PURPOSE: Load the feature data as a data frame

load_data <- function(filename,sheetname){
  #load the excel sheet of interest (col_names= T means labels in first row)
  featureData = read_excel(filename, sheet = sheetname, col_names = T)
  # store featureData as a data frame (minus the first row of labels)
  feature_data <- data.frame(featureData[,(c(2:length(featureData)))])
  return(list(featureData,feature_data))
}

#------------------------------------Global PCA----------------------------------
#PURPOSE: Perform PCA on the features to reduce dimensionality 
global_pca <- function(featureDF,significance_level){
  # featureDF: data frame of features (col) and cases (row)
  # significance_level: level of variance accounted for by PCs
  
  #perform PCA on featureDF (scale = TRUE scales to unit variance)
  features_pca<-prcomp(featureDF, center= TRUE, scale= TRUE)
  
  #plot the principal components that account for most of the varaince:
  #autoplot(features_pca)#, data = featuresDF, colour = 'Species')
  
  # figure out how much of the variance is explained by each component
  pca_summary <- summary(features_pca)
  for (i in 1:length(pca_summary$importance)){
    if (pca_summary$importance[3,i] > significance_level){
      pcs=i
      break}
  }
  pca_score = as.data.frame(features_pca$x)
#pick out only the significantly contributing PCs (up to significance_level)
  pca_score <-data.frame(pca_score[,(c(1:pcs))])
  return(list(features_pca,pca_score, pca_summary))

}


#-------------------------------Contributing Features----------------------------
#PURPOSE: Pick out features that contribute above the significance level in globa_pca()
contributing_features <-function(featurePC,pca_score,featureDF,file){
  # featuresPC: output from global_pca function, first argument in returned list
  # featureDF: data frame of features, outputted from load_data function
  # file: file name to which to write the output, PC_contrib
  
  # loop through PCA scores of most significant features to get those features
  for(i in 1:length(pca_score)) {
    summary <-facto_summarize(featurePC, element="var", axes=i)
    #pick out only the features which contribute more than the expected amount
    #it's like you have 100%, and each feature should contribute at least equally?
    #the contrib column of each PC sums to 100%!
    contrib_features <- summary[which(summary$contrib>100/length(featureDF)),]
    contrib_features <- contrib_features[order(-contrib_features$contrib),]
    keeps <- c("name", "contrib")
    PC_contrib <- contrib_features[keeps]
    
    #commented out because difficulty loading appropriate packages
    #fviz_contrib(features_pca, choice= "var", axes=i)
    
    #write.xlsx(PC_contrib, file, sheetName = paste("PC",i, sep = ""),
    #           col.names = TRUE, row.names = FALSE, append = TRUE)
    
    # write.xlsx() wasn't working so wrote to csv instead
    write.csv(PC_contrib, paste("PCA_contributingFeatures_", as.character(i),".csv"),row.names = TRUE)} 
  return(list(PC_contrib, summary))
}

#-----------------------------------Cluster Features-----------------
#PURPOSE: Cluster the features so that similar features are grouped together
cluster_features <-function(featureDF, clustAlg, nrep, maxClust,pval){
  # featureDF: data frame of features, outputted from load_data
  # clustAlg: which clustering algorithm for the ConcensusClusterPlus() func.
  # nrep: number of reps for the ConcensusClusterPlus() func.
  # maxClust: how many k clusters for ConcensusClusterPlus() func.
  # pval: pmax in ConsensusClusterPlus() func.
  
  featuresMatrix<-as.matrix(featureDF)  
  feature_distance = "pearson"
  feature_linkage = "complete" #options: single, average,complete, ward.d, ward.d2
  result_features = ConsensusClusterPlus(featuresMatrix, maxK=maxClust, reps=nrep, pItem= pval, clusterAlg= clustAlg, distance=feature_distance, innerLinkage=feature_linkage, finalLinkage=feature_linkage)
  #now this one is Definitely 4 (subjectively looks more like 3 to me)
  return(result_features)
}


#---------------------------------PCA within clusters----------------------
#PURPOSE: To purform PCA within the clusters identified by cluster_features()
fancy_pca <- function(featureDF, clustered_features, significance_level,num_clusters,dowrite,file){
  # featureDF: data frame of features, outputted from load_data
  # clustered_features: output from cluster_features()
  # significance_level: level of variance accounted for by PCs
  # num_clusters: number of clusters identified by cluster_features()
  # dowrite: flag for writing the results to file
  
  feat_class <- clustered_features[[num_clusters]]$consensusClass
  feat_class <- as.matrix(feat_class)
  
  
  clusterlist=list() 
  count = 0
  for (i in 1:num_clusters){
    cluster_comp = prcomp(featureDF[,which(feat_class==i)], center= TRUE, scale= TRUE)
    cluster_summ = summary(cluster_comp)
    cluster_score = as.data.frame(cluster_comp$x)
    
    #now get all the relevant scores
    for (k in 1:length(cluster_summ$importance)){
      count=count+1
      summary <-facto_summarize(cluster_comp, element="var", axes=k)
      contrib_features <- summary[which(summary$contrib>100/length(featureDF)),]
      contrib_features <- contrib_features[order(-contrib_features$contrib),]
      keeps <- c("name","contrib")
      PC_contrib_within_cluster <- contrib_features[keeps]
      
      if (dowrite==1){
        #write.xlsx(PC_contrib_within_cluster, file, sheetName = paste("Cluster ",i,", PC", k, sep = ""),
        #           col.names = TRUE, row.names = FALSE, append = TRUE)
        write.csv(PC_contrib_within_cluster, file, row.names = TRUE)
      }
      
      clusterlist[[count]]<-cluster_score[k]
      if (cluster_summ$importance[3,k] > significance_level){
         pcs_cluster=k
        break}
    }
    
     cluster_score <-data.frame(cluster_score[,(c(1:pcs_cluster))])
    
    
    new_pca_score = do.call(cbind,clusterlist)}
  return(new_pca_score)
}
# commented out because the required class wasn't installing properly
# fviz_contrib(cluster_comp, choice= "var", axes=3)


# #--------------------------Dendrograms!--------------------------------
# 
#PURPOSE: To plot dendrogram
dendrogram <- function(colCluster_1,colCluster_2){
  dend1<- as.dendrogram(colCluster_1)
  dend2<- as.dendrogram(colCluster_2)
  dendlist(dend2, dend1) %>%    ## This provides a visualization of how the clustering assignment, and resulting dendrograms change between each experiment
    untangle(method = "step1side") %>% # Find the best alignment layout
    tanglegram() %>%
    entanglement()
}


#-------------------- Consensus Clustering & SigClust ---------------------
#PURPOSE: To do consensus clustering on cases and test significance of results
cluster_data <- function(featureDF,numClusters,linkage_method,distance_metric, fiSaveDir, figSaveName,assertClust){
  
  # data_to_cluster<-featureData
  # numClusters<-5
  # linkage_method<-"ward.D2"
  # distance_metric<-"euclidean"
  # if assertClust is not 0, then instead of taking the results from significance testing,...
  # ... assert whatever the value of assertClust is as the number of clusters
  
  # First we cluster the cases
  caseCluster = as.dist(dist(featureDF, method = distance_metric));
  #cluster distance based on Ward's method
  colCluster = hclust(caseCluster, method = linkage_method);
  # calcaute distance between features (rows)
  featureCluster = t(featureDF)
  # now cluster the features too
  featureCluster = as.dist(dist(featureCluster, method= distance_metric));
  rowCluster = hclust(featureCluster, method = linkage_method)
  
  
  result = ConsensusClusterPlus(t(featureDF), maxK=numClusters, reps=50, pItem= 0.8, pFeature=1, clusterAlg= "hc", distance=distance_metric, innerLinkage=linkage_method, finalLinkage = linkage_method)
  #test = calcICL(result,title="untitled_consensus_cluster")

  #test significance of each cluster split
  signif = 0
  leadingi = 0
    for (i in numClusters:2) {
      colCluster.assignment <- cutree(colCluster,k=i)
      Cluster = as.data.frame(as.matrix(colCluster.assignment)) #if this is the first test, it should have the original assignments
      count = 0;
      #cat("Cluster$V1: ", Cluster$V1)
      for (k in i:1){
        Cluster$split[Cluster$V1 != k]=2
        Cluster$split[Cluster$V1 == k]=1
        sig<-sigclust(featureDF, 1000, labflag = 1, label = Cluster$split, icovest=3)
        
        cat("Cluster",i,": p =", sig@pval,fill=TRUE)
        if (sig@pval < 0.05)
        {count = count+1}
      }
      if (count/i > leadingi)
      {leadingi = count/i
      signif = i}
      cat("count/i = ",count/i, "signif = ", signif, fill=TRUE)
    }

  if(assertClust!=0)
    {signif = assertClust}
  end
  
  out_cluster <- as.data.frame(cutree(colCluster,signif)) 
  colors=c("red","green","blue","yellow","purple")
  ColumnCluster_col<-colByValue(as.matrix(out_cluster),
                                col=colors[1:signif])
  # In Andrew's version of this code, the information needs to be "colored by value"
  # in this same way that ColumnCluster_col is here...
  # Assign colors to each of those values, for example age, survival time, etc
  # Use cbind() and put all those colored vectors together
  
  col_combos.train<-ColumnCluster_col 
  colnames(col_combos.train)<-c("Phenotypes")
  
  ColSideColors=col_combos.trainColSideColors=col_combos.train
  # #run this if you're getting the invalid graphics state error on the heatmap
  # dev.off()
  ## This is how you generate the heatmap
  out_heatmap<-heatmap3(t(as.matrix(featureDF)), margin = c(6,6),  # Make sure the feature data is input as a transpose (cases as columns, features as rows)
                       balanceColor=TRUE,
                       Rowv = as.dendrogram(rowCluster),  # This is the feature dendrogram made from clustering features above
                       Colv = as.dendrogram(colCluster),  # This is the patient/case dendrogram made from clustering cases above
                       col = colorpanel(500, 'blue', 'white', 'red4'),                # This is the color pallete made above
                       scale = "none",
                       labCol = "",
                       ColSideColors = col_combos.train,# This is the matrix of colors representing the clinical features, where you assign additional coolumns of information
                       RowSideLabs = "",
                       revC = T,
                       #breaks= break_vals,
                       cexRow = 1,    # Text sizes
                       cexCol = 1)    # Text sizes

  return(list(out_cluster,out_heatmap, result, sig))
}


################### RUN THE COMMANDS! User Input: #############################

### Define all inputs to the functions here:

# define your path and feature file and sheet within feature file
workingDir = "/Users/kalina/Documents/CBIG/Project_DCIS_R01/Final_results" #where code lives and where R will save things
setwd(workingDir) 
feature_file = "first_postcontrast_features_allResults_n297_local_resampv2.xlsx"
feature_sheet = "processed_harm_features2_ex"
#covar_file = "first_postcontrast_features_covars_n297_resamp_experiments.xls"
#covar_sheet = "clinical_covars_numeric"


# define some file names for saving things
figSaveDir = workingDir
withinClust_pca_file = "withinClust_PCA.csv" #name of file for saving fancy_pca() output
contrib_pca_file = "contrib_75.csv" # file for saving output from contributing_features(), but currently obsolete.
clustResults_file = "out_cluster.csv"
pcaResults_file = "out_PCA.csv"
figSaveName = "heatmap_dend_feat.tiff"


# define the parameters for clustering, PCA, and file saving
significance_level=0.85 # for getting top contributing PCs
pval = 0.8 #pItem in clustering in cluster_features()
clustAlg = "hc" # define clust algorithm ("hc" for heirarchical, kmeans, etc.)
nk = 8 # max num of k clusters to consider in cluster_features()
nkdata = 5 # max num of k clusters to consider in cluster_data
assertClust = 2
nclust = 3 # input into fancy_pca()
nrep = 50 # num reps in clustering in cluster_features()
dowrite = 0 # flag to write file in fancy_pca()

  
### Run the functions with the parameters defined above!

# load the data
featureLoaded <- load_data(feature_file,feature_sheet) #file name, sheet name
featureData = featureLoaded[[2]] # data frame
feature_data = featureLoaded[[1]] # raw loaded data (not as a data frame)

# run PCA on the feature data using global_pca()
{pca <- global_pca(featureData,significance_level)
  features_pca<-pca[[1]] 
  pca_score<-pca[[2]]
  global_pca_summary<-pca[[3]]}

result_features <- cluster_features(featureData, clustAlg, nrep, nk, significance_level)
new_pca_score <- fancy_pca(featureData, result_features, pval, nclust,dowrite,fancy_pca_file)

# cluster the cases into radiomic phenotypes using cluster_data()
#data.frame(global_pca_summary$x) to cluster on PCs
{cluster_result<-cluster_data(featureData,nkdata,"ward.D2","euclidean", figSaveDir,figSaveName, assertClust) #build heatmap; relies on consensus_sig
  out_cluster<-as.data.frame(cluster_result[[1]])
  out_heatmap<-cluster_result[[2]]
  consClustResult = cluster_result[[3]]
  sig_result = cluster_result[[4]]}

feature_index = out_heatmap$rowInd
feature_names = colnames(featureData)
cluster_table = cbind(feature_names,feature_index)

### Write the clustering results to a csv
write.csv(out_cluster, clustResults_file, row.names = TRUE)
write.csv(pca_score, pcaResults_file, row.names = TRUE)

#extras
contrib_results=contributing_features(features_pca, pca_score,featureData, contrib_pca_file) #look at what features are contributing to 
PC_contrib = contrib_results[[1]]
PC_contrib_summary = contrib_results[[2]]
colCluster_pruned = cluster_result[[1]]
colCluster_pca = new_pca_score[1]

#dendrogram(colCluster_pruned,colCluster_pca)
# Plot the data w.r.t the top two PCs, colored by radiomic phenotype
dtp <- data.frame('Cluster' = out_cluster$cutree, features_pca$x[,1:3]) # the first two componets are selected (NB: you can also select 3 for 3D plottings or 3+)
ggplot(data = dtp) + 
  geom_point(aes(x = PC1, y = PC2, col = Cluster)) + 
  theme(legend.title = element_text(color = "black", size = 15),
  legend.text = element_text(color = "black")) + 
  scale_color_gradient2(low="red", mid="blue",high="green",midpoint=1.5) + theme_light()
ggsave(path = workingDir, filename = "pca_clust_scatterPlot_pureDCIS_n295.tiff", width = 6, height = 4, device='tiff', dpi=700)

# Okay let's try clustering on the PCs
#featureData, pca_score, or new_pca_score

#{cluster_result_pca<-cluster_data(pca_score,2,"ward.D2","euclidean", figSaveDir,figSaveName, assertClust) #build heatmap; relies on consensus_sig
#  out_cluster_pca<-as.data.frame(cluster_result_pca[[1]])
#  out_heatmap_pca<-cluster_result_pca[[2]]
#  consClustResult_pca = cluster_result_pca[[3]]
#  sig_result_pca = cluster_result_pca[[4]]}
