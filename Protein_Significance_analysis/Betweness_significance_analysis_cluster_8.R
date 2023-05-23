################## Analyzing the Betweness ratios
library(data.table)

intermediate_stats_df <- read.table("intermediate_stats_df.txt",sep="\t",header=T)
intermediate_stats_df <- intermediate_stats_df[,1:4]
metaloop_clusters <- read.table("Custom_CTCF_lncap_ENCFF155SPQ_8_all-nodes_overlap-metaloop_clusters.tsv",sep="\t",header=T)
metaloop_clusters <- metaloop_clusters[,c("Gene","cluster","CTCF","GWAS_Cat_prostate.carcinoma")]
colnames(metaloop_clusters)[1] <- "network"
colnames(intermediate_stats_df)[4] <- "network"

size_clusters <- table(metaloop_clusters$cluster)

merge_tab <- merge(intermediate_stats_df,metaloop_clusters,by = "network")

all_genes <- sort(unique(intermediate_stats_df$Gene))

# You need a table -per cluster- where each gene has, mean ratio in /out,

list_cluster_tabs <- list()

for(i in 1:8){
  red_cluster <- merge_tab[merge_tab$cluster==i,]
  non_red_cluster <- merge_tab[merge_tab$cluster!=i,]
  cuantas_in <- length(unique(red_cluster$network))
  mean_bet_in <- c()
  mean_bet_out <- c()
  mean_deg_in <- c()
  mean_deg_out <- c()
  ratio_bet <- c()
  ratio_deg <- c()
  
  Sys.time()
  for(j in 1:length(all_genes)){
    cuales <- red_cluster[red_cluster$Gene==all_genes[j],]
    non_cuales <- non_red_cluster[non_red_cluster$Gene==all_genes[j],]
    cuantas_veces <- nrow(cuales)
    non_cuantas_veces <- nrow(non_cuales)
    mean_betweeneess_in <- mean(cuales$ll_btw)
    if(cuantas_veces==0){
      mean_betweeneess_in <- 0
    }
    mean_bet_in <- c(mean_bet_in,mean_betweeneess_in)
    mean_degree_in <- mean(cuales$ll_deg)
    if(cuantas_veces==0){
      mean_degree_in <- 0
    }
    mean_deg_in <- c(mean_deg_in,mean_degree_in)
    mean_betweeneess_out <- mean(non_cuales$ll_btw)
    if(non_cuantas_veces==0){
      mean_betweeneess_out <- 0
    }
    mean_bet_out <- c(mean_bet_out,mean_betweeneess_out)
    mean_degree_out <- mean(non_cuales$ll_deg)
    if(non_cuantas_veces==0){
      mean_degree_out <- 0
    }
    mean_deg_out <- c(mean_deg_out,mean_degree_out)
    ratio_betweeness <- mean_betweeneess_in / mean_betweeneess_out
    ratio_bet <- c(ratio_bet,ratio_betweeness)
    ratio_degree <- mean_degree_in / mean_degree_out
    ratio_deg <- c(ratio_deg,ratio_degree)
  }
  Sys.time()
  
  cluster_tab <- data.frame(all_genes,mean_bet_in,mean_bet_out,ratio_bet,mean_deg_in,mean_deg_out,ratio_deg)
  
  list_cluster_tabs[[i]] <- cluster_tab
  #
}

# Remove NA rows in all tabs
 list_cluster_tabs_clean <- lapply(list_cluster_tabs,na.omit)

#Round up function for plotting
round.off <- function (x, digits=0)
{
  posneg = sign(x)
  z = trunc(abs(x) * 10 ^ (digits + 1)) / 10
  z = floor(z * posneg + 0.5) / 10 ^ digits
  return(z)
}

list_plots <- list()
for(q in 1:length(list_cluster_tabs_clean)){
  top_5 <- list_cluster_tabs_clean[[q]][list_cluster_tabs_clean[[q]][,4] > quantile(list_cluster_tabs_clean[[q]][,4],0.95),]
  vec_fin <- top_5$ratio_bet
  names(vec_fin) <- top_5$all_genes
  vec_fin <- sort(vec_fin)
  wolf <- round.off(max(vec_fin))
  
  name_plot <- paste0("betweeness_barplot_cluster_",q,".pdf")
  cluster_name <- paste0("Cluster ",q)
  pdf(name_plot,height = 7,width = 10)
  plotito <- barplot(vec_fin,horiz = F,col = "cyan",las=2,ylim= c(0,wolf),cex.names = 0.6,main=paste0("Top 5% Proteins with highest Betweeness ratio for ",cluster_name))
  dev.off()
}


############### Randomizations

cluster_to_rand <- 8
metaloop_clusters_2 <- metaloop_clusters
list_randomized <- list()

set.seed(2022)
for(k in 1:1000){
  print(Sys.time())
  print(k)
  red_cluster <- merge_tab[merge_tab$cluster==cluster_to_rand,]
  non_red_cluster <- merge_tab[merge_tab$cluster!=cluster_to_rand,]
  cuantas_in <- length(unique(red_cluster$network))
  metaloop_clusters_2[,2] <- 5
  subsetillo <- sample(1:nrow(metaloop_clusters),cuantas_in,replace=F)
  metaloop_clusters_2[subsetillo,2] <- cluster_to_rand
  new_merge_tab <- merge(intermediate_stats_df,metaloop_clusters_2,by = "network")
  red_cluster <- new_merge_tab[new_merge_tab$cluster==cluster_to_rand,]
  non_red_cluster <- new_merge_tab[new_merge_tab$cluster!=cluster_to_rand,]
  
  mean_bet_in <- c()
  mean_bet_out <- c()
  mean_deg_in <- c()
  mean_deg_out <- c()
  ratio_bet <- c()
  ratio_deg <- c()
  
  
  for(u in 1:length(all_genes)){
    cuales <- red_cluster[red_cluster$Gene==all_genes[u],]
    non_cuales <- non_red_cluster[non_red_cluster$Gene==all_genes[u],]
    cuantas_veces <- nrow(cuales)
    non_cuantas_veces <- nrow(non_cuales)
    mean_betweeneess_in <- mean(cuales$ll_btw)
    if(cuantas_veces==0){
      mean_betweeneess_in <- 0
    }
    mean_bet_in <- c(mean_bet_in,mean_betweeneess_in)
    mean_degree_in <- mean(cuales$ll_deg)
    if(cuantas_veces==0){
      mean_degree_in <- 0
    }
    mean_deg_in <- c(mean_deg_in,mean_degree_in)
    mean_betweeneess_out <- mean(non_cuales$ll_btw)
    if(non_cuantas_veces==0){
      mean_betweeneess_out <- 0
    }
    mean_bet_out <- c(mean_bet_out,mean_betweeneess_out)
    mean_degree_out <- mean(non_cuales$ll_deg)
    if(non_cuantas_veces==0){
      mean_degree_out <- 0
    }
    mean_deg_out <- c(mean_deg_out,mean_degree_out)
    ratio_betweeness <- mean_betweeneess_in / mean_betweeneess_out
    ratio_bet <- c(ratio_bet,ratio_betweeness)
    ratio_degree <- mean_degree_in / mean_degree_out
    ratio_deg <- c(ratio_deg,ratio_degree)
  }
  cluster_tab <- data.frame(all_genes,mean_bet_in,mean_bet_out,ratio_bet,mean_deg_in,mean_deg_out,ratio_deg)
  
  list_randomized[[k]] <- cluster_tab
}

list_rand_binary_ratio <- lapply(list_randomized,function(x) as.numeric(list_cluster_tabs[[cluster_to_rand]][,4]) <= x[,4])
list_rand_binary_tab_ratio <- do.call(rbind,list_rand_binary_ratio)
pvalue_ratio <- colSums(list_rand_binary_tab_ratio) / 1000

list_rand_binary_degree <- lapply(list_randomized,function(x) as.numeric(list_cluster_tabs[[cluster_to_rand]][,7]) <= x[,7])
list_rand_binary_tab_degree <- do.call(rbind,list_rand_binary_degree)
pvalue_degree <- colSums(list_rand_binary_tab_degree) / 1000

list_cluster_tabs[[8]] <- cbind(list_cluster_tabs[[8]],pvalue_ratio,pvalue_degree)
double_sig <- list_cluster_tabs[[8]][intersect(which(pvalue_ratio < 0.05),which(pvalue_degree<0.05)),]
write.table(x=double_sig,"significant_proteins_cluster_8.txt",sep = "\t",row.names = F)
write.table(x=list_cluster_tabs[[8]],"full_betweness_proteins_cluster_8.txt",sep = "\t",row.names = F)

  
# Remove NA rows in all tabs
#list_randomized_clean <- lapply(list_cluster_tabs,na.omit)
