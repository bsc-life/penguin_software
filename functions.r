
find_network <- function(ppi_sp_local, e_dbp_local, p_dbp_local)
{
  e_dbp_local <- unique(e_dbp_local)
  p_dbp_local <- unique(p_dbp_local)
  bin <- unique(e_dbp_local$position)
  
  if(length(c(e_dbp_local$TF,p_dbp_local$TF)) == 0){return(data.frame())}
  
  ppi_spsb = graph.union(lapply(ppi_sp_local, function(x) induced.subgraph(graph = ppi_net, x)))
  
  #here delete the connections in the same side.
  
  p_dbp_only <- V(ppi_spsb)$name[V(ppi_spsb)$name %in% p_dbp_local$TF[p_dbp_local$TF %nin% e_dbp_local$TF]]
  e_dbp_only <- V(ppi_spsb)$name[V(ppi_spsb)$name %in% e_dbp_local$TF[e_dbp_local$TF %nin% p_dbp_local$TF]]
  intermediates <- V(ppi_spsb)$name[V(ppi_spsb)$name %nin% union(e_dbp_local$TF, p_dbp_local$TF)]
  both <- V(ppi_spsb)$name[V(ppi_spsb)$name %in% intersect(p_dbp_local$TF, e_dbp_local$TF)]
  
  
  ####################
  ##########  E and E Connections to intermediate nodes
  p_network <- induced.subgraph(ppi_spsb, c(intermediates, both, p_dbp_only)) # this includes only nodes in the P and in the intermediate
  e_network <- induced.subgraph(ppi_spsb, c(intermediates, both, e_dbp_only)) # this includes only nodes in the E and in the intermediate
  
  ####################
  ##########  Direct E-P Connections for Common nodes.
  #First extract the direct network of connections
  direct_network <- induced.subgraph(ppi_spsb, both)
  # Create network that contains the direct connections from common nodes. Remember we are trying to duplicate the common nodes in the graph
  direct_net_df <- data.frame()
  for(i in E(direct_network)){
    from <- rbind(paste0("PROM_", ends(direct_network, E(direct_network)[i])[1]),
                  paste0("PROM_", ends(direct_network, E(direct_network)[i])[2]))
    to <- rbind(paste0("ENHA_", bin,"_", ends(direct_network, E(direct_network)[i])[2]),
                paste0("ENHA_", bin,"_", ends(direct_network, E(direct_network)[i])[1]))
    
    direct_net_df <- rbind(direct_net_df,
                           data.frame(cbind(from, to)))
    
  }
  ####################
  ##########  Direct E-P Connections for non-Common nodes.
  #First extract the netwrork of connections
  direct_network2 <- induced.subgraph(ppi_spsb, c(p_dbp_only,e_dbp_only))
  V(direct_network2)$name[V(direct_network2)$name %in% p_dbp_only] = paste0("PROM_", V(direct_network2)$name[V(direct_network2)$name %in% p_dbp_only])
  V(direct_network2)$name[V(direct_network2)$name %in% e_dbp_only] = paste0("ENHA_", bin,"_",  V(direct_network2)$name[V(direct_network2)$name %in% e_dbp_only])
  ####################
  
  # Change the names of the P and E network
  V(p_network)$name[V(p_network)$name %in% c(p_dbp_only, both)] = paste0("PROM_", V(p_network)$name[V(p_network)$name %in% c(p_dbp_only, both)])
  V(e_network)$name[V(e_network)$name %in% c(e_dbp_only, both)] = paste0("ENHA_", bin,"_", V(e_network)$name[V(e_network)$name %in% c(e_dbp_only, both)])
  
  if(dim(direct_net_df)[1] > 0) # if there is direct network
  {
    direct_network <- simplify(graph_from_data_frame(direct_net_df, directed=F))
    g_merge <- graph.union(p_network, e_network, direct_network, direct_network2) # make union on the symbolic names
  }else
  {
    g_merge <- graph.union(p_network, e_network, direct_network2) # make union on the symbolic names
  }
  
  g_merge <- delete.vertices( g_merge, degree(g_merge)==0)
  return(g_merge)
}


filter_fimo <- function(my_fimo, width)
{
  
  step_promoter <- data.frame()
  
  for(window in seq(1, width, 25))
  {
    n_t = my_fimo %>% filter(start >= window & start <= window + sliding_window & stop >= window & stop <= window + sliding_window) %>%
      distinct(motif_alt_id) %>% dim()
    step_promoter <- rbind(step_promoter,
                           data.frame(start_window=window,
                                      n_tfs = n_t[1] ))
  }
  step_promoter <- step_promoter %>% filter(n_tfs >= mean(n_tfs) + sd(n_tfs))
  filter_fimo_df <- as.data.frame(sqldf("select * from step_promoter
                                        left join my_fimo on
                                        step_promoter.start_window <= my_fimo.start and
                                        step_promoter.start_window + 100 >=  my_fimo.start and
                                        step_promoter.start_window <= my_fimo.stop and
                                        step_promoter.start_window + 100 >=  my_fimo.stop
                                        ")) %>% select(-start_window, -n_tfs) %>% unique
  return(filter_fimo_df)
  
}


file_name <- function(x){
  elm <- strsplit(x,"_")[[1]]
  s <- paste0("chr",elm[1],"_",elm[2],"-",elm[3],".tsv")
  return(s)
}


read_fimo <- function(fimo_dir , anchor, step )
{
  if(!file.exists(paste0(fimo_dir,file_name(anchor))) ){return(data.frame())} 
  if(length(readLines(paste0(fimo_dir,file_name(anchor)))) == 4){return(data.frame())} 
  my_df <- read.table(paste0(fimo_dir,file_name(anchor)), header=T,sep='\t', stringsAsFactors = F)
  
  my_df <- my_df %>%
    filter(start>=step ) %>%
    filter(  p.value <= FIMO_pvalue ) %>%
    filter(  q.value <= 0.05 ) %>%
    filter(strand == "+") %>%
    mutate(DNA_segment2 = sequence_name) %>%
    separate(DNA_segment2, c("chr_DNA", "start_DNA_stop_DNA"), ":") %>%
    separate(start_DNA_stop_DNA, c("start_DNA", "stop_DNA"), "-") %>%
    mutate(start_DNA = as.numeric(start_DNA),
           stop_DNA = as.numeric(stop_DNA),
           anchor = anchor,
           start = start +start_DNA,
           stop = stop + start_DNA)
  
  #I filter to make the merge faster
  paintor_gwas_tmp = paintor_gwas %>% filter(CHR == unique(my_df$chr_DNA) )
  
  my_df1 <- as.data.frame(sqldf("select * from my_df
                                left join paintor_gwas_tmp on
                                paintor_gwas_tmp.CHR == my_df.chr_DNA and
                                paintor_gwas_tmp.SNP_start >= my_df.start and
                                paintor_gwas_tmp.SNP_stop <= my_df.stop+ 1")) %>%
    select(-c(CHR, SNP_stop, Pvalue)) %>%
    mutate(SNP_in_Motif = ifelse(!is.na(rsid),1,0))
  
  return(my_df1)
}



retrieve_DBP_after_comparison <- function(myfimo, p)
{
  mydf <- myfimo %>%
    group_by(motif_alt_id) %>%
    summarise(n_motifs = n_distinct(start),
              SNP_in_Motif = max(SNP_in_Motif, na.rm = T) ) %>%
    inner_join(TF_stats, by=c("motif_alt_id")) %>%
    group_by(motif_alt_id) %>%
    mutate(l = length(unlist(l_motifs)), q = quantile(unlist(l_motifs), probs = p)) %>%
    filter( (n_motifs > q | SNP_in_Motif == 1) & motif_alt_id %in% V(ppi_net)$name )
  
  return(unique(mydf$motif_alt_id))
}



## return the X,Y coordinates for the vertical/horizontal layout
return_lo_linear <- function(i, x_cor, my_net)
{
  GroupV = which(V(my_net)$Group1 == i)
  names <- V(my_net)$gene_name[GroupV]
  if(i == "red") # make horizontal the common
  {
    l_o <- cbind(rescale(matrix(1:length(GroupV)), to = c(0,2)), matrix(rep(-1, length(GroupV))) , GroupV)
  }else if(i == "lightblue"){ # make vertical in the middle the intermediates
    x_noisy <- c(x_cor - 0.3, x_cor - 0.15, x_cor, x_cor+ 0.15, x_cor+0.3 )
    if(length(GroupV)%%length(x_noisy) == 0)
    {
      l_o <- cbind(matrix(rep(x_noisy,  length(GroupV)/length(x_noisy)) ) , rescale(matrix(1:length(GroupV)), to = c(1,40)), GroupV)
    }else
    {
      l_o <- cbind(matrix(c(rep(x_noisy,  length(GroupV)/length(x_noisy)) , x_noisy[1:(length(GroupV)%%length(x_noisy))])) , rescale(matrix(1:length(GroupV)), to = c(1,40)), GroupV)
    }
  }else { # make vertical the EP binders
    l_o <- cbind(matrix(rep(x_cor, length(GroupV))) , rescale(matrix(1:length(GroupV)), to = c(5,35)), GroupV)
  }
  return(cbind(l_o))
}
