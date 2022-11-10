#libraries install and load
source("packages_installer.R")
# helper functions
source("functions.r")

####################
#### To run:
#### Rscript reconstruct_epins.R \
#### src/EPIN_reconstruction_data/LNCaP/E-P_loops_refined-regions.tsv \
#### src/EPIN_reconstruction_data/LNCaP/fimo_promoters \
#### src/EPIN_reconstruction_data/LNCaP/fimo_enhancers
####################


args <- commandArgs(trailingOnly = TRUE)


ppi_filter = "filtered"
## ppi_filter: filtered , filtered Protein Proteins interactions to be LNCaP-specific and expressed in LNCaP cell line )
## ppi_filter: unfiltered , unfiltered Protein Proteins interactions include all
cell_line = "LNCaP" 
data_folder = file.path("src/EPIN_reconstruction_data", cell_line)
general_data_folder = file.path("src/EPIN_reconstruction_data")
output_folder = paste0("outputs_", cell_line, "_EPINS")
dir.create(output_folder,  showWarnings = FALSE)
dir.create(file.path(output_folder, "tables"),  showWarnings = FALSE)
dir.create(file.path(output_folder, "tables", ppi_filter),  showWarnings = FALSE)
dir.create(file.path(output_folder, "tables", ppi_filter ,"EP_graph_edges"), showWarnings = FALSE)

loop_file_path = args[1]
promoters_fimo_file_path = args[2]
enhancers_fimo_file_path = args[3]

## variables
FPKM_threshold = 0.003 # both replicates to be above
sliding_window = 100  # sliding window for filtering fimo
FIMO_pvalue = 1e-4
number_intermediate_nodes = 1

##### READ read counts 
expression <- as.data.frame(read.table(file.path(data_folder, "Cuff_Gene_Counts.txt"), header = T)) %>%
  filter(LNCaP_1 > FPKM_threshold & LNCaP_2 > FPKM_threshold)

#### READ SNPS
paintor_gwas_SNPs <- as.data.frame(read.table( file.path(data_folder, "paintor_1causals.txt"), header = TRUE, stringsAsFactors = F)) %>%
  select(-c(BP, A0, A1, Z, BETA, region, index_rsid, n, SE, N_CONTROLS, N_CASES, ID, Posterior_Prob, ch_pos,N)) %>%
  rename(SNP_start = start, SNP_stop = stop)

#### READ PPIs. Downloaded PPIS from IID database
ppi_list <- as.data.frame(read.delim(file.path(general_data_folder, "human_annotated_PPIs_IID.txt"),  header = TRUE, stringsAsFactors = F)) %>% 
  select(uniprot1, uniprot2, symbol1, symbol2, prostate, prostate.carcinoma, prostate.cancer, nucleus, evidence.type, drug.targets, methods) %>%
  filter(grepl(pattern = "exp", evidence.type)) %>%
  filter(nucleus == 1 ) %>%
  mutate(num_methods = str_count(methods, ";") + 1) %>%
  filter(num_methods >= 2) %>%
  rename(to = symbol2, from = symbol1 )

if(ppi_filter == "filtered")
{
  ppi_list <- ppi_list %>%
    filter(prostate.carcinoma == 1 | prostate.cancer == 1 | prostate == 1 ) %>%
    filter(to %in% expression$Gene_ID & from %in% expression$Gene_ID ) %>%
    select(-methods) %>%
    select(from , to)
}else{
  ppi_list <- ppi_list %>%
    select(-methods) %>%
    select(from , to)
}


ppi_net <- simplify(graph_from_data_frame(ppi_list %>% select(from, to), directed=F))

#REEAD PROMOTER regions
promoters <- as.data.frame(read.table(file.path(data_folder, "genepos.txt"), sep = "\t", header = F)) %>%
  rename(seqnames = V1, start = V2, end = V3, gene = V4) %>%
  mutate(seqnames = paste("chr", seqnames, sep = ""))


## READ the refined enhancer anchors
#loop_file = file.path(data_folder,  "E-P_loops_refined-regions.tsv")
loop_file = loop_file_path
loops <- as_tibble(read.table(loop_file, comment.char="", header = T,  sep = "\t", stringsAsFactors = F)) %>%
  mutate(enhancer_anchor_id = paste(chromosome.Enhancer.bin, start.Enhancer.bin, end.Enhancer.bin, sep = "_"),
         promoter_anchor_id = paste(chromosome.Promoter, start.Promoter, end.Promoter, sep = "_"),
         loopid = paste(enhancer_anchor_id, promoter_anchor_id, sep =  "_"),
         promoter_gene = gene.name) %>%
  distinct(gene.name, promoter_gene, loopid, enhancer_anchor_id, promoter_anchor_id, start.Enhancer.bin, chromosome.Enhancer.bin, end.Enhancer.bin, chromosome.Promoter , start.Promoter, end.Promoter) 



##### READ CTCF binding sites in LCNaP cell line
chipseq_CTCF <- as.data.frame(read.table(file.path(data_folder, "CTCF_lncap_ENCFF155SPQ.bed")))



make_graph_per_gene <- function(loop, number_intermediate_nodes)
{
  ptm=proc.time()
  my_gene <- unique(loop$promoter_gene)

  ppi_sp <- c()
  p_dbp <- data.frame()
  e_dbp <- data.frame()
  anchor_binding_sites_wSNPs <- data.frame()
  
  my_enhancers <- loop %>% filter(gene.name == my_gene)
  if(dim(my_enhancers)[1] == 0){return(data.frame())}
  
  ##### PROMOTERS
  my_promoter <- promoters %>% filter(gene == my_gene)
  if(dim(my_promoter)[1] == 0){return(data.frame())}
  
  
  chromo_promo = str_replace(my_promoter$seqnames, "chr", "")
  #p_fimo_file <- file.path(data_folder, "fimo_promoters/")
  p_fimo_file <- promoters_fimo_file_path
  promoter = paste(chromo_promo,  my_promoter$start, my_promoter$end, sep = "_")
  p_fimo <- read_fimo(p_fimo_file, promoter, 0)
  
  ## add the chipseq_info of CTCF
  add_ctcf_chipseq_info <- function(fimo , chr , start, end)
  { 
    chipseq <- chipseq_CTCF %>% filter(V1 == chr &  V2 >= start &  V3 <= end) %>% 
      select(V1, V2, V3) %>% mutate(Chipseq_CTCF = 1) %>%
      rename(ctcf_start = V2, ctcf_chr = V1, ctcf_end = V3) %>% 
      as.data.table()
    
    if(dim(chipseq)[1] == 0){return(fimo %>% mutate(Chipseq_CTCF = NA))}
    setkey(chipseq, ctcf_chr, ctcf_start, ctcf_end)
    p_fimo <- as.data.frame(foverlaps(as.data.table(p_fimo), chipseq, by.x=c("chr_DNA", "start", "stop"), type="any")) %>%
      select(-ctcf_start,-ctcf_end )
    return(p_fimo)
  }
  
  p_fimo <- add_ctcf_chipseq_info(p_fimo, unique(my_promoter$seqnames), min(my_promoter$start),  min(my_promoter$end))

  p_dbp_current <- data.frame(TF = unique(p_fimo$motif_alt_id))
  p_dbp_current <- p_dbp_current %>% filter(TF %in% V(ppi_net)$name)
  if(dim(p_dbp_current)[1] == 0){return(data.frame())}
  p_dbp <- rbind(p_dbp, p_dbp_current)
  
  if(dim(p_fimo)[1] == 0){return(data.frame())}
  
  #keep the binding sites in P that cary SNPs..
  anchor_binding_sites_wSNPs <- rbindlist(list(anchor_binding_sites_wSNPs,
                                               p_fimo  %>%
                                                 filter(!is.na(rsid)) %>%
                                                 mutate(anchor = promoter,
                                                        type = "PROM") %>%
                                                 select(anchor, motif_alt_id, rsid)))
  
  
  ##### ENHANCER(S)
  
  for(enhancer in unique(my_enhancers$enhancer_anchor_id))
  {
    print(paste("enhancer", enhancer))
    
    #e_fimo_file <- file.path(data_folder, "fimo_enhancers/")
    e_fimo_file <- enhancers_fimo_file_path
    e_fimo <- read_fimo(e_fimo_file, enhancer, 0)
    
    if(dim(e_fimo)[1] == 0){next}
    
    #keep the binding sites in E that cary SNPs..
    anchor_binding_sites_wSNPs <- rbindlist(list(anchor_binding_sites_wSNPs,
                                                 e_fimo  %>%
                                                   filter(!is.na(rsid)) %>%
                                                   mutate(anchor = enhancer,
                                                          type = "ENHA") %>%
                                                   select(anchor, motif_alt_id, rsid)))
                                                   
                                          
    
    
    
    e_fimo <- add_ctcf_chipseq_info(e_fimo, str_split(enhancer, "_")[[1]][1], str_split(enhancer, "_")[[1]][2], str_split(enhancer, "_")[[1]][3])
    

    e_dbp_current <- data.frame(TF = unique(e_fimo$motif_alt_id), position = enhancer)
    e_dbp_current <- e_dbp_current %>% filter(TF %in% V(ppi_net)$name)
    if(dim(e_dbp_current)[1] == 0){next}
    e_dbp <- rbind(e_dbp, e_dbp_current)
    
    ppi_sp <- c()
    for(i in p_dbp_current$TF){ppi_sp <- c(ppi_sp, all_shortest_paths(ppi_net, from=i, to=as.character(e_dbp_current$TF))$res)}
    ppi_sp <- ppi_sp[lapply(ppi_sp, length) <= number_intermediate_nodes + 2 ]
    
    if(length(ppi_sp)  > 0 )
    {
      if(exists("net")) { net <- graph.union(net, find_network(ppi_sp, e_dbp_current, p_dbp)) }
      if(!exists("net")) { net <-  find_network(ppi_sp, e_dbp_current, p_dbp) }
    }
  }

  
  anchor_binding_sites_wSNPs = data.frame(anchor_binding_sites_wSNPs)
  
  if(dim(e_dbp)[1] == 0){return(data.frame())}
  
  #### Delete an intermediate if it is just connecting two DBPs of the same anchor
  ValidIntermediate <- function(i){
    my_neighbors <- paste(names(i),collapse ="_" )
    if( str_count(my_neighbors, "ENHA") < length(i)  && str_count(my_neighbors, "PROM") < length(i))
    {
      return(TRUE) 
    }else{
      return(FALSE)
    }
  }
  intermediates <- V(net)$name[!grepl("ENHA_", V(net)$name) & !grepl("PROM_", V(net)$name) ]
  
  if(length(intermediates) > 0){
  
    neig <- adjacent_vertices(net, intermediates)
    valid_intermediates <- intermediates[which(sapply(neig, ValidIntermediate))]
    remove_intermediates <- intermediates[intermediates %nin% valid_intermediates ]
  
    net <- delete.vertices(net, remove_intermediates)
  
  }
  net <- delete.vertices(net, degree(net)==0)
  

  
  as_long_data_frame(net) %>% 
    select(3,4) %>%
    write.table(file.path(output_folder,  "tables", ppi_filter, "EP_graph_edges", paste(my_gene, "txt", sep = ".")), quote = F, row.names = F, sep = "\t", col.names = F)
  
  return(anchor_binding_sites_wSNPs)
  print(proc.time() - ptm)
  
}


all_anchors_binding_sites_wSNPs = data.frame()
#for(gene in unique(loops$promoter_gene))
for(gene in c("MYC"))
{
  print(gene)
  loop <- loops %>% filter(promoter_gene == gene)
  all_anchors_binding_sites_wSNPs = rbindlist(list(all_anchors_binding_sites_wSNPs, 
                                                    make_graph_per_gene(loop, number_intermediate_nodes)))
}

data.frame(all_anchors_binding_sites_wSNPs ) %>% 
  filter(type == "ENHA") %>%
  write.table(file.path(output_folder, "tables", ppi_filter, "binding_sites_snp_enha.txt"), quote = F, row.names = F, sep = "\t")

data.frame(all_anchors_binding_sites_wSNPs ) %>% 
  filter(type == "PROM") %>%
  write.table(file.path(output_folder, "tables", ppi_filter, "binding_sites_snp_prom.txt"), quote = F, row.names = F, sep = "\t")



