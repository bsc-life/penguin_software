README File for node importance significance analysis

### Input files (Produced by Alex)
- Custom_CTCF_lncap_ENCFF155SPQ_8_all-nodes_overlap-metaloop_clusters.tsv
- intermediate_stats_df.txt

### Scripts (For running at a SLURM-based cluster machine)
- Base randomization scripts: Betweness_significance_analysis_cluster_*.R . 1 per cluster
- Job text files: Betw_sig_analysis_cluster_*.cmd . 1 per cluster. Expected runtime ~ 4h. Each job asks for 1 cpu of 1 node.

### Outputs
- Full results per cluster: full_betweness_proteins_cluster_*.cmd 1 per cluster
- Double Significant rows from the full results (significant for both degree and betweenness -pvalue < 0.05): significant_proteins_cluster_*.txt
- PLEASE NOTE THAT AS INPUT FOR THE ENRICHMENT ANALYSIS WITH G:PROFILER WE CONSIDER DOUBLE SIGNIFICANT ROWS WITH P-VALUE < 0.01.

### Columns from output files (for cluster output i) (for a given row j) :
-all_genes : The node protein j
-Mean_bet_in: Mean betweenneess of the protein j inside networks belonging to Cluster i
-Mean_bet_out: Mean betweenneess of the protein j inside networks NOT belonging to Cluster i
-ratio_bet: (Mean_bet_in + 1) / (Mean_bet_out + 1)
-Mean_deg_in: Mean interaction degree of the protein j inside networks belonging to Cluster i
-Mean_deg_out: Mean interaction degree of the protein j inside networks NOT belonging to Cluster i
-ratio_deg: (Mean_deg_in + 1) / (Mean_deg_out + 1)
-pvalue_ratio: Probability of finding an equal or higher BETWEENNESS ratio for protein j than the real value (ratio_bet) out of 1000 randomly generated clusters of the same size as Cluster i
-pvalue_degree: Probability of finding an equal or higher DEGREE ratio for protein j than the real value (ratio_deg) out of 1000 randomly generated clusters of the same size as Cluster i
