# GSEA_Pathway_analysis

## Methods

#### Parameters
The data from the RNA seq analysis on the three models of interest 
    Res_allgenes_highvslow_model8_Deseq2_RSEM.csv
    
    
The genes are classified according to their CFL. We then get a list of genes with each a CFL value associated.
Genes sets Hallmark (species 'hs' category 'H') composed of 50 pathways recognized in homo sapiens. 
Fgsea : a cluster will have a minimum of 15 genes and a maximum of 500 genes.

#### Collapse list of enriched pathways to independent ones
This is a method that reduces similar paths and therefore returns a list of enriched main paths 

#### Clustering similar gene sets
Sometimes, even after grouping gene sets, there may still be a fairly large number of pathways to decipher. To facilitate the interpretation of a large number of biological pathways, we can group similar gene sets together (based on the overlap of the genes they contain). This will allow us to see if many sets of genes related to the same function are similarly affected and provide a more comprehensive view of pathway deregulation.

Use of distance Jaccards (set1 union set2)
Cluster 0 = cluster poubelle

## Results
HIGH : worst pronostic
![image](https://github.com/dinaOuahbi/GSEA_Pathway_analysis/blob/main/KM_high_low.png)
#### Total pathways
![image](https://github.com/dinaOuahbi/GSEA_Pathway_analysis/blob/main/total_paths.png)
#### Selection of independent channels
![image](https://github.com/dinaOuahbi/GSEA_Pathway_analysis/blob/main/select_indep_channels.png)
#### Clustering similar gene sets - Evaluation
![image](https://github.com/dinaOuahbi/GSEA_Pathway_analysis/blob/main/down_reg_paths.png)
![image](https://github.com/dinaOuahbi/GSEA_Pathway_analysis/blob/main/up_reg_paths.png)




####### Thanks for following (OD)
