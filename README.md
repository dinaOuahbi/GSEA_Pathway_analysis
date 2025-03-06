# **GSEA_Pathway_Analysis**  
**Gene Set Enrichment Analysis (GSEA) of Pathways in High vs. Low Prognosis Models**  

## **Methods**  

### **Parameters**  

- **Data File**: RNA-seq data for the three models of interest:  
  - `Res_allgenes_highvslow_model8_Deseq2_RSEM.csv`  

- **Gene Classification**: Genes are classified based on their **CFL** value. A list of genes is generated, each associated with a CFL score.  
- **Gene Sets**: The **Hallmark gene sets** for Homo sapiens (species 'hs', category 'H') consist of **50 pathways**.  
- **Fgsea**: A cluster is composed of a minimum of 15 genes and a maximum of 500 genes.

### **Methods Overview**  

- **Collapse List of Enriched Pathways**: This step reduces similar pathways, returning a list of independent enriched pathways.  
- **Clustering Similar Gene Sets**: To interpret a large number of biological pathways, we group similar gene sets based on gene overlap. This approach reveals if gene sets related to the same function are similarly affected, providing a clearer understanding of pathway deregulation.  
  - **Distance Measure**: Jaccard distance is used to measure overlap between gene sets (Set1 ∪ Set2).  
  - **Cluster 0**: Represents a "trash" cluster that is discarded.

---

## **Results**  

### **Prognosis Groups**  
- **HIGH**: Worst prognosis  
![Kaplan-Meier Plot](https://github.com/dinaOuahbi/GSEA_Pathway_analysis/blob/main/KM_high_low.png)

### **Total Pathways**  
![Total Pathways](https://github.com/dinaOuahbi/GSEA_Pathway_analysis/blob/main/total_paths.png)

### **Selection of Independent Pathways**  
![Selection of Independent Pathways](https://github.com/dinaOuahbi/GSEA_Pathway_analysis/blob/main/select_indep_channels.png)

### **Clustering of Similar Gene Sets**  
- **Downregulated Pathways**  
![Downregulated Pathways](https://github.com/dinaOuahbi/GSEA_Pathway_analysis/blob/main/down_reg_paths.png)

- **Upregulated Pathways**  
![Upregulated Pathways](https://github.com/dinaOuahbi/GSEA_Pathway_analysis/blob/main/up_reg_paths.png)
