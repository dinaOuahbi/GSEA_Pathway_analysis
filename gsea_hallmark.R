

# GESA hallmarks
library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(msigdbr)
library(ggrepel)
library(stringr)

# functions
GSEA = function(gene_list, GO_file, pval) {
  set.seed(54321)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = GO_file#fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = gene_list,
                        minSize=15, ## minimum gene set size
                        maxSize=500, ## maximum gene set size
                        nperm =10000) %>% 
    as.data.frame() %>%
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGO,
                                      stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  
  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("Total pathways: Up=", total_up,", Down=",    total_down)
  
  colos = setNames(c("firebrick2", "dodgerblue2"),
                   c("Up-regulated", "Down-regulated"))
  
  g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=header)
  
  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}

plot_geneset_clusters = function( gs_results, GO_file, min.sz = 4, main="GSEA clusters"){
  
  myGO = GO_file#fgsea::gmtPathways(GO_file)
  df = matrix(nrow=nrow(gs_results), ncol = nrow(gs_results), data = 0)
  rownames(df) = colnames(df) = gs_results$pathway
  print(gs_results$pathway[i])
  
  for ( i in 1:nrow(gs_results)) {
    genesI =  unlist(myGO[names(myGO) == gs_results$pathway[i] ])
    for (j in 1:nrow(gs_results)) {
      genesJ = unlist(myGO[names(myGO) == gs_results$pathway[j] ])
      ## Jaccards distance  1 - (intersection / union )
      overlap = sum(!is.na(match(genesI, genesJ )))
      jaccards = overlap / length(unique(c(genesI, genesJ) ))
      df[i,j] = 1-jaccards
    }
  }
  
  ## Cluster nodes using dynamic tree cut, for colors
  distMat = as.dist(df)
  dendro = hclust(distMat, method = "average" )
  clust = dynamicTreeCut::cutreeDynamicTree( dendro, minModuleSize = min.sz )
  ## Note: in dynamicTreeCut, cluster 0, is a garbage cluster for things that dont cluster, so we remove it
  
  gs_results$Cluster = clust
  gs_results = gs_results[gs_results$Cluster != 0, ]
  
  ## select gene sets to label for each clusters
  bests = gs_results %>%  
    group_by( Cluster ) %>% 
    top_n(wt = abs(size), n = 1) %>% 
    .$pathway
  ## determine cluster order for plotting
  clust_ords = gs_results %>% 
    group_by( Cluster ) %>% 
    summarise("Average" = NES ) %>% 
    arrange(desc(Average)) %>% 
    .$Cluster %>% 
    unique
  
  gs_results$Cluster = factor(gs_results$Cluster, levels = clust_ords)
  
  gs_results$Label = ""
  gs_results$Label[gs_results$pathway %in% bests ] = gs_results$pathway[gs_results$pathway %in% bests ]
  gs_results$Label = str_remove(gs_results$Label, "GO_")
  gs_results$Label = tolower(gs_results$Label)
  
  g1 = ggplot(gs_results, aes(x = Cluster, y = NES, label = Label )) +
    geom_jitter( aes(color = Cluster,  size = size), alpha = 0.8, height = 0, width = 0.2 ) +
    scale_size_continuous(range = c(0.5,5)) +
    geom_text_repel( force = 2, max.overlaps = Inf) +
    ggtitle(main)
  
  
  return(g1)
}



setwd('/work/shared/ptbc/CNN_Pancreas_V2/RNAseq_TCGA_analysis/results/')
S4table = read.csv2('model8/Res_allgenes_highvslow_model8_Deseq2_RSEM.csv',header=TRUE, sep=';') %>%
  filter(X != "")

#S4table = read.csv('model8/Res_allgenes_highvslow_model8_Deseq2.csv',header=TRUE) %>%
#  filter(X != "")


#Create ranks
gene_list <- S4table$log2FoldChange
names(gene_list) = S4table$X
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(names(gene_list))]
head(gene_list)

#HALLMARK
h_gene_sets = msigdbr(species = "human", category = "H")
head(h_gene_sets)

# fixing format to work with fgsea
pathwaysH = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

# FGSEA
fgseaRes <- fgsea(pathways = pathwaysH, 
                  stats    = gene_list,
                  minSize  = 15,
                  maxSize  = 500)


# The resulting table contains enrichment scores and p-values
head(fgseaRes[order(padj), ], n=10)


#One can make an enrichment plot for a pathway
plotEnrichment(pathwaysH[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
               gene_list) + labs(title="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")


#Or make a table plot for a bunch of selected pathways:
topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(padj), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathwaysH[topPathways], gene_list, fgseaRes, 
              gseaParam=0.5) #Il s'agit de la puissance à laquelle toutes les statistiques de niveau génique d'entrée sont élevées


#From the plot above one can see that there are very similar pathways in the table
# select only independent pathways

collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], 
                                      pathwaysH, gene_list)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(pathwaysH[mainPathways], gene_list, fgseaRes, 
              gseaParam = 0.5) 


# plot les voies enrichie avec le collapse
res = GSEA(gene_list, pathwaysH, pval = 0.05)
res$Plot

# To save the results in a text format data:table::fwrite function can be used
fwrite(fgseaRes, file="/work/shared/ptbc/CNN_Pancreas_V2/RNAseq_TCGA_analysis/results/model10/fgseaRes_model10_Hallmark.txt", sep="\t", sep2=c("", " ", ""))



#To make leading edge more human-readable it can be converted using mapIdsList 


library(org.Mm.eg.db)
fgseaResMain <- fgseaRes[match(mainPathways, pathway)]
fgseaResMain[, leadingEdge := mapIdsList(
  x=org.Mm.eg.db, 
  keys=leadingEdge,
  keytype="SYMBOL", 
  column="SYMBOL")]
fwrite(fgseaResMain, file="/work/shared/ptbc/CNN_Pancreas_V2/RNAseq_TCGA_analysis/results/model10/fgseaResMain_model10_Hallmark..txt", sep="\t", sep2=c("", " ", ""))


res = GSEA(gene_list, pathwaysH, pval = 0.05)
res$Plot


# Clustering similar gene setsres$Plot
barplot(sort(gene_list, decreasing = T))

plot_geneset_clusters( gs_results = res$Results[res$Results$NES > 0, ], 
                       main = "Up-regulated GSEA clusters",
                       GO_file = pathwaysH,
                       min.sz = 4 )

plot_geneset_clusters( gs_results = res$Results[res$Results$NES < 0, ], 
                       main = "Down-regulated GSEA clusters",
                       GO_file = pathwaysH,
                       min.sz = 4 )






















