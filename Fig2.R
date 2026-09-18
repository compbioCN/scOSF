library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggstatsplot)
library(gground)
library(ggprism)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(ggpubr)
library(ggvenn)

load("~/OSF/fib.Rdata")

dp <- cbind(fib@meta.data, fib@reductions$umap@cell.embeddings)

if(any(!c("UMAP_1", "UMAP_2") %in% colnames(dp))){
  colnames(dp)[colnames(dp) %in% c("umap_1", "umap_2")] <- c("UMAP_1", "UMAP_2")
}
unique(dp$GENE1_group)

colors <- c("#9BC6B5","#8F5F95")

C <- c("#726BAE", "#C7AED5" "#F7DBF0" "#EA945A" "#60A897" "#F5BC6E" "#86A667"
"#276D9F" "#6488B9" "#ACD48A" "#A5AA99" "#C7AED5" "#E05F48" "#5D88BF"
"#78C4D4" "#8AB6D6" "#E7C5DB" "#FCEDDC" "#FCE8E3" "#CCD7DD")
centroid <- ddply(dp, .(GENE1_group), summarise, label_x = median(UMAP_1), label_y = median(UMAP_2))
colnames(dp) <- make.unique(colnames(dp))
ggplot(dp, aes(x=UMAP_1, y=UMAP_2, color=GENE1_group)) +
  geom_point(aes(fill = GENE1_group, color = GENE1_group), 
             size = 0.5, 
             shape = 21,   
             alpha = 0.8,     
             stroke = 0.4) +  
  scale_fill_manual(values = colors) +         
  scale_color_manual(values = colors) +
  geom_text_repel(data = centroid, aes(x=label_x, y=label_y, label=GENE1_group), color='black',size=5)+  
  theme_bw()+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(family = "Arial"))+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  labs(color = "")+
  theme(
    panel.background = element_rect(fill = "transparent"),
    panel.grid = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    title = element_text(size = 16,family = "Arial"),
    text = element_text(size = 15,family = "Arial"),
    plot.title = element_text(size = 16,family = "Arial"),
    axis.title = element_text(size = 16,family = "Arial"),
    axis.text = element_text(size = 15,family = "Arial"),
    axis.text.x = element_text(size = 15,angle = 0, family = "Arial",hjust = 1,vjust = 1),
    axis.text.y = element_text(size = 15,family = "Arial"),
    legend.text = element_text(size = 15,family = "Arial"),
    legend.title = element_text(size = 16,family = "Arial"),
    legend.position = 'right',
    plot.margin = margin(t = 20, r = 10, b = 20, l = 10))
ggsave(filename = paste0("~/OSF/Fig/260504/UMAP.pdf"),device = cairo_pdf, family = "Arial",width = 7,height = 5,dpi = 300)

load("~/OSF/DEG.Rdata")

DEGs_osf <- DEGs_osf[order(DEGs_osf$p_val,DEGs_osf$avg_log2FC,decreasing = c(FALSE,TRUE)),]
DEGs_diff_osf <- subset(DEGs_osf,abs(DEGs_osf$avg_log2FC)>0.5&DEGs_osf$p_val<0.05)
DEGs_nodiff_osf <- subset(DEGs_osf,!(abs(DEGs_osf$avg_log2FC)>0.5&DEGs_osf$p_val<0.05))
DEGs_osf <- bind_rows(DEGs_diff_osf,DEGs_nodiff_osf)
DEGs_up_osf <- DEGs_osf[which(DEGs_osf$avg_log2FC >= 0.5 & DEGs_osf$p_val<0.05),]
DEGs_down_osf <- DEGs_osf[which(DEGs_osf$avg_log2FC <= -0.5 & DEGs_osf$p_val<0.05),]
DEGs_up_osf <- DEGs_up_osf[order(DEGs_up_osf$avg_log2FC,decreasing = TRUE),]
DEGs_down_osf <- DEGs_down_osf[order(DEGs_down_osf$avg_log2FC,decreasing = FALSE),]
Upvals <- rownames(DEGs_up_osf)[1:5]
Downvals <- rownames(DEGs_down_osf)[1:5]
vals <- c(Upvals,Downvals)
group<-ifelse(
  DEGs_osf$avg_log2FC<(-0.5)&DEGs_osf$p_val<0.05,'#4D4398',
  ifelse(DEGs_osf$avg_log2FC>(0.5)&DEGs_osf$p_val<0.05,'indianred1',
         '#b5b5b5'))
max <- max(abs(na.omit(DEGs_osf[,c("avg_log2FC")])))
x_max <- max*(11/10)
x_min <- -x_max
group[is.na(group)]<-'#b5b5b5'
names(group)[group=='indianred1']<-'Up'
names(group)[group=='#b5b5b5']<-'No_sig'
names(group)[group=='#4D4398']<-'Down'

pdf(
  file = paste0("~/OSF/Fig/260504/Volcano.pdf"),
  width = 5,
  height = 4.5,
  bg = "white"
)
EnhancedVolcano(DEGs_osf,
                x = "avg_log2FC",
                y = "p_val",
                lab = rownames(DEGs_osf),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize=c(ifelse(rownames(DEGs_osf) %in% vals,1.5,1)),
                labSize = 5,
                xlim = c(x_min,x_max),
                colCustom = group,
                title = paste0("OMD_Negative (UP)","          "," OMD_Positive(UP)"),
                subtitle = NULL,
                caption = NULL,
                legendPosition = "right",
                selectLab = c(Upvals,Downvals),
                xlab = bquote(~Log[2]~'fold change'),
                legendLabSize = 12,
                legendIconSize = 6,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.8,
                max.overlaps = 30
                
)+theme(
  axis.title = element_text(family = "Arial",size = 14),
  axis.text = element_text(family = "Arial"),
  legend.text = element_text(family = "Arial"),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
  plot.title = element_text(
    size = 15,
    face = "bold",
    hjust = 0.5,
    vjust = 2,
    family = "Arial",
    margin = margin(b = 10)  
  ),
  legend.title = element_blank(),                                       
  text = element_text()                        
)
ggsave("~/OSF/Fig/260504/Volcano.pdf",width = 6.5,height = 5,dpi = 300,device = cairo_pdf, family = "Arial")

##GO enrichment analysis
GO_BP_UP <- read_csv("~/OSF/GO_BP_UP.csv")
GO_data <- GO_BP_UP
GO_data <- GO_data[order(GO_data$p.adjust,decreasing = c(F)),]
GO_data <- GO_data[c(1:12),]
GO_data <- GO_data[-c(7,11),]
GO_data$Description <- str_to_sentence(GO_data$Description)
GO_data$category <- "GO:BP"
plot_df <- GO_data %>%
  dplyr::group_by(category) %>%
  dplyr::top_n(20, wt = -as.numeric(p.adjust)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(
    p.adjust = as.numeric(p.adjust),
    Count = as.numeric(Count)
  ) %>%
  dplyr::arrange(category, dplyr::desc(p.adjust)) %>% 
  dplyr::mutate(Description = factor(Description, levels = unique(Description))) %>%
  tibble::rowid_to_column('index')
width <- 0.6
plot_df <- plot_df %>%
  mutate(p.adjust = as.numeric(p.adjust),
         Count = as.numeric(Count))
rect_data <- plot_df %>%
  dplyr::group_by(Description) %>%  
  dplyr::summarise(n = n()) %>%
  dplyr::ungroup() %>%
  mutate(
    ymax = cumsum(n) + 0.4,
    ymin = lag(ymax, default = 0.4) + 0.2,
    xmin = -7, 
    xmax = -6
  )
xaxis_max <- max(-log10(plot_df$p.adjust), na.rm = TRUE) + 2
target_terms <- c("Myeloid cell differentiation") 
plot_df <- plot_df %>%
  mutate(is_target = ifelse(Description %in% target_terms, "Yes", "No"))
custom_colors1 <- c(
  "#726BAE", "#C7AED5", "#F7DBF0", "#EA945A", "#60A897", "#F5BC6E", "#86A667", 
  "#276D9F", "#6488B9", "#ACD48A" 
)
plot_df <- plot_df %>%
  mutate(
    logP = -log10(p.adjust)
  )
xaxis_max <- max(plot_df$FoldEnrichment, na.rm = TRUE) * 1.05
plot_df <- plot_df %>%
  mutate(
    Description = fct_reorder(Description, FoldEnrichment, .desc = FALSE)
  )
p_x <- -1.5
count_x <- -0.5
custom_colorsrev <- rev(custom_colors)
plot_df %>%
  ggplot(aes(x = FoldEnrichment, y = Description)) +
  geom_point(
    aes(x = p_x, fill = logP),
    shape = 21,
    size = 7,
    color = "black",
    stroke = 0.6
  ) +
  geom_text(
    aes(x = p_x, label = round(logP, 1)),
    size = 3.2,
    color = "black"
  ) +
  scale_fill_gradientn(
    colours = c("#F7FCFD", "#D0EAE7", "#A8D5CF", "#5AAEAE", "#0B525B"),
    name = expression(-log[10](p.adjust))
  ) +
  ggnewscale::new_scale_fill() +
  geom_round_col(
    aes(fill = Description),
    width = 0.6,
    alpha = 0.8
  ) +
  geom_point(
    aes(x = count_x, size = Count, fill = Description),
    shape = 21,
    color = "black",
    stroke = 0.6
  ) +
  geom_text(
    aes(x = count_x, label = Count),
    size = 3.5,
    color = "black"
  ) +
  scale_fill_manual(
    name = "Category",
    values = rev(custom_colors1)
  ) +
  scale_size_continuous(
    name = "Count",
    range = c(5, 16)
  ) +
  geom_text(
    aes(
      x = 0.05,
      label = Description,
      colour = ifelse(Description %in% target_terms, "red", "black")
    ),
    hjust = 0,
    size = 5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_colour_identity() +
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(
    limits = c(-2, xaxis_max),
    breaks = seq(0, xaxis_max, 2),
    expand = expansion(c(0, 0))
  ) +
  scale_y_discrete(
    expand = expansion(add = c(0.9, 0.5))
  )+
  labs(
    x = "Fold Enrichment",
    y = "Description"
  ) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )+
  annotate(
    "text",
    x = p_x,
    y = 0.45,
    label = "-log10(p.adjust)",
    angle = 90,
    fontface = "bold",
    size = 4,
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = count_x,
    y = 0.45,
    label = "Count",
    angle = 90,
    fontface = "bold",
    size = 4,
    hjust = 0.5
  )+
  coord_cartesian(clip = "off")
ggsave(filename = "~/OSF/Fig/260504/GO_BP_UP.pdf",width = 9.5,height = 8,dpi = 300,device = cairo_pdf, family = "Arial")

##GSEA-----
df <- bitr(rownames(DEGs_osf), 
           fromType = "SYMBOL",
           toType =  "ENTREZID",
           OrgDb = "org.Hs.eg.db") 
DEGs_osf$SYMBOL <- rownames(DEGs_osf)
need_DEG <- merge(DEGs_osf, df, by='SYMBOL')  
geneList <- need_DEG$avg_log2FC
names(geneList) <- need_DEG$ENTREZID
geneList <- sort(geneList, decreasing = T)   
GO_kk_entrez <- gseGO(geneList     = geneList,
                      ont          = "BP",  
                      OrgDb        = "org.Hs.eg.db",
                      keyType      = "ENTREZID",
                      pvalueCutoff = 0.05)   
GO_kk <- DOSE::setReadable(GO_kk_entrez, 
                           OrgDb="org.Hs.eg.db",
                           keyType='ENTREZID')#转化id
kk_gse <- GO_kk
kk_gse_entrez <- GO_kk_entrez
kk_gse_cut <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$p.adjust<0.25 & abs(kk_gse$NES)>1]
kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0,]
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]
down_gsea <- kk_gse_cut_down[tail(order(kk_gse_cut_down$NES,decreasing = T),10),]
up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES,decreasing = T),10),]
diff_gsea <- kk_gse_cut[head(order(abs(kk_gse_cut$NES),decreasing = T),10),]
up_gsea$Description
i=10  


plot_gsea_pub <- function(gsea_obj, geneSetID,
                          line_col = "#D73027",
                          title_width = 70,
                          base_size = 12) {
  gsdata <- enrichplot:::gsInfo(gsea_obj, geneSetID)
  res <- as.data.frame(gsea_obj)
  res_i <- res[res$ID == geneSetID, ]
  
  plot_title <- str_wrap(res_i$Description[1], width = title_width)
  subtitle_text <- paste0(
    "NES = ", round(res_i$NES[1], 3),
    "    p = ", signif(res_i$pvalue[1], 3),
    "    FDR = ", signif(res_i$p.adjust[1], 3)
  )
  hit_df <- gsdata %>%
    dplyr::filter(position == 1)
  pub_theme <- theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(color = "black", size = base_size + 2, face = "bold"),
      plot.title = element_text(size = base_size + 3, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "grey25", hjust = 0),
      plot.margin = margin(4, 8, 2, 8),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.major.y = element_blank()
    )
  p_es <- ggplot(gsdata, aes(x = x, y = runningScore)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey45") +
    geom_line(color = line_col, linewidth = 1.1) +
    labs(
      title = plot_title,
      subtitle = subtitle_text,
      x = NULL,
      y = "Running Enrichment Score"
    ) +
    pub_theme +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  p_hits <- ggplot(gsdata, aes(x = x)) +
    geom_tile(
      aes(y = 0, fill = geneList),
      height = 0.35
    ) +
    geom_linerange(
      data = hit_df,
      aes(x = x, ymin = 0.25, ymax = 1.25),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.35
    ) +
    scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      name = "Ranked metric"
    ) +
    scale_y_continuous(limits = c(-0.25, 1.35), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = base_size) +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 8, 0, 8)
    )
  p_metric <- ggplot(gsdata, aes(x = x, y = geneList)) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50") +
    geom_col(width = 1, fill = "grey70", color = "grey70") +
    labs(
      x = "Rank in Ordered Dataset",
      y = "Ranked List Metric"
    ) +
    pub_theme +
    theme(
      plot.margin = margin(2, 8, 4, 8)
    )
  p_es / p_hits / p_metric +
    plot_layout(heights = c(1.6, 0.35, 1))
}
plot_gsea_pub(
  gsea_obj = kk_gse,
  geneSetID = up_gsea$ID[i],
  line_col = "#D73027",
  base_size = 13
)
ggsave(
  filename = "~/OSF/Fig/260504/GSEA_TGF_beta_receptor_pathway.pdf",
  width = 8,
  height = 5.5,
  device = cairo_pdf, family = "Arial",dpi=300
)
{
  entrez_ids <- names(geneList)
  id_map <- bitr(entrez_ids, 
                 fromType = "ENTREZID",
                 toType = "SYMBOL",
                 OrgDb = "org.Hs.eg.db")
  aligned_symbols <- id_map$SYMBOL[match(entrez_ids, id_map$ENTREZID)]
  aligned_symbols[is.na(aligned_symbols)] <- entrez_ids[is.na(aligned_symbols)]
  names(geneList) <- aligned_symbols
  head(geneList, 50)
  top_gene_names <- names(geneList)
  print(head(top_gene_names,50))
  
  genes <- unique(c(top_gene_names[1:20], "COL1A1", "COL1A2"))
  cell_exp <- GetAssayData(fib, assay = "RNA", slot = "data")[genes, ]
  cell_exp_scaled <- t(scale(t(as.matrix(cell_exp))))
  group_info <- fib$GENE1_group
  avg_exp_scaled <- t(apply(cell_exp_scaled, 1, function(x) {
    tapply(x, group_info, mean)
  }))
  avg_exp_scaled[is.na(avg_exp_scaled)] <- 0
  avg_exp_scaled[avg_exp_scaled > 2] <- 2
  avg_exp_scaled[avg_exp_scaled < -2] <- -2
  col_fun <- colorRamp2(
    c(-1, 0, 1),
    c("#4575B4", "white", "#D73027")
  )
  pdf("~/OSF/Fig/260504/heatmap_GENE1_group.pdf", width = 2.5, height = 8)
  Heatmap(
    avg_exp_scaled,
    name = "Z-score",
    col = col_fun,
    cluster_rows = F,
    cluster_columns = FALSE,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    row_names_gp = gpar(fontsize = 11, fontface = "italic"),
    column_names_gp = gpar(fontsize = 12, fontface = "bold"),
    column_names_rot =90,
    rect_gp = gpar(col = "white", lwd = 1.5),
    heatmap_legend_param = list(
      title = "Exp(zScore)",
      title_gp = gpar(fontsize = 11, fontface = "bold"),
      labels_gp = gpar(fontsize = 10),
      legend_height = unit(4, "cm")
    ),
    row_title = NULL,
    column_title = "",
    column_title_gp = gpar(fontsize = 14, fontface = "bold")
  )
  dev.off()

}

##Gene Correlation---------
Gene1 <- "OMD"
Gene2 <- "COL1A2"
data <- data.frame(
  Gene1_exp = as.numeric(fib@assays$RNA$data[Gene1, ]),
  Gene2_exp = as.numeric(fib@assays$RNA$data[Gene2, ])
)
colnames(data) <- c(Gene1, Gene2)
p <- ggscatter(
  data,
  x = Gene1,
  y = Gene2,
  color = "#2C3E50",
  alpha = 0.55,
  size = 0.8,
  add = "reg.line",
  add.params = list(
    color = "#D73027",
    fill = "#F4A6A6",
    size = 1.1
  ),
  conf.int = TRUE,
  cor.coef = TRUE,
  cor.method = "pearson",
  cor.coef.size = 5,
  cor.coeff.args = list(
    label.x.npc = "left",
    label.y.npc = "top"
  )
) +
  labs(
    x = paste0("Expression of ", Gene1),
    y = paste0("Expression of ", Gene2)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 13, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.18, "cm"),    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),    
    plot.margin = margin(8, 8, 8, 8)
  )
p
ggsave(
  filename = "~/OSF/Fig/260504/OMD_COL1A2_pub.pdf",
  plot = p,
  width = 5.5,
  height = 5.5,
  dpi = 300,
  device = cairo_pdf, 
  family = "Arial"
)

##Venn
{
  load("~/OSF/gene_list.Rdata")
  intersect_genes <- intersect(intersect(gene_list[[1]], gene_list[[2]]), gene_list[[3]])
  print(intersect_genes)
  p <- ggvenn(
    data = gene_list,                        
    columns = c("TGF-beta receptor superfamily\nsignaling pathway", 
                "Response to\ntransforming growth factor beta", 
                "Cellular response to\nTGF-beta stimulus"), 
    fill_color = c("#6baed6", "#fec44f", "#78c679"), 
    fill_alpha = 0.5,                
    stroke_color = "white",          
    stroke_alpha = 1, 
    stroke_size = 0.8,               
    stroke_linetype = "solid", 
    set_name_color = "black",        
    set_name_size = 4.5,             
    text_color = "black",            
    text_size = 5,                   
    show_elements = FALSE,           
    show_percentage = FALSE          
  ) +
    theme_void() +                   
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.margin = margin(20, 20, 20, 20)
    )
  print(p)
  ggsave("~/OSF/Fig/260504/TGF_beta_Venn_diagram.pdf", plot = p, width = 7, height = 6, dpi = 300)  
}
