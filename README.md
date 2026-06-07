# scRNA-seq-R-scripts-Cell-type-identification
#==============================使用SingleR工具进行注释=========================
setwd('/Users/apple/Desktop/data_file/数据和代码 3')
seu_obj <- readRDS('seu_obj.rds')

# 安装和加载注释工具
# remotes::install_github("immunogenomics/singleR")
library(SingleR)
library(celldex)
library(Seurat)
library(dplyr)
library(ggplot2)

# 使用参考数据库进行自动注释
# 加载人类原代细胞图谱数据作为参考数据集
# 这个数据集包含了多种人类原代细胞的基因表达谱和对应的细胞类型标签
ref <- HumanPrimaryCellAtlasData()

# 如果是分析小鼠数据，应使用小鼠参考数据库：
# ref <- MouseRNAseqData()

# 使用SingleR算法进行细胞类型注释
# SingleR通过比较测试数据集与参考数据集的基因表达相似性来预测细胞类型
pred <- SingleR(
  test = GetAssayData(seu_obj, assay = "RNA", layer = "data"),  # 获取待注释的单细胞数据
  ref = ref,  # 使用上面加载的参考数据集
  labels = ref$label.main  # 指定参考数据集中的主要细胞类型标签
)

# 将注释结果添加到原始单细胞对象的metadata中
# 这样可以将预测的细胞类型与每个细胞关联起来，便于后续分析
seu_obj$singleR_labels <- pred$labels


# 查看注释结果
DimPlot(seu_obj, reduction = "umap", group.by = "singleR_labels", label = TRUE)
DimPlot(seu_obj, reduction = "tsne", group.by = "singleR_labels", label = TRUE)

saveRDS(seu_obj, 'seu_obj.rds')

# 计算每种细胞类型的数量和比例
cell_type_summary <- seu_obj@meta.data %>%
  group_by(singleR_labels) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  arrange(desc(count))

# 绘制条形图
ggplot(cell_type_summary, aes(x = reorder(singleR_labels, -count), y = count, fill = singleR_labels)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            vjust = -0.5, size = 3) +
  labs(title = "Cell Type Composition", 
       x = "Cell Type", 
       y = "Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")



#=============================手动注释==================================
# 找出每个cluster的差异表达基因
cluster_markers <- FindAllMarkers(
  object = seu_obj,
  only.pos = TRUE,          # 只保留上调的基因
  min.pct = 0.25,           # 在至少25%的细胞中表达
  logfc.threshold = 0.25    # log2 fold change阈值
)

#FindMarkers()   找两组数据之间的marker进行对比

# 查看每个cluster的前几个marker基因
top_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

print(top_markers)

DimPlot(seu_obj, reduction = 'umap', label = TRUE)
DimPlot(seu_obj, reduction = 'tsne', label = TRUE)


#===========================方法一====================================
# http://117.50.127.228/CellMarker/


# 定义一些常见的细胞类型marker基因
feature_genes <- c(
  "CD3D", "CD3E", "CD4", "CD8A", # T细胞
  "CD19", "CD79A", "MS4A1",      # B细胞
  "CD14", "CD68", "FCGR3A",      # 单核/巨噬细胞
  "NKG7", "GNLY",                # NK细胞
  "PPBP", "PF4",                 # 血小板
  "COL1A1", "COL1A2",            # 成纤维细胞
  "PECAM1", "VWF",               # 内皮细胞
  "EPCAM", "KRT8", "KRT18"       # 上皮细胞
)

# 可视化这些marker基因的表达
DotPlot(seu_obj, features = feature_genes) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 或者在UMAP上查看特定marker的表达
p1 <- FeaturePlot(seu_obj, features = c("CD3D", "CD3E", "CD4", "CD8A"))
p2 <- DimPlot(seu_obj, reduction = 'umap', label = TRUE)

p1 | p2

# 添加celltype
seu_obj@meta.data$cell_type <- 'unknown'

seu_obj@meta.data$cell_type[seu_obj@meta.data$seurat_clusters %in% c(29,27,5,8,24)] <- 'T_cells'

DimPlot(seu_obj,group.by = 'cell_type', reduction = 'umap')
