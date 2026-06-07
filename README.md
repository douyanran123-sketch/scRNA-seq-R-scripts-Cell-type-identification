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



#===================================方法二===============================

# 定义不同细胞类型的标记基因集
# 这里创建了一个列表，每个元素代表一种细胞类型和其对应的特征基因
marker_gene_sets <- list(
  T_cell = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B"),        # T细胞标记基因
  B_cell = c("CD19", "CD79A", "MS4A1", "CD79B"),            # B细胞标记基因
  Monocyte = c("CD14", "CD68", "FCGR3A", "LYZ"),            # 单核细胞标记基因
  NK_cell = c("NKG7", "GNLY", "KLRD1", "KLRF1"),            # NK细胞标记基因
  Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5"),        # 内皮细胞标记基因
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM"),         # 成纤维细胞标记基因
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19")         # 上皮细胞标记基因
)

# 为每种细胞类型计算得分
# 使用循环为每种细胞类型计算模块得分
for(cell_type in names(marker_gene_sets)) {
  seu_obj <- AddModuleScore(
    object = seu_obj,                                      # 目标Seurat对象
    features = list(marker_gene_sets[[cell_type]]),        # 该细胞类型的标记基因列表
    name = paste0(cell_type, "_score"),                    # 新列的名称前缀
    ctrl = 100                                             # 控制基因的数量（用于背景校正）
  )
}
# AddModuleScore会为每个细胞计算一个得分，表示该细胞表达这些标记基因的程度
# 控制基因用于校正技术偏差和批次效应

# 查看新添加的列
head(seu_obj@meta.data)
# 查看元数据的前几行，确认新的得分列已添加（如T_cell_score1, B_cell_score1等）
# 注意：函数会自动在名称后添加"1"，所以实际列名为"T_cell_score1"等

# 方法1: 小提琴图查看每个cluster的得分分布
VlnPlot(seu_obj, features = c("T_cell_score1", "B_cell_score1", "Monocyte_score1"), 
        pt.size = 0, ncol = 3)
# 绘制小提琴图，显示每个cluster中三种细胞类型得分的分布情况
# pt.size = 0 不显示个体数据点，使图形更清晰
# ncol = 3 将三个图排列成3列

# 方法2: 特征图查看得分在UMAP上的分布
FeaturePlot(seu_obj, features = c("T_cell_score1", "B_cell_score1", "Monocyte_score1"),
            cols = c("lightgrey", "red"), ncol = 3)
# 在UMAP降维图上可视化得分分布
# 灰色表示低得分，红色表示高得分
# 可以直观地看到哪些区域高表达特定细胞类型的标记基因

# 方法3: 点图比较不同cluster的得分
DotPlot(seu_obj, features = c("T_cell_score1", "B_cell_score1", "Monocyte_score1", 
                              "NK_cell_score1", "Endothelial_score1"),
        group.by = "seurat_clusters") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 创建点图，每个点代表一个cluster的平均得分
# 点的大小表示表达该特征的细胞比例，颜色深浅表示平均表达水平
# 可以快速比较不同cluster对各种细胞类型的倾向性



# 加载必要的包
library(dplyr)    # 用于数据操作和转换
library(tidyr)    # 用于数据重塑和整理

# 计算每个cluster对各种细胞类型的平均得分
score_data <- seu_obj@meta.data %>%              # 从元数据开始
  group_by(seurat_clusters) %>%                  # 按cluster分组
  summarise(across(ends_with("score1"), mean)) %>%  # 计算每个得分列的平均值
  pivot_longer(cols = -seurat_clusters,          # 将宽数据转换为长数据
               names_to = "cell_type",           # 将列名转换为cell_type列
               values_to = "mean_score")         # 将值转换为mean_score列

# 找出每个cluster得分最高的细胞类型
best_annotations <- score_data %>%
  group_by(seurat_clusters) %>%                  # 按cluster分组
  slice_max(mean_score, n = 1) %>%               # 选择每个组中得分最高的行
  mutate(cell_type = gsub("_score1", "", cell_type))  # 清理cell_type名称（去掉"_score1"）

print(best_annotations)
# 输出结果，显示每个cluster最可能的细胞类型注释

# 完全自动注释
auto_annotation <- best_annotations$cell_type[match(seu_obj@meta.data$seurat_clusters, 
                                                    best_annotations$seurat_clusters)]
# 使用match函数将每个细胞的cluster编号映射到对应的最佳注释
# match函数返回位置索引，然后用这些索引提取对应的细胞类型名称

seu_obj$auto_celltype <- auto_annotation
# 将自动注释结果添加到Seurat对象的元数据中，列名为"auto_celltype"

DimPlot(seu_obj, group.by = 'auto_celltype', reduction = 'umap', label = TRUE)


# 计算每种细胞类型的数量和比例
cell_type_summary <- seu_obj@meta.data %>%
  group_by(auto_celltype) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  arrange(desc(count))

# 绘制条形图
ggplot(cell_type_summary, aes(x = reorder(auto_celltype, -count), y = count, fill = auto_celltype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            vjust = -0.5, size = 3) +
  labs(title = "Cell Type Composition", 
       x = "Cell Type", 
       y = "Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
