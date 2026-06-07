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
