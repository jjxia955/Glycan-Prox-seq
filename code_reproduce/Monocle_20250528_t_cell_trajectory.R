library(Seurat)
library(dplyr)
library(monocle3)
library(ggplot2)
library(scattermore)

# Load dataset
data_dir <- "/Users/junjie/Desktop/OneDrive - The University of Chicago (1)/Python/GPS/20250521/Ctrl_filtered_feature_bc_matrix"  # folder with matrix.mtx, barcodes.tsv, features.tsv
seurat_obj <- Read10X(data.dir = data_dir) %>%
  CreateSeuratObject(project = "T cell differentiation trajectory")

# load cell metadata
annotation <- read.csv("/Users/junjie/Desktop/OneDrive - The University of Chicago (1)/Python/GPS/20250521/ctrl_annotations_0801.csv", row.names = 1)

#Add annotation to Seurat object metadata (only for matched barcodes)
matched_barcodes <- intersect(rownames(seurat_obj@meta.data), rownames(annotation))
seurat_obj <- AddMetaData(seurat_obj, metadata = annotation[matched_barcodes, , drop = FALSE])

#Subset to CD4+ T cell subtypes
cd4_subtypes <- c("CD4+ Naive T cells", "CD4+ Central Memory", "CD4+ Effector Memory")
seurat_cd4 <- subset(seurat_obj, subset = manual_annotation_reclustered %in% cd4_subtypes)

#expression matrix
expr_matrix <- GetAssayData(seurat_cd4, slot = "counts")
#cell_metadata
cell_metadata <- seurat_cd4@meta.data
#gene_matadata
gene_metadata <- data.frame(
  gene_short_name = rownames(seurat_cd4),
  biotype = NA,       
  gc_content = NA      
)
rownames(gene_metadata) <- gene_metadata$gene_short_name


# Create the Monocle3 cell_data_set object
cds_cd4 <- new_cell_data_set(
  as(expr_matrix, "sparseMatrix"),
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

## Step 1: Normalize and pre-process the data
cds_cd4 <- preprocess_cds(cds_cd4, num_dim = 13)

plot_pc_variance_explained(cds_cd4)

plot_cells(cds_cd4, reduction_method = "PCA",
           color_cells_by = "manual_annotation_reclustered", group_label_size = 3.5,
           label_groups_by_cluster = FALSE)

#UMAP reduction
cds_cd4 <- reduce_dimension(cds_cd4) # umap.min_dist = 0.3,umap.n_neighbors = 30)

plot_cells(cds_cd4, color_cells_by = 'manual_annotation_reclustered')

cds_cd4 <- cluster_cells(cds_cd4, resolution=1e-2)
plot_cells(cds_cd4)

#filter noisy cluster
cds_cd4$monocle_cluster <- clusters(cds_cd4)
cds_cd4<- cds_cd4[, !cds_cd4$monocle_cluster %in% c('19','23','18','10','7','9','22','3')]

#update
## Step 1: Normalize and pre-process the data
cds_cd4 <- preprocess_cds(cds_cd4, num_dim = 20)
plot_pc_variance_explained(cds_cd4)
plot_cells(cds_cd4, reduction_method = "PCA",
           color_cells_by = "manual_annotation_reclustered", group_label_size = 3.5,
           label_groups_by_cluster = FALSE)

#UMAP reduction
cds_cd4 <- reduce_dimension(cds_cd4) # umap.min_dist = 0.3,umap.n_neighbors = 30)

plot_cells(cds_cd4, color_cells_by = 'manual_annotation_reclustered')

cds_cd4 <- cluster_cells(cds_cd4, resolution=1e-2)
plot_cells(cds_cd4)
#filter noisy cluster
cds_cd4$monocle_cluster <- clusters(cds_cd4)
cds_cd4<- cds_cd4[, !cds_cd4$monocle_cluster %in% c('2')]

#update
cds_cd4 <- preprocess_cds(cds_cd4, num_dim = 50)
plot_pc_variance_explained(cds_cd4)
plot_cells(cds_cd4, reduction_method = "PCA",
           color_cells_by = "manual_annotation_reclustered", group_label_size = 3.5,
           label_groups_by_cluster = FALSE)

#UMAP reduction
cds_cd4 <- reduce_dimension(cds_cd4) # umap.min_dist = 0.3,umap.n_neighbors = 30)

plot_cells(cds_cd4, color_cells_by = 'manual_annotation_reclustered')

#recluster
cds_cd4 <- cluster_cells(cds_cd4, resolution=1e-4)
plot_cells(cds_cd4)
plot_cells(cds_cd4, color_cells_by="partition", group_cells_by="partition")

cds_cd4 <- learn_graph(cds_cd4,use_partition = TRUE)
plot_cells(cds_cd4,
           color_cells_by = "manual_annotation_reclustered",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)


cds_cd4 <- order_cells(cds_cd4)

p1 <- plot_cells(
  cds_cd4,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
) +
  scale_color_viridis_c(option = "inferno")

# Define genes in the desired order
de_genes <- c("CCR7", "SELL", "IL7R", "CCL5")

# Get CD4 cells and their annotations
cd4_cells <- colnames(cds_cd4)[grepl("^CD4", colData(cds_cd4)$manual_annotation_reclustered)]
cd4_annotations <- colData(cds_cd4)[cd4_cells, "manual_annotation_reclustered"]

# Define desired CD4 order
cd4_order <- c("CD4+ Naive T cells", "CD4+ Central Memory", "CD4+ Effector Memory")

# Sort CD4 cells by annotation order
cd4_cells_ordered <- cd4_cells[order(match(cd4_annotations, cd4_order))]

# Subset to selected genes and ordered CD4 cells
cd4_lineage_cds <- cds_cd4[rowData(cds_cd4)$gene_short_name %in% de_genes, cd4_cells_ordered]

# Set gene plotting order
rowData(cd4_lineage_cds)$gene_short_name <- factor(
  rowData(cd4_lineage_cds)$gene_short_name,
  levels = de_genes
)

# Ensure cells are ordered by pseudotime (required by Monocle)
cd4_lineage_cds <- order_cells(cd4_lineage_cds)

# Plot in pseudotime with custom gene order and cell grouping
plot_genes_in_pseudotime(
  cd4_lineage_cds,
  color_cells_by = "manual_annotation_reclustered",
  min_expr = 0.5
)


####================================================###

#DE genes
marker_test_res <- top_markers(cds, group_cells_by="cluster", 
                               reference_cells=500, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)

plot_cells(cds, genes=c("SELL", "LEF1", "TCF7", "IL7R", "CCR7","GZMB"))


# Extract pseudotime and UMAP coordinates
plot_data <- data.frame(
  reducedDims(cds)$UMAP,
  pseudotime = monocle3::pseudotime(cds),
  annotation = cds@colData$manual_annotation
)

colnames(plot_data)[1:2] <- c("UMAP_1", "UMAP_2")

# Create scattermore plot with rasterized points
p <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_scattermore(pointsize = 2) +
  scale_color_viridis_c(option = "inferno") +
  labs(title = "Trajectory colored by Pseudotime", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()


ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

plot_cells(cds, genes=c("CCR7", "SELL", "IL7R", "CCL5"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE) 


######====================CD8 T cell trajectory=================================

#Subset to CD8+ T cell subtypes
cd8_subtypes <- c("CD8+ Naive T cells", "CD8+ Central Memory", "CD8+ Effector Memory", "CD8+ Effector Memory RA")
seurat_cd8 <- subset(seurat_obj, subset = manual_annotation_reclustered %in% cd8_subtypes)

#expression matrix
expr_matrix <- GetAssayData(seurat_cd8, slot = "counts")
#cell_metadata
cell_metadata <- seurat_cd8@meta.data
#gene_matadata
gene_metadata <- data.frame(
  gene_short_name = rownames(seurat_cd8),
  biotype = NA,       
  gc_content = NA      
)
rownames(gene_metadata) <- gene_metadata$gene_short_name


# Create the Monocle3 cell_data_set object
cds <- new_cell_data_set(
  as(expr_matrix, "sparseMatrix"),
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)

#UMAP reduction
cds <- reduce_dimension(cds)
plot_cells(cds, color_cells_by = 'manual_annotation_reclustered')

cds <- cluster_cells(cds, resolution=1e-2)
plot_cells(cds)

#filter noisy cluster
cds$monocle_cluster <- clusters(cds)
cds<- cds[, !cds$monocle_cluster %in% c('2','4','19','16')]

#update
## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)

#UMAP reduction
cds <- reduce_dimension(cds)
plot_cells(cds, color_cells_by = 'manual_annotation_reclustered')

#recluster
cds <- cluster_cells(cds, resolution=1e-4)
plot_cells(cds)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")

cds <- learn_graph(cds,use_partition = TRUE)
plot_cells(cds,
           color_cells_by = "manual_annotation_reclustered",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

# Choose a cell from Naive T group
cds <- order_cells(cds)

plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
) +
  scale_color_viridis_c(option = "inferno")

de_genes <- c("CCR7", 'LEF1',"IL7R",'NKG7')
# Check actual gene symbols
unique(rowData(cds)$gene_short_name)

# Subset only CD8-related cells (use substring matching if needed)
cd8_cells <- colnames(cds)[grepl("^CD8", colData(cds)$manual_annotation_reclustered)]

# Subset CDS to those genes and CD8 cells
lineage_cds <- cds[rowData(cds)$gene_short_name %in% de_genes, cd8_cells]

# Order cells along pseudotime (if trajectory is learned)
lineage_cds <- order_cells(lineage_cds)

# Plot genes in pseudotime
plot_genes_in_pseudotime(
  lineage_cds,
  color_cells_by = "manual_annotation_reclustered",
  min_expr = 0.5
)

###=====================================================####


#DE genes
marker_test_res <- top_markers(cds, group_cells_by="cluster", 
                               reference_cells=500, cores=8)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)

plot_cells(cds, genes=c("SELL", "LEF1", "TCF7", "IL7R", "CCR7","GZMB"))



ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

plot_cells(cds, genes=c("ISG15", "TNFRSF4", "TNFRSF18", "TNFRSF14", "TNFRSF25",'B3GALT6'),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)



