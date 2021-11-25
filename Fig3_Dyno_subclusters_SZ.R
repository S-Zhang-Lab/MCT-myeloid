load("MCT.singlet_small_analysis_manuscript.rda")

setwd("/Volumes/GoogleDrive/Shared drives/Zhang_Lab_MS_Bhavana/Paper/Data/Bioinfo/Dyno_SZ")

library(Seurat)
library(dyno)
library(tidyverse)
library(Matrix)
library(ggplot2)

l <- c(ls())
l
l <- l[-72]
rm(list = (l)) # remove all other clusters. Only keep Dynosubclusters

###  Seurat plots
# get the metadata
Meta <- Dynosubclusters@meta.data

# confirm Seurat UMAP and feature violin plots
DimPlot(Dynosubclusters, split.by = "HTO_classification_group")
table(Dynosubclusters@meta.data$HTO_classification_group)
table(Dynosubclusters@meta.data$seurat_clusters)

table(Idents(Dynosubclusters), Dynosubclusters@meta.data$HTO_classification_group) # cell distribution group

# stacked bar chart
ggplot(Dynosubclusters@meta.data, aes(x=HTO_classification_group, fill=seurat_clusters)) + geom_bar(position = "fill")

ggplot(Dynosubclusters@meta.data, aes(x = HTO_classification_group, fill=seurat_clusters)) +  
  geom_bar(aes(y = (..count..)/sum(..count..)))

ggplot(Dynosubclusters@meta.data, aes(seurat_clusters, group = HTO_classification_group)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") + 
  scale_y_continuous(labels=scales::percent) +
  ylab("relative frequencies") +
  facet_grid(~HTO_classification_group)



# use blow to draw % barchart. 
ggplot(Dynosubclusters@meta.data, aes(seurat_clusters, group = HTO_classification_group)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
  geom_text(aes(label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="seurat_clusters") +
  facet_grid(~HTO_classification_group) +
  scale_y_continuous(labels = scales::percent)

# look for markers
All.markers <- FindAllMarkers(Dynosubclusters, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- All.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC) 
top5

cluster.markers <- c("S100a9","mt-Co3","Arg1","H2-Ab1","Cd74", "Mrc1","Ccl8","Ly6c2")
VlnPlot(Dynosubclusters, features =cluster.markers, group.by = "HTO_classification_group")
VlnPlot(Dynosubclusters, features =cluster.markers, split.by = "HTO_classification_group")

### Dyno
plot <- FeaturePlot(Dynosubclusters, features = "mt-Co3")
HoverLocator(plot = plot, information = FetchData(Dynosubclusters, vars = c("ident", "HTO_classification_group", "orig.ident")))
select.cells.start <- CellSelector(plot = plot)

plot <- FeaturePlot(Dynosubclusters, features = "mt-Co3")
HoverLocator(plot = plot, information = FetchData(Dynosubclusters, vars = c("ident", "HTO_classification_group", "orig.ident")))
select.cells.start <- CellSelector(plot = plot)


plot <- FeaturePlot(Dynosubclusters, features = "S100a9")
HoverLocator(plot = plot, information = FetchData(Dynosubclusters, vars = c("ident", "HTO_classification_group", "orig.ident")))
select.cells.end.1 <- CellSelector(plot = plot)

select.cells.end <- c(select.cells.end.1)

# get metadata
meta.Dynosubclusters <- Dynosubclusters@meta.data

Dynosubclusters.dataset <- list(
  cell_ids = row.names(Dynosubclusters@meta.data),
  expression = t(Dynosubclusters@assays$RNA@data), 
  counts = t(Dynosubclusters@assays$RNA@counts),
  group_ids = Dynosubclusters@meta.data$HTO_classification_group,
  group_ssn = Dynosubclusters@meta.data$seurat_clusters
)

# establish Dyno dataset
dataset_Dynosubclusters <- wrap_expression(
  expression = Dynosubclusters.dataset$expression, 
  counts = Dynosubclusters.dataset$counts,
  group_ids = Dynosubclusters.dataset$group_ids,
  group_ssn = Dynosubclusters@meta.data$seurat_clusters
)

dataset <- dataset_Dynosubclusters

# add staring ID
dataset <- add_prior_information(dataset, start_id = select.cells.start, end_id = select.cells.end)

#select methods
guidelines <- guidelines_shiny(dataset)


# Reproduces the guidelines as created in the shiny app
answers <- dynguidelines::answer_questions(
  multiple_disconnected = FALSE, 
  expect_topology = FALSE, 
  expected_topology = NULL, 
  expect_cycles = FALSE, 
  expect_complex_tree = TRUE, 
  n_cells = 3263, 
  n_features = 31053, 
  memory = "10GB", 
  prior_information = c("start_id", "end_id"), 
  docker = FALSE
)
guidelines <- dynguidelines::guidelines(answers = answers) 

# select methods
methods_selected <- guidelines$methods_selected


# infer model with the first method selected
# model <- infer_trajectory(dataset, first(methods_selected))

set.seed(1)
model <- infer_trajectory(dataset, "paga_tree")


# or provide more parameters. In this case, using umap for embedding instead of force altals. 
set.seed(1)
model <- infer_trajectory(dataset, ti_paga_tree(
  n_neighbors = 15L, 
  n_comps = 50L, 
  n_dcs = 15L, 
  resolution = 1L, 
  embedding_type = "umap"))


## combined plot
patchwork::wrap_plots(
  plot_dimred(model) + ggtitle("Cell ordering"),
  plot_dimred(model, grouping = group_onto_nearest_milestones(model)) + ggtitle("Cell grouping"),
  plot_dimred(model, feature_oi = "Cd274", expression_source = dataset) + ggtitle("Feature expression"),
  plot_dimred(model, "pseudotime", pseudotime = calculate_pseudotime(model)) + ggtitle("Pseudotime")
)

## Plotting the trajectory
# plot by group
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  grouping = dataset$group_ids
)


# plot by UMAP ssn 
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  grouping = dataset$group_ssn
)


# expression of a certain gene. Coloring by expression
plot_dimred(
  model, 
  expression_source = dataset$expression,
  feature_oi = "mt-Co3"
)

plot_dimred(
  model, 
  expression_source = dataset$expression,
  feature_oi = "S100a9"
)


plot_dimred(
  model, 
  expression_source = dataset$expression,
  feature_oi = "Cd274"
)

plot_dimred(
  model, 
  expression_source = dataset$expression,
  feature_oi = "Stat1"
)

plot_dimred(
  model, 
  expression_source = dataset$expression,
  feature_oi = "Mrc1"
)

plot_dimred(
  model, 
  expression_source = dataset$expression,
  feature_oi = "Arg1"
)

plot_dimred(
  model, 
  expression_source = dataset$expression,
  feature_oi = "H2-Ab1"
)


# visualised gene using a background color




## -- starting --- 
# cluster 1 marker  
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_cells = "feature",
  feature_oi = "mt-Co3",
  color_density = "grouping",
  grouping = dataset$group_ssn,
  label_milestones = FALSE
)

# cluster 8 marker 
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_cells = "feature",
  feature_oi = "Eps8",
  color_density = "grouping",
  grouping = dataset$group_ssn,
  label_milestones = FALSE
)

#3 -- ending ---
# cluster 0 marker 
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_cells = "feature",
  feature_oi = "S100a9",
  color_density = "grouping",
  grouping = dataset$group_ssn,
  label_milestones = FALSE
)

plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_cells = "feature",
  feature_oi = "Ccr1",
  color_density = "grouping",
  grouping = dataset$group_ssn,
  label_milestones = FALSE
)


plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_cells = "feature",
  feature_oi = "Cd274",
  color_density = "grouping",
  grouping = dataset$group_ssn,
  label_milestones = FALSE
)


plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_cells = "feature",
  feature_oi = "Stat1",
  color_density = "grouping",
  grouping = dataset$group_ssn,
  label_milestones = FALSE
)


# cluster 4 marker 
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_cells = "feature",
  feature_oi = "Arg1",
  color_density = "grouping",
  grouping = dataset$group_ssn,
  label_milestones = FALSE
)

# cluster 3 marker 
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_cells = "feature",
  feature_oi = "H2-Ab1",
  color_density = "grouping",
  grouping = dataset$group_ssn,
  label_milestones = FALSE
)

# cluster 7 marker 
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_cells = "feature",
  feature_oi = "Mrc1",
  color_density = "grouping",
  grouping = dataset$group_ssn,
  label_milestones = FALSE
)

# cluster 14 marker 
plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_cells = "feature",
  feature_oi = "Siglech",
  color_density = "grouping",
  grouping = dataset$group_ssn,
  label_milestones = FALSE
)


plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_cells = "feature",
  feature_oi = "Ccr1",
  color_density = "grouping",
  grouping = dataset$group_ssn,
  label_milestones = FALSE
)
### A global overview of the most predictive genes
# overall_feature_importances
overall_feature_importances <- dynfeature::calculate_overall_feature_importance(model, expression_source = dataset$expression)
features <- overall_feature_importances %>% 
  top_n(40, importance) %>% 
  pull(feature_id)

# Lineage/branch markers
branch_feature_importance <- calculate_branch_feature_importance(model)
features <- branch_feature_importance %>% 
  filter(to == "5") %>% 
  top_n(20, importance) %>% 
  pull(feature_id)

# Genes important at bifurcation points
branching_milestone <- "5"
branch_feature_importance <- calculate_branching_point_feature_importance(model, milestones_oi = branching_milestone)

features <- branch_feature_importance %>% top_n(20, importance) %>% pull(feature_id)

## Heatmap
# No features of interest provided, selecting the top 50 features automatically
# You can define to us dynfeature for selecting the top 10 features
# Coloring by grouping

plot_heatmap(
  model,
  expression_source = dataset$expression,
  grouping = dataset$group_ssn,
  features_oi = 50
)

# change rooting after referring to RNA velocity
model_rooted <- model %>% add_root(root_milestone_id = "5")
plot_heatmap(
  model_rooted,
  expression_source = dataset$expression,
  grouping = dataset$group_ssn,
  features_oi = 50
)
