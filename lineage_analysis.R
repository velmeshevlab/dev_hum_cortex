#import functions
source(file = "functions.R")

##############################
##############################
#Analysis of specific lineages
##############################
##############################

#you can get the matrices and metadata for each lineage group (excitatory neurons, interneurons etc) at https://cells.ucsc.edu/?ds=pre-postnatal-cortex (snRNA-seq)
##############################################
#excitatory neuron analysis
##############################################

#select excitatory neurons and recluster
meta = data@'meta.data'
cells = rownames(meta[meta$seurat_clusters %in% c("34", "3", "16", "37", "32", "11", "21", "18", "0", "23", "7", "33", "28", "20", "30", "24", "13", "14", "17", "36", "22", "12", "29") & !(meta$region_broad %in% c("GE", "MGE", "CGE", "LGE")),])
data = data[,cells]
N = 19
data <-FindVariableFeatures(data,selection.method='mean.var.plot')
data <-RunPCA(data,features = VariableFeatures(data),verbose=F)
data@reductions$pca2 <- data@reductions$pca
data@reductions$pca2@"cell.embeddings" <- data@reductions$pca2@"cell.embeddings"[,1:N]
data<- RunHarmony(data, group.by.vars = "chemistry", theta = 2, max.iter.harmony = 20, reduction = 'pca2', dims.use = 1:N)
data<- RunUMAP(data, dims = 1:N, reduction = 'harmony', return.model=TRUE)
data<- FindNeighbors(data, reduction = "harmony", dims = 1:N) %>% FindClusters()
#subcluster select clusters
sel.cluster = "3"
data = FindSubCluster(data, sel.cluster, graph.name = "RNA_snn", resolution = 0.1)
Idents(data.sub) <- "sub.cluster"
sel.cluster = "10"
data = FindSubCluster(data, sel.cluster, graph.name = "RNA_snn", resolution = 0.1)
Idents(data.sub) <- "sub.cluster"
sel.cluster = "11"
data = FindSubCluster(data, sel.cluster, graph.name = "RNA_snn", resolution = 0.1)
Idents(data.sub) <- "sub.cluster"

#run monocle 3
d = GetAssayData(object = data, assay = "RNA", slot = "counts")
gene_annotation = as.data.frame(rownames(d))
rownames(gene_annotation) = rownames(d)
colnames(gene_annotation) = 'gene_short_name'
meta = data@"meta.data"
cds = new_cell_data_set(d, cell_metadata = meta, gene_metadata = gene_annotation)
s.umap <- data@"reductions"$umap[[]]
s.umap = s.umap[colnames(cds),]
reducedDims(cds)$"UMAP" <- s.umap
s.clusters = as.character(Idents(data))
names(s.clusters) <- names(Idents(data))
s.clusters = s.clusters[colnames(cds)]
cds@clusters$"UMAP"$"clusters" <- s.clusters
cds@clusters$UMAP$partitions <- cds@clusters$UMAP$clusters
cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions != "1"] <- "1"
graph_control = setNames(list(3, 1, 5, FALSE, TRUE, FALSE, NULL, 10, 1e-5, 0.05, 0.01), c("euclidean_distance_ratio","geodesic_distance_ratio","minimal_branch_len","orthogonal_proj_tip","prune_graph","scale","rann.k","maxiter","eps","L1.gamma","L1.sigma"))
cds <- learn_graph(cds, use_partition = F, learn_graph_control = graph_control)

#import monocle object and clean up trajectories
cds = import_monocle(cds)
node_plot(cds)
#connect nodes
cds = connect_nodes(cds, "Y_2103", "Y_1803", add_node = T)
cds = connect_nodes(cds, "Y_646", "Y_283", add_node = T)
cds = connect_nodes(cds, "Y_2212", "Y_55", add_node = T)
cds = connect_nodes(cds, "Y_1853", "Y_426", add_node = T)
cds = connect_nodes(cds, "Y_2234", "Y_2443", add_node = T)
cds = connect_nodes(cds, "Y_2443", "Y_275", add_node = T)
cds = connect_nodes(cds, "Y_302", "Y_201", add_node = T)
cds = connect_nodes(cds, "Y_775", "Y_654", add_node = T)

######################
###isolate lineages###
######################
#L2-3
lineage = "L2_3"
start = 1334
end = 2344
inc.node = c("Y_1853", "Y_572", "Y_2177")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("1", "28", "0", "9", "13", "8", "2", "27", "29", "39", "5")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 5)

#L4
lineage = "L4"
start = 1334
end = 93
inc.node = c("Y_1225", "Y_957", "Y_201")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("1", "4", "0", "40", "14", "6", "38", "20", "26", "16", "43", "15", "22")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 5)

#L5
lineage = "L5"
start = 1320
end = 2506
inc.node = c("Y_2103", "Y_775", "Y_654")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("24", "36", "11", "25", "3", "35")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 5)

#L6
lineage = "L6"
start = 1320
end = 2514
inc.node = c("Y_2103")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("24", "36", "11", "25", "23")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 5)

#L5-6-IT
lineage = "L5_6_IT"
start = 1334
end = 37
inc.node = c("Y_1853", "Y_2234")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("1", "28", "0", "9", "10", "7", "45", "26", "41")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 5)

#SP
lineage = "SP"
start = 1334
end = 762
inc.node = c("Y_1177")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("1", "4", "33", "18", "17", "12", "46", "37")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 5)

#add the connection between the two graphs
lineage = "RG"
start = 1334
end = 1320
cds <- isolate_graph(cds, start, end, lineage)
sel.cluster = c("1", "24")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 5)

#combine lineages and calculate pseudotime
cds_new = combine_lineages(cds, c(1334, 1320))
cds_new = order_cells(cds_new, root_pr_nodes = c("Y_1334", "Y_1320"))

#compress into metacells
cds_new <- compress_lineage_v2(cds_new, "L5_6_IT", start = 1334)
cds_new <- compress_lineage_v2(cds_new, "SP", start = 1334)
cds_new <- compress_lineage_v2(cds_new, "L2_3", start = 1334)
cds_new <- compress_lineage_v2(cds_new, "L4", start = 1334)
cds_new <- compress_lineage_v2(cds_new, "L5", start = 1320)
cds_new <- compress_lineage_v2(cds_new, "L6", start = 1320)

#identify lineage-specific genes
#modify monocle 3 functions to enable linear model
trace('graph_test', edit = T, where = asNamespace("monocle3"))
#replace with "https://github.com/velmeshevlab/dev_hum_cortex/blob/main/graph_test_lm.R"
trace('my.moran.test', edit = T, where = asNamespace("monocle3"))
#replace with "https://github.com/velmeshevlab/dev_hum_cortex/blob/main/my.moran.test_lm.R"

lineages = c("L2_3", "L4", "L5_6_IT", "SP")
for(lineage in lineages){
factor = 0.05
N = 10000
print(lineage)
cds.sub = get_lineage_object(cds_new, lineage, 1334, N = N)
meta.sub = meta[colnames(cds.sub),]
colData(cds.sub)$sex <- meta.sub$sex
colData(cds.sub)$region_broad <- meta.sub$region_broad
colData(cds.sub)$chemistry <- meta.sub$chemistry
colData(cds.sub)$PMI <- meta.sub$PMI
colData(cds.sub)$UMI <- meta.sub$nCount_RNA
data = counts(cds.sub)
rows=rownames(data)[rowSums(data> 0) > factor*ncol(data)]
cds.sel = cds.sub[rows]
cds_pr_test_res <- monocle3:::graph_test(cds.sel, neighbor_graph="principal_graph", verbose = T, cores = 24)
save(cds_pr_test_res, file = paste("Moran_", lineage,".R", sep = ""))
write.table(subset(cds_pr_test_res, q_value < 0.05), paste("pt_DGE_", lineage,".txt", sep = ""), sep ="\t", quote = F)
}

lineages = c("L5", "L6")
for(lineage in lineages){
factor = 0.05
N = 10000
print(lineage)
cds.sub = get_lineage_object(cds_new, lineage, 1320, N = N)
meta.sub = meta[colnames(cds.sub),]
colData(cds.sub)$sex <- meta.sub$sex
colData(cds.sub)$region_broad <- meta.sub$region_broad
colData(cds.sub)$chemistry <- meta.sub$chemistry
colData(cds.sub)$PMI <- meta.sub$PMI
colData(cds.sub)$UMI <- meta.sub$nCount_RNA
data = counts(cds.sub)
rows=rownames(data)[rowSums(data> 0) > factor*ncol(data)]
cds.sel = cds.sub[rows]
cds_pr_test_res <- monocle3:::graph_test(cds.sel, neighbor_graph="principal_graph", verbose = T, cores = 24)
save(cds_pr_test_res, file = paste("Moran_", lineage,".R", sep = ""))
write.table(subset(cds_pr_test_res, q_value < 0.05), paste("pt_DGE_", lineage,".txt", sep = ""), sep ="\t", quote = F)
}

##############################################
#analysis of interneurons
##############################################
#select interneurons and recluster
meta = data@'meta.data'
cells1 = rownames(meta[meta$seurat_clusters %in% c("27", "8", "31", "2", "9", "10"),])
cells2 = rownames(meta[meta$seurat_clusters == "15" & meta$region_broad %in% c("GE", "MGE", "CGE", "LGE"),])
data = data[,c(cells1, cells2)]
data <-FindVariableFeatures(data,selection.method='mean.var.plot')
data <-RunPCA(data,features = VariableFeatures(data),verbose=F)
data@reductions$pca2 <- data@reductions$pca
N = 14
data@reductions$pca2@"cell.embeddings" <- data@reductions$pca2@"cell.embeddings"[,1:N]
data<- RunHarmony(data, group.by.vars = "chemistry", theta = 2, max.iter.harmony = 20, reduction = 'pca2', dims.use = 1:N)
data<- RunUMAP(data, dims = 1:N, reduction = 'harmony', return.model=TRUE)
data<- FindNeighbors(data, reduction = "harmony", dims = 1:N) %>% FindClusters()
#subcluster select clusters
sel.cluster = "3"
data.sub = FindSubCluster(data, sel.cluster, graph.name = "RNA_snn", resolution = 0.1)
Idents(data.sub) <- "sub.cluster"
sel.cluster = "10"
data.sub = FindSubCluster(data, sel.cluster, graph.name = "RNA_snn", resolution = 0.1)
Idents(data.sub) <- "sub.cluster"
sel.cluster = "11"
data.sub = FindSubCluster(data, sel.cluster, graph.name = "RNA_snn", resolution = 0.1)
Idents(data.sub) <- "sub.cluster"

#run monocle 3
d = GetAssayData(object = data, assay = "RNA", slot = "counts")
gene_annotation = as.data.frame(rownames(d))
rownames(gene_annotation) = rownames(d)
colnames(gene_annotation) = 'gene_short_name'
meta = data@"meta.data"
cds = new_cell_data_set(d, cell_metadata = meta, gene_metadata = gene_annotation)
s.umap <- data@"reductions"$umap[[]]
s.umap = s.umap[colnames(cds),]
reducedDims(cds)$"UMAP" <- s.umap
s.clusters = as.character(Idents(data))
names(s.clusters) <- names(Idents(data))
s.clusters = s.clusters[colnames(cds)]
cds@clusters$"UMAP"$"clusters" <- s.clusters
cds@clusters$UMAP$partitions <- cds@clusters$UMAP$clusters
cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions != "1"] <- "1"
graph_control = setNames(list(3, 1, 1, FALSE, TRUE, FALSE, NULL, 10, 1e-5, 0.05, 0.01), c("euclidean_distance_ratio","geodesic_distance_ratio","minimal_branch_len","orthogonal_proj_tip","prune_graph","scale","rann.k","maxiter","eps","L1.gamma","L1.sigma"))
cds <- learn_graph(cds, use_partition = F, learn_graph_control = graph_control)

#import monocle object and clean up trajectories
cds = import_monocle(cds)
node_plot(cds)
ggsave(file="graph.png",width = 7, height = 6, units = "in",  dpi = 1200)
cds = connect_nodes(cds, "Y_19", "Y_414", add_node = T)
cds = connect_nodes(cds, "Y_381", "Y_106", add_node = T)
cds = connect_nodes(cds, "Y_457", "Y_459", add_node = T)
cds = connect_nodes(cds, "Y_108", "Y_117", add_node = T)


######################
###isolate lineages###
######################
#VIP
lineage = "VIP"
start = 203
end = 36
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("16", "1", "5", "3_1")
starting.clusters = c("19", "4", "0")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 2, start_regions = c("CGE", "GE"), starting_clusters = starting.clusters)

#CALB2
lineage = "CALB2"
start = 203
end = 386
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("16", "1", "5", "3_0")
starting.clusters = c("19", "4", "0")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 2, start_regions = c("CGE", "GE"), starting_clusters = starting.clusters)

#CCK
lineage = "CCK"
start = 203
end = 17
inc.node = c("Y_72", "Y_364")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("16", "1", "5", "10_0")
starting.clusters = c("19", "4", "0")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 2, start_regions = c("CGE", "GE"), starting_clusters = starting.clusters)

#RELN
lineage = "RELN"
start = 203
end = 62
inc.node = c("Y_72", "Y_364")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("16", "1", "5", "10_1")
starting.clusters = c("19", "4", "0")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 2, start_regions = c("CGE", "GE"), starting_clusters = starting.clusters)

#SV2C
lineage = "SV2C"
start = 203
end = 429
inc.node = c("Y_171", "Y_85")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("16", "1", "13", "11_2", "11_0")
starting.clusters = c("19", "4", "0")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 2, start_regions = c("CGE", "GE"), starting_clusters = starting.clusters)

#NOS
lineage = "NOS"
start = 203
end = 61
inc.node = c("Y_171", "Y_85")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("16", "1", "13", "11_2", "11_1")
starting.clusters = c("19", "4", "0")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 2, start_regions = c("CGE", "GE"), starting_clusters = starting.clusters)

#PV_MP
lineage = "PV_MP"
start = 203
end = 409
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("2", "15", "9", "8", "17")
starting.clusters = c("19", "4", "0")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 2, start_regions = c("MGE", "GE"), starting_clusters = starting.clusters)

#PV
lineage = "PV"
start = 203
end = 482
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("2", "15", "9", "8")
starting.clusters = c("19", "4", "0")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 2, start_regions = c("MGE", "GE"), starting_clusters = starting.clusters)
#SST
lineage = "SST"
start = 203
end = 12
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("2", "14", "6", "12", "7")
starting.clusters = c("19", "4", "0")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 2, start_regions = c("MGE", "GE"), starting_clusters = starting.clusters)

#SST_RELN
lineage = "SST_RELN"
start = 203
end = 465
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("2", "14", "6", "12")
starting.clusters = c("19", "4", "0")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 2, start_regions = c("MGE", "GE"), starting_clusters = starting.clusters)

combine lineages and calculate pseudotime
cds = combine_lineages(cds, 203)
cds_new = order_cells(cds_new, root_pr_nodes = "Y_203")

#compress into metacells
cds_new <- compress_lineages_v2(cds_new, start = 203)

#identify lineage-specific genes
meta = cds_new@colData
lineages = names(cds_new@lineages)
for(lineage in lineages){
  factor = 0.05
  N = 10000
  print(lineage)
  cds.sub = get_lineage_object(cds_new, lineage, 203, N = N)
  meta.sub = meta[colnames(cds.sub),]
  colData(cds.sub)$sex <- meta.sub$sex
  colData(cds.sub)$region_broad <- meta.sub$region_broad
  colData(cds.sub)$chemistry <- meta.sub$chemistry
  colData(cds.sub)$PMI <- meta.sub$PMI
  colData(cds.sub)$UMI <- meta.sub$nCount_RNA
  data = counts(cds.sub)
  rows=rownames(data)[rowSums(data> 0) > factor*ncol(data)]
  cds.sel = cds.sub[rows]
  cds_pr_test_res <- monocle3:::graph_test(cds.sel, neighbor_graph="principal_graph", verbose = T, cores = 24)
  save(cds_pr_test_res, file = paste("Moran_", lineage,".R", sep = ""))
  write.table(subset(cds_pr_test_res, q_value < 0.05), paste("pt_DGE_", lineage,".txt", sep = ""), sep ="\t", quote = F)
}

#############################
#analysis of macroglial cells (astrocytes, oligos)
#############################

######################
#analysis of microglia
######################

###########################
#analysis of vascular cells
###########################

