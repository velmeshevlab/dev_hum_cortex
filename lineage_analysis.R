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
#select ventral progenitors and interneurons
sel.regions = c("GE", "MGE", "LGE", "CGE")
cells1 = rownames(data@'meta.data'[data@'meta.data'$seurat_cluster == "5" & data@'meta.data'$region %in% sel.regions,])
#sel.clusters = c("16", "8", "14", "31", "32", "17", "27")
sel.clusters = c("16", "8", "14", "31", "32", "17")
cells2 = rownames(data@'meta.data'[data@'meta.data'$seurat_cluster %in% sel.clusters,])
cells = c(cells1, cells2)
data = subset(data, cells = cells)
data <- FindVariableFeatures(object = data)
data <- ScaleData(object = data)
data <- RunPCA(object = data)
pca=data@reductions$pca
write.table(pca@"stdev","sdev.txt",sep="\t")
N = 9
data <- FindNeighbors(object = data, dims = 1:N)
data <- FindClusters(object = data)
data <- RunUMAP(object = data, dims = 1:N)
#subcluster cluster 13
sel.cluster = "13"
data = FindSubCluster(data, sel.cluster, graph.name = "RNA_snn", resolution = 0.15)
Idents(data) <- "sub.cluster"
#create monocle object and construct trajectories
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
names(s.clusters) <- rownames(data@"meta.data")
s.clusters = s.clusters[colnames(cds)]
cds@clusters$"UMAP"$"clusters" <- s.clusters
cds@clusters$UMAP$partitions <- cds@clusters$UMAP$clusters
cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions != "1"] <- "1"
graph_control = setNames(list(3, 1, 1, FALSE, TRUE, FALSE, NULL, 10, 1e-5, 0.01, 0.01), c("euclidean_distance_ratio","geodesic_distance_ratio","minimal_branch_len","orthogonal_proj_tip","prune_graph","scale","rann.k","maxiter","eps","L1.gamma","L1.sigma"))
cds <- learn_graph(cds, use_partition = F, learn_graph_control = graph_control)

#import monocle object and join nodes
cds = import_monocle(cds)
node_plot(cds)
#CR+RELN
cds = connect_nodes(cds, "Y_492", "Y_718")
#SV2C+NOS
cds = connect_nodes(cds, "Y_1199", "Y_856")
#PV+PV_MP
cds = connect_nodes(cds, "Y_842", "Y_712")
#SST-MAF
cds = connect_nodes(cds, "Y_657", "Y_573")

#isolate lineages
#VIP
lineage = "VIP"
start = 93
end = 421
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("1", "10", "8", "3", "0", "11", "14", "7")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
#CALB2
lineage = "CALB2"
end = 502
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("1", "10", "8", "3", "0", "15", "14", "13_0")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
CALB2 = get_lineage_object(cds, "CALB2", 93)
#RELN
lineage = "RELN"
end = 1296
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("1", "10", "8", "3", "0", "15", "14", "13_1")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
RELN = get_lineage_object(cds, "RELN", 93)
#SV2C
lineage = "SV2C"
end = 1315
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("1", "10", "8", "3", "0", "21", "20", "17")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
SV2C = get_lineage_object(cds, "SV2C", 93)
#NOS
lineage = "NOS"
end = 354
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("1", "10", "8", "3", "0", "21", "20", "23")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
NOS = get_lineage_object(cds, "NOS", 93)
#PV_MP
lineage = "PV_MP"
end = 1293
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("1", "10", "8", "3", "2", "4", "18")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
PV_MP = get_lineage_object(cds, "PV_MP", 93)
#PV
lineage = "PV"
end = 1266
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("1", "10", "8", "3", "2", "4", "12", "22", "6")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
PV = get_lineage_object(cds, "PV", 93)
#SST_MAF
lineage = "SST_MAF"
end = 1324
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("1", "10", "8", "3", "2", "4", "16", "9", "24")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
SST_MAF = get_lineage_object(cds, "SST_MAF", 93)
#SST
lineage = "SST"
end = 400
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("1", "10", "8", "3", "2", "4", "16", "5")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
SST = get_lineage_object(cds, "SST", 93)

#combine trajectories
cds = combine_lineages(cds, 93)

#Moran's test and lineage-specific gene analysis are performed using the same command as for excitatory neurons

#branch-specific gene analysis
#(1)
res = get_branch_genes(cds, c("VIP", "CALB2", "RELN", "SV2C", "NOS"), names(cds@lineages))
write.table(res, "IN_1.txt", sep = "\t", quote = F)

#(2)
res = get_branch_genes(cds, c("PV_MP", "PV", "SST_MAF", "SST"), names(cds@lineages))
write.table(res, "IN_2.txt", sep = "\t", quote = F)

#(3)
res = get_branch_genes(cds, c("VIP", "CALB2", "RELN"), names(cds@lineages))
write.table(res, "IN_3.txt", sep = "\t", quote = F)

#(4)
res = get_branch_genes(cds, c("CALB2", "RELN"), names(cds@lineages))
write.table(res, "IN_4.txt", sep = "\t", quote = F)

#(5)
res = get_branch_genes(cds, c("SV2C", "NOS"), names(cds@lineages))
write.table(res, "IN_5.txt", sep = "\t", quote = F)

#(6)
res = get_branch_genes(cds, c("PV_MP", "PV"), names(cds@lineages))
write.table(res, "IN_6.txt", sep = "\t", quote = F)

#(7)
res = get_branch_genes(cds, c("SST_MAF", "SST"), names(cds@lineages))
write.table(res, "IN_7.txt", sep = "\t", quote = F)

#determine dynamic expression classification and age of max expression for lineage genes
#lineage-specific analysis is performed using the same command as for excitatory neurons
#classification of branch-specific genes
branch_genes = read.table("branch_genes.txt", sep = "\t", header = T)
genes = branch_genes[branch_genes$branch == "IN_1", 1]
res = get_peak_age_branches(cds, genes, c("VIP", "CALB2", "RELN", "SV2C", "NOS"), meta, 93)
write.table(res, "IN1_age.txt", quote = F, sep = "\t")

genes = branch_genes[branch_genes$branch == "IN_2", 1]
res = get_peak_age_branches(cds, genes, c("PV_MP", "PV", "SST_MAF", "SST"), meta, 93)
write.table(res, "IN2_age.txt", quote = F, sep = "\t")

genes = branch_genes[branch_genes$branch == "IN_3", 1]
res = get_peak_age_branches(cds, genes, c("VIP", "CALB2", "RELN"), meta, 93)
write.table(res, "IN3_age.txt", quote = F, sep = "\t")

genes = branch_genes[branch_genes$branch == "IN_4", 1]
res = get_peak_age_branches(cds, genes, c("CALB2", "RELN"), meta, 93)
write.table(res, "IN4_age.txt", quote = F, sep = "\t")

branch_genes = read.table("branch_genes.txt", sep = "\t", header = T)
genes = branch_genes[branch_genes$branch == "IN_5", 1]
res = get_peak_age_branches(cds, genes, c("SV2C", "NOS"), meta, 93)
write.table(res, "IN5_age.txt", quote = F, sep = "\t")

genes = branch_genes[branch_genes$branch == "IN_6", 1]
res = get_peak_age_branches(cds, genes, c("PV_MP", "PV"), meta, 93)
write.table(res, "IN6_age.txt", quote = F, sep = "\t")

genes = branch_genes[branch_genes$branch == "IN_7", 1]
res = get_peak_age_branches(cds, genes, c("SST_MAF", "SST"), meta, 93)
write.table(res, "IN7_age.txt", quote = F, sep = "\t")

#sex-specific expression
#male
cells = rownames(colData(cds)[colData(cds)$sex == "Male",])
cds.male = cds[,cells]
cds.male = get_lineage_object(cds.male, start = 1041)
#female
cells = rownames(colData(cds)[colData(cds)$sex == "Female",])
cds.female = cds[,cells]
cds.female = get_lineage_object(cds.female, start = 1041)
#compress lineages
cds.male = compress_lineages(cds.male, 1041, N = 500, cores = 6)
#compress lineages
cds.female = compress_lineages(cds.female, 1041, N = 500, cores = 6)
#combine
cds_new = combine_objects(cds.male, cds.female, "Male", "Female")

#downstream sex-specific analysis is performed using the same command as for excitatory neurons

#analysis of region-specific genes is performed using the same command as for excitatory neurons

#############################
#analysis of oligodendrocytes
#############################
#subset the clusters with ventral/dorsal progenitors, OPCs, and oligodendrocytes
subset2 <- subset(data, idents = c("18", "11", "9", "5"))
subset2 <- FindVariableFeatures(object = subset2)
#run PCA. Select significant PCs based on a scree plot. Look for the last point before the plot becomes flat
subset2 <- RunPCA(object = subset2)
ElbowPlot(subset2, ndims = 30)
#cluster data
subset2 <- FindNeighbors(object = subset2, dims = 1:N)
subset2 <- FindClusters(object = subset2)
#embed in UMAP coordinates
subset2 <- RunUMAP(object = subset2, dims = 1:N)
#monocle3 trajectory
#data is a Seurat object after clustering and UMAP
#get gene counts, gene names and metadata
d = GetAssayData(object = subset2, assay = "RNA", slot = "counts")
gene_annotation = as.data.frame(rownames(d))
rownames(gene_annotation) = rownames(d)
colnames(gene_annotation) = 'gene_short_name'
meta = subset2@"meta.data"
#create monocle object
cds = new_cell_data_set(d, cell_metadata = meta, gene_metadata = gene_annotation)
#import Seurat umap to the monocle object
s.umap <- subset2@"reductions"$umap[[]]
s.umap = s.umap[colnames(cds),]
reducedDims(cds)$"UMAP" <- s.umap
#import Seurat clusters to the monocle object
s.clusters = as.character(Idents(subset2))
names(s.clusters) <- names(Idents(subset2))
s.clusters = s.clusters[colnames(cds)]
cds@clusters$"UMAP"$"clusters" <- s.clusters
cds@clusters$UMAP$partitions <- cds@clusters$UMAP$clusters
cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions != "1"] <- "1"
#learn trajectory graph
cds <- learn_graph(cds, use_partition = F)
#Select the lineage trajectory
lineage = "Oligo"
cds <- isolate_graph(cds, 224, 109, lineage)
#analysis of lineage-specific, sex- and region-enriched genes is performed using the code as for excitatory neurons

######################
#analysis of microglia
######################
#select microglia and recluster
cells = rownames(data@'meta.data'[data@'meta.data'$seurat_cluster == c("25"),])
data = subset(data, cells = cells)
meta = read.data("meta_ext.txt")
meta = meta[colnames(data),]
data <- AddMetaData(data, meta)
data <- FindVariableFeatures(object = data)
data <- ScaleData(object = data)
data <- RunPCA(object = data)
pca=data@reductions$pca
write.table(pca@"stdev","sdev.txt",sep="\t")
N = 6
data <- FindNeighbors(object = data, dims = 1:N)
data <- FindClusters(object = data)
data <- RunUMAP(object = data, dims = 1:N)
#import to monocle and reconstruct trajectories
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
graph_control = setNames(list(1, 1/3, 10, FALSE, TRUE, FALSE, 25, 10, 1e-5, 1, 0.01), c("euclidean_distance_ratio","geodesic_distance_ratio","minimal_branch_len","orthogonal_proj_tip","prune_graph","scale","rann.k","maxiter","eps","L1.gamma","L1.sigma"))
cds <- learn_graph(cds, use_partition = F, learn_graph_control = graph_control)
#isolate lineages
cds = import_monocle(cds)
node_plot(cds, filter = F)
ggsave(file="graph.png",width = 14, height = 12, units = "in",  dpi = 600)
#MG_1
start = 266
end = 289
lineage = "MG_1"
inc.node = c("Y_360", "Y_163", "Y_135")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
cds <- isolate_lineage(cds, lineage, cl = 4, N = 5)
#MG_2
start = 266
end = 289
lineage = "MG_2"
inc.node = c("Y_306", "Y_91", "Y_157")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
cds <- isolate_lineage(cds, lineage, cl = 4, N = 5)
cds = combine_lineages(cds, 266)

#analysis of lineage-specific, sex- and region-enriched genes is performed using the code as for excitatory neurons

###########################
#analysis of vascular cells
###########################
#subset cells and cluster
cells = rownames(data@'meta.data'[data@'meta.data'$seurat_cluster == c("20"),])
data = subset(data, cells = cells)
meta = read.data("meta_ext.txt")
meta = meta[colnames(data),]
data <- AddMetaData(data, meta)
data <- FindVariableFeatures(object = data)
data <- ScaleData(object = data)
data <- RunPCA(object = data)
pca=data@reductions$pca
write.table(pca@"stdev","sdev.txt",sep="\t")
N = 7
data <- FindNeighbors(object = data, dims = 1:N)
data <- FindClusters(object = data)
data <- RunUMAP(object = data, dims = 1:N)
#subset endothelial cells and pericytes 
cells = names(Idents(data)[Idents(data) %in% c("0", "1", "2", "3", "5", "6", "7", "8", "10")])
data = subset(data, cells = cells)
#import to monocle and reconstruct trajectories
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
graph_control = setNames(list(2, 1/3, 10, FALSE, TRUE, FALSE, 25, 10, 1e-5, 1, 0.01), c("euclidean_distance_ratio","geodesic_distance_ratio","minimal_branch_len","orthogonal_proj_tip","prune_graph","scale","rann.k","maxiter","eps","L1.gamma","L1.sigma"))
cds <- learn_graph(cds, use_partition = F, learn_graph_control = graph_control)

cds = import_monocle(cds)
node_plot(cds, filter = F)
ggsave(file="graph.png",width = 14, height = 12, units = "in",  dpi = 600)
#endothelial cells
start = 386
end = 77
lineage = "END_1"
inc.node = c("Y_461", "Y_343", "Y_61")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
cds <- isolate_lineage(cds, lineage, cl = 4, N = 5)
#pericytes
start = 397
end = 144
lineage = "PER"
inc.node = c("Y_287")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
cds <- isolate_lineage(cds, lineage, cl = 4, N = 5)
#activated endothelial cells
start = 428
end = 77
lineage = "END_2"
cds<- isolate_graph(cds, start, end, lineage)
cds <- isolate_lineage(cds, lineage, cl = 4, N = 5)
cds = combine_lineages(cds, c(386, 397, 428))

#analysis of lineage-specific, sex- and region-enriched genes is performed using the code as for excitatory neurons
