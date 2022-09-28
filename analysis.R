#import functions
source(file = "functions.R")

#UMAP ebmedding and clustering using standard Seurat pipeline
data <- CreateSeuratObject(counts = d)
meta = read.table("meta_ext.txt", row.names=1, header=T, sep="\t", check.names = F)
data <- AddMetaData(data, meta, col.name = NULL)
data <- NormalizeData(object = data)
data <- FindVariableFeatures(object = data)
data <- ScaleData(object = data)
data <- RunPCA(object = data)
#use elbow plot to select significant PCs
ElbowPlot(object = data, ndims = 40)
data <- FindNeighbors(object = data, dims = 1:N)
data <- FindClusters(object = data)
data <- RunUMAP(object = data, dims = 1:N)

#sex determination
agg = AggregateExpression(data, group.by = "individual")
agg = agg[[1]][c("XIST", "TSIX", "DDX3Y", "KDM5D", "USP9Y", "ZFY", "EIF1AY", "UTY"),]
agg = t(agg)
rownames(agg)[60] <- "GW19-1-17-20"
sex = unique(cbind(data@"meta.data"$individual, data@"meta.data"$sex_original))
sex = as.data.frame(sex)
rownames(sex) <- sex[,1]
sex = sex[rownames(agg),]
sex.exp = cbind(agg, sex)
colnames(sex.exp)[9] <- "individual"
colnames(sex.exp)[10] <- "sex"
#plot for each sex-specific gene
ggplot(data = sex.exp, aes(x = sex, y = XIST, fill = sex, label=individual)) + geom_violin(scale = "width", trim = FALSE) + geom_jitter(shape=16, position=position_jitter(0.2)) + geom_text(hjust=0, vjust=0)

##############################
##############################
#Analysis of specific lineages
##############################
##############################

#get the Seurat object (data) and metadata at https://cells.ucsc.edu/?ds=cortical-dev
##############################################
#analysis of excitatory neurons and astrocytes
##############################################

#select clusters and brain regions
sel.regions = c("BA9","BA24","PFC","BA46","FIC","FC","Cing","INS","cortex","BA13","BA22","cing","STG","temp","ACC")
sel.clusters = c("0","3","4","6","7","10","12","13","15","18","19","21","22","23","24","26")
cells = rownames(data@'meta.data'[data@'meta.data'$seurat_cluster %in% sel.clusters & data@'meta.data'$region %in% sel.regions,])
data.sub = subset(data, cells = cells)
#import Surat object and UMAP coordinates to monocle3
d = GetAssayData(object = data.sub, assay = "RNA", slot = "counts")
gene_annotation = as.data.frame(rownames(d))
rownames(gene_annotation) = rownames(d)
colnames(gene_annotation) = 'gene_short_name'
meta = data.sub@"meta.data"
cds = new_cell_data_set(d, cell_metadata = meta, gene_metadata = gene_annotation)
s.umap <- data.sub@"reductions"$umap[[]]
s.umap = s.umap[colnames(cds),]
reducedDims(cds)$"UMAP" <- s.umap
s.clusters = as.character(Idents(data.sub))
names(s.clusters) <- names(Idents(data.sub))
s.clusters = s.clusters[colnames(cds)]
cds@clusters$"UMAP"$"clusters" <- s.clusters
cds@clusters$UMAP$partitions <- cds@clusters$UMAP$clusters
cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions != "1"] <- "1"
#reconstruct trajectory graph
graph_control = setNames(list(3, 1, 1, FALSE, TRUE, FALSE, NULL, 10, 1e-5, 0.005, 0.01), c("euclidean_distance_ratio","geodesic_distance_ratio","minimal_branch_len","orthogonal_proj_tip","prune_graph","scale","rann.k","maxiter","eps","L1.gamma","L1.sigma"))
cds <- learn_graph(cds, use_partition = F, learn_graph_control = graph_control)
#import the monocle object to add new slots
cds = import_monocle(cds)
#connect nodes for trajectories not joined in the original graph
#plot nodes
Y <- cds@principal_graph_aux[["UMAP"]]$dp_mst
d = as.data.frame(t(Y))
ggplot(data=d, aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.01) + geom_text_repel(data=d, aes(x=UMAP_1, y=UMAP_2), label=rownames(d), size=0.3, hjust = 2, color = "red", max.overlaps = Inf, segment.size = 0.1) + monocle_theme_opts()
#connect nodes
#layers 6-IT/3-5/4
cds = connect_nodes(cds, "Y_2026", "Y_645")
#layers 6-IT/3-5
cds = connect_nodes(cds, "Y_1475", "Y_697")
#layer 4
cds = connect_nodes(cds, "Y_2032", "Y_653")
#layer 6-IT
cds = connect_nodes(cds, "Y_122", "Y_554")
#layer 6b
cds = connect_nodes(cds, "Y_602", "Y_585")
#layer 5-6-IT
cds = connect_nodes(cds, "Y_353", "Y_2045")

#select individual lineages
lineage = "L2_3"
start = 2203
end = 204
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("2", "1", "0", "5", "6", "24", "8", "23", "13", "29", "3", "4")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
lineage = "L4"
end = 163
inc.node = c("Y_2032")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("2", "1", "0", "5", "6", "24", "8", "23", "17", "28", "18")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
lineage = "L6_IT"
end = 133
inc.node = c("Y_1475", "Y_1996")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("2", "1", "0", "5", "6", "24", "8", "23", "17", "3", "14", "20")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
lineage = "L3_5"
end = 192
inc.node = c("Y_1475")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("2", "1", "0", "5", "6", "24", "8", "23", "17", "3", "14")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
lineage = "L5_6_IT"
start = 2203
end = 550
inc.node = c("Y_353", "Y_1304")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("2", "1", "0", "5", "6", "24", "8", "15", "30")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
lineage = "AST"
end = 56
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("2", "0", "10", "7", "11", "31")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
lineage = "L5"
end = 252
inc.node = c("Y_1459")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("2", "1", "0", "5", "6", "21", "19", "9", "12", "16_0", "26")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
lineage = "L6_CT"
end = 75
inc.node = c("Y_602")
cds<- isolate_graph(cds, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("2", "1", "0", "5", "6", "21", "19", "9", "27", "16_0", "16_1")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
lineage = "L6b"
end = 57
cds<- isolate_graph(cds, start, end, lineage)
sel.cluster = c("2", "1", "0", "5", "6", "21", "19", "9", "27", "22")
cds <- isolate_lineage(cds, lineage, sel_clusters = sel.cluster, cl = 4, N = 10)
cds = combine_lineages(cds, start)

#perform Moran's I test
#modify functions
trace('graph_test', edit = T, where = asNamespace("monocle3"))
#replace with "graph_test_lm.R"
trace('my.moran.test', edit = T, where = asNamespace("monocle3"))
#replace with "my.moran.test_lm.R"
#run the test
lineages = names(cds@lineages)
factor = 0.05
N = 10000
for(lineage in lineages){
  print(lineage)
  cds.sub = get_lineage_object(cds, lineage, 2203, N = N)
  data = counts(cds.sub)
  rows=rownames(data)[rowSums(data> 0) > factor*ncol(data)]
  cds.sel = cds.sub[rows]
  print(paste0("Testing ", ncol(cds.sel)," cells and ", nrow(cds.sel), " genes"))
  cds_pr_test_res <- monocle3:::graph_test(cds.sel, neighbor_graph="principal_graph", verbose = T, cores = 8)
  save(cds_pr_test_res, file = paste("Moran_", lineage,".R", sep = ""))
  write.table(subset(cds_pr_test_res, q_value < 0.05), paste("pt_DGE_", lineage,".txt", sep = ""), sep ="\t", quote = F)
}

#compress data
lineages = names(cds@lineages)
for(lineage in lineages){
  print(lineage)
  cds = compress_lineage(cds, lineage, 2203, N = 500, cores = 8)
}

#find lineage-specific genes
for(lineage in lineages){
  print(lineage)
  d = get_pt_exp(cds, lineage, I = 0)
  names = gsub(paste0(lineage, "__"), "", rownames(d))
  input = paste0("fit = ","cds","@expectation$", lineage)
  eval(parse(text=input))
  names = names[names %in% colnames(fit)]
  res = sapply(names, AUC_window_sub, cds = cds, lineage = lineage, comp_lineages = lineages[lineages != lineage], factor = 1.2, window_ratio = 0.01)
  res = res[res != 0]
  res = cbind(as.data.frame(res), rep(lineage, length(res)))
  write.table(res, paste0(lineage, "_spec.txt"), quote = F, sep = "\t")
}

#find branch-specific genes
#1
d = get_pt_exp(cds, "L2_3", I = 0)
names = gsub(paste0("L2_3", "__"), "", rownames(d))
fit = cds@expectation$L2_3
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L2_4", I = 0)
names = gsub(paste0("L2_4", "__"), "", rownames(d))
fit = cds@expectation$L2_4
names2 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L4", I = 0)
names = gsub(paste0("L4", "__"), "", rownames(d))
fit = cds@expectation$L4
names3 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L6_IT", I = 0)
names = gsub(paste0("L6_IT", "__"), "", rownames(d))
fit = cds@expectation$L6_IT
names4 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L5_6_IT", I = 0)
names = gsub(paste0("L5_6_IT", "__"), "", rownames(d))
fit = cds@expectation$L5_6_IT
names5 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L5", I = 0)
names = gsub(paste0("L5", "__"), "", rownames(d))
fit = cds@expectation$L5
names6 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L6b", I = 0)
names = gsub(paste0("L6b", "__"), "", rownames(d))
fit = cds@expectation$L6b
names7 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L6_CT", I = 0)
names = gsub(paste0("L6_CT", "__"), "", rownames(d))
fit = cds@expectation$L6_CT
names8 = names[names %in% colnames(fit)]
names = intersect(intersect(intersect(intersect(intersect(intersect(intersect(names1, names2),names3),names4),names5),names6),names7),names8)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "L2_3", comp_lineages = c("AST"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "L2_4", comp_lineages = c("AST"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res3 = sapply(names, AUC_window_sub, cds = cds, lineage = "L4", comp_lineages = c("AST"), factor = 1.2, window_ratio = 0.01)
res3 = res3[res3 != 0]
res4 = sapply(names, AUC_window_sub, cds = cds, lineage = "L6_IT", comp_lineages = c("AST"), factor = 1.2, window_ratio = 0.01)
res4 = res4[res4 != 0]
res5 = sapply(names, AUC_window_sub, cds = cds, lineage = "L5_6_IT", comp_lineages = c("AST"), factor = 1.2, window_ratio = 0.01)
res5 = res5[res5 != 0]
res6 = sapply(names, AUC_window_sub, cds = cds, lineage = "L5", comp_lineages = c("AST"), factor = 1.2, window_ratio = 0.01)
res6 = res6[res6 != 0]
res7 = sapply(names, AUC_window_sub, cds = cds, lineage = "L6b", comp_lineages = c("AST"), factor = 1.2, window_ratio = 0.01)
res7 = res7[res7 != 0]
res8 = sapply(names, AUC_window_sub, cds = cds, lineage = "L6_CT", comp_lineages = c("AST"), factor = 1.2, window_ratio = 0.01)
res8 = res8[res8 != 0]
res = intersect(intersect(intersect(intersect(intersect(intersect(intersect(names(res1), names(res2)),names(res3)),names(res4)),names(res5)),names(res6)),names(res7)),names(res8))
res = cbind(res1[res], res2[res], res3[res], res4[res], res5[res], res6[res], res7[res], res8[res])
colnames(res) <- c("L2_3", "L2_4", "L4", "L6_IT", "L5_6_IT", "L5", "L6b", "L6_CT")
write.table(res, "1_specific.txt", sep = "\t", quote = F)

#2
d = get_pt_exp(cds, "L6b", I = 0)
names = gsub(paste0("L6b", "__"), "", rownames(d))
fit = cds@expectation$L6b
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L6_CT", I = 0)
names = gsub(paste0("L6_CT", "__"), "", rownames(d))
fit = cds@expectation$L6_CT
names2 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L5", I = 0)
names = gsub(paste0("L5", "__"), "", rownames(d))
fit = cds@expectation$L5
names3 = names[names %in% colnames(fit)]
names = intersect(intersect(names1, names2),names3)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "L6b", comp_lineages = c("AST", "L2_3", "L2_4", "L4", "L6_IT",  "L5_6_IT"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "L5", comp_lineages = c("AST", "L2_3", "L2_4", "L4", "L6_IT",  "L5_6_IT"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res3 = sapply(names, AUC_window_sub, cds = cds, lineage = "L6_CT", comp_lineages = c("AST", "L2_3", "L2_4", "L4", "L6_IT",  "L5_6_IT"), factor = 1.2, window_ratio = 0.01)
res3 = res3[res3 != 0]
res = intersect(intersect(names(res1), names(res2)),names(res3))
res = cbind(res1[res], res2[res], res3[res])
colnames(res) <- c("L6b", "L6_CT", "L5")
write.table(res, "2.txt", sep = "\t", quote = F)

#3
d = get_pt_exp(cds, "L2_3", I = 0)
names = gsub(paste0("L2_3", "__"), "", rownames(d))
fit = cds@expectation$L2_3
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L2_4", I = 0)
names = gsub(paste0("L2_4", "__"), "", rownames(d))
fit = cds@expectation$L2_4
names2 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L4", I = 0)
names = gsub(paste0("L4", "__"), "", rownames(d))
fit = cds@expectation$L4
names3 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L6_IT", I = 0)
names = gsub(paste0("L6_IT", "__"), "", rownames(d))
fit = cds@expectation$L6_IT
names4 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L5_6_IT", I = 0)
names = gsub(paste0("L5_6_IT", "__"), "", rownames(d))
fit = cds@expectation$L5_6_IT
names5 = names[names %in% colnames(fit)]
names = intersect(intersect(intersect(intersect(names1, names2),names3),names4),names5)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "L2_3", comp_lineages = c("AST", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "L2_4", comp_lineages = c("AST", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res3 = sapply(names, AUC_window_sub, cds = cds, lineage = "L4", comp_lineages = c("AST", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res3 = res3[res3 != 0]
res4 = sapply(names, AUC_window_sub, cds = cds, lineage = "L6_IT", comp_lineages = c("AST", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res4 = res4[res4 != 0]
res5 = sapply(names, AUC_window_sub, cds = cds, lineage = "L5_6_IT", comp_lineages = c("AST", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res5 = res5[res5 != 0]
res = intersect(intersect(intersect(intersect(names(res1), names(res2)),names(res3)),names(res4)),names(res5))
res = cbind(res1[res], res2[res], res3[res], res4[res], res5[res])
colnames(res) <- c("L2_3", "L2_4", "L4", "L6_IT", "L5_6_IT")
write.table(res, "3_specific.txt", sep = "\t", quote = F)

#4
d = get_pt_exp(cds, "L2_3", I = 0)
names = gsub(paste0("L2_3", "__"), "", rownames(d))
fit = cds@expectation$L2_3
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L2_4", I = 0)
names = gsub(paste0("L2_4", "__"), "", rownames(d))
fit = cds@expectation$L2_4
names2 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L4", I = 0)
names = gsub(paste0("L4", "__"), "", rownames(d))
fit = cds@expectation$L4
names3 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L6_IT", I = 0)
names = gsub(paste0("L6_IT", "__"), "", rownames(d))
fit = cds@expectation$L6_IT
names4 = names[names %in% colnames(fit)]
names = intersect(intersect(intersect(names1, names2),names3),names4)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "L2_3", comp_lineages = c("AST", "L5_6_IT", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "L2_4", comp_lineages = c("AST", "L5_6_IT", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res3 = sapply(names, AUC_window_sub, cds = cds, lineage = "L4", comp_lineages = c("AST", "L5_6_IT", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res3 = res3[res3 != 0]
res4 = sapply(names, AUC_window_sub, cds = cds, lineage = "L6_IT", comp_lineages = c("AST", "L5_6_IT", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res4 = res4[res4 != 0]
res = intersect(intersect(intersect(names(res1), names(res2)),names(res3)),names(res4))
res = cbind(res1[res], res2[res], res3[res], res4[res])
colnames(res) <- c("L2_3", "L2_4", "L4", "L6_IT")
write.table(res, "4_specific.txt", sep = "\t", quote = F)

#5
d = get_pt_exp(cds, "L4", I = 0)
names = gsub(paste0("L4", "__"), "", rownames(d))
fit = cds@expectation$L4
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L2_4", I = 0)
names = gsub(paste0("L2_4", "__"), "", rownames(d))
fit = cds@expectation$L2_4
names2 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L6_IT", I = 0)
names = gsub(paste0("L6_IT", "__"), "", rownames(d))
fit = cds@expectation$L6_IT
names3 = names[names %in% colnames(fit)]
names = intersect(intersect(names1, names2),names3)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "L4", comp_lineages = c("AST", "L2_3", "L5_6_IT", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "L2_4", comp_lineages = c("AST", "L2_3", "L5_6_IT", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res3 = sapply(names, AUC_window_sub, cds = cds, lineage = "L6_IT", comp_lineages = c("AST", "L2_3", "L5_6_IT", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res3 = res2[res2 != 0]
res = intersect(intersect(names(res1), names(res2)), names(res3))
res = cbind(res1[res], res2[res], res3[res])
colnames(res) <- c("L4", "L2_4", "L6_IT")
write.table(res, "5_specific.txt", sep = "\t", quote = F)

#6
d = get_pt_exp(cds, "L2_4", I = 0)
names = gsub(paste0("L2_4", "__"), "", rownames(d))
fit = cds@expectation$L2_4
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L6_IT", I = 0)
names = gsub(paste0("L6_IT", "__"), "", rownames(d))
fit = cds@expectation$L6_IT
names2 = names[names %in% colnames(fit)]
names = intersect(names1, names2)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "L2_4", comp_lineages = c("AST", "L2_3", "L4", "L5_6_IT", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "L6_IT", comp_lineages = c("AST", "L2_3", "L4", "L5_6_IT", "L5", "L6b", "L6_CT"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res = intersect(names(res1), names(res2))
res = cbind(res1[res], res2[res])
colnames(res) <- c("L2_4", "L6_IT")
write.table(res, "6_specific.txt", sep = "\t", quote = F)

#7
d = get_pt_exp(cds, "L6b", I = 0)
names = gsub(paste0("L6b", "__"), "", rownames(d))
fit = cds@expectation$L6b
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "L6_CT", I = 0)
names = gsub(paste0("L6_CT", "__"), "", rownames(d))
fit = cds@expectation$L6_CT
names2 = names[names %in% colnames(fit)]
names = intersect(names1, names2)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "L6b", comp_lineages = c("AST", "L2_3", "L2_4", "L4", "L5", "L6_IT", "L5_6_IT"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "L6_CT", comp_lineages = c("AST", "L2_3", "L2_4", "L4", "L5", "L6_IT", "L5_6_IT"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res = intersect(names(res1), names(res2))
res = cbind(res1[res], res2[res])
colnames(res) <- c("L6b", "L6_CT")
write.table(res, "7_specific.txt", sep = "\t", quote = F)

#determine dynamic expression classification and age of max expression for lineage genes
meta = read.data("meta_ext.txt")
#lineage_genes is a text file with three columns: gene, score from AUC_window_sub function and lineage specificity
lin_genes = read.table("lineage_genes.txt", sep = "\t", header = T)
#specific lineages
lineages = unique(lin_genes$lineage)
out = matrix(nrow = 0, ncol = 5,0)
for(lineage in lineages){
  print(lineage)
  genes = lin_genes[lin_genes$lineage == lineage, 1]
  res = get_max_age_v2(cds, meta, genes, lineage, 2203)
  res = cbind(t(res), rep(lineage, ncol(res)))
  res = as.data.frame(res)
  res$gene <- rownames(res)
  out = rbind(out, res)
}
write.table(out, "lineages_age.txt", quote = F, sep = "\t")
#branches
#lineage_genes is a text file with three columns: gene, score from AUC_window_sub function and branch specificity
branch_genes = read.table("branch_genes.txt", sep = "\t", header = T)
genes = branch_genes[branch_genes$branch == "Ex_1", 1]
res = get_peak_age_branches(cds, genes, c("L2_3","L2_4","L4","L6_IT","L5_6_IT","L5","L6b","L6_CT"), meta, 2203)
write.table(res, "Ex1_age.txt", quote = F, sep = "\t")

genes = branch_genes[branch_genes$branch == "Ex_2", 1]
res = get_peak_age_branches(cds, genes, c("L6b","L6_CT","L5"), meta, 2203)
write.table(res, "Ex2_age.txt", quote = F, sep = "\t")

genes = branch_genes[branch_genes$branch == "Ex_3", 1]
res = get_peak_age_branches(cds, genes, c("L2_3","L2_4","L4","L6_IT","L5_6_IT"), meta, 2203)
write.table(res, "Ex3_age.txt", quote = F, sep = "\t")

genes = branch_genes[branch_genes$branch == "Ex_4", 1]
res = get_peak_age_branches(cds, genes, c("L2_3","L2_4","L4","L6_IT"), meta, 2203)
write.table(res, "Ex4_age.txt", quote = F, sep = "\t")

genes = branch_genes[branch_genes$branch == "Ex_5", 1]
res = get_peak_age_branches(cds, genes, c("L2_4","L4","L6_IT"), meta, 2203)
write.table(res, "Ex5_age.txt", quote = F, sep = "\t")

genes = branch_genes[branch_genes$branch == "Ex_6", 1]
res = get_peak_age_branches(cds, genes, c("L2_4", "L6_IT"), meta, 2203)
write.table(res, "Ex6_age.txt", quote = F, sep = "\t")

genes = branch_genes[branch_genes$branch == "Ex_7", 1]
res = get_peak_age_branches(cds, genes, c("L6b","L6_CT"), meta, 2203)
write.table(res, "Ex7_age.txt", quote = F, sep = "\t")

#sex-specific analysis
#select male cells
cells = rownames(colData(cds)[colData(cds)$sex == "Male",])
cds.male = cds[,cells]
cds.male = get_lineage_object(cds.male, start = 2203)
#select female cells
cells = rownames(colData(cds)[colData(cds)$sex == "Female",])
cds.female = cds[,cells]
cds.female = get_lineage_object(cds.female, start = 2203)
#compress lineages
cds.male = compress_lineages(cds.male, 2203, N = 500, cores = 6)
#compress lineages
cds.female = compress_lineages(cds.female, 2203, N = 500, cores = 6)
#combine
cds_new = combine_objects(cds.male, cds.female, "_Male", "_Female")
#Moran's I test
#modify regression formula in my.moran.test_lm.R to exp ~ region
lineages = names(cds.male@lineages)
factor = 0.05
N = 10000
options(mc.cores = 1)
for(lineage in lineages){
  print(lineage)
  cds.sub = get_lineage_object(cds.male, lineage, 2203, N = N)
  data = counts(cds.sub)
  rows=rownames(data)[rowSums(data> 0) > factor*ncol(data)]
  cds.sel = cds.sub[rows]
  print(paste0("Testing ", ncol(cds.sel)," cells and ", nrow(cds.sel), " genes"))
  cds_pr_test_res <- monocle3:::graph_test(cds.sel, neighbor_graph="principal_graph", verbose = T, cores = 8)
  save(cds_pr_test_res, file = paste("Moran_", lineage,"_Male.R", sep = ""))
  write.table(subset(cds_pr_test_res, q_value < 0.05), paste("pt_DGE_", lineage,"_Male.txt", sep = ""), sep ="\t", quote = F)
}
#find sex-specific genes
lineages = lineages = names(cds.male@lineages)
for(lin in lineages){
  print(lin)
  #male
  lineage = paste0(lin,"_Male")
  comp.lineage = paste0(lin,"_Female")
  d = get_pt_exp(cds_new, lineage, I = 0.1)
  names = gsub(paste0(lineage, "__"), "", rownames(d))
  cds_name = deparse(substitute(cds_new))
  input = paste0("fit = ",cds_name,"@expectation$", lineage)
  eval(parse(text=input))
  names1 = names[names %in% colnames(fit)]
  d = get_pt_exp(cds_new, comp.lineage, I = 0.1)
  names = gsub(paste0(comp.lineage, "__"), "", rownames(d))
  input = paste0("fit = ",cds_name,"@expectation$", comp.lineage)
  eval(parse(text=input))
  names2 = names[names %in% colnames(fit)]
  names = union(names1, names2)
  res = sapply(names, AUC_window_sub, cds = cds_new, lineage = lineage, comp_lineages = comp.lineage, factor = 1.2, window_ratio = 0.01)
  res = res[res != 0]
  res = cbind(as.data.frame(res), rep(lineage, length(res)))
  write.table(res, paste0(lineage, "_spec.txt"), quote = F, sep = "\t")
  #female
  lineage = paste0(lin,"_Female")
  comp.lineage = paste0(lin,"_Male")
  d = get_pt_exp(cds_new, lineage, I = 0.1)
  names = gsub(paste0(lineage, "__"), "", rownames(d))
  cds_name = deparse(substitute(cds_new))
  input = paste0("fit = ",cds_name,"@expectation$", lineage)
  eval(parse(text=input))
  names1 = names[names %in% colnames(fit)]
  d = get_pt_exp(cds_new, comp.lineage, I = 0.1)
  names = gsub(paste0(comp.lineage, "__"), "", rownames(d))
  input = paste0("fit = ",cds_name,"@expectation$", comp.lineage)
  eval(parse(text=input))
  names2 = names[names %in% colnames(fit)]
  names = union(names1, names2)
  res = sapply(names, AUC_window_sub, cds = cds_new, lineage = lineage, comp_lineages = comp.lineage, factor = 1.2, window_ratio = 0.01)
  res = res[res != 0]
  res = cbind(as.data.frame(res), rep(lineage, length(res)))
  write.table(res, paste0(lineage, "_spec.txt"), quote = F, sep = "\t")
}

#region-specific analysis
#select cells from specific regions
#frontal
cells = rownames(colData(cds)[colData(cds)$region %in% c("BA9", "PFC", "FC", "BA46"),])
cds.fc = cds[,cells]
cds.fc = get_lineage_object(cds.fc, start = 2203)
#cingulate
cells = rownames(colData(cds)[colData(cds)$region %in% c("BA24", "Cing", "cing", "ACC"),])
cds.cc = cds[,cells]
cds.cc = get_lineage_object(cds.cc, start = 2203)
#compress lineages
cds.male = compress_lineages(cds.fc, 2203, N = 500, cores = 6)
#compress lineages
cds.female = compress_lineages(cds.cc, 2203, N = 500, cores = 6)
#combine
cds_new = combine_objects(cds.fc, cds.cc, "_FC", "_CC")
#Moran's I test
#modify regression formula in my.moran.test_lm.R to exp ~ sex
lineages = names(cds.fc@lineages)
factor = 0.05
N = 10000
for(lineage in lineages){
  print(lineage)
  cds.sub = get_lineage_object(cds.fc, lineage, 2203, N = N)
  data = counts(cds.sub)
  rows=rownames(data)[rowSums(data> 0) > factor*ncol(data)]
  cds.sel = cds.sub[rows]
  print(paste0("Testing ", ncol(cds.sel)," cells and ", nrow(cds.sel), " genes"))
  cds_pr_test_res <- monocle3:::graph_test(cds.sel, neighbor_graph="principal_graph", verbose = T, cores = 8)
  save(cds_pr_test_res, file = paste("Moran_", lineage,".R", sep = ""))
  write.table(subset(cds_pr_test_res, q_value < 0.05), paste("pt_DGE_FC_", lineage,".txt", sep = ""), sep ="\t", quote = F)
}

#region-specific genes
lineages = names(cds.fc@lineages)
for(lin in lineages){
  print(lin)
  #frontal
  lineage = paste0(lin,"_FC")
  comp.lineage = paste0(lin,"_CC")
  d = get_pt_exp(cds, lineage, I = 0.1)
  names = gsub(paste0(lineage, "__"), "", rownames(d))
  cds_name = deparse(substitute(cds))
  input = paste0("fit = ",cds_name,"@expectation$", lineage)
  eval(parse(text=input))
  names1 = names[names %in% colnames(fit)]
  d = get_pt_exp(cds, comp.lineage, I = 0.1)
  names = gsub(paste0(comp.lineage, "__"), "", rownames(d))
  input = paste0("fit = ",cds_name,"@expectation$", comp.lineage)
  eval(parse(text=input))
  names2 = names[names %in% colnames(fit)]
  names = union(names1, names2)
  res = sapply(names, AUC_window_sub, cds = cds, lineage = lineage, comp_lineages = comp.lineage, factor = 1.2, window_ratio = 0.01)
  res = res[res != 0]
  res = cbind(as.data.frame(res), rep(lineage, length(res)))
  write.table(res, paste0(lineage, "_spec.txt"), quote = F, sep = "\t")
  #cingulate
  lineage = paste0(lin,"_CC")
  comp.lineage = paste0(lin,"_FC")
  d = get_pt_exp(cds, lineage, I = 0.1)
  names = gsub(paste0(lineage, "__"), "", rownames(d))
  cds_name = deparse(substitute(cds))
  input = paste0("fit = ",cds_name,"@expectation$", lineage)
  eval(parse(text=input))
  names1 = names[names %in% colnames(fit)]
  d = get_pt_exp(cds, comp.lineage, I = 0.1)
  names = gsub(paste0(comp.lineage, "__"), "", rownames(d))
  input = paste0("fit = ",cds_name,"@expectation$", comp.lineage)
  eval(parse(text=input))
  names2 = names[names %in% colnames(fit)]
  names = union(names1, names2)
  res = sapply(names, AUC_window_sub, cds = cds, lineage = lineage, comp_lineages = comp.lineage, factor = 1.2, window_ratio = 0.01)
  res = res[res != 0]
  res = cbind(as.data.frame(res), rep(lineage, length(res)))
  write.table(res, paste0(lineage, "_spec.txt"), quote = F, sep = "\t")
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
d = get_pt_exp(cds, "VIP", I = 0)
names = gsub(paste0("VIP", "__"), "", rownames(d))
fit = cds@expectation$VIP
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "CALB2", I = 0)
names = gsub(paste0("CALB2", "__"), "", rownames(d))
fit = cds@expectation$CALB2
names2 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "RELN", I = 0)
names = gsub(paste0("RELN", "__"), "", rownames(d))
fit = cds@expectation$RELN
names3 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "SV2C", I = 0)
names = gsub(paste0("SV2C", "__"), "", rownames(d))
fit = cds@expectation$SV2C
names4 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "NOS", I = 0)
names = gsub(paste0("NOS", "__"), "", rownames(d))
fit = cds@expectation$NOS
names5 = names[names %in% colnames(fit)]
names = intersect(intersect(intersect(intersect(names1, names2),names3),names4),names5)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "VIP", comp_lineages = c("PV_MP", "PV", "SST_MAF", "SST"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "CALB2", comp_lineages = c("PV_MP", "PV", "SST_MAF", "SST"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res3 = sapply(names, AUC_window_sub, cds = cds, lineage = "RELN", comp_lineages = c("PV_MP", "PV", "SST_MAF", "SST"), factor = 1.2, window_ratio = 0.01)
res3 = res3[res3 != 0]
res4 = sapply(names, AUC_window_sub, cds = cds, lineage = "SV2C", comp_lineages = c("PV_MP", "PV", "SST_MAF", "SST"), factor = 1.2, window_ratio = 0.01)
res4 = res3[res3 != 0]
res5 = sapply(names, AUC_window_sub, cds = cds, lineage = "NOS", comp_lineages = c("PV_MP", "PV", "SST_MAF", "SST"), factor = 1.2, window_ratio = 0.01)
res5 = res3[res3 != 0]
res = intersect(intersect(intersect(intersect(names(res1), names(res2)),names(res3)),names(res4)),names(res5))
res = cbind(res1[res], res2[res], res3[res], res4[res], res5[res])
colnames(res) <- c("VIP", "CALB2", "RELN", "SV2C", "NOS")
write.table(res, "CGE_specific.txt", sep = "\t", quote = F)

#(2)
d = get_pt_exp(cds, "PV_MP", I = 0)
names = gsub(paste0("PV_MP", "__"), "", rownames(d))
fit = cds@expectation$PV_MP
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "PV", I = 0)
names = gsub(paste0("PV", "__"), "", rownames(d))
fit = cds@expectation$PV
names2 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "SST_MAF", I = 0)
names = gsub(paste0("SST_MAF", "__"), "", rownames(d))
fit = cds@expectation$SST_MAF
names3 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "SST", I = 0)
names = gsub(paste0("SST", "__"), "", rownames(d))
fit = cds@expectation$SST
names4 = names[names %in% colnames(fit)]
names = intersect(intersect(intersect(names1, names2),names3),names4)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "PV_MP", comp_lineages = c("VIP", "CALB2", "RELN", "SV2C", "NOS"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "PV", comp_lineages = c("VIP", "CALB2", "RELN", "SV2C", "NOS"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res3 = sapply(names, AUC_window_sub, cds = cds, lineage = "SST_MAF", comp_lineages = c("VIP", "CALB2", "RELN", "SV2C", "NOS"), factor = 1.2, window_ratio = 0.01)
res3 = res3[res3 != 0]
res4 = sapply(names, AUC_window_sub, cds = cds, lineage = "SST", comp_lineages = c("VIP", "CALB2", "RELN", "SV2C", "NOS"), factor = 1.2, window_ratio = 0.01)
res4 = res3[res3 != 0]
res = intersect(intersect(intersect(names(res1), names(res2)),names(res3)),names(res4))
res = cbind(res1[res], res2[res], res3[res], res4[res])
colnames(res) <- c("PV_MP", "PV", "SST_MAF", "SST")
write.table(res, "MGE_specific.txt", sep = "\t", quote = F)

#(3)
d = get_pt_exp(cds, "VIP", I = 0)
names = gsub(paste0("VIP", "__"), "", rownames(d))
fit = cds@expectation$VIP
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "CALB2", I = 0)
names = gsub(paste0("CALB2", "__"), "", rownames(d))
fit = cds@expectation$CALB2
names2 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "RELN", I = 0)
names = gsub(paste0("RELN", "__"), "", rownames(d))
fit = cds@expectation$RELN
names3 = names[names %in% colnames(fit)]
names = intersect(intersect(names1, names2),names3)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "VIP", comp_lineages = c("PV_MP", "PV", "SST_MAF", "SST", "SV2C", "NOS"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "CALB2", comp_lineages = c("PV_MP", "PV", "SST_MAF", "SST", "SV2C", "NOS"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res3 = sapply(names, AUC_window_sub, cds = cds, lineage = "RELN", comp_lineages = c("PV_MP", "PV", "SST_MAF", "SST", "SV2C", "NOS"), factor = 1.2, window_ratio = 0.01)
res3 = res3[res3 != 0]
res = intersect(intersect(names(res1), names(res2)),names(res3))
res = cbind(res1[res], res2[res], res3[res])
colnames(res) <- c("VIP", "CALB2", "RELN")
write.table(res, "3_specific.txt", sep = "\t", quote = F)

#(4)
d = get_pt_exp(cds, "CALB2", I = 0)
names = gsub(paste0("CALB2", "__"), "", rownames(d))
fit = cds@expectation$CALB2
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "RELN", I = 0)
names = gsub(paste0("RELN", "__"), "", rownames(d))
fit = cds@expectation$RELN
names2 = names[names %in% colnames(fit)]
names = intersect(names1, names2)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "CALB2", comp_lineages = c("PV_MP", "PV", "SST_MAF", "SST", "SV2C", "NOS", "VIP"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "RELN", comp_lineages = c("PV_MP", "PV", "SST_MAF", "SST", "SV2C", "NOS", "VIP"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res = intersect(names(res1), names(res2))
res = cbind(res1[res], res2[res])
colnames(res) <- c("CALB2", "RELN")
write.table(res, "4_specific.txt", sep = "\t", quote = F)

#(5)
d = get_pt_exp(cds, "SV2C", I = 0)
names = gsub(paste0("SV2C", "__"), "", rownames(d))
fit = cds@expectation$SV2C
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "NOS", I = 0)
names = gsub(paste0("NOS", "__"), "", rownames(d))
fit = cds@expectation$NOS
names2 = names[names %in% colnames(fit)]
names = intersect(names1, names2)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "SV2C", comp_lineages = c("PV_MP", "PV", "SST_MAF", "SST", "CALB2", "RELN", "VIP"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "NOS", comp_lineages = c("PV_MP", "PV", "SST_MAF", "SST", "CALB2", "RELN", "VIP"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res = intersect(names(res1), names(res2))
res = cbind(res1[res], res2[res])
colnames(res) <- c("SV2C", "NOS")
write.table(res, "5_specific.txt", sep = "\t", quote = F)

#(6)
d = get_pt_exp(cds, "PV_MP", I = 0)
names = gsub(paste0("PV_MP", "__"), "", rownames(d))
fit = cds@expectation$PV_MP
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "PV", I = 0)
names = gsub(paste0("PV", "__"), "", rownames(d))
fit = cds@expectation$PV
names2 = names[names %in% colnames(fit)]
names = intersect(names1, names2)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "PV_MP", comp_lineages = c("VIP", "CALB2", "RELN", "SV2C", "NOS", "SST_MAF", "SST"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "PV", comp_lineages = c("VIP", "CALB2", "RELN", "SV2C", "NOS", "SST_MAF", "SST"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res = intersect(names(res1), names(res2))
res = cbind(res1[res], res2[res])
colnames(res) <- c("PV_MP", "PV")
write.table(res, "6_specific.txt", sep = "\t", quote = F)

#(7)
d = get_pt_exp(cds, "SST_MAF", I = 0)
names = gsub(paste0("SST_MAF", "__"), "", rownames(d))
fit = cds@expectation$SST_MAF
names1 = names[names %in% colnames(fit)]
d = get_pt_exp(cds, "SST", I = 0)
names = gsub(paste0("SST", "__"), "", rownames(d))
fit = cds@expectation$SST
names2 = names[names %in% colnames(fit)]
names = intersect(names1, names2)
res1 = sapply(names, AUC_window_sub, cds = cds, lineage = "SST_MAF", comp_lineages = c("VIP", "CALB2", "RELN", "SV2C", "NOS", "PV_MP", "PV"), factor = 1.2, window_ratio = 0.01)
res1 = res1[res1 != 0]
res2 = sapply(names, AUC_window_sub, cds = cds, lineage = "SST", comp_lineages = c("VIP", "CALB2", "RELN", "SV2C", "NOS", "PV_MP", "PV"), factor = 1.2, window_ratio = 0.01)
res2 = res2[res2 != 0]
res = intersect(names(res1), names(res2))
res = cbind(res1[res], res2[res])
colnames(res) <- c("SST_MAF", "SST")
write.table(res, "7_specific.txt", sep = "\t", quote = F)

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
