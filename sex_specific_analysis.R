#######################################
Excitatry neurons
#######################################
ExNeu = readRDS(file = "ExNeu_AL.R"DS)
#isolate male cells and compress lineages
cells = rownames(colData(ExNeu)[colData(ExNeu)$sex == "Male",])
cds.male = ExNeu[,cells]
cds.male = get_lineage_object(cds.male, start = c(1334, 1320))
cds.male <- compress_lineage_v2(cds.male, "L5_6_IT", start = 1334, cores = 12)
cds.male <- compress_lineage_v2(cds.male, "SP", start = 1334, cores = 12)
cds.male <- compress_lineage_v2(cds.male, "L2_3", start = 1334, cores = 12)
cds.male <- compress_lineage_v2(cds.male, "L4", start = 1334, cores = 12)
cds.male <- compress_lineage_v2(cds.male, "L5", start = 1320, cores = 12)
cds.male <- compress_lineage_v2(cds.male, "L6", start = 1320, cores = 12)
#isolate female cells and compress lineages
cells = rownames(colData(cds_new)[colData(cds_new)$sex == "Female",])
cds.female = cds_new[,cells]
cds.female = get_lineage_object(cds.female, start = c(1334, 1320))
cds.female <- compress_lineage_v2(cds.female, "L5_6_IT", start = 1334, cores = 12)
cds.female <- compress_lineage_v2(cds.female, "SP", start = 1334, cores = 12)
cds.female <- compress_lineage_v2(cds.female, "L2_3", start = 1334, cores = 12)
cds.female <- compress_lineage_v2(cds.female, "L4", start = 1334, cores = 12)
cds.female <- compress_lineage_v2(cds.female, "L5", start = 1320, cores = 12)
cds.female <- compress_lineage_v2(cds.female, "L6", start = 1320, cores = 12)

#Run Moran's I test
#male
factor = 0.05
N = 10000
lineages = c("L2_3", "L4", "L5_6_IT", "SP")
for(lineage in lineages){
options(mc.cores = 1)
print(lineage)
cds.sub = get_lineage_object(cds.male, lineage, 1334, N = N)
data = counts(cds.sub)
rows=rownames(data)[rowSums(data> 0) > factor*ncol(data)]
cds.sel = cds.sub[rows]
print(paste0("Testing ", ncol(cds.sel)," cells and ", nrow(cds.sel), " genes"))
cds_pr_test_res <- monocle3:::graph_test(cds.sel, neighbor_graph="principal_graph", verbose = T, cores = 24)
save(cds_pr_test_res, file = paste("Moran_", lineage,"_Male.R", sep = ""))
write.table(subset(cds_pr_test_res, q_value < 0.05), paste("pt_DGE_", lineage,"_Male.txt", sep = ""), sep ="\t", quote = F)
}
lineages = c("L5", "L6")
for(lineage in lineages){
options(mc.cores = 1)
print(lineage)
cds.sub = get_lineage_object(cds.male, lineage, 1320, N = N)
data = counts(cds.sub)
rows=rownames(data)[rowSums(data> 0) > factor*ncol(data)]
cds.sel = cds.sub[rows]
print(paste0("Testing ", ncol(cds.sel)," cells and ", nrow(cds.sel), " genes"))
cds_pr_test_res <- monocle3:::graph_test(cds.sel, neighbor_graph="principal_graph", verbose = T, cores = 24)
save(cds_pr_test_res, file = paste("Moran_", lineage,"_Male.R", sep = ""))
write.table(subset(cds_pr_test_res, q_value < 0.05), paste("pt_DGE_", lineage,"_Male.txt", sep = ""), sep ="\t", quote = F)
}

#female
lineages = c("L2_3", "L4", "L5_6_IT", "SP")
for(lineage in lineages){
options(mc.cores = 1)
print(lineage)
cds.sub = get_lineage_object(cds.female, lineage, 1334, N = N)
data = counts(cds.sub)
rows=rownames(data)[rowSums(data> 0) > factor*ncol(data)]
cds.sel = cds.sub[rows]
print(paste0("Testing ", ncol(cds.sel)," cells and ", nrow(cds.sel), " genes"))
cds_pr_test_res <- monocle3:::graph_test(cds.sel, neighbor_graph="principal_graph", verbose = T, cores = 24)
save(cds_pr_test_res, file = paste("Moran_", lineage,"_Female.R", sep = ""))
write.table(subset(cds_pr_test_res, q_value < 0.05), paste("pt_DGE_", lineage,"_Female.txt", sep = ""), sep ="\t", quote = F)
}
lineages = c("L5", "L6")
for(lineage in lineages){
options(mc.cores = 1)
print(lineage)
cds.sub = get_lineage_object(cds.female, lineage, 1320, N = N)
data = counts(cds.sub)
rows=rownames(data)[rowSums(data> 0) > factor*ncol(data)]
cds.sel = cds.sub[rows]
print(paste0("Testing ", ncol(cds.sel)," cells and ", nrow(cds.sel), " genes"))
cds_pr_test_res <- monocle3:::graph_test(cds.sel, neighbor_graph="principal_graph", verbose = T, cores = 24)
save(cds_pr_test_res, file = paste("Moran_", lineage,"_Female.R", sep = ""))
write.table(subset(cds_pr_test_res, q_value < 0.05), paste("pt_DGE_", lineage,"_Female.txt", sep = ""), sep ="\t", quote = F)
}

#get sex-enriched lineage genes
exneu_sex = combine_objects(cds.male, cds.female, "_Male", "_Female")
cds_new = combine_objects(exneu_sex, in_sex, "", "")
get_sex_enriched_genes(cds_new, factor = 1.2, window_ratio = 0.01, I = 0.1)
                       
#classify genes
meta = read.table("meta.txt", row.names=1, header=T, sep="\t", check.names = F)
lineages = c("L2_3", "L4", "L5_6_IT", "SP")
for(lineage in lineages){
print(lineage)
res = get_max_age_v2(cds_new, meta = meta, lineage = lineage, start = 1334)
write.table(t(res), paste0("mod_", lineage,".txt"), sep = "\t", quote = F)
}
lineages = c("L5", "L6")
for(lineage in lineages){
print(lineage)
res = get_max_age_v2(cds_new, meta = meta, lineage = lineage, start = 1320)
write.table(t(res), paste0("mod_", lineage,".txt"), sep = "\t", quote = F)
}
