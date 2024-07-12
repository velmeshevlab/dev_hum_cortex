###########################
#combine all lineages
###########################
ExNeu = readRDS(file = "ExNeu_AL.RDS")
IN = readRDS(file = "IN_AL.RDS")
cds_comb = combine_objects(ExNeu, IN, "", "")
Glia = readRDS(file = "glia_AL.RDS")
cds_comb = combine_objects(cds_comb, Glia, "", "")
MG = readRDS(file = "MG_AL.RDS")
cds_comb = combine_objects(cds_comb, MG, "", "")
END = readRDS(file = "END_AL.RDS")
cds_comb = combine_objects(cds_comb, END, "", "")
PER = readRDS(file = "PER_AL.RDS")
cds_comb = combine_objects(cds_comb, PER, "", "")
saveRDS(cds_comb, file = "ALL_AL.RDS")

###########################
#find lineage-specific genes
###########################
lineages = names(cds_comb@lineages)
lineages = lineages[lineages!="SP"]
for(lineage in lineages){
print(lineage)
d = get_pt_exp(cds_comb, lineage, I = 0)
names = gsub(paste0(lineage, "__"), "", rownames(d))
input = paste0("fit = ","cds_comb","@expectation$", lineage)
eval(parse(text=input))
names = names[names %in% colnames(fit)]
res = sapply(names, AUC_window_sub, cds = cds_comb, lineage = lineage, comp_lineages = lineages[lineages != lineage], factor = 1.2, window_ratio = 0.01)
res = res[res != 0]
res = cbind(as.data.frame(res), rep(lineage, length(res)))
write.table(res, paste0(lineage, "_spec.txt"), quote = F, sep = "\t")
}

###########################
#classify genes based on expression trend and assign to biological age.
#final formatting
###########################
meta = read.data("meta.txt")
#excitatory neurons
for(lineage in lineages){
  print(lineage)
  res = get_max_age_v2(ExNeu, meta = meta, lineage = lineage, start = 1334)
  write.table(t(res), paste0("mode_", lineage,".txt"), sep = "\t", quote = F)
}
lineages = c("L5", "L6")
for(lineage in lineages){
  print(lineage)
  res = get_max_age_v2(ExNeu, meta = meta, lineage = lineage, start = 1320)
  write.table(t(res), paste0("mode_", lineage,".txt"), sep = "\t", quote = F)
}
out = format_lineage_out(c("L2_3", "L5_6_IT", "L4", "L5", "L6"))
write.table(out, "ExNeu_lineage.txt", sep = "\t", quote = F)

#interneurons
lineages = names(IN@lineages)
for(lineage in lineages){
print(lineage)
res = get_max_age_v2(IN, meta = meta, lineage = lineage, start = 203)
write.table(t(res), paste0("mode_", lineage,".txt"), sep = "\t", quote = F)
}
out = format_lineage_out(lineages)
write.table(out, "IN_lineage.txt", sep = "\t", quote = F)

#macroglial cells
lineages = names(Glia@lineages)
for(lineage in lineages){
print(lineage)
res = get_max_age_v2(Glia, meta = meta, lineage = lineage, start = 203)
write.table(t(res), paste0("mode_", lineage,".txt"), sep = "\t", quote = F)
}
out = format_lineage_out(lineages)
write.table(out, "glia_lineage.txt", sep = "\t", quote = F)

###########################
#find branch-specific genes
###########################
#Excitatory neurons
#ExNeu1
lineages = names(cds_comb@lineages)
lineages = lineages[lineages!="SP"]
target.lins = c("L2_3", "L4", "L5_6_IT")
res = get_branch_genes(cds_comb, target.lins, lineages)
write.table(res, "Ex1.txt", sep = "\t", quote = F)
target.lins = c("L2_3", "L4", "L5_6_IT")
genes = as.character(rownames(read.table("Ex1_spec.txt", sep = "\t")))
res = get_peak_age_branches(ExNeu, genes, target.lins, meta, 1334)
write.table(res, "mode_Ex1.txt", quote = F, sep = "\t")

#ExNeu2
lineages = names(cds_comb@lineages)
lineages = lineages[lineages!="SP"]
target.lins = c("L5", "L6")
res = get_branch_genes(cds_comb, target.lins, lineages)
write.table(res, "Ex2.txt", sep = "\t", quote = F)
target.lins = c("L5", "L6")
genes = as.character(rownames(read.table("Ex2_spec.txt", sep = "\t")))
res = get_peak_age_branches(ExNeu, genes, target.lins, meta, 1320)
write.table(res, "mode_Ex2.txt", quote = F, sep = "\t")

#ExNeu3
lineages = names(cds_comb@lineages)
lineages = lineages[lineages!="SP"]
target.lins = c("L2_3", "L5_6_IT")
res = get_branch_genes(cds_comb, target.lins, lineages)
write.table(res, "Ex3.txt", sep = "\t", quote = F)
target.lins = c("L2_3", "L5_6_IT")
genes = as.character(rownames(read.table("Ex3_spec.txt", sep = "\t")))
res = get_peak_age_branches(ExNeu, genes, target.lins, meta, 1334)
write.table(res, "mode_Ex3.txt", quote = F, sep = "\t")

out = format_branches_out(c("Ex1", "Ex2", "Ex3"))
write.table(out, "Ex_branches.txt", sep = "\t", quote = F)

#interneurons
#IN1
lineages = names(cds_comb@lineages)
lineages = lineages[lineages!="SP"]
target.lins = c("PV_MP","PV","SST","SST_RELN")
res = get_branch_genes(cds_comb, target.lins, lineages)
write.table(res, "IN1_spec.txt", sep = "\t", quote = F)
genes = as.character(rownames(read.table("IN1_spec.txt", sep = "\t")))
res = get_peak_age_branches(IN, genes, target.lins, meta, 203)
write.table(res, "mode_IN1.txt", quote = F, sep = "\t")

#IN2
lineages = names(cds_comb@lineages)
lineages = lineages[lineages!="SP"]
target.lins = c("VIP","CALB2","CCK","RELN","SV2C","NOS")
res = get_branch_genes(cds_comb, target.lins, lineages)
write.table(res, "IN2_spec.txt", sep = "\t", quote = F)
genes = as.character(rownames(read.table("IN2_spec.txt", sep = "\t")))
res = get_peak_age_branches(IN, genes, target.lins, meta, 203)
write.table(res, "mode_IN2.txt", quote = F, sep = "\t")

#IN3
lineages = names(cds_comb@lineages)
lineages = lineages[lineages!="SP"]
target.lins = c("SST","SST_RELN")
res = get_branch_genes(cds_comb, target.lins, lineages)
write.table(res, "IN3_spec.txt", sep = "\t", quote = F)
genes = as.character(rownames(read.table("IN3_spec.txt", sep = "\t")))
res = get_peak_age_branches(IN, genes, target.lins, meta, 203)
write.table(res, "mode_IN3.txt", quote = F, sep = "\t")

#IN4
lineages = names(cds_comb@lineages)
lineages = lineages[lineages!="SP"]
target.lins = c("PV_MP","PV")
res = get_branch_genes(cds_comb, target.lins, lineages)
write.table(res, "IN4_spec.txt", sep = "\t", quote = F)
genes = as.character(rownames(read.table("IN4_spec.txt", sep = "\t")))
res = get_peak_age_branches(IN, genes, target.lins, meta, 203)
write.table(res, "mode_IN4.txt", quote = F, sep = "\t")

#IN5
lineages = names(cds_comb@lineages)
lineages = lineages[lineages!="SP"]
target.lins = c("VIP","CALB2","CCK","RELN")
res = get_branch_genes(cds_comb, target.lins, lineages)
write.table(res, "IN5_spec.txt", sep = "\t", quote = F)
genes = as.character(rownames(read.table("IN5_spec.txt", sep = "\t")))
res = get_peak_age_branches(IN, genes, target.lins, meta, 203)
write.table(res, "mode_IN5.txt", quote = F, sep = "\t")

#IN6
lineages = names(cds_comb@lineages)
lineages = lineages[lineages!="SP"]
target.lins = c("SV2C","NOS")
res = get_branch_genes(cds_comb, target.lins, lineages)
write.table(res, "IN6_spec.txt", sep = "\t", quote = F)
genes = as.character(rownames(read.table("IN6_spec.txt", sep = "\t")))
res = get_peak_age_branches(IN, genes, target.lins, meta, 203)
write.table(res, "mode_IN6.txt", quote = F, sep = "\t")

#IN7
lineages = names(cds_comb@lineages)
lineages = lineages[lineages!="SP"]
target.lins = c("VIP","CALB2")
res = get_branch_genes(cds_comb, target.lins, lineages)
write.table(res, "IN7_spec.txt", sep = "\t", quote = F)
genes = as.character(rownames(read.table("IN7_spec.txt", sep = "\t")))
res = get_peak_age_branches(IN, genes, target.lins, meta, 203)
write.table(res, "mode_IN7.txt", quote = F, sep = "\t")

#IN8
lineages = names(cds_comb@lineages)
lineages = lineages[lineages!="SP"]
target.lins = c("CCK","RELN")
res = get_branch_genes(IN, target.lins, lineages)
write.table(res, "IN8_spec.txt", sep = "\t", quote = F)
genes = as.character(rownames(read.table("IN8_spec.txt", sep = "\t")))
res = get_peak_age_branches(IN, genes, target.lins, meta, 203)
write.table(res, "mode_IN8.txt", quote = F, sep = "\t")

out = format_branches_out(c("IN1", "IN2", "IN3", "IN4", "IN5", "IN6", "IN7", "IN8"))
write.table(out, "IN_branches.txt", sep = "\t", quote = F)

