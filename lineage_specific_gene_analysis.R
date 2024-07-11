###########################
#combine all lineages
###########################
l = load(file = "IN_AL.R")
IN = cds_new
l = load(file = "ExNeu_AL.R")
ExNeu = cds_new
cds_comb = combine_objects(ExNeu, IN, "", "")
l = load(file = "glia_AL.R")
Glia = cds_new
cds_comb = combine_objects(cds_comb, Glia, "", "")
l = load(file = "MG_AL.R")
MG = cds_new
cds_comb = combine_objects(cds_comb, MG, "", "")
l = load(file = "END_AL.R")
END = cds_new
cds_comb = combine_objects(cds_comb, END, "", "")
l = load(file = "PER_AL.R")
PER = cds_new
cds_comb = combine_objects(cds_comb, PER, "", "")

###########################
#find lineage-specific genes
###########################
l = load(file = "ALL_AL.R")
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
