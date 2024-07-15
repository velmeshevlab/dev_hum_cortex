#' @export
setClass("cell_data_set_ext", contains = "cell_data_set", slots=c(graphs = "list", lineages="list", expression="list", expectation="list", pseudotime="list")) -> cell_data_set_ext

#' @export
import_monocle <-function(cds){
cds <- as(cds,"cell_data_set_ext")
return(cds)
}

monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(plot.title = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(axis.line.y = element_line(size=1, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

#' @export
combine_objects <- function(obj1, obj2, name1, name2){
  cds_new = new("cell_data_set_ext")
  #cds_new@'preprocess_aux'<-obj1@'preprocess_aux'
  cds_new@'reduce_dim_aux'<-obj1@'reduce_dim_aux'
  cds_new@'principal_graph_aux'<-obj1@'principal_graph_aux'
  cds_new@'principal_graph'<-obj1@'principal_graph'
  cds_new@'clusters'<-obj1@'clusters'
  cds_new@'int_elementMetadata'<-obj1@'int_elementMetadata'
  cds_new@'int_colData'<-obj1@'int_colData'
  cds_new@'int_metadata'<-obj1@'int_metadata'
  cds_new@'rowRanges'<-obj1@'rowRanges'
  cds_new@'colData'<-obj1@'colData'
  cds_new@'assays'<-obj1@'assays'
  cds_new@'NAMES'<-obj1@'NAMES'
  cds_new@'elementMetadata'<-obj1@'elementMetadata'
  cds_new@'metadata'<-obj1@'metadata'
  cds_new@'graphs'<-c(obj1@'graphs', obj2@'graphs')
  cds_new@'lineages'<-c(obj1@'lineages', obj2@'lineages')
  cds_new@'expression'<-c(obj1@'expression', obj2@'expression')
  cds_new@'expectation'<-c(obj1@'expectation', obj2@'expectation')
  cds_new@'pseudotime'<-c(obj1@'pseudotime', obj2@'pseudotime')
  names(cds_new@'graphs') <- c(paste(names(obj1@'graphs'), name1, sep = ""), paste(names(obj2@'graphs'), name2, sep = ""))
  names(cds_new@'lineages') <- c(paste(names(obj1@'lineages'), name1, sep = ""), paste(names(obj2@'lineages'), name2, sep = ""))
  names(cds_new@'expression') <- c(paste(names(obj1@'expression'), name1, sep = ""), paste(names(obj2@'expression'), name2, sep = ""))
  names(cds_new@'expectation') <- c(paste(names(obj1@'expectation'), name1, sep = ""), paste(names(obj2@'expectation'), name2, sep = ""))
  names(cds_new@'pseudotime') <- c(paste(names(obj1@'pseudotime'), name1, sep = ""), paste(names(obj2@'pseudotime'), name2, sep = ""))
  cds_new
  }

#' @export
combine_lineages <- function(cds, start){
  cds_name = deparse(substitute(cds))
  lineage = names(cds@lineages)[1]
  input = paste0(cds_name, "@graphs$", lineage)
  if(length(names(cds@lineages)) > 1){
    for(lineage in names(cds@lineages)[2:length(names(cds@lineages))]){
      input = paste0(input, ",", cds_name,"@graphs$", lineage)
    }
    input = paste0("igraph::union(", input, ")")
  }
  g = eval(parse(text=input))
  nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
  principal_graph(cds)[["UMAP"]] <- g
  cds@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP[,names(V(g))]
  cells_UMAP = as.data.frame(reducedDims(cds)["UMAP"])
  closest_vertex = apply(cells_UMAP[,c("UMAP_1", "UMAP_2")], 1, calculate_closest_vertex, nodes = as.matrix(nodes_UMAP[,names(V(g))]))
  closest_vertex = as.data.frame(closest_vertex)
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex
  source_url("https://raw.githubusercontent.com/cole-trapnell-lab/monocle3/master/R/learn_graph.R")
  cds <- project2MST(cds, project_point_to_line_segment, F, T, "UMAP", nodes_UMAP[,names(V(g))])
  cds <- order_cells(cds, root_pr_nodes = as.character(paste0("Y_",start)))
  return(cds)
}

#' @export
node_plot <- function(cds, filter = F, N = 50, size = 0.5){
Y <- cds@principal_graph_aux[["UMAP"]]$dp_mst
d = as.data.frame(t(Y))
if(filter == T){
g = principal_graph(cds)[["UMAP"]]
dd = degree(g)
names1 = names(dd[dd > 2 | dd == 1])
names2 = names(dd[dd == 2])
names2 = sample(names2, length(names2)/N, replace = F)
d.f = d[c(names1, names2),]
ggplot(data=d, aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.01) + geom_text_repel(data=d.f, aes(x=UMAP_1, y=UMAP_2), label=rownames(d.f), size=size, hjust = 2, color = "red", max.overlaps = Inf, segment.size = 0.1) + monocle_theme_opts()
}
else{
ggplot(data=d, aes(x=UMAP_1, y=UMAP_2)) + geom_point(size=0.01) + geom_text(data=d, aes(x=UMAP_1, y=UMAP_2), label=rownames(d), size=size, hjust = 1, color = "red") + monocle_theme_opts()
}
}

#' @export
path.distance <- function(path){
dists=c()
for(i in 2:nrow(path)){
x1 = path[i-1,1]
y1 = path[i-1,2]
x2 = path[i,1]
y2 = path[i,2]
d.x = x2 - x1
d.y = y2 - y1
dist = sqrt(d.x*d.x + d.y*d.y)
dists = append(dist,dists)
}
return(mean(dists))
}

cell.selector_sub2 <- function(cell, coords, r){
x2 = cell[1]
y2 = cell[2]
d.x = x2 - coords[1]
d.y = y2 - coords[2]
dist = sqrt(d.x*d.x + d.y*d.y)
if(dist <= r){
return(TRUE)
}
else{
return(FALSE)
}
}

#' @export
selector_sub <- function(node, cells, r){
x1 = node[1]
y1 = node[2]
res = apply(cells, 1, cell.selector_sub2, coords = c(x1, y1), r = r, simplify = T)
res = names(res[res == TRUE])
return(res)
}

#' @export
cell.selector <- function(path, cells, r, cl){
sel.cells = c()
sel.cells = pbapply(path, 1, selector_sub, cells = cells, r = r, cl = cl, simplify = T)
return(unique(unlist(sel.cells)))
}

#' @export
make_graph <- function(sub.graph){
edges = names(sub.graph)
start.edges = c()
end.edges = c()
for(i in 1:(length(edges)-1)){
start.edges = append(start.edges, edges[i])
end.edges = append(end.edges, edges[i+1])
}
d = cbind(start.edges, end.edges)
g = graph_from_data_frame(d, directed = F)
return(g)
}

#' @export
included <- function(graph, include_nodes){
all(include_nodes %in% names(graph))
}

isolate_graph_sub <- function(cds, start, end, lineage, include_nodes = NULL){
#get lineage graph
reduction_method = "UMAP"
graph = cds@principal_graph[[reduction_method]]
#select cells that are 1) progenitor cells from the region of interest (MGE, CGE) or 2) lineage-committed cells
sub.graph = all_simple_paths(graph, paste0("Y_", start), paste0("Y_", end))
if(length(include_nodes) > 0){
sub.graph = sub.graph[sapply(sub.graph, included, include_nodes = include_nodes)]
}
lengths = lengths(sub.graph)
#get the shortest path
n = which(lengths==min(lengths))[1]
sub.graph = sub.graph[[n]]
}

#' @export
isolate_graph <- function(cds, start, end, lineage, include_nodes = NULL){
#get lineage graph
cds_name = deparse(substitute(cds))
sub.graph = isolate_graph_sub(cds, start, end, lineage, include_nodes = include_nodes)
input = paste0(cds_name, "@graphs$", lineage, " <- make_graph(sub.graph)")
eval(parse(text=input))
eval(parse(text=paste0("return(", cds_name, ")")))
}

get_lineage_object <- function(cds, lineage = FALSE, start, N = FALSE)
{
cds_name = deparse(substitute(cds))
if(lineage != FALSE){
input = paste0("sub.graph = ",cds_name,"@graphs$", lineage)
eval(parse(text=input))
input = paste0("sel.cells = ",cds_name,"@lineages$", lineage)
eval(parse(text=input))
}
else{
sel.cells = colnames(cds)
}
sel.cells = sel.cells[sel.cells %in% colnames(cds)]
nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
if(N != FALSE){
if(N < length(sel.cells)){
sel.cells = sample(sel.cells, N)
}
}
#subset the moncole object
cds_subset = cds[,sel.cells]
#set the graph, node and cell UMAP coordinates
if(lineage == FALSE){
sub.graph = principal_graph(cds_subset)[["UMAP"]]
}
principal_graph(cds_subset)[["UMAP"]] <- sub.graph
cds_subset@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP[,names(V(sub.graph))]
cds_subset@clusters[["UMAP"]]$partitions <- cds_subset@clusters[["UMAP"]]$partitions[colnames(cds_subset)]
#recalculate closest vertex for the selected cells
cells_UMAP = as.data.frame(reducedDims(cds_subset)["UMAP"])
closest_vertex = apply(cells_UMAP[,c("UMAP_1", "UMAP_2")], 1, calculate_closest_vertex, nodes = as.matrix(nodes_UMAP[,names(V(sub.graph))]))
closest_vertex = as.data.frame(closest_vertex)
cds_subset@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex <- closest_vertex
source_url("https://raw.githubusercontent.com/cole-trapnell-lab/monocle3/master/R/learn_graph.R")
cds_subset <- project2MST(cds_subset, project_point_to_line_segment, F, T, "UMAP", nodes_UMAP[,names(V(sub.graph))])
cds_subset <- order_cells(cds_subset, root_pr_nodes = c(paste0("Y_", as.character(start))))
return(cds_subset)
}

isolate_lineage_sub <- function(cds, lineage, sel_clusters = NULL, start_regions = NULL, starting_clusters = NULL, subset = FALSE, N = 5, cl = 1){
  cds_name = deparse(substitute(cds))
  input = paste0("sub.graph = ",cds_name,"@graphs$", lineage)
  eval(parse(text=input))
  nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
  if(subset == F){
    nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names(V(sub.graph))]))
  }
  else{
    g = principal_graph(cds)[["UMAP"]]
    dd = degree(g)
    names1 = names(dd[dd > 2 | dd == 1])
    names2 = names(dd[dd == 2])
    names2 = sample(names2, length(names2)/subset, replace = F)
    names = c(names1, names2)
    names = intersect(names(V(sub.graph)), names)
    nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names]))
  }
  #select cells along the graph
  mean.dist = path.distance(nodes_UMAP.sub)
  r = mean.dist*N
  cells_UMAP = as.data.frame(reducedDims(cds)["UMAP"])
  cells_UMAP = cells_UMAP[,c("UMAP_1", "UMAP_2")]
  sel.cells = cell.selector(nodes_UMAP.sub, cells_UMAP, r, cl = cl)
  #only keep cells in the progenitor and lineage-specific clusters
  sel.cells1 = c()
  sel.cells2 = sel.cells
  if(length(starting_clusters) > 0){
    sel.cells1 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% starting_clusters])
  }
  if(length(start_regions) > 0){
    sel.cells1 = sel.cells1[sel.cells1 %in% rownames(cds@colData[cds@colData$region %in% start_regions,])]
  }
  if(length(sel_clusters) > 0){
    sel.cells2 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% sel_clusters])
  }
  cells = unique(c(sel.cells1, sel.cells2))
  sel.cells = sel.cells[sel.cells %in% cells]
  return(sel.cells)
}

#' @export
isolate_lineage <- function(cds, lineage, sel_clusters = NULL, start_regions = F, starting_clusters = F, subset = FALSE, N = 5, cl = 1){
cds_name = deparse(substitute(cds))
sel.cells = isolate_lineage_sub(cds, lineage, sel_clusters = sel_clusters, start_regions = start_regions, starting_clusters = starting_clusters, subset = subset, N = N, cl = cl)
input = paste0(cds_name, "@lineages$", lineage, " <- sel.cells")
eval(parse(text=input))
return(cds)
}

#' @export
isolate_fork <- function(cds, start, end1, end2, lineage, include_nodes = F, sel_clusters = F, start_regions = F, starting_clusters = F, subset = FALSE, N = 5, cl = 1){
#get lineage graphs
graph1 = isolate_graph_sub(cds, start, end1, paste0(lineage, "_fork_1"), include_nodes = include_nodes)
graph1 = make_graph(graph1)
input = paste0("cds@graphs$", paste0(lineage, "_fork_1"), " <- graph1")
eval(parse(text=input))
graph2 = isolate_graph_sub(cds, start, end2, paste0(lineage, "_fork_2"), include_nodes = include_nodes)
graph2 = make_graph(graph2)
input = paste0("cds@graphs$", paste0(lineage, "_fork_2"), " <- graph2")
eval(parse(text=input))
graph = igraph::union(graph1, graph2)
input = paste0("cds@graphs$", paste0(lineage, "_fork"), " <- graph")
eval(parse(text=input))
cells1 = isolate_lineage_sub(cds, paste0(lineage, "_fork_1"), sel_clusters = sel_clusters, start_regions = start_regions, starting_clusters = starting_clusters, subset = subset, N = N, cl = cl)
input = paste0("cds@lineages$", paste0(lineage, "_fork_1"), " <- cells1")
eval(parse(text=input))
cells2 = isolate_lineage_sub(cds, paste0(lineage, "_fork_2"), sel_clusters = sel_clusters, start_regions = start_regions, starting_clusters = starting_clusters, subset = subset, N = N, cl = cl)
input = paste0("cds@lineages$", paste0(lineage, "_fork_2"), " <- cells2")
eval(parse(text=input))
cells = unique(c(cells1, cells2))
input = paste0("cds@lineages$", paste0(lineage, "_fork"), " <- cells")
eval(parse(text=input))
return(cds)
}

#' @export
calculate_closest_vertex <- function(cells, nodes){
new.pos = as.numeric(cells)
nearest.idx <- which.min(colSums((nodes - new.pos)^2))
out = as.integer(gsub("Y_", "", names(nearest.idx)))
}

#' @export
connect_nodes <- function(cds, node1, node2, add_node = F){
graph.old = cds@principal_graph[["UMAP"]]
if(add_node == F){
graph.new <- add_edges(graph.old, c(node1, node2))
}
else{
node_coords = cds@principal_graph_aux[["UMAP"]]$dp_mst
node_X = (node_coords[1,node1] + node_coords[1,node2])/2
node_Y = (node_coords[2,node1] + node_coords[2,node2])/2
new_name = paste0("Y_", as.character(length(names(V(graph.old)))+1))
node_coords = as.data.frame(c(node_X, node_Y))
colnames(node_coords) = new_name
rownames(node_coords) = c("UMAP_1", "UMAP_2")
cds@principal_graph_aux[["UMAP"]]$dp_mst <- cbind(cds@principal_graph_aux[["UMAP"]]$dp_mst, node_coords)
graph.new <- add_vertices(graph.old, 1,attr = list(name = new_name))
graph.new <- add_edges(graph.new, c(node1, new_name))
graph.new <- add_edges(graph.new, c(new_name, node2))
}
cds@principal_graph[["UMAP"]] <- graph.new
return(cds)
}

#' @export
compress_UMAP <- function(cds, lineage, start, N = 2500){
  cds_name = deparse(substitute(cds))
  input = paste0("sel.cells = ",cds_name,"@lineages$", lineage)
  eval(parse(text=input))
  RNA_umap = reducedDims(cds)$"UMAP"[sel.cells,]
  window = nrow(RNA_umap)/N
  step = ((nrow(RNA_umap)-window)/N)
  #use sliding window to compress expression values and pseudotime
  print(paste0("Window: ", window))
  print(paste0("Step: ", step))
  umap_X = SlidingWindow("mean", RNA_umap[,1], window, step)
  umap_Y = SlidingWindow("mean", RNA_umap[,2], window, step)
  umap = as.data.frame(cbind(umap_X, umap_Y))
  umap
}

#' @export
common <- function(df){
  names(which.max(table(df)))
}

compress_meta <- function(meta, N){
  window = nrow(meta)/N
  step = ((nrow(meta)-window)/N)
  cluster = meta$cluster
  age = meta$age
  cluster.comp = SlidingWindow(common, cluster, window, step)
  age.comp = SlidingWindow(common, age, window, step)
  meta.comp = cbind(age.comp, cluster.comp)
  meta.comp = as.data.frame(meta.comp)
  colnames(meta.comp) <- c("age", "cluster")
  meta.comp
}

#' @export
isolate_lineage_ATAC <- function(cds, ATAC.integrated, lineage, sel_clusters = NULL, starting_clusters = NULL, subset = FALSE, N = 5, cl = 1){
  cds_name = deparse(substitute(cds))
  input = paste0("sub.graph = ",cds_name,"@graphs$", lineage)
  eval(parse(text=input))
  nodes_UMAP = cds@principal_graph_aux[["UMAP"]]$dp_mst
  if(subset == F){
    nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names(V(sub.graph))]))
  }
  else{
    g = principal_graph(cds)[["UMAP"]]
    dd = degree(g)
    names1 = names(dd[dd > 2 | dd == 1])
    names2 = names(dd[dd == 2])
    names2 = sample(names2, length(names2)/subset, replace = F)
    names = c(names1, names2)
    names = intersect(names(V(sub.graph)), names)
    nodes_UMAP.sub = as.data.frame(t(nodes_UMAP[,names]))
  }
  #select cells along the graph
  mean.dist = path.distance(nodes_UMAP.sub)
  r = mean.dist*N
  cells_UMAP = as.data.frame(ATAC.integrated[['UMAP']][[]])
  colnames(cells_UMAP) = c("UMAP_1", "UMAP_2")
  sel.cells = cell.selector(nodes_UMAP.sub, cells_UMAP, r, cl = cl)
  #only keep cells in the progenitor and lineage-specific clusters
  sel.cells1 = c()
  sel.cells2 = sel.cells
  if(length(starting_clusters) > 0){
    sel.cells1 = names(Idents(ATAC.integrated)[Idents(ATAC.integrated) %in% starting_clusters])
  }
  if(length(sel_clusters) > 0){
    sel.cells2 = names(Idents(ATAC.integrated)[Idents(ATAC.integrated) %in% sel_clusters])
  }
  cells = unique(c(sel.cells1, sel.cells2))
  sel.cells = sel.cells[sel.cells %in% cells]
  return(sel.cells)
}

#' @export
import_ATAC <- function(cds_ref, ATAC_counts, UMAP, ATAC_meta, lineage){
  cds_name = deparse(substitute(cds_ref))
  input = paste0("sub.graph = ",cds_name,"@graphs$", lineage)
  eval(parse(text=input))
  nodes_UMAP = cds_ref@principal_graph_aux[["UMAP"]]$dp_mst
  #create monocle object with ATAC counts or gene activity counts
  gene_annotation = as.data.frame(rownames(ATAC_counts))
  rownames(gene_annotation) = rownames(ATAC_counts)
  colnames(gene_annotation) = 'gene_short_name'
  cds = new_cell_data_set(ATAC_counts, cell_metadata = ATAC_meta, gene_metadata = gene_annotation)
  reducedDims(cds)$"UMAP" <- UMAP
  s.clusters = as.character(ATAC_meta$cluster)
  names(s.clusters) <- rownames(ATAC_meta)
  cds@clusters$"UMAP"$"clusters" <- s.clusters
  cds@clusters$UMAP$partitions <- cds@clusters$UMAP$clusters
  cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions != "1"] <- "1"
  principal_graph(cds)[["UMAP"]] <- sub.graph
  cds@principal_graph_aux[["UMAP"]]$dp_mst <- nodes_UMAP[,names(V(sub.graph))]
  cds = import_monocle(cds)
  ref.umap = reducedDims(cds_ref)$"UMAP"
  atac.umap = reducedDims(cds)$"UMAP"
  n_neighbors = 10
  nn <- get.knnx(ref.umap, atac.umap, k=n_neighbors)
  pt <- cds_ref@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pt_atac = apply(nn$"nn.index", 1, function(x) mean(pt[x]))
  cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]] <- pt_atac
  cds@principal_graph_aux[["UMAP"]]$pseudotime <- pt_atac
  names(cds@principal_graph_aux[["UMAP"]]$pseudotime) <- row.names(colData(cds))
  sel.cells = colnames(ATAC_counts)
  cds_name = "cds"
  input = paste0(cds_name, "@graphs$", lineage, " <- sub.graph")
  eval(parse(text=input))
  input = paste0(cds_name, "@lineages$", lineage, " <- sel.cells")
  eval(parse(text=input))
  cds
}


#' @export
compress <- function(df, window = 3, step){
df.comp = SlidingWindow("mean", df, window, step)
}

#' @export
fit.m3 <- function(exp.sel, pt, max.pt, model = "expression ~ splines::ns(pseudotime, df=3)", N = 500){
require(speedglm)
family = stats::quasipoisson()
exp_data.sel = cbind(pt, exp.sel)
colnames(exp_data.sel) <- c("pseudotime","expression")
exp_data.sel = as.data.frame(exp_data.sel)
exp_data.sel$pseudotime <- as.numeric(as.character(exp_data.sel$pseudotime))
exp_data.sel$expression <- as.numeric(as.character(exp_data.sel$expression))
tryCatch({fit = speedglm(model, data = exp_data.sel, family = family, acc=1e-3, model=FALSE, y=FALSE)
d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
colnames(d) <- c("pseudotime")
fit = stats::predict(fit, newdata=d, type="response")
return(fit)
}, error=function(cond) {return(rep("NA", N))})
}

#' @export
as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
    
  for (i in seq_along(val)){
      tmp[row_pos[i],col_pos[i]] <- val[i]
  }
    
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

#' @export
compress2 <- function(df, window, step){
  df.comp = SlidingWindow("mean", df, window, step)
}

#' @export
pull <- function(df, window, step){
  df.comp = SlidingWindow("sum", df, window, step)
}

#' @export
compress_lineages_v2 <- function(cds, start, window = F, N = 500, cores = F){
  lineages = names(cds@lineages)
  for(lineage in lineages){
    print(lineage)
    cds = compress_lineage_v2(cds, lineage = lineage, start = start, window = window, gene = FALSE, N = N, cores = cores)
    gc()
  }
  return(cds)
}

#' @export
compress_lineage_v2 <- function(cds, lineage, start, window = F, gene = FALSE, N = 500, cores = F, cells = FALSE){
  cds_name = deparse(substitute(cds))
  if(gene == FALSE){
    input = paste0("compress_expression_v2(",cds_name,", lineage = '", lineage, "', start = ", start, ", window = ", window, ", gene = ", gene, ", N = ", N, ", cores = ", cores, ")")
  }
  else{
    input = paste0("compress_expression_v2(",cds_name,", lineage = '", lineage, "', start = ", start, ", window = ", window, ", gene = '", gene, "', N = ", N, ", cores = ", cores, ")")
  }
  exp = eval(parse(text=input))
  input = paste0(cds_name, "@expression$", lineage, " <- exp$expression")
  eval(parse(text=input))
  input = paste0(cds_name, "@expectation$", lineage, " <- exp$expectation")
  eval(parse(text=input))
  input = paste0(cds_name, "@pseudotime$", lineage, " <- exp$pseudotime")
  eval(parse(text=input))
  eval(parse(text=paste0("return(",cds_name, ")")))
}

#' @export
compress_expression_v2 <- function(cds, lineage, start, window = F, gene = FALSE, N = 500, cores = F){
  cds_name = deparse(substitute(cds))
  if(cores != F){
    cl <- makeCluster(cores)
    clusterEvalQ(cl, c(library(evobiR)))
  }
  input = paste0("get_lineage_object(",cds_name,", lineage = '", lineage, "', start = ", start, ")")
  cds_subset = eval(parse(text=input))
  family = stats::quasipoisson()
  model = "expression ~ splines::ns(pseudotime, df=3)"
  names(cds_subset) <- rowData(cds_subset)$gene_short_name
  exp = as_matrix(exprs(cds_subset))
  exp = t(exp) /  pData(cds_subset)[, 'Size_Factor']
  pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pt = pt[order(pt)]
  exp = exp[names(pt),]
  if(window == FALSE){
      window = nrow(exp)/N
  }
  step = ((nrow(exp)-window)/N)
  #use sliding window to compress expression values and pseudotime
  print(paste0("Window: ", window))
  print(paste0("Step: ", step))
  pt.comp = SlidingWindow("mean", pt, window, step)
  max.pt = max(pt.comp)
  if(gene != F){
    exp.comp = compress2(exp[,gene], window = window, step = step)
  }
  else{
    print(paste0("Compressing lineage ", lineage, " and fitting curves"))
    step = ((nrow(exp)-window)/N)
    if(cores != F){
    exp.comp = pbapply(exp, 2, compress2, window = window, step = step, cl = cl)
    }
    else{
    exp.comp = pbapply(exp, 2, compress2, window = window, step = step)
    }
  }
  if(gene != F){
    exp_data.sel = cbind(pt.comp, exp.comp)
    exp_data.sel = as.data.frame(exp_data.sel)
    colnames(exp_data.sel) <- c("pseudotime", "expression")
    exp_data.sel$pseudotime <- as.numeric(as.character(exp_data.sel$pseudotime))
    exp_data.sel$expression <- as.numeric(as.character(exp_data.sel$expression))
    d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
    colnames(d) <- c("pseudotime")
    tryCatch({fit = speedglm(model, data = exp_data.sel, family = family, acc=1e-3, model=FALSE, y=FALSE)
    fit = predict(fit, newdata = d, type='response')
    }, error=function(cond) {fit = as.data.frame(rep(0, N))})
    exp = as.data.frame(cbind(exp.comp, fit, d))
    colnames(exp) <- c("expression", "expectation", "pseudotime")
    exp$expression[exp$expression < 0] <- 0
    exp$expectation[exp$expectation < 0] <- 0
  }
  else{
    d = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
    if(cores != F){
    fit = pbapply(exp.comp, 2, fit.m3, pt = d, max.pt = max(d), N = N, cl = cl)
    }
    else{
    fit = pbapply(exp.comp, 2, fit.m3, pt = d, max.pt = max(d), N = N)
    }
    fit = apply(fit, 2, as.numeric)
    return(list("expression" = exp.comp, "expectation" = fit, "pseudotime" = d))
  }
  exp$expression[exp$expression < 0] <- 0
  exp$expectation[exp$expectation < 0] <- 0
  if(cores != F){
    stopCluster(cl)
  }
  return(exp)
}

#' @export
compress_lineages <- function(cds, start, N = 500, cores = F){
lineages = names(cds@lineages)
for(lineage in lineages){
print(lineage)
cds = compress_lineage(cds, lineage, start, gene = FALSE, N = N, cores = cores)
}
return(cds)
}

#' @export
compress_ATAC_expression <- function(cds, lineage, start, genes = NULL, N = 500, cores = F){
  cds_name = deparse(substitute(cds))
  if(cores != F){
    cl <- makeCluster(cores)
    clusterEvalQ(cl, c(library(evobiR)))
  }
  input = paste0("get_lineage_object(",cds_name,", lineage = '", lineage, "', start = ", start, ")")
  cds_subset = eval(parse(text=input))
  family = stats::quasipoisson()
  model = "expression ~ splines::ns(pseudotime, df=3)"
  names(cds_subset) <- rowData(cds_subset)$gene_short_name
  exp = as_matrix(exprs(cds_subset))
  if(length(genes) > 0){
    exp = exp[genes,]
  }
  pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pt = pt[order(pt)]
  exp = exp[,names(pt)]
  window = ncol(exp)/N
  step = ((ncol(exp)-window)/N)
  #use sliding window to compress expression values and pseudotime
  print(paste0("Window: ", window))
  print(paste0("Step: ", step))
  pt.comp = SlidingWindow("mean", pt, window, step)
  max.pt = max(pt.comp)
  print(paste0("Compressing lineage ", lineage, " and fitting curves"))
  step = ((ncol(exp)-window)/N)
  print("pulling counts into metacells")
  if(cores > 1){
    print("multicore processing")
    exp.comp = pbapply(exp, 1, pull, window = window, step = step, cl = cl)
  }
  else{
    exp.comp = pbapply(exp, 1, pull, window = window, step = step)
  }
  exp.comp = round(t(exp.comp))
  return(list(exp.comp, pt.comp))
}

#' @export
get_pt <- function(lineage){
load(file = paste0(lineage, ".R"))
pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pt = pt[order(pt)]
return(pt)
}

#' @export
common <- function(df){
  names(which.max(table(df)))
}

#' @export
compress_meta <- function(meta, N){
  window = nrow(meta)/N
  step = ((nrow(meta)-window)/N)
  cluster = meta$cluster
  age = meta$age
  cluster.comp = SlidingWindow(common, cluster, window, step)
  age.comp = SlidingWindow(common, age, window, step)
  meta.comp = cbind(age.comp, cluster.comp)
  meta.comp = as.data.frame(meta.comp)
  colnames(meta.comp) <- c("age", "cluster")
  meta.comp
}

#' @export
plot_multiple <- function(cds, gene, lineages, meta = NULL, points = T, age.scale = T, scale.lineage = NULL, age.points = c("3rd trimester", "0-1 years", "2-4 years", "4-10 years"), breaks.labels = c("2nd", "3rd", "birth", "1y", "4y"), point_size = 0.1, line_size = 1, text.size = 14, plot.title.size = 36, legend.key.size = 0.5, legend.text.size = 10, colors = c("red", "blue", "green", "cyan", "magenta", "purple", "orange", "black", "yellow", "tan"), N = 500, legend_position = "none"){
  cds_name = deparse(substitute(cds))
  input = paste0(cds_name,"@expression$", lineages[1])
  N = nrow(eval(parse(text = input)))
  pts = c()
  if(length(scale.lineage) == 0){
  for(lineage in lineages){
    input = paste0(cds_name,"@pseudotime$", lineage)
    pt = eval(parse(text = input))[,1]
    pts = c(pts, pt)
  }
  max.pt = max(pts)
  }
  else{
  input = paste0(cds_name,"@pseudotime$", scale.lineage)
  pt = eval(parse(text = input))[,1]
  max.pt = max(pt)
  }
  print(max.pt)
  if(points == T){
    dd = as.data.frame(seq(from=0, to=max.pt, by = max.pt/(N-1)))
    cols = c("pseudotime")
    fits = c()
    exps = c()
    for(lineage in lineages){
      input = paste0("exp = ",cds_name,"@expression$", lineage)
      eval(parse(text=input))
      if(gene %in% colnames(exp)){
        input = paste0("exp = ",cds_name,"@expression$", lineage,"[,'",gene,"']")
        eval(parse(text=input))
        input = paste0("fit = ",cds_name,"@expectation$", lineage,"[,'",gene,"']")
        eval(parse(text=input))
      }
      else{
        exp = rep(0, N)
        fit = rep(0, N)
      }
      dd = cbind(dd, exp, fit)
      cols = append(cols, paste0("exp_", lineage))
      cols = append(cols, paste0("fit_", lineage))
      fits = c(fits, fit)
      exps = c(exps, exp)
    }
    colnames(dd) <- cols
    ymax = max(fits)
  }
  else{
    fits = c()
    dd = matrix(ncol = 3, nrow = 0,)
    for(lineage in lineages){
      input = paste0("exp = ",cds_name,"@expression$", lineage)
      eval(parse(text=input))
      if(gene %in% colnames(exp)){
        input = paste0("fit = ",cds_name,"@expectation$", lineage,"[,'",gene,"']")
        eval(parse(text=input))
      }
      else{
        fit = rep(0, N)
      }
      fits = c(fits, fit)
      dd = rbind(dd, cbind(seq(from=0, to=max.pt, by = max.pt/(N-1)), fit, rep(lineage, length(fit))))
    }
    ymax = max(fits)
    colnames(dd) <- c("pseudotime", "fit", "lineage")
    dd = as.data.frame(dd)
    dd$pseudotime <- as.numeric(dd$pseudotime)
    dd$fit <- as.numeric(dd$fit)
    dd$lineage <- factor(dd$lineage, levels = lineages)
  }
  q <- ggplot(data = dd)
  if(points == T){
    for(M in 1:length(lineages)){
      loop_input1 = paste0("geom_point(aes_string(x='pseudotime',y = '", paste0('exp_', lineages[M]), "',color='pseudotime'), size=I(", point_size, "))")
      loop_input2 = paste0("scale_color_gradient2(lineages[M],low='grey', ", "high='",colors[M],"')")
      loop_input3 = "new_scale_color()"
      loop_input4 = paste0("geom_line(aes_string(x='pseudotime', y = '", paste0('fit_', lineages[M]), "',size = I(", line_size, ")), color = '", colors[M],"')")
      q <- q + eval(parse(text=loop_input1)) + eval(parse(text=loop_input2)) + eval(parse(text=loop_input3)) + eval(parse(text=loop_input4))
    }
  }
  else{
    q <- q + geom_line(aes(x = pseudotime, y = fit, color = lineage), size = I(line_size)) + scale_color_manual(values = colors)
  }
  q <- q + scale_y_log10()
  if(age.scale == T){
    if(length(scale.lineage) == 1){
    input = paste0(cds_name,"@lineages$", scale.lineage)
    cells = eval(parse(text = input))
    age = meta[cells,c("age_num", "age_range")]
    }
    else{
    age = meta[,c("age_num", "age_range")]
    }
    age = age[order(age$age_num),]
    window = nrow(age)/N
    step = ((nrow(age)-window)/N)
    age.comp = SlidingWindow("mean", age$age_num, window, step)
    d = seq(from=0, to=max.pt, by = max.pt/(N-1))
    d = cbind(as.data.frame(d), age.comp) 
    breaks.list = c(0)
    for(age.point in age.points){
      age.break = quantile(age[age$age_range == age.point,]$age_num, 0.95)
      age.break = d[which.min(abs(d[,2]-age.break)),1]
      breaks.list = append(breaks.list, age.break)
    }
    q <- q + scale_x_continuous(breaks = breaks.list, labels = breaks.labels)
  }
  q <- q + ylim(y = c(0,ymax))
  q <- q + monocle_theme_opts() + ylab("Expression") + xlab("Pseudotime") + ggtitle(gene) + theme(legend.key.size = unit(legend.key.size, 'cm'), plot.title = element_text(size = plot.title.size, face="bold", hjust = 0.5), axis.text=element_text(size=text.size), axis.text.x=element_text(angle = 60, hjust=1), axis.title=element_blank(), legend.text=element_text(size=legend.text.size), legend.title=element_text(size=text.size, face = "bold"), legend.position = legend_position)
  q
}

#' @export
get_max_age_v2 <- function(cds, meta, genes = NULL, lineage, start){
  if(length(genes) == 0){
    genes = as.character(rownames(read.table(paste0(lineage, "_spec.txt"), sep = "\t")))
  }
  cds_name = deparse(substitute(cds))
  input = paste0("fit = ",cds_name,"@expectation$", lineage)
  eval(parse(text=input))
  N = nrow(fit)
  age = meta[,c("age_range", "age_num")]
  input = paste0("get_lineage_object(",cds_name,", '", lineage, "',", start, ")")
  cds_subset = eval(parse(text=input))
  pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pt = pt[order(pt)]
  age_sel = age[names(pt), 2]
  window = length(age_sel)/N
  step = ((length(age_sel)-window)/N)
  age.comp = SlidingWindow("mean", age_sel, window, step)
  res = pbsapply(genes, phase_sub_v2, fit = fit, age = age, age.comp = age.comp)
  res
}

#' @export
phase_sub_v2 <- function(gene, fit, age, age.comp, factor = 0.2, factor2 = 0.2){
  age_max = get_max_age_sub(gene, fit, age, age.comp)
  age_mid = get_max_age_sub(gene, fit, age, age.comp, age_factor = 0.5)
  fit = fit[,gene]
  locmin = rollapply(fit, 3, function(x) which.min(x)==2)
  locmin = which(locmin == TRUE)
  if(length(locmin)>1){
  locmin = locmin[1]
  }
  locmax = rollapply(fit, 3, function(x) which.max(x)==2)
  locmax = which(locmax == TRUE)
  if(length(locmax)>1){
    locmax = locmax[1]
  }
  class = c()
  if(min(which(fit == max(fit))) < min(which(fit == min(fit)))){
    direction = "descrease"
  }
  else{
    direction = "increase"
  }
  #one local min and one local max
  if(length(locmin) == 1 & length(locmax) == 1 & length(locmin) > 0){
  #plateau classification
    if(abs(fit[locmax]-fit[length(fit)])/fit[locmax] < factor){
      mode = "plateau"
      class = append(class, mode)
    }
  #biphasic 1 classification
    if(locmin<locmax & (max(fit[1:locmin]) - fit[locmin])/max(fit) > factor2 & (max(fit) - fit[locmin])/max(fit) > factor2){
      mode = "biphasic"
      class = append(class, mode)
    }
  #biphasic 2 classification
    if(locmin>locmax & (max(fit[locmin:length(fit)]) - fit[locmin])/max(fit) > factor2 & (fit[locmax] - fit[locmin])/max(fit) > factor2){
      mode = "biphasic"
      class = append(class, mode)
    }
  #biphasic 3 classification
    if(locmin>locmax & (max(fit[locmin:length(fit)]) - fit[locmin])/max(fit) > factor2 & (fit[locmax]-min(fit[1:locmax]))/max(fit) > factor2){
      mode = "biphasic"
      class = append(class, mode)
    }
  #burst classification
    if(locmax<locmin & (max(fit[locmin:length(fit)]) - fit[locmin])/max(fit) > factor2){
      mode = "burst"
      class = append(class, mode)
    }
  }
  #transient classification
  if((max(fit) - fit[1]) > (max(fit)*factor2) & (max(fit) - fit[length(fit)]) > (max(fit)*factor2)){
    mode = "transient"
    class = append(class, mode)
  }
  #no localmin
  if(length(locmin) == 0){
    change1 = abs((fit[as.integer(length(fit))/2] - fit[1]))
    change2 = abs((fit[length(fit)] - fit[as.integer(length(fit))/2]))
    if(change1/change2 > factor2){
      mode = "steady"
      class = append(class, mode)
    }
    else{
      mode = "burst"
      class = append(class, mode)
    }
  }
  #at least one localmin
  if(length(locmin) > 0){
    #burst classification
    if((fit[length(fit)-fit[locmin]])/fit[locmin] > factor2){
      mode = "burst"
      class = append(class, mode)
    }
  }
  final_class = "unclassified"
  class = unique(class)
  if("plateau" %in% class){
    final_class = "plateau"
  }
  if("transient" %in% class){
    final_class = "transient"
  }
  if("biphasic" %in% class){
    final_class = "biphasic"
  }
  if(length(class) == 1){
    if(class == "burst" & direction == "descrease"){
      final_class = "drop"
    }
    if(class == "burst" & direction == "increase"){
      final_class = "burst"
    }
    if(class == "steady"){
      final_class = "steady"
    }
    if(class == "plateau"){
      final_class = "plateau"
    }
    if("biphasic" %in% class){
      final_class = "biphasic"
    }
    if("transient" %in% class){
      final_class = "transient"
    }
    }
  res = c(final_class,age_mid,age_max,direction)
  names(res) <- c("pattern", "age_mid", "age_max", "direction")
  res
}
                     
#' @export
format_lineage_out <- function(lineages){
  out = matrix(nrow = 0, ncol = 8,) 
  for(lineage in lineages){
    spec = read.data(paste0(lineage, "_spec.txt"))
    class = read.data(paste0("mod_", lineage, ".txt"))
    moran = read.data(paste0("pt_DGE_", lineage, ".txt"))
    moran = moran[rownames(spec),]
    res = cbind(spec, class, moran[,4])
    names = paste0(rownames(res), "_", lineage)
    res$gene <- rownames(res)
    rownames(res) <- names
    out = rbind(out, res)
  }
  colnames(out) <- c("diff", "lineage", "pattern", "age_mid", "age_max", "direction", "I", "gene")
  out
}

#' @export
read.data=function(x){
  data=read.table(x,sep="\t",header=TRUE,row.names=1,check.names=FALSE);
  return(data);
}

#' @export
phase_sub <- function(gene, fit, age, age.comp, factor = 0.2, factor2 = 0.5, age_factor = 0.9){
  fit = fit[,gene]
  locmin = rollapply(fit, 3, function(x) which.min(x)==2)
  locmin = which(locmin == TRUE)
  locmax = rollapply(fit, 3, function(x) which.max(x)==2)
  locmax = which(locmax == TRUE)
  inc = min(which((fit-min(fit))>((max(fit)-min(fit))*age_factor)))
  age_num = age.comp[inc]
  target_age = min(age[age$age_num > age_num,2])
  age_range = unique(age[age$age_num == target_age,1])
  if(min(which(fit == max(fit))) < min(which(fit == min(fit)))){
    mode = "descrease"
  }
  else{
    mode = "increase"
  }
  #if there are no local minima or maxima - steady increase
  if(length(locmax) == 0 & length(locmin) == 0){
    c("steady",age_range,mode)
  }
  #if there are no local maxima or but there is a minimum - burst increase
  else if(length(locmin) > 0 & length(locmax) == 0){
    c("burst",age_range,mode)
  }
  #if there are local maxima, there are several situations
  else if(length(locmax) == 1){
    age_max = max(age.comp[locmax])
    #if there is a local min, can be burst, plateau or biphasic
    if(length(locmin) > 0){
      #if there is a local min after local max
      if(locmin > locmax){
        #if there is at least 20% increase of exp after local min
        if((max(fit[locmin:length(fit)]) - fit[locmin])/max(fit) > factor){
          #if there are two significant bends, it is a biphasic curve
          if((fit[locmax] - min(fit)) > max(fit)*factor){
            c("biphasic",age_range,mode)
          }
          #otherwise, it's a burst curve
          else if(max(fit) == fit[length(fit)]){
            c("burst",age_range,mode)
          }
          else{
            c("plateau",age_range,mode)
          }
        }
        else{c("plateau",age_range,mode)}
      }
      else if((fit[locmax] - min(fit[locmax:length(fit)]))/fit[locmax] > factor2){
        c("transient",age_range,mode)
      }
      else if(((max(fit[1:locmin]) - fit[locmin])/max(fit)) > factor & ((fit[locmax] - fit[locmin])/fit[locmax]) > factor){
        c("biphasic",age_range,mode)
      }
      else if(max(fit) == fit[length(fit)]){
        c("burst",age_range,mode)
      }
      else{c("plateau",age_range,mode)}
    }
    else{if((fit[locmax] - min(fit[locmax:length(fit)]))/fit[locmax] > factor2){
      c("transient",age_range,mode)
    }
      else{c("plateau",age_range,mode)}
    }
  }
  else if(length(locmax) > 1){
    c("biphasic",age_range,mode)
  }
  else{
    c("other","Adult",mode)
  }
}

#' @export
get_peak_age_branches <- function(cds, genes, lineages, meta, start){
  age_max_mat = matrix(nrow = length(genes), ncol = 0,)
  age_mid_mat = matrix(nrow = length(genes), ncol = 0,)
  mode_mat = matrix(nrow = length(genes), ncol = 0,)
  direction_mat = matrix(nrow = length(genes), ncol = 0,)
  lineage_mat = matrix(nrow = length(genes), ncol = 0,)
  for(lin in lineages){
  print(lin)
  cds_name = deparse(substitute(cds))
  input = paste0("fit = ",cds_name,"@expectation$", lin)
  eval(parse(text=input))
  N = nrow(fit)
  age = meta[,c("age_range", "age_num")]
  input = paste0("get_lineage_object(",cds_name,", '", lin, "',", start, ")")
  cds_subset = eval(parse(text=input))
  pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  pt = pt[order(pt)]
  age_sel = age[names(pt), 2]
  window = length(age_sel)/N
  step = ((length(age_sel)-window)/N)
  age.comp = SlidingWindow("mean", age_sel, window, step)
  res = pbsapply(genes, phase_sub_v2, fit = fit, age = age, age.comp = age.comp)
  res = cbind(t(res), rep(lin, ncol(res)))
  res = as.data.frame(res)
  res$gene <- rownames(res)
  colnames(res) <- c("mode", "age_max", "age_mid", "direction", "lineage", "gene")
  age_max_mat = cbind(age_max_mat, res$age_max)
  age_mid_mat = cbind(age_mid_mat, res$age_mid)
  mode_mat = cbind(mode_mat, res$mode)
  direction_mat = cbind(direction_mat, res$direction)
  }
  ages = c("2nd trimester","0-1 years","1-2 years","2-4 years","4-10 years", "10-20 years", "Adult")
  rownames(age_max_mat) <- genes
  colnames(age_max_mat) <- lineages
  rownames(age_mid_mat) <- genes
  colnames(age_mid_mat) <- lineages
  rownames(mode_mat) <- genes
  colnames(mode_mat) <- lineages
  rownames(direction_mat) <- genes
  colnames(direction_mat) <- lineages
  common_age <- function(x) if(length(x) != length(unique(x))) {names(which.max(table(x)))} else{ages[ages%in%x][1]}
  common <- function(x) if(length(x) != length(unique(x))) {names(which.max(table(x)))} else{paste(x, collapse = "/")}
  age_max = apply(age_max_mat, 1, common_age)
  age_mid = apply(age_mid_mat, 1, common_age)
  pattern = apply(mode_mat, 1, common)
  direction = apply(direction_mat, 1, common)
  out = cbind(pattern, age_max, age_mid, direction)
  rownames(out) <- genes
  out
}
    
#' @export
get_peak_age <- function(cds, genes, lineage, meta, age_factor = 1){
  cds_name = deparse(substitute(cds))
  input = paste0("fit = ",cds_name,"@expectation$", lineage)
  eval(parse(text=input))
  N = nrow(fit)
  age = meta[,c("age_range", "age_num")]
  input = paste0("cells = ",cds_name,"@lineages$", lineage)
  eval(parse(text=input))
  pt <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
  cells = cells[cells %in% names(pt)]
  pt = pt[cells]
  pt = pt[order(pt)]
  age_sel = age[names(pt), 2]
  window = length(age_sel)/N
  step = ((length(age_sel)-window)/N)
  age.comp = SlidingWindow("mean", age_sel, window, step)
  res = pbsapply(genes, get_max_age_sub, age = age, fit = fit, age.comp = age.comp, age_factor = age_factor)
  res
}
                     
#' @export
get_max_age_sub  <- function(gene, fit, age, age.comp, age_factor = 1, shift_factor = 0.1, max_age = "Adult", min_age = "2nd trimester"){
  fit = fit[,gene]
  age = age[order(age$age_num),]
  inc = min(which((fit-min(fit))>=((max(fit)-min(fit))*age_factor)))
  if(inc == length(fit)){
    age_range = max_age
  }
  else if (inc == 1){
    age_range = min_age
  }
  else{
  if(inc+(inc*shift_factor) <= length(fit)){
    mmin = age.comp[inc+(inc*shift_factor)]
  }
  else{
    mmin = age.comp[inc]
  }
  if(inc-(inc*shift_factor) > 0){
    mmax = age.comp[inc-(inc*shift_factor)]
  }
  else{
    mmax = age.comp[inc]
  }
  age_num_min = unique(age$age_num[which(abs(age$age_num - min(mmax, mmin)) == min(abs(age$age_num - min(mmax, mmin))))])
  age_num_max = unique(age$age_num[which(abs(age$age_num - max(mmax, mmin)) == min(abs(age$age_num - max(mmax, mmin))))])
  target_age = age[age$age_num >= age_num_min & age$age_num <= age_num_max,1]
  age_range = names(sort(table(target_age),decreasing=TRUE)[1])
  }
  age_range
}

#' @export
get_max_age <- function(cds, meta, genes = NULL, lineage, start, age_factor = 0.9){
if(length(genes) == 0){
genes = as.character(rownames(read.table(paste0(lineage, "_spec.txt"), sep = "\t")))
}
cds_name = deparse(substitute(cds))
input = paste0("fit = ",cds_name,"@expectation$", lineage)
eval(parse(text=input))
N = nrow(fit)
age = meta[,c("age_range", "age_num")]
input = paste0("get_lineage_object(",cds_name,", '", lineage, "',", start, ")")
cds_subset = eval(parse(text=input))
pt <- cds_subset@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
pt = pt[order(pt)]
step = ((length(pt)-3)/N)
pt.comp = SlidingWindow("mean", pt, 3, step)
age_sel = age[names(pt), 2]
age.comp = SlidingWindow("mean", age_sel, length(age_sel)/N, step)
#res = pbsapply(genes, get_max_age_sub, age = age, fit = fit, pt.comp = pt.comp, age.comp = age.comp)
res = pbsapply(genes, phase_sub, fit = fit, age = age, age.comp = age.comp, age_factor = age_factor)
res
}

#' @export
annotate_gene_peaks_sub <- function(gene, d, cells, age){
exp_age = cbind(d[gene,cells], age[cells,])
colnames(exp_age)[1] <- "exp"
res = exp_age[exp_age$exp >= quantile(exp_age$exp, 0.99),]
res = c(Mode(res$age_range), mean(res$age_num))
names(res) <- c("age_range", "age_num")
res
}

#' @export
annotate_gene_peaks <- function(data, cds, genes, lineage, age){
d = GetAssayData(object = data, assay = "RNA", slot = "data")
cds_name = deparse(substitute(cds))
input = paste0("cells = ",cds_name,"@lineages$", lineage)
eval(parse(text=input))
res = pbsapply(genes, annotate_gene_peaks_sub, d = d, cells = cells, age = age)
res
}

#' @export
AUC_window_sub <- function(gene, cds, lineage, comp_lineages, factor, window_ratio){
  cds_name = deparse(substitute(cds))
  input = paste0("fit = ",cds_name,"@expectation$", lineage)
  eval(parse(text=input))
  fit.sel = fit[,gene]
  if(any(is.na(fit.sel))){
  return(FALSE)
  }
  N = length(fit.sel)
  window = N*window_ratio
  out = c()
  top = 0
  mat = matrix(nrow = N, ncol = 0, )
  for(lin in comp_lineages){
    cds_name = deparse(substitute(cds))
    input = paste0("fit = ",cds_name,"@expectation$", lin)
    eval(parse(text=input))
    if(gene %in% colnames(fit) & !(any(is.na(fit[,gene])))){
    fit = fit[,gene]
    mat = cbind(mat, fit)
    top = max(fit, top)
    }
  }
      res = c()
      score = c()
      for(i in 1:N){
          if(fit.sel[i] > top*factor){
            res = append(res, 1)
          }
          else{
            res = append(res, 0)
          }
          score = append(score, fit.sel[i]-max(mat[i,]))
        }
      target = paste(res,collapse="")
      res = rle(res)$lengths[rle(res)$values==1]
      if(any(res > window))
      {
        search = paste(rep(1, max(rle(res)$values)),collapse="")
        index = gregexpr(pattern = search, target)[[1]][1]
        score = sum(score[index:(index+max(rle(res)$values))-1])/top
        score
      }
      else{
        FALSE
      }
    }

#' @export
AUC_window_sub_complex <- function(gene, cds, lineage, comp_lineages, factor, window_ratio){
  cds_name = deparse(substitute(cds))
  input = paste0("fit = ",cds_name,"@expectation$", lineage)
  eval(parse(text=input))
  fit.sel = fit[,gene]
  N = length(fit.sel)
  window = N*window_ratio
  out = c()
  res = c()
      for(i in 1:N){
        top = 0
        for(lin in comp_lineages){
          cds_name = deparse(substitute(cds))
          input = paste0("fit = ",cds_name,"@expectation$", lin)
          eval(parse(text=input))
          if(gene %in% colnames(fit) & !(is.na(fit[,gene]))){
            fit = fit[,gene]
            top = max(fit[i], top)
          }
        }
          if(fit.sel[i] > top*factor){
            res = append(res, TRUE)
          }
          else{
            res = append(res, FALSE)
          }
        }
      res = rle(res)$lengths[rle(res)$values==TRUE]
      if(any(res > window))
      {
        TRUE
      }
      else{
        FALSE
      }
    }

#' @export
get_dist <- function(genes, lineage, lineages, dist)
out = sapply(genes, get_dist_sub, lineage = lineage, lineages = lineages, dist = dist)

get_dist_sub <- function(gene, lineage, lineages, dist){
gene.name = paste0(lineage, "__", gene)
dist.sel = dist[gene.name,]
dist.sel = dist.sel[names(dist.sel) %in% paste0(lineages, "__", gene)]
names(dist.sel) <- gsub(paste0("__", gene), "", names(dist.sel))
if(length(dist.sel) != length(lineages)){
dist.sel = rep(0, length(lineages))
names(dist.sel) <- lineages
}
dist.sel
}

#' @export
get_Is <-  function(res, lineages, action){
out = matrix(nrow = length(res), ncol = 0,)
for(lin in lineages){
load(file = paste("Moran_", lin,".R", sep = ""))
Is = cds_pr_test_res$morans_I
names(Is) <- rownames(cds_pr_test_res)
I = sapply(res, lookup_I, I = Is, action = 0)
out = cbind(out, I)
}
colnames(out) <- lineages
rownames(out) <- res
out
}

#' @export
get_lineage_genes <- function(lineage, lineages, exlude_lineages = F, factor = 4, high_I = 0.1){
load(file = paste("Moran_", lineage,".R", sep = ""))
lin.genes = cds_pr_test_res[cds_pr_test_res$morans_I >= high_I & cds_pr_test_res$q_value < 0.05, ]
lin.genes = lin.genes[,c("morans_I", "q_value")]
test_lineages = lineages[!(lineages == lineage)]
if(exlude_lineages != F){
test_lineages = test_lineages[!(test_lineages %in% exlude_lineages)]
}
Is = get_Is(rownames(lin.genes), test_lineages, action = 0)
out = c()
for(i in 1:nrow(lin.genes)){
test = lin.genes$morans_I[i]
d = Is[i,]
res = (d* factor) < test
res = sum(res == TRUE)
if(res == length(test_lineages)){
out = append(out, rownames(Is)[i])
}
}
out
}

#' @export
get_lineage_genes_multiple <- function(test_lineages, lineages, factor = 4, high_I = 0.1){
res = get_lineage_genes(test_lineages[1], lineages, exlude_lineages = test_lineages, high_I = high_I, factor = factor)
for(test_lineage in test_lineages[2:length(test_lineages)]){
local_res = get_lineage_genes(test_lineage, lineages, exlude_lineages = test_lineages, high_I = high_I, factor = factor)
res = intersect(res, local_res)
}
out = get_Is(res, lineages, action = 0)
out
}

prep_mat <- function(gene, mat){
name = paste0(lineage, "__", gene)
cols = colnames(mat)[grepl(paste0("__", gene), colnames(mat))]
res = mat[name,cols]
unlist(res)
}

#' @export
get_distance <- function(lineage, genes){
d = read.table("adult_dist_mat.txt", sep="\t",header=TRUE,row.names=1,check.names=FALSE)
res = sapply(genes, prep_mat, mat = d)
res = t(res)
colnames(res) <- gsub(paste0("__", genes[1]), "", colnames(res))
res
}

lookup_I <- function(gene, I, action = "NA"){
if(gene %in% names(I)){
res = I[gene]
}
else{
res = action
}
res
}

#' @export
get_Moran_I <- function(lineages, genes){
out = matrix(nrow = length(genes), ncol = 0)
for(lineage in lineages){
load(file = paste("Moran_", lineage,".R", sep = ""))
I = cds_pr_test_res$morans_I
names(I) <- rownames(cds_pr_test_res)
res = sapply(genes, lookup_I, I = I)
out = cbind(out, res)
}
colnames(out) <- lineages
out
}

#' @export
calc_dist <- function(lineage1, lineage2, pt_genes = FALSE, dist_method = "CORT", q = 0.05, I = 0.15, adjust = F, cores = 1){
if(pt_genes == FALSE){
pt_genes = c()
for(lineage in c(lineage1, lineage2)){
d = read.table(paste0("pt_DGE_", lineage,".txt"), sep = "\t", header = T)
d = d[d$q_value < q & d$morans_I >= I,]
if(length(pt_genes) == 0){
pt_genes = toupper(d$gene_short_name)
}
else{
pt_genes = union(pt_genes, toupper(d$gene_short_name))
}
}
}
input = paste0("fit = cds@expectation$", lineage1)
eval(parse(text=input))
colnames(fit) <- toupper(colnames(fit))
fit1 = fit[,colnames(fit)%in%pt_genes]
input = paste0("fit = cds@expectation$", lineage2)
eval(parse(text=input))
colnames(fit) <- toupper(colnames(fit))
fit2 = fit[,colnames(fit)%in%pt_genes]
pt_genes = intersect(pt_genes, intersect(colnames(fit1), colnames(fit2)))
fit1 = fit1[,pt_genes]
fit2 = fit2[,pt_genes]
fit1 = t(fit1)
fit2 = t(fit2)
fit1 = log10(fit1+1)
fit2 = log10(fit2+1)
fit1.list <- split(fit1, seq(nrow(fit1)))
names(fit1.list) <- rownames(fit1)
fit2.list <- split(fit2, seq(nrow(fit2)))
names(fit2.list) <- rownames(fit2)
if(dist_method == "CORT"){
out = pbmapply(CortDistance, fit1.list, fit2.list, deltamethod="DTW")
}
if(dist_method == "ACF"){
out = pbmapply(ACFDistance, fit1.list, fit2.list)
}
if(dist_method == "DTW"){
out = pbmapply(ACFDistance, fit1.list, fit2.list)
}
return(out)
}

#' @export
calc_dist_lineages <- function(lineages, pt_genes = FALSE, dist_method = "CORT", suffix = "", q = 0.05, I = 0.15){
if(pt_genes == FALSE){
pt_genes = c()
for(lineage in lineages){
d = read.table(paste0("pt_DGE_", lineage,".txt"), sep = "\t", header = T)
d = d[d$q_value < q & d$morans_I >= I,]
if(length(pt_genes) == 0){
pt_genes = d$gene_short_name
}
else{
pt_genes = union(pt_genes, d$gene_short_name)
}
}
}
input = paste0("fit = cds@expectation$", lineages[1])
eval(parse(text=input))
fit.comb = matrix(ncol = 0, nrow = nrow(fit), 0)
for(lineage in lineages){
input = paste0("fit = cds@expectation$", lineage)
eval(parse(text=input))
colnames(fit) <- colnames(fit)
fit = fit[,colnames(fit)%in%pt_genes]
colnames(fit) <- paste(lineage, colnames(fit), sep = "_")
fit.comb = cbind(fit.comb, fit)
}
fit.comb = apply(fit.comb, 2, as.numeric)
fit.comb = as.data.frame(fit.comb)
d = log10(t(fit.comb)+1)
if(dist_method == "CORT"){
dis = TSclust::diss(d, dist_method, deltamethod="DTW")
}
else{
dis = TSclust::diss(d, dist_method)
}
return(list("dis" = dis, "d" = d))
}

#' @export
get_pt_exp <-function(cds, lineages, pt_genes = FALSE, q = 0.05, I = 0.15){
if(pt_genes == FALSE){
pt_genes = c()
for(lineage in lineages){
d = read.table(paste0("pt_DGE_", lineage,".txt"), sep = "\t", header = T)
d = d[d$q_value < q & d$morans_I >= I,]
if(length(pt_genes) == 0){
pt_genes = d$gene_short_name
}
else{
pt_genes = union(pt_genes, d$gene_short_name)
}
}
}
cds_name = deparse(substitute(cds))
input = paste0("fit = ",cds_name,"@expectation$", lineages[1])
eval(parse(text=input))
fit.comb = matrix(ncol = 0, nrow = nrow(fit), 0)
for(lineage in lineages){
input = paste0("fit = ",cds_name,"@expectation$", lineage)
eval(parse(text=input))
colnames(fit) <- colnames(fit)
fit = fit[,colnames(fit)%in%pt_genes]
colnames(fit) <- paste(lineage, colnames(fit), sep = "__")
fit.comb = cbind(fit.comb, fit)
}
fit.comb = apply(fit.comb, 2, as.numeric)
fit.comb = as.data.frame(fit.comb)
d = log10(t(fit.comb)+1)
}

#' @export
clust_lineages <- function(d, dis, k, lineage_names, colors, suffix = "", myCol = c("pink1", "violet", "mediumpurple1", "slateblue1", "purple", "purple3",  "turquoise2", "skyblue", "steelblue", "blue2", "navyblue",  "orange", "tomato", "coral2", "khaki1", "violetred", "red2",  "springgreen2", "yellowgreen", "palegreen4",  "wheat2", "tan", "tan2", "tan3", "brown",  "grey70", "grey50", "grey30", "aquamarine", "bisque3", "cornflowerblue", "darkseagreen1", "darkred", "lightgreen","hotpink"), method = "ward.D2"){
tree = hclust(dis, method = method)
png(filename = paste0("tree_",suffix,".png"), width = 3500, height = 2080, bg = "white",  res = 300);
plot(tree)
rect.hclust(tree, k = k, border = 2:10)
dev.off();
clust = cutree(tree, k = k)
out = matrix(nrow = length(clust)/length(lineage_names),ncol = 0, 0)
for(lineage_name in lineage_names){
clust.lin = clust[grepl(paste0(lineage_name,"_"), names(clust))]
out = cbind(out, as.character(clust.lin))
}
rownames(out) = gsub(paste0(lineage_name, "_"), "", names(clust.lin))
colnames(out) <- lineage_names
look = cbind(unique(clust[tree$order]), myCol[1:k])
new <- out
new[] <- look[,2][match(unlist(out), look[,1])]
write.table(new, paste0("dynamic_clust_",suffix,".txt"), sep = "\t", quote = F)
hcd = as.dendrogram(tree)
cols = rownames(d)
N = 1
for(lineage_name in lineage_names){
cols = gsub(paste0(lineage_name, "_.+"), colors[N], cols)
N = N+1
}
png(filename = paste0("dendo_L_",suffix,".png"), width = 3500, height = 2080, bg = "white",  res = 300);
hcd = as.dendrogram(tree)
labels(hcd) <- c(1:nrow(d))
hcd %>% set("branches_lwd", 5) %>% plot(horiz = F, leaflab = "none")
colored_bars(cols, hcd, horiz = F, y_shift = -40, y_scale = 200)
hcd %>% rect.dendrogram(cluster = clust, k = k, horiz = F, lwd = 5, lower_rect = -280, text = unique(clust[tree$order]), border = myCol[1:k])
dev.off()
png(filename = paste0("dendo_",suffix,".png"), width = 3500, height = 2080, bg = "white",  res = 300);
hcd = as.dendrogram(tree)
labels(hcd) <- c(1:nrow(d))
hcd %>% set("branches_lwd", 5) %>% plot(horiz = F, leaflab = "none")
colored_bars(cols, hcd, horiz = F, y_shift = -40, y_scale = 200)
hcd %>% rect.dendrogram(cluster = clust, k = k, horiz = F, lwd = 5, lower_rect = -280, border = myCol[1:k])
dev.off()
}

get_sex_enriched_genes <-function(cds_new, factor = 1.2, window_ratio = 0.01, I = 0.1){                     
lineages = gsub("_Male", "", names(cds_new@lineages))
lineages = gsub("_Female", "", lineages)
lineages = unique(lineages)
lineages = lineages[lineages != "MG_3"]
for(lin in lineages){
  print(lin)
  #male
  lineage = paste0(lin,"_Male")
  comp.lineage = paste0(lin,"_Female")
  d = get_pt_exp(cds_new, lineage, I = I)
  names = gsub(paste0(lineage, "__"), "", rownames(d))
  cds_name = deparse(substitute(cds_new))
  input = paste0("fit = ",cds_name,"@expectation$", lineage)
  eval(parse(text=input))
  names1 = names[names %in% colnames(fit)]
  d = get_pt_exp(cds_new, comp.lineage, I = 0)
  names = gsub(paste0(comp.lineage, "__"), "", rownames(d))
  input = paste0("fit = ",cds_name,"@expectation$", comp.lineage)
  eval(parse(text=input))
  names2 = names[names %in% colnames(fit)]
  names = union(names1, names2)
  res = sapply(names, AUC_window_sub, cds = cds_new, lineage = lineage, comp_lineages = comp.lineage, factor = factor, window_ratio = window_ratio)
  res = res[res != 0]
  res = cbind(as.data.frame(res), rep(lineage, length(res)))
  write.table(res, paste0(lineage, "_spec.txt"), quote = F, sep = "\t")
  #female
  lineage = paste0(lin,"_Female")
  comp.lineage = paste0(lin,"_Male")
  d = get_pt_exp(cds_new, lineage, I = 0)
  names = gsub(paste0(lineage, "__"), "", rownames(d))
  cds_name = deparse(substitute(cds_new))
  input = paste0("fit = ",cds_name,"@expectation$", lineage)
  eval(parse(text=input))
  names1 = names[names %in% colnames(fit)]
  d = get_pt_exp(cds_new, comp.lineage, I = I)
  names = gsub(paste0(comp.lineage, "__"), "", rownames(d))
  input = paste0("fit = ",cds_name,"@expectation$", comp.lineage)
  eval(parse(text=input))
  names2 = names[names %in% colnames(fit)]
  names = union(names1, names2)
  res = sapply(names, AUC_window_sub, cds = cds_new, lineage = lineage, comp_lineages = comp.lineage, factor = factor, window_ratio = window_ratio)
  res = res[res != 0]
  res = cbind(as.data.frame(res), rep(lineage, length(res)))
  write.table(res, paste0(lineage, "_spec.txt"), quote = F, sep = "\t")
}
  }
