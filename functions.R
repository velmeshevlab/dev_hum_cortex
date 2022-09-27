connect_nodes <- function(cds, node1, node2){
  graph.old = cds@principal_graph[["UMAP"]]
  graph.new <- add_edges(graph.old, c(node1, node2))
  cds@principal_graph[["UMAP"]] <- graph.new
  return(cds)
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

isolate_graph <- function(cds, start, end, lineage, include_nodes = NULL){
  #get lineage graph
  cds_name = deparse(substitute(cds))
  sub.graph = isolate_graph_sub(cds, start, end, lineage, include_nodes = include_nodes)
  input = paste0(cds_name, "@graphs$", lineage, " <- make_graph(sub.graph)")
  eval(parse(text=input))
  eval(parse(text=paste0("return(", cds_name, ")")))
}

isolate_lineage_sub <- function(cds, lineage, sel_clusters = F, start_regions = F, starting_clusters = NULL, subset = FALSE, N = 5, cl = 1){
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
  if(starting_clusters != F){
    sel.cells1 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% starting_clusters])
  }
  if(start_regions != F){
    sel.cells1 = sel.cells1[sel.cells1 %in% rownames(cds@colData[cds@colData$region %in% start_regions,])]
  }
  if(length(sel_clusters) > 0){
    sel.cells2 = names(cds@"clusters"[["UMAP"]]$clusters[cds@"clusters"[["UMAP"]]$clusters %in% sel_clusters])
  }
  cells = unique(c(sel.cells1, sel.cells2))
  sel.cells = sel.cells[sel.cells %in% cells]
  return(sel.cells)
}

isolate_lineage <- function(cds, lineage, sel_clusters = NULL, start_regions = F, starting_clusters = F, subset = FALSE, N = 5, cl = 1){
  cds_name = deparse(substitute(cds))
  sel.cells = isolate_lineage_sub(cds, lineage, sel_clusters = sel_clusters, start_regions = start_regions, starting_clusters = starting_clusters, subset = subset, N = N, cl = cl)
  input = paste0(cds_name, "@lineages$", lineage, " <- sel.cells")
  eval(parse(text=input))
  return(cds)
}

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

