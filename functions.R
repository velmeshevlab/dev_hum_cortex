setClass("cell_data_set_ext", contains = "cell_data_set", slots=c(graphs = "list", lineages="list", expression="list", expectation="list", pseudotime="list")) -> cell_data_set_ext

import_monocle <-function(cds){
  cds <- as(cds,"cell_data_set_ext")
  return(cds)
}

read.data=function(x){
  data=read.table(x,sep="\t",header=TRUE,row.names=1,check.names=FALSE);
  return(data);
}

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

compress_lineages <- function(cds, start, window = F, N = 500, cores = F){
  lineages = names(cds@lineages)
  for(lineage in lineages){
    print(lineage)
    cds = compress_lineage(cds, lineage = lineage, start = start, window = window, gene = FALSE, N = N, cores = cores)
    gc()
  }
  return(cds)
}

#' @export
compress_lineage <- function(cds, lineage, start, window = F, gene = FALSE, N = 500, cores = F, cells = FALSE){
  cds_name = deparse(substitute(cds))
  if(gene == FALSE){
    input = paste0("compress_expression(",cds_name,", lineage = '", lineage, "', start = ", start, ", window = ", window, ", gene = ", gene, ", N = ", N, ", cores = ", cores, ")")
  }
  else{
    input = paste0("compress_expression(",cds_name,", lineage = '", lineage, "', start = ", start, ", window = ", window, ", gene = '", gene, "', N = ", N, ", cores = ", cores, ")")
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
compress_expression <- function(cds, lineage, start, window = F, gene = FALSE, N = 500, cores = F){
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
    exp.comp = compress(exp[,gene], window = window, step = step)
  }
  else{
    print(paste0("Compressing lineage ", lineage, " and fitting curves"))
    exp.comp = pbapply(exp, 2, compress, window = window, step = step)
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
    fit = pbapply(exp.comp, 2, fit.m3, pt = d, max.pt = max(d), N = N)
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

compress <- function(df, window, step){
  df.comp = SlidingWindow("mean", df, window, step)
}

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
phase_sub_v2 <- function(gene, fit, age, age.comp, factor = 0.2, factor2 = 0.4){
  age_max = get_max_age_sub(gene, fit, age, age.comp)
  fit = fit[,gene]
  locmin = rollapply(fit, 3, function(x) which.min(x)==2)
  locmin = which(locmin == TRUE)
  locmax = rollapply(fit, 3, function(x) which.max(x)==2)
  locmax = which(locmax == TRUE)
  if(min(which(fit == max(fit))) < min(which(fit == min(fit)))){
    direction = "descrease"
  }
  else{
    direction = "increase"
  }
  mode = "unclassified"
  if(length(locmax) > 0){
    if(length(locmin) > 0){
      if(abs(fit[locmax]-fit[length(fit)])/fit[locmax] < factor){
        mode = "plateau"
      }
      if(locmax<locmin & (max(fit[locmin:length(fit)]) - fit[locmin])/max(fit) > factor2){
        mode = "burst"
      }
      if(locmin>locmax & (max(fit[locmin:length(fit)]) - fit[locmin])/max(fit) > factor2 & (fit[locmax] - fit[locmin])/max(fit) > factor2){
        mode = "biphasic"
      }
      if(locmin>locmax & (max(fit[locmin:length(fit)]) - fit[locmin])/max(fit) > factor2 & (fit[locmax]-min(fit[1:locmax]))/max(fit) > factor2){
        mode = "biphasic"
      }
      if(locmin<locmax & (max(fit[1:locmin]) - fit[locmin])/max(fit) > factor2 & (max(fit) - fit[locmin])/max(fit) > factor2){
        mode = "biphasic"
      }
      if(length(locmin) > 0 & (max(fit)-fit[locmax])/max(fit) < 0.05 & (max(fit) - fit[1]) > (max(fit)*factor2) & (max(fit) - fit[length(fit)]) > (max(fit)*factor2)){
        mode = "transient"
      }
    }
    else if(length(locmin) == 0 & (max(fit)-fit[locmax])/max(fit) < 0.05 & (max(fit) - fit[1]) > (max(fit)*factor2) & (max(fit) - fit[length(fit)]) > (max(fit)*factor2)){
      mode = "transient"
    }
    else if(abs(fit[locmax]-fit[length(fit)])/fit[locmax] < factor){
      mode = "plateau" 
    }
  }
  else{
    if(length(locmin) > 0){
      if((fit[length(fit)-fit[locmin]])/fit[locmin] > factor2){
        mode = "burst"
      }
      if((max(fit[1:locmin]) - fit[locmin])/max(fit) > factor2 & (fit[length(fit)] - fit[locmin])/max(fit) > factor2){
        mode = "biphasic"
      }
    }
    if(length(locmin) == 0){
      change1 = abs((fit[as.integer(length(fit))/2] - fit[1]))
      change2 = abs((fit[length(fit)] - fit[as.integer(length(fit))/2]))
      if(change1/change2 > factor2)
        mode = "steady"
      else{
        mode = "burst"
      }
    }
  }
  if(direction == "descrease"){
    if(mode == "burst" | mode == "steady")
      mode = "drop"
  }
  c(age_max, mode, direction)
}
