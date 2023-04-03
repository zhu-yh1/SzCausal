graphSampling <- function(data, LV, FV, data.use, name, group = NULL, metadata = NULL, n = 20, lambda1 = 0.1) {
	graph_list = list()
	set.seed(42)
	for (i in 1:n) {
		data.sub = unique(sample(data.use, replace=T))
		nt.data = cbind(t(data[data.sub,]), t(LV), FV)
		if(!is.null(group)) {
			nt.data = cbind(nt.data, group)
		}
		if (!is.null(metadata)) {
			nt.data = cbind(nt.data, metadata)
		}
		colnames(nt.data) = c(data.sub, name)
		ntres = notearsInterceptMultiLoss(nt.data, lambda1 = lambda1)
		graph = ntres$graph
		graph_list[[i]]=graph 
	}
	graph_list
}

graphConstruction <- function (graph_list, name, thres = 0.2, n = 20, n_cutoff = 10, symmetric = F) {
	l = length(name)
	newgraph = matrix(data=c(0.0), nrow = l, ncol = l)
	rownames(newgraph) = name
	colnames(newgraph) = rownames(newgraph)
	if (symmetric) {
		for (i in 1:n) {
			graph_list[[i]] = graph_list[[i]] + t(graph_list[[i]])
		}
	}
	for (i in 1:l) {
		for (j in 1:l) {
			cnt = 0;
			sum = 0.0;
			row_name=name[i]
			col_name=name[j]
			for (k in 1:n) {
				namek = rownames(graph_list[[k]])
				if (row_name %in% namek & col_name %in% namek) {
					if(abs(graph_list[[k]][row_name, col_name]) >= thres) {
						cnt = cnt+1;
						sum = sum + graph_list[[k]][row_name, col_name]
					}
				}
			}
			if (cnt >= n_cutoff) {
				newgraph[i,j] = 1.0*sum/cnt
			}
		}
	}
	newgraph
}