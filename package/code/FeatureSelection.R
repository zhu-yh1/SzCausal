featureSelection = function (data, focalVariable, contr = NULL, limmaCutoff = NULL, glmnetCutoff = NULL, class = NULL, nfolds = 3, bootstrap = T, alpha = 0.9, ...) {
	# limma
	if (class == "binomial") {
		limRes = fitLimma(data, as.factor(focalVariable), contr = contr)
		print(topTable(limRes))
	}
	else
	{
		limRes = fitLimma(data, focalVariable)
	}
	
	if(!is.null(limmaCutoff)) {
		limSelected = rownames(topTable(limRes, n=limmaCutoff))
	}
	else {
		limSelected = rownames(data)
	}
	
	
	# elastic net
	resGN = runGnetMulti(x = t(data), y = focalVariable, nfolds = nfolds, bootstrap = bootstrap, alpha = alpha)

	if (!is.null(glmnetCutoff)) {
		cutoff = sort(resGN$features, decreasing = T)[glmnetCutoff]
		resGNSelected = rownames(data)[which(resGN$features >= cutoff)]
	}
	else {
		resGNSelected = rownames(data)
	}
	selected = unique(c(limSelected, resGNSelected))
	results = list(limRes = limRes, limSelected = limSelected, resGN = resGN, resGNSlected = resGNSelected, selected = selected)
	class(results) = "featureSelectionResults"
	results
}

plotSelection = function(data, res) {
	col = rep("Unselected", nrow(data))
	col[rownames(data) %in% res$selected] = c("Selected")
	# y = topTable(res$limRes, n=length(res$limSelected))
	# y = min(abs(y$t))
	y = topTable(res$limRes, n=length(res$limSelected))
	y = -log(max(y$P.Value))
	x = sort(res$resGN$features, decreasing = T)[length(res$resGNSlected)]
	p=myggscatter(res$resGN$features, -log(res$limRes$p.value), line = F,
				  label = rownames(data),
				  col = col,
				  palette = c("#ff0000", "#bebebe"),
				  xlab="Chosen frequency in elastic-net regressions",
				  ylab="-log(p) from Limma",
				  label.select = res$selected,
				  font.label = c(10, "black"),
				  title = paste0("models > ", x, " or -log(p-value) > ", round(y, digits=2))) + 
		geom_hline(yintercept = round(y, digits=2), linetype="dashed", color = "red", size=1) +
		geom_vline(xintercept = x, linetype="dashed", color = "red", size=1)

	p
}