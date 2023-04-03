library(rsvd)

rotateSVD = function(svdres) {
	upos = svdres$u
	uneg = svdres$u
	upos[upos < 0] = 0
	uneg[uneg >= 0] = 0
	uneg = -uneg
	sumposu = colSums(upos)
	sumnegu = colSums(uneg)
	
	for(i in 1:ncol(svdres$u)){
		if(sumnegu[i] > sumposu[i]){
			svdres$u[,i] = -svdres$u[,i] 
			svdres$v[,i] = -svdres$v[,i] 
		}
	}
	svdres
}

simpleDecomp = function(Y, k, svdres = NULL, L1 = NULL, L2 = NULL, Zpos = T, max.iter = 200,
						tol = 5e-6, trace = F, rseed = NULL, B = NULL, scale = 1,
						pos.adj=3, adaptive.frac=0.05, adaptive.iter=30,  cutoff=0){
	ng = nrow(Y)
	ns = ncol(Y)
	
	Bdiff = Inf
	BdiffTrace = double()
	BdiffCount = 0
	message("****")
	
	
	if(is.null(svdres)) {
		message("Computing SVD")
		set.seed(123)
		svdres = rsvd(Y, k = k) 
		
		svdres=rotateSVD(svdres)
		# show(svdres$d[k])
	}
	
	
	if(is.null(L1)) {
		L1 = svdres$d[k] * scale
		if(!is.null(pos.adj)) {
			L1 = L1 / pos.adj
		}
	}
	
	
	if(is.null(L2)) {
		L2 = svdres$d[k] * scale
	}
	# L1=svdres$d[k]/2*scale
	print(paste0("L1 is set to ", L1))
	print(paste0("L2 is set to ", L2))
	
	
	if(is.null(B)) {
		#initialize B with svd
		message("Init")
		B = t(svdres$v[1:ncol(Y), 1:k] %*% diag(sqrt(svdres$d[1:k])))
		# B=t(svdres$v[1:ncol(Y), 1:k]%*%diag((svdres$d[1:k])))
		# B=t(svdres$v[1:ncol(Y), 1:k])
	}
	else {
		message("B given")
	}
	
	
	if (!is.null(rseed)) {
		message("using random start")
		set.seed(rseed)
		B = t(apply(B, 1, sample))
	}
	
	
	round2 = function(x) { signif(x,4) }
	
	getT = function(x) { -quantile(x[x<0], adaptive.frac) }
	
	
	for (i in 1:max.iter) {
		#main loop    
		Zraw = Z = (Y %*% t(B)) %*% solve(tcrossprod(B) + L1 * diag(k))
		
		if(i >= adaptive.iter && adaptive.frac > 0) {
			cutoffs = apply(Zraw, 2, getT)
			
			for(j in 1:ncol(Z)) {
				Z[Z[,j] < cutoffs[j], j] = 0
			}
		}
		
		else if(Zpos) {
			Z[Z < cutoff]=0
		}
		
		oldB = B
		
		B = solve(t(Z) %*% Z + L2 * diag(k)) %*% (t(Z) %*% Y)
		
		#update error
		Bdiff = sum((B - oldB) ^ 2) / sum(B ^ 2)
		BdiffTrace = c(BdiffTrace, Bdiff)
		err0 = sum((Y - Z %*% B) ^ 2) + sum((Z) ^ 2) * L1 + sum(B ^ 2) * L2
		if(trace) {
			message(paste0("iter", i, "errorY = ",erry <- round2(mean((Y - Z %*% B) ^ 2)), ", Bdiff = ", round2(Bdiff), ", Bkappa=", round2(kappa(B))))
		}
		
		#check for convergence
		if(i > 52 && Bdiff > BdiffTrace[i - 50]) {
			BdiffCount = BdiffCount + 1
		}
		else if(BdiffCount > 1) {
			BdiffCount = BdiffCount - 1
		}
		
		if(Bdiff < tol && i > 40) {
			message(paste0("converged at iteration ", i))
			break
		}
		if(BdiffCount > 5 && i > 40) {
			message(paste0("stopped at  iteration ", i, " Bdiff is not decreasing"))
			break
		}
	}
	rownames(B) = colnames(Z) = paste("LV",1:k)
	return(list(B = B, Z = Z, Zraw = Zraw, L1 = L1, L2 = L2))
}