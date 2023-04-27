library(pROC)
library(epiR)

#################################################################################
#
# pltROC_scores <- function(data, type, scores, title, legends, cols, pngfile=NULL)
#
# A simple function to plot multiple ROC curves in a single plot.
#
# data: data frame containing the data 
# type: the variable defining the disease status: 1, disease; 0, non-disease
# scores: a list of predictors whose ROC curves will be plotted
# title: the title for the polt
# legends: a list of legends defining the predictors
# cols: a lost of colors used for plotting the predictors
# pngfile: if provided, the png file used to save the plot
#
pltROC_scores <- function(data, type, scores, title, legends, cols, pngfile=NULL,
                          fnt=1.2, adj=0.5, cex=1, grid=1, sgn=1){
	# ROC & AUC
	rocs <- list()
	cis <- c()
	nsmpl <- c()
	for (i in 1:length(scores)){	
		# dataset
		set <- which(!is.na(data[,scores[i]]) & !is.na(data[,type]))
		nsmpl <- c(nsmpl, length(set))

		# ROC		
		roc <- roc(data[set, type], sgn*data[set,scores[i]], direction="<", quiet = T)
		ci <- round(as.numeric(ci.auc(roc, method="delong")), 3)
		
		rocs[[i]] <- roc
		cis <- rbind(cis, ci)
	}

	# ROC plot
	if(!is.null(pngfile)){
		png(pngfile)
	}

	fnt <- fnt
	omar <- par()$mar
	opty <- par()$pty
	omgp <- par()$mgp
	par(mar=c(3,3,1,1)+0.1, pty="s", mgp=c(2,1,0))
	
	plot(1-rocs[[1]]$specificities, rocs[[1]]$sensitivities, type="l", 
	     xlab="1-Specificity", ylab="Sensitivity", col=cols[1], lwd=2, 
	     cex.lab=fnt, cex.axis=fnt, asp=1, xlim=c(0,1), ylim=c(0,1))
	if (i > 1){
	  for (i in 2:length(scores)){	
	    lines(1-rocs[[i]]$specificities, rocs[[i]]$sensitivities, col=cols[i], lwd=2)
	  }
	}
	if(grid){
  	abline(v=seq(0,1,0.1),h=seq(0,1,0.1),col="lightgrey",lty=2)
	}
	abline(0,1,lty=2)
	title(title,adj=adj)
	legend("bottomright", 
		paste0(legends, " (n=", nsmpl, "): ", 
			sprintf("%.3f",cis[,2]), " (", 
			sprintf("%.3f",cis[,1]), "-", sprintf("%.3f",cis[,3]), ")"), 
			col=cols, lwd=2, cex=cex, title="AUROC (95% CI)", bty='n')
	par(mar=omar, pty=opty, mgp=omgp)

	if(!is.null(pngfile)){
		dev.off()
	}
}


################################################################################
#
# pltROC_sets <- function(data, type, score, datasets, title, legends, cols, pngFile=NULL, cex=1)
# A simple function to plot ROC curves of subsets of data in a single plot.
#
# data: data frame containing the data 
# type: the variable defining the disease status: 1, disease; 0, non-disease
# score: the predictor defining performance in disease prediction
# datasets: a list of subsets of data whose ROC curves will be plotted
# title: the title for the plot
# legends: a list of legends defining the subsets of data
# cols: a lost of colors used for plotting the predictors
# pngfile: if provided, the png file used to save the plot
# cex: font size in legends
#
pltROC_sets <- function(data, type, score, datasets, title="", legends, cols, pngfile=NULL, 
                        fnt=1.2, adj=0.5, cex=1, grid=1, sgn=1){
  # ROC & AUC
  rocs <- list()
  cis <- c()
  nsmpl <- c()
  for (i in 1:length(datasets)){	
    # dataset
    dset <- datasets[[i]]
    set <- which(!is.na(data[dset,type]) & !is.na(data[dset,score]))
    nsmpl <- c(nsmpl, length(set))
    
    # ROC		
    roc <- roc(data[dset, type], sgn*data[dset,score], direction="<", quiet=T)
    ci <- round(as.numeric(ci.auc(roc, method="delong")), 3)

    rocs[[i]] <- roc
    cis <- rbind(cis, ci)
  }
  
  # ROC plot
  if(!is.null(pngfile)){
    png(pngfile)
  }
  
  fnt <- fnt
  omar <- par()$mar
  opty <- par()$pty
  omgp <- par()$mgp
  par(mar=c(3,3,1,1)+0.1, pty="s", mgp=c(2,1,0))
  
  plot(1-rocs[[1]]$specificities, rocs[[1]]$sensitivities, type="l", 
       xlab="1-Specificity", ylab="Sensitivity", col=cols[1], lwd=2, 
       cex.lab=fnt, cex.axis=fnt, asp=1, xlim=c(0,1), ylim=c(0,1))
  if (i > 1){
    for (i in 2:length(datasets)){	
      lines(1-rocs[[i]]$specificities, rocs[[i]]$sensitivities, col=cols[i], lwd=2)
    }
  }
  if(grid){
    abline(v=seq(0,1,0.1),h=seq(0,1,0.1),col="lightgrey",lty=2)
  }
  abline(0,1,lty=2)
  title(title,adj=adj)
  legend("bottomright", 
         paste0(legends, " (n=", nsmpl, "): ", 
                sprintf("%.3f",cis[,2]), " (", 
                sprintf("%.3f",cis[,1]), "-", sprintf("%.3f",cis[,3]), ")"), 
         col=cols, lwd=2, cex=cex, title="AUROC (95% CI)", bty='n')
  par(mar=omar, pty=opty, mgp=omgp)
  
  if(!is.null(pngfile)){
    dev.off()
  }
}


# ###################################################################################################
# #
# # performance
# #
# 
# performanceCI <- function(data, type, score, thresholds){
#   
#   neg <- length(which(data[,type] == 0))
#   pos <- length(which(data[,type] == 1))
#   
#   perf <- c()	
#   for (thres in thresholds){
#     tn <- length(which(data[,type] == 0 & data[,score] < thres))
#     fn <- length(which(data[,type] == 1 & data[,score] < thres))
#     tp <- pos - fn
#     fp <- neg - tn
#     
#     # contigency table
#     m <- matrix(c(tp, fn, fp, tn), nrow=2)
#     
#     # performance
#     p <- epi.tests(m)
#     
#     # sensitivity
#     val <- 100*as.double(p$rval$se)
#     osen <- sprintf("%.1f (%.1f-%.1f)", val[1], val[2], val[3])
# 
#     # specificity
#     val <- 100*as.double(p$rval$sp)
#     ospe <- sprintf("%.1f (%.1f-%.1f)", val[1], val[2], val[3])
#     
#     # accuracy
#     val <- 100*as.double(p$rval$diag.acc)
#     oacc <- sprintf("%.1f (%.1f-%.1f)", val[1], val[2], val[3])
#     
#     # PPV
#     val <- 100*as.double(p$rval$ppv)
#     oppv <- sprintf("%.1f (%.1f-%.1f)", val[1], val[2], val[3])    
#     
#     # NPV
#     val <- 100*as.double(p$rval$npv)
#     onpv <- sprintf("%.1f (%.1f-%.1f)", val[1], val[2], val[3])
#     
#     # odds ratio
#     val <- as.double(p$rval$diag.or)
#     oor <- sprintf("%.3f (%.3f-%.3f)", val[1], val[2], val[3])
#     
#     # PLR
#     val <- as.double(p$rval$plr)
#     oplr <- sprintf("%.3f (%.3f-%.3f)", val[1], val[2], val[3])
#     
#     # NLR
#     val <- as.double(p$rval$nlr)
#     onlr <- sprintf("%.3f (%.3f-%.3f)", val[1], val[2], val[3])
#     
#     perf <- rbind(perf, data.frame(Index = thres, TP=tp, TN=tn, FP=fp, FN=fn, 
#                                    Sensitivity = osen, Specificity = ospe,
#                                    Accuracy = oacc, PPV = oppv, NPV = onpv, 
#                                    Odds_Ratio = oor, PLR = oplr, NLR = onlr))
#   }
#   
#   
#   return(perf)
# }
# 
# 
# 
# 
