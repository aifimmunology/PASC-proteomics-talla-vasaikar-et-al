library(pROC)
library(dplyr)
library(ggplot2)
library(doMC)
library(foreach)
library(doRNG)


xvGroups <- function(nsmpl, nfold, rev){
  groups <- rep(nsmpl%/%nfold, nfold)
  mod <- nsmpl%%nfold
  if (mod){
    groups[1:mod] <- groups[1:mod]+1
  }
  if(rev%%2){
    groups <- rev(groups)
  }
  return(groups)
}

xvSplit <- function(nsmpl, nfold){
  groups <- rep(nsmpl%/%nfold, nfold)
  mod <- nsmpl%%nfold
  if (mod){
    groups[1:mod] <- groups[1:mod]+1
  }
  return(groups)
}

xvGrouping <- function(data, type, nfold){
  tgrps <- tapply(1:nrow(data), data[,type], function(x)xvSplit(length(x), nfold)) %>% do.call(rbind, .)

  grps <- rbind(rev(tgrps[1,]), tgrps[2,])
  ngrps <- nrow(tgrps)     
  if(ngrps > 2){
    for (i in 3:ngrps){
      grps <- rbind(grps[,order(apply(grps, 2, sum))], tgrps[i,])
    }      
  }
  row.names(grps) <- row.names(tgrps)
  
  return (grps)
}




comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

xvBinaryMC <- function(data, type, biomarkers, nboots, nfold){
  ntype <- 2
  
  # XV grouping & samples
  scounts <- c(0, 0)
  grouping <- c()
  samples <- c()
  for (i in 1:ntype){
    set <- which(data[, type] == i-1)
    scounts[i] <- length(set)
    grouping <- c(grouping, xvGroups(scounts[i], nfold, i)) 
    samples <- c(samples, set) 
  }
  grouping <- matrix(grouping, byrow=T, nrow=ntype)
  

  
  # Formula
  formula <- as.formula(paste0(type, "~", paste0(biomarkers, collapse="+")))
  
  # Bootstraping & Cross Validation
  xvResult <- foreach (i=1:nboots, .combine='comb', .inorder=F, .multicombine=TRUE, .options.RNG=123) %dorng% {

    # make a copy
    bdata <- data

    # bootstraping
    rsamples <- samples
    indx <- 0
    for (j in 1:ntype){
      rsamples[indx+1:scounts[j]] <- sample(rsamples[indx+1:scounts[j]])
      indx <- indx + scounts[j]
    }
    
    # cross validation
    tmpdata <- data.frame(Type=bdata[,type], Score=1)
    tcoeffs <- c()
    for (k in 1:nfold){
      test <- c()
      indx <- 0
      for (j in 1:ntype){
        sel <- 1:grouping[j,k]
        if (k > 1){
          sel <- sel + sum(grouping[j,1:(k-1)])
        }
        test <- c(test, rsamples[indx+sel])
        indx <- indx + scounts[j]
      }
      lgMod <- glm(formula, family=binomial(logit), data=bdata[-test,])
      tcoeffs <- rbind(tcoeffs, lgMod$coeff)
      
      tmpdata[test,"Score"] <- as.numeric(predict(lgMod, newdata=bdata[test,], type="response"))
    }
    

    # AUC
    tauc <- auc(tmpdata$Type, tmpdata$Score, quiet=T)
    if(!is.na(tauc) & all(!is.na(tcoeffs))){
      return(list(tcoeffs, data.frame(AUC=tauc), tmpdata))
    } else {
      return(list(NULL, NULL, NULL))
    }
  }
  
  # output
  coeffs <- xvResult[[1]]
  aucs <- as.double(unlist(xvResult[[2]]))
  xvData <- xvResult[[3]]
  
  coeffTable <- cbind(apply(coeffs, 2, summary, digits=7), summary(aucs))
  colnames(coeffTable)[ncol(coeffTable)] <- "AUC"
  
  coeffTable <- rbind(coeffTable, 
                      round(c(apply(coeffs, 2, sd), sd(aucs)),6))
  rownames(coeffTable)[nrow(coeffTable)] <- "SD"
  
  coeffTable <- rbind(coeffTable, 
                      round(c(abs(apply(coeffs, 2, sd)/apply(coeffs, 2, mean)), 
                              sd(aucs)/mean(aucs)), 6))
  rownames(coeffTable)[nrow(coeffTable)] <- "CV"
  
  return(coeffTable)
}



print.XV <- function(xvResults){
  z <- 1/max(xvResults["CV", 2:(ncol(xvResults)-1)])
  pval <- formatC(2*(1-pnorm(z)), format = "g", digits = 4)
  print(paste("Biomarkers:", ncol(xvResults)-2, 
              "   Median AUC:", round(xvResults["Median", "AUC"], 4), 
              "   Max p:", pval))
  flush.console()
}


xvStable <- function(data, type, biomarkers, nboots, nfold, MC_cores=4){
  xvResults <- list()
    
  xvResults[[1]] <- xvBinaryMC(data, type, biomarkers, nboots, nfold)
  print.XV(xvResults[[1]])
  
  # remove any insignificant, unstable biomarkers
  mxcv <- 1/qnorm(0.975)
  nmks <- length(biomarkers)
  CVs <- xvResults[[1]]["CV",1:nmks+1]
  if (max(CVs) > mxcv){
    imx <- which.max(CVs)		
    xvResults2 <- xvStable(data, type, biomarkers[-imx], nboots, nfold, MC_cores)
    xvResults[1:length(xvResults2)+1] <- xvResults2
  } 
  
  return(xvResults)
}




xvSelection <- function(data, type, biomarkers, nboots, nfold, MC_cores=4){
  if(MC_cores > 1){
    registerDoMC(cores=MC_cores)
  }
  
  # remove any insignificant, unstable biomarkers
  xvResults <- xvStable(data, type, biomarkers, nboots, nfold, MC_cores)
  
  # remaining biomarkers
  n <- length(xvResults)
  nmks <- ncol(xvResults[[n]]) - 2
  submarkers <- colnames(xvResults[[n]])[1:nmks+1]
  
  # remove the biomarker whose removal results in the best AUC & stable coefficients	
  mxcv <- 1/qnorm(0.975)
  
  # lower bound of AUC 
  AUC1SD <- xvResults[[n]]["Median","AUC"] - xvResults[[n]]["SD","AUC"]
  while(nmks > 1){
    tmpResults <- list()
    CVs <- numeric(nmks)
    AUCs <- numeric(nmks)
    for (i in 1:nmks){        
      tmpResults[[i]] <- xvBinaryMC(data, type, submarkers[-i], nboots, nfold)

      CVs[i] <- max(tmpResults[[i]]["CV", 2:nmks])
      AUCs[i] <- tmpResults[[i]]["Median", "AUC"]
    }
    
    # identify the biomarker whose removal results in the best AUC & stable coefficients
    imx <- which.max(AUCs)
    if(AUCs[imx] > 1 | CVs[imx] > mxcv){ # unstable
      set <- which(AUCs <= 1 & CVs <= mxcv) # stable panels
      if(length(set)){ # there are stable panels
        indx <- order(AUCs, decreasing=T)
        for (i in 1:nmks){
          if (AUCs[indx[i]] <= 1 & CVs[indx[i]] <= mxcv){ # first one with stable coefficients
            imx <- indx[i]
            break
          }
        }
      } else { # no stable panel
        # most stable panel
        imx <- which.min(CVs)
        n <- n + 1
        xvResults[[n]] <- tmpResults[[imx]]
        print.XV(xvResults[[n]])
        
        # remove the least stable biomarker
        x <- colnames(xvResults[[n]])
        bmks <- x[2:(length(x)-1)]
        ii <- which.max(xvResults[[n]]["CV", 2:(length(x)-1)])
        bmks <- bmks[-ii]
        
        # get remaining panels
        tmpResults <- xvSelection(data, type, bmks, nboots, nfold, MC_cores)
        xvResults[1:length(tmpResults)+n] <- tmpResults
        
        return (xvResults)
      }
    } 
    
    # store
    n <- n + 1
    xvResults[[n]] <- tmpResults[[imx]]
    nmks <- nmks - 1
    submarkers <- submarkers[-imx]
    
    #print(xvResults[[n]])
    #flush.console()
    print.XV(xvResults[[n]])
    
  }
  
  return(xvResults)
}



write.XV <- function(xvResults, file){
  if(is.list(xvResults)){
    write("[[1]]", file)
    suppressWarnings(write.table(xvResults[[1]], file, sep=",", append=T, col.names=NA))
    
    nbmks <- length(xvResults)
    if (nbmks > 1){
      for (i in 2:nbmks){
        write(paste0(c("\n[[", i, "]]"),collapse=""), file, append=T)
        suppressWarnings(write.table(xvResults[[i]], file, sep=",", append=T, col.names=NA))
      }
    }
  } else {
    write.table(xvResults, file, sep=",", col.names=NA)
  }
}

read.XV <- function(file){
  xvData <- strsplit(readLines(con = file, n = -1L), ",")
  
  set <- grep("[[[:digit:]+]]", xvData)
  
  xvResults <- list()
  for (i in set){
    id <- as.integer(sub("\\[\\[", "", sub("\\]\\]", "", xvData[[i]])))
    
    n <- length(xvData[[i+1]])
    colnames <- as.character(strsplit(gsub("\"", "", xvData[[i+1]]), ",")[-1])
    
    result <- matrix(nrow=8, ncol=n-1)
    rownames <- c()
    for (j in 1:8){
      x <- strsplit(gsub("\"", "", xvData[[j+i+1]]), ",")
      rownames <- c(rownames, as.character(x[1]))
      result[j,] <- as.numeric(x[-1])
    }
    result <- as.data.frame(result)
    
    names(result) <- colnames
    row.names(result) <- rownames
    
    xvResults[[id]] <- result
  }
  
  return(xvResults)
}

summary.XV <- function(xvResults){
  mxcv <- 1/qnorm(0.975)
  nmods <- length(xvResults)
  
  # collect data
  results <- data.frame()
  cnt <- 0 
  for (i in 1:nmods){
    if(xvResults[[i]]["Median", "AUC"] < 0 
       | xvResults[[i]]["Median", "AUC"] > 1){
      next
    }
    
    nbmks <- ncol(xvResults[[i]]) - 2
    cv <- max(xvResults[[i]]["CV", 1:nbmks+1])
    stable <- 1
    if(cv > mxcv){
      stable <- 0
    }
    results <- rbind(results, 
                     data.frame(Model=i, Biomarkers=nbmks, Stable=stable,
                                AUC=xvResults[[i]]["Median", "AUC"],
                                AUC_SE=xvResults[[i]]["SD", "AUC"],
                                Max_CV=cv,
                                Int_CV=xvResults[[i]]["CV", 1])) 
    cnt <- cnt + 1
  }
  
  results <- results[order(-results$Stable, -results$AUC, results$Max_CV),]  
  rownames(results) <- 1:cnt
  # return
  return(results)
}


summary.XVCoeffs <- function(model){
  tmod <- model[,colnames(model) != "AUC"]
  coeffs <- data.frame(Biomarker=colnames(tmod), Coefficient=tmod["Median",], PValue=2*(1-pnorm(1/tmod["CV",])))
  
  return(coeffs)
}


best.XV <- function(xvResults){
  # order models by AUC
  models <- summary.XV(xvResults)
  
  set <- which(models$Stable==1)
  if(length(set)){
    return (models[set[1], "Model"])
  } else {
    print("No stable model found!")
    flush.console()
    return(0)
  }
}



AUC.XV <- function(xvResults){
  nbmks <- length(xvResults)
  
  # collect data
  AUCs <- data.frame()
  for (i in 1:nbmks){
    if(xvResults[[i]]["Median", "AUC"] < 0 
       | xvResults[[i]]["Median", "AUC"] > 1){
      next
    }
    
    AUCs <- rbind(AUCs, data.frame(biomarkers=ncol(xvResults[[i]])-2, 
                                   AUC=xvResults[[i]]["Median", "AUC"],
                                   SE=xvResults[[i]]["SD", "AUC"])) 
  }
  
  # return
  return(AUCs)
  
}

plot.XV <- function(xvResults, pngfile=NULL){
  AUCs <- AUC.XV(xvResults)
  
  # plot
  plt <- ggplot(AUCs, aes(x=biomarkers, y=AUC, colour="red")) + 
           geom_errorbar(aes(ymin=pmax(0,AUC-SE), ymax=pmin(1,AUC+SE)), width=.1, show.legend = F) +
           geom_line(colour="black") +
           geom_point(show.legend = F)

  if(is.null(pngfile)){
    plt
  } else {
    ggsave(pngfile, plot=plt)
  }
}


addModelScore <- function(data, model=NULL, coeffFile=NULL){
  if(is.null(model) & !is.null(coeffFile)){
    model <- read.csv(coeffFile, header=TRUE, stringsAsFactors=F)
  }
  
  bmks <- names(model)[-1]
  
  # Score
  scores <- as.double(model[1]) + data.matrix(data[,bmks]) %*% data.matrix(as.double(model[bmks]))
  scores <- as.double(1/(1+exp(-scores)))
  
  return(scores)
}

