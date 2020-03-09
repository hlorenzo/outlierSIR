#' Check for outliers in SIR models.
#'
#' This function takes a matrix \eqn{X} defining the same \eqn{n} individuals and a vector \eqn{y} defining also those individuals.
#'
#' @param X A matrix nxp describinf n observations through p varaibles.
#' @param y A vector of \strong{n} observations corresponding to the n rows of \emph{X}.
#' @param H A positive real value. The number of slices in the SIR model. Default is 5.
#' @param method Which method to be performed. One of \emph{BOOT} for bootstrap, \emph{TTR} for train-test replication and \emph{MONO} for mono-model testing.
#' @param kernel The kernel to be used, one of "normal" or "box". Default is "normal". Can be abbreviated. See \code{\link{ksmooth}}.
#' @param bandwidth A real indicating the width of the kernel. If NULL, automatic kernel width tuning is used. Default is "NULL".
#' @param nb.replications An interger giving the number of repliccations for \emph{BOOT} and \emph{TTR} methods. Default is "1000".
#' @param pourcent A real between 0 and 1. The proportion of observations to keep in the test dataset. Default is "0.1".
#' @param plotBandwidthEst An integer. Wether or not to plot the EDR direction with the kernel estimation. Default is "FALSE".
#' @param NCORES An integer. Fiwing the number of cores to be used in the sampling computations. Only BOOT method has been parallelized due to its computation time.
#'
#' @return A list containing indices of outliers in the case of the 3 methods and also the extreme observations in the case of \emph{BOOT}.
#'
#' @export
#' @importFrom stats ksmooth median na.omit
#' @importFrom grDevices boxplot.stats
#' @importFrom graphics abline lines mtext par plot points text title legend
#' @importFrom changepoint cpt.var
#' @importFrom dr dr dr.basis
#' @import foreach
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#'
#' @references{
#'   \insertAllCited{}
#' }
#'
#' @examples
#' data(ozone)
#' # MONO method
#' output <- outlierSIR(ozone[,-c(1,12,13)],ozone[,1],method="MONO")
#' # TTR method
#' # output <- outlierSIR(ozone[,-c(1,12,13)],ozone[,1],method="TTR",
#' # plotBandwidthEst = FALSE,kernel="normal",nb.replications = 1000)
#' # BOOT method, using 5 cores of computation
#' # output <- outlierSIR(ozone[,-c(1,12,13)],ozone[,1],method="BOOT",
#' # plotBandwidthEst = TRUE,kernel="normal",nb.replications = 5000,NCORES=5)
outlierSIR <- function(X,y,method=c("BOOT","TTR","MONO"),
                    kernel="normal",bandwidth=NULL,H=5,
                    nb.replications=1000,
                    pourcent=0.1,plotBandwidthEst=F,NCORES=1){

  ## Reset personnal plot par() settings
  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))
  ## ----------------------------------
  ##########
  ## TEST ##
  ##########
  ## Unitarian test
  test_vec <- function(vec){
    a <- boxplot.stats(vec)$out
    tete <- which(a<median(vec))
    if(length(tete)>0){
      a <- a[-tete]
    }
    which(vec %in% a)
  }
  # BOOT Version
  detectOutliers <- function(BOOTS){
    #-------------------
    n <- nrow(BOOTS[[1]])
    er <- BOOTS[[3]]
    nb <- BOOTS[[1]]
    mean_moins_sd_1 <- rep(NA,n)
    for(i in 1:n){
      df <- data.frame(x=nb[i,],y=er[i,])
      mean_moins_sd_1[i] <- mean(df$y[-which(df$x %in% 0)],na.rm = T)
    }
    ## Find outliers
    outli <- test_vec(log(mean_moins_sd_1))
    ## Find outliers and extremes
    extre <- test_vec(mean_moins_sd_1)
    if(length(outli)>0 & length(extre)>0){
      extre <- extre[-which(extre %in% outli)]
    }
    outli <- list(outliers = outli,extreme = extre)
  }

  bd_edr_est <- function(DONNEES,kernel="normal",H=5,plotBandwidthEst=F){
    CV.prog <- function(x,y,hmin=(max(x)-min(x))/20,hmax=(max(x)-min(x))/2,nbh=25,
                        kernel="normal",plotBandwidthEst=F){
      n<-length(x)
      vecth<-seq(from=hmin,to=hmax,length=nbh)
      matCV<-cbind(vecth,rep(0,nbh))
      for (h in 1:nbh){
        ypred<-rep(0,n)
        for (i in 1:n){
          ypred[i]<-ksmooth(x[-i], y[-i], kernel=kernel, bandwidth=vecth[h],x.points=x[i])$y
        }
        matCV[h,2]<-sum((y-as.matrix(ypred,ncol=1))^2)
      }
      matCV <- na.omit(matCV)
      hopt<-matCV[which(matCV[,2]==min(matCV[,2])),1]
      if (plotBandwidthEst){
        par(mfrow=c(1,3))
        plot(matCV,type="l",xlab="h",ylab="CV(h)")
        abline(v=hopt,col=2,lwd=2)
        title(paste("hopt=",round(hopt,digits=3)," with kernel '",kernel,"'",sep=""))
        plot(x,y,main=paste("Regression line with kernel '",kernel,"'",sep=""))
        lines(ksmooth(x, y, bandwidth=hopt,n.points=50,kernel=kernel))
      }
      list(matCV=matCV,hopt=hopt)
    }
    ## ----------
    res.est <- dr(y~., method="sir", nslices = H, data=DONNEES)
    B.est <- dr.basis(res.est)[,1,drop=FALSE]
    indice.est <- as.vector(as.matrix(DONNEES[,-1])%*%B.est)
    list(bwd=CV.prog(indice.est,
                     DONNEES$y,
                     plotBandwidthEst=plotBandwidthEst,kernel = kernel)$hopt,B.est=B.est)
  }

  ###########
  ### TTR ###
  ###########
  method_TTR <- function(DONNEES,nb.replications=1000,
                         pourcent=0.35,H=5,bandwidth=NULL,
                         plotBandwidthEst=F){

    # Bandwith estimation
    res <- bd_edr_est(DONNEES,H=H,plotBandwidthEst=plotBandwidthEst,kernel=kernel)
    B.est <- res$B.est
    if(is.null(bandwidth)){
      bandwidth <- res$bwd
    }

    matrice.distances <- matrix(NA,nrow=nrow(DONNEES),ncol=nb.replications)
    for (r in 1 :nb.replications){
      indice.train <- sample(1:nrow(DONNEES),replace=FALSE,size=round(nrow(DONNEES)*(1-pourcent),digits = 0))
      donnees.train <- DONNEES[indice.train,]
      donnees.test <- as.matrix(DONNEES[-indice.train,])
      res.r <- dr(y~., method="sir", nslices = H, data=donnees.train)
      B.r <- dr.basis(res.r)[,1,drop=FALSE]
      donnees.train <-as.matrix(donnees.train)
      indice.test <- c(1:nrow(DONNEES))[-indice.train]
      for (i in 1:length(indice.test)){
        matrice.distances[indice.test[i],r] <- abs(
          ksmooth(x = as.vector(donnees.train[,-1]%*%B.r),y=donnees.train[,1],
                  bandwidth = bandwidth,kernel="normal",x.points = sum(donnees.test[i,-1]*B.r))$y
          - donnees.test[i,1])
      }
    }
    vecteur.distances <- apply(matrice.distances,1,function(l){mean(na.omit(l))})
    # Automatic determination of the number of outliers
    #--------------------------------------------------
    indice.outlier.NaN <- which(is.na(vecteur.distances))
    if (length(indice.outlier.NaN)>0){
      err_mean <- apply(matrice.distances[-indice.outlier.NaN,],1,function(rr){mean(na.omit(rr))})
    }
    if (length(indice.outlier.NaN)==0){
      err_mean <- apply(matrice.distances,1,function(rr){mean(na.omit(rr))})
    }
    ord_mean <- rev(order(err_mean))
    seg1=cpt.var(err_mean[ord_mean],method="BinSeg",Q=1)
    nbre.choisi.outliers <- seg1@cpts[1]
    indice.outlier <- rev(order(vecteur.distances))[1:(nbre.choisi.outliers+length(indice.outlier.NaN))]   #
    nb_outliers_ttr <- length(indice.outlier)
    indice.outlier
  }

  ############
  ### BOOT ###
  ############
  method_BOOT <- function(DONNEES,#zerosAndOnes=F,
                          nb.replications=1000,
                          H=5,bandwidth=NULL,
                          kernel="normal",
                          plotBandwidthEst=F){
    matrice.distances=matrice.realizations=
      matrice.nb_input <- matrix(0,nrow=nrow(DONNEES),ncol=nb.replications)

    # Bandwith estimation
    res <- bd_edr_est(DONNEES,H=H,plotBandwidthEst=plotBandwidthEst,kernel=kernel)
    B.est <- res$B.est
    if(is.null(bandwidth)){
      bandwidth <- res$bwd
    }


    `%my_do%` <- ifelse(NCORES!=1,
                        {
                          out<-`%dopar%`
                          cl <- makeCluster(NCORES)
                          registerDoParallel(cl)
                          out},
                        {
                          out <- `%do%`
                          out})
    MATRICES <- foreach(r=1:nb.replications,
                        .combine = rbind,.packages = "dr") %my_do% {
                          # for (r in 1:nb.replications){
                          # if(!zerosAndOnes){
                          indice.train <- sample(1:nrow(DONNEES),replace=TRUE,size=nrow(DONNEES))
                          # }else{
                          # indice.train_01 <- sample(c(0,1),replace=TRUE,size=nrow(DONNEES))
                          # indice.train <- (1:nrow(DONNEES))[which(indice.train_01==1)]
                          # }
                          matrice.distances=matrice.realizations=matrice.nb_input <- rep(0,nrow(DONNEES))
                          id_train <- as.numeric(names(table(indice.train)))
                          matrice.nb_input[id_train] <- table(indice.train)
                          donnees.train <- DONNEES[indice.train,]
                          donnees.test <- as.matrix(DONNEES)
                          res.r <- dr(y~., method="sir", nslices = H, data=donnees.train)
                          B.r <- dr.basis(res.r)[,1,drop=FALSE]
                          donnees.train <-as.matrix(donnees.train)
                          indice.test <- 1:nrow(DONNEES)
                          x.points <- as.vector(donnees.test[indice.test,-1,drop=F]%*%B.r)
                          realisations <- ksmooth(x = as.vector(donnees.train[,-1]%*%B.r),
                                                  y=donnees.train[,1],
                                                  bandwidth = bandwidth,
                                                  x.points = x.points,
                                                  kernel=kernel)
                          index_reals <- unlist(lapply(x.points,function(xo){which(realisations$x==xo)}))
                          matrice.realizations <- realisations$y[index_reals]
                          matrice.distances <-
                            abs(
                              matrice.realizations
                              - donnees.test[,1]
                            )
                          c(matrice.nb_input,matrice.realizations,matrice.distances)
                        }
    if(NCORES!=1){
      stopCluster(cl)
    }
    pp <- ncol(MATRICES)
    matrice.nb_input <- t(MATRICES[,1:(pp/3),drop=F])
    matrice.realizations <- t(MATRICES[,pp/3+1:(pp/3),drop=F])
    matrice.distances <- t(MATRICES[,2*pp/3+1:(pp/3),drop=F])


      list(
        matrice.nb_input,
        matrice.realizations,
        matrice.distances,
        B.est)
  }

  ############
  ### MONO ###
  ############
  method_MONO <- function(DONNEES,
                          H=5,bandwidth=NULL,
                          kernel="normal",plotBandwidthEst=F){
    nb.replications <- nrow(DONNEES)
    matrice.distances=
      matrice.realizations=
      matrice.nb_input <- matrix(0,nrow=nrow(DONNEES),ncol=1)

    # Bandwith estimation
    res <- bd_edr_est(DONNEES,H=H,plotBandwidthEst=plotBandwidthEst,kernel=kernel)
    B.est <- res$B.est
    if(is.null(bandwidth)){
      bandwidth <- res$bwd
    }

    r <- 1
    indice.train <- 1:nrow(DONNEES)
    matrice.nb_input[as.numeric(names(table(indice.train))),r] <-
      table(indice.train)
    donnees.train <- DONNEES[indice.train,]
    donnees.test <- as.matrix(DONNEES)
    res.r <- dr(y~., method="sir", nslices = H, data=donnees.train)
    B.r <- dr.basis(res.r)[,1,drop=FALSE]
    donnees.train <-as.matrix(donnees.train)
    indice.test <- 1:nrow(DONNEES)
    for(j in indice.test){
      matrice.realizations[j,r] <-
        ksmooth(x = as.vector(donnees.train[,-1]%*%B.r),
                y=donnees.train[,1],
                bandwidth = bandwidth,
                x.points = as.vector(donnees.test[j,-1,drop=F]%*%B.r),
                kernel=kernel
        )$y
    }
    matrice.distances[,r] <-
      abs(
        matrice.realizations[,r]
        - donnees.test[,1]
      )


    list(
      matrice.nb_input,
      matrice.realizations,
      matrice.distances,
      B.est)
  }

  ############
  ### PLOT ###
  ############
  plot_solution_no_TRUE <- function(DONNEES,id_in,H=5){
    if(is.list(id_in)){
      indice.outlier <- id_in$outliers
      extr <- id_in$extreme
    }else{
      indice.outlier <- id_in
    }
    # Estimation de la direction EDR sur les donnees completes
    #---------------------------------------------------------
    ids <- if(length(indice.outlier)>0)-indice.outlier else 1:nrow(DONNEES)
    res.est <- dr(y~., method="sir", nslices = H, data=DONNEES[ids,])
    B.est <- dr.basis(res.est)[,1,drop=FALSE]
    plot(as.matrix(DONNEES[,-1])%*%B.est,DONNEES$y, pch=16,ylab="y",
         xlab="",main=method)
    mtext(text = "Estimated index",side = 1,line = 2)
    if(length(indice.outlier)>0){
      yy <- DONNEES$y[indice.outlier]
      id_left <- which(yy<median(yy))
      id_right <- which(yy>=median(yy))
      points(as.matrix(DONNEES[indice.outlier,-1])%*%B.est, DONNEES$y[indice.outlier],
             col="blue",pch=1,cex=2)
      text(as.matrix(DONNEES[indice.outlier[id_left],-1])%*%B.est,
           DONNEES$y[indice.outlier][id_left],
           as.character(indice.outlier)[id_left],cex=1.5,col="blue",pos = 3)
      text(as.matrix(DONNEES[indice.outlier[id_right],-1])%*%B.est,
           DONNEES$y[indice.outlier][id_right],
           as.character(indice.outlier)[id_right],cex=1.5,col="blue",pos = 1)
    }
    if(is.list(id_in)){
      legend("topleft",legend = c("Normal observation","Extreme observation", "Outlier"),
             fill=c("black","darkorange","blue"),bty="n")
      if(length(extr)>0){
        yy <- DONNEES$y[extr]
        id_left <- which(yy<median(yy))
        id_right <- which(yy>=median(yy))
        points(as.matrix(DONNEES[extr,-1])%*%B.est, DONNEES$y[extr],
               col="darkorange",pch=1,cex=2)
        text(as.matrix(DONNEES[extr[id_left],-1])%*%B.est,
             DONNEES$y[extr][id_left],
             as.character(extr)[id_left],cex=1.5,col="darkorange",pos = 3)
        text(as.matrix(DONNEES[extr[id_right],-1])%*%B.est,
             DONNEES$y[extr][id_right],
             as.character(extr)[id_right],cex=1.5,col="darkorange",pos = 1)
      }
    }else{
      legend("topleft",legend = c("Normal observation", "Outlier"),
             fill= c("black","blue"),bty="n")
      }
  }

  ## ----------- ##
  #################
  ### MAIN PART ###
  #################
  ## ----------- ##

  DONNEES <- data.frame(list(y=y,X))

  time <- system.time({
    if(method=="MONO"){
      MONOS <- method_MONO(DONNEES,H = H, kernel=kernel,bandwidth=bandwidth,
                           plotBandwidthEst=plotBandwidthEst)
      output <- test_vec(MONOS[[3]])
    }else if(method=="TTR"){
      output <- method_TTR(DONNEES,nb.replications=nb.replications,
                           pourcent=pourcent,H=H,bandwidth=bandwidth,
                           plotBandwidthEst=plotBandwidthEst)
    }else{
      BOOTS <- method_BOOT(DONNEES,
                           nb.replications = nb.replications,
                           H = H, kernel=kernel,bandwidth=bandwidth,
                           plotBandwidthEst=plotBandwidthEst)
      output <- detectOutliers(BOOTS)
    }
  })[3]
  cat(paste("Time elapsed",round(time,3),"\n"))
  plot_solution_no_TRUE(DONNEES,output)
  output
}
