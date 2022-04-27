######### random minipatch
MPCC <-function(d=NULL,
                K = NULL,
                reps=300, ## number of iterations
                pItem=0.25, ## minipatch size for observations
                pFeature=0.1, ## minipatch size for features
                innerLinkage="ward.D", ## internal HC linkage
                distance="manhattan",## internal HC distance
                h=0.95, ## cut internal tree at h quantile
                finalAlgorithm='hclust',
                finalLinkage='ward.D',
                early_stop=TRUE, ## whether perform early stop criteria5
                num_unchange = 5,
                eps = 0.00001,
                verbose=TRUE){
    ###########
    ## The MPCC function perform MPCC (MiniPatch Consensus Clustering) with uniform sampling scheme on both observations and features.
    ## It returns final N by N consensus matrix
    #######################
    ## check validity of data input
    if ( !is( d )[1] %in% c( "data.frame","matrix" ) ) {
        stop("d must be a data.frame or matrix")
    }

    ## check distance input
    acceptable.distance <- c( "euclidean", "maximum", "manhattan", "canberra", "binary","minkowski",
                              "pearson", "spearman" )
    ## check distance function
    if ( inherits( distance, "character" ) ) {
        if ( ! distance %in%  acceptable.distance  &  ( is(try(get(distance),silent=TRUE))!="function") ){
            stop("unsupported distance.")}
    }else{stop("unsupported distance specified.")}


    d <-data.frame(scalematrix(as.matrix(d)))

    #######################################
    #######function to update consensus matrix
    #######################################

    update_MPCC  <- function(){
        if(verbose){
            message(paste("  i =",i))
        }
        #################
        ## create MP
        #################
        sample_x <-sample_MiniPatch(d, pItem, pFeature,pi_feature=1,pi_item=1)
        this_assignment <-cluster_algo(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage)
        ##mCount stores number of times a sample pair was sampled together.
        mCount <<- connectivityMatrix(rep(1,length(sample_x$subcols)),
                                      mCount,
                                      sample_x$subcols)

        ##ML stores number of times a sample pair was CLUSTERED together.
        ml <<- connectivityMatrix(this_assignment,
                                  ml,
                                  sample_x$subcols)
        this_Co <-ml / mCount
        this_Co[mCount==0]=0
        CoAsso<<-this_Co
    }

    ##############################
    ############ initialize
    n=ncol(d)
    CoAsso=mCount=ml=matrix(c(0),ncol=n,nrow=n)

    if (early_stop == FALSE){
        ###############################################################
        ### if not perform early stopping, run given number of iterations
        ################################################################
        for (i in seq_len(reps)){
            ## update consensus matrix
            update_MPCC()
        }
    }else{
        conf_q <-NULL ## record quantile of confusions
        continue <-TRUE
        i <-1
        while (continue==TRUE & i<reps){
            update_MPCC() ## update css matrix
            conf_q <-c(conf_q,quantile(rowMeans(CoAsso*(1-CoAsso)),0.9)) ## record 90% quantile of the confusion
            ########################
            ### stop the iteration if change of 90% quantile is less than 0.00001 for past 5 iterations
            #######################
            if (i>num_unchange){
                ## find changes of quantile
                cm <-vapply(c(i:(i-(num_unchange-1))), function(x) abs(conf_q[x]-conf_q[x-1]), numeric(1))
                if (max(cm)<eps){
                    continue <-FALSE
                    message(paste0('Stop at iteration ',i))
                }else{
                    continue <-TRUE
                }
            }
            i <-i+1
        }
    }
    message('Done')

    labels <-IMPACC_cluster(css=CoAsso,K=K,finalAlgorithm=finalAlgorithm,finalLinkage=finalLinkage)

    return(list(ConsensusMatrix = CoAsso,labels = labels,nIter = i))
}







######## IMPACC
IMPACC <-function(d=NULL,
                  K=NULL,
                  adaptiveFeature = TRUE,
                  reps=300, ## number of iterations
                  pItem=0.25, ## minipatch size for observations
                  pFeature=0.1, ## minipatch size for features
                  innerLinkage="ward.D", ## internal HC linkage
                  distance="manhattan",## internal HC distance
                  h=0.95, ## cut internal tree at h quantile
                  E= 3, ## number of epochs in burn-in stage
                  qI=0.95, ## observations with greater than (qI) quantile of weights are classified to the set of high uncertain observations in adaptive sampling
                  qF=1, ## features with weights greater than (mean+qF*std) are classified to the set of high importance genes in adaptive sampling
                  alpha_I=0.5, ## learning rate for updating observation weights
                  alpha_F=0.5, ## learning rate for updating feature weights
                  pp=0.05, ## feature support threshold, a feature will be add to feature support if its p-value in ANOVA is smaller than pp
                  finalAlgorithm='hclust',
                  finalLinkage='ward.D',
                  early_stop=TRUE, ## whether perform early stop criteria
                  num_unchange = 5,
                  eps = 0.00001,
                  feature_evaluation = 'ANOVA',
                  verbose=TRUE){
    ###########
    ## The IMPACC function perform IMPACC (Interpretable MiniPatch Adaptive Consensus Clustering) with adaptive sampling scheme on both observations and features.
    ## It returns (1)consensus: final N by N consensus matrix; (2):feature_importance: feature importance scores
    ############

    ######################
    ## check validity of data input
    if ( !is( d )[1] %in% c( "data.frame","matrix" ) ) {
        stop("d must be a data.frame or matrix")
    }

    ## check distance input
    acceptable.distance <- c( "euclidean", "maximum", "manhattan", "canberra", "binary","minkowski",
                              "pearson", "spearman" )
    ## check distance function
    if ( inherits( distance, "character" ) ) {
        if ( ! distance %in%  acceptable.distance  &  ( is(try(get(distance),silent=TRUE))!="function") ){
            stop("unsupported distance.")}
    }else{stop("unsupported distance specified.")}

    if (feature_evaluation == 'ANOVA'){
        get_pv = pv_anova
    }else if(feature_evaluation == 'rankANOVA'){
        get_pv = pv_anova_rank
    }else if(feature_evaluation == 'multinomial'){
        get_pv= pv_multinom
    }else{stop("unsupported feature evaluation method")}
        d <-data.frame(scalematrix(as.matrix(d)))

    update_IMPACC <-function(){
        if(verbose){
            message(paste("  i =",i))
        }
        ### update observation weights
        confusion <-rowMeans(CoAsso*(1-CoAsso))
        ww <-(i+n_burnin)/subsample_o*confusion
        if (sum(ww)!=0){
            wi <<- alpha_I*wi+(1-alpha_I)*(ww/sum(ww))
        }

        if (adaptiveFeature==TRUE){
            wi_p<<- alpha_F*wi_p+(1-alpha_F)*feature_score
        }else{
            wi_p<<-NULL
        }
        ## subsample
        sample_x <-sample_MiniPatch(d, pItem, pFeature,pi_item=pi_item[i],pi_feature=pi_feature[i],
                                    weightsItem = wi, weightsFeature = wi_p,qI=qI,qF=qF)
        subsample_o<<- subsample_o+colnames(d)%in%colnames(d)[sample_x$subcols]
        subsample_f<<- subsample_f+rownames(d)%in%rownames(d)[sample_x$subrows]
        ############################
        ## clustering
        ########################
        this_assignment <-cluster_algo(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage)

        #######################################
        ######### for feature importance
        #######################################
        if (adaptiveFeature==TRUE){
            pvalu <-get_pv(as.matrix(sample_x$submat),as.factor(this_assignment))
            pv <-quantile(na.omit(pvalu),pp)
            pvalue <-pvalu<=pv
            feature_support[rownames(sample_x$submat)[which(pvalue)]]<<- 1+ feature_support[rownames(sample_x$submat)[which(pvalue)]]
            feature_score<<- feature_support/subsample_f
        }else{
            feature_support <-feature_score <<- NULL
        }

        ##########################################
        ##mCount stores number of times a sample pair was sampled together.
        mCount<<- connectivityMatrix(rep(1,length(sample_x$subcols)),
                                     mCount,
                                     sample_x$subcols)

        ##ml stores number of times a sample pair was sampled together.
        ml<<- connectivityMatrix(this_assignment,
                                 ml,
                                 sample_x$subcols)

        this_coA <-ml / mCount
        this_coA[mCount==0]=0
        CoAsso  <<- this_coA

    }
    ## initialize parameters
    n <-ncol(d)
    if (adaptiveFeature==TRUE){
        pi_item <- pi_feature <- seq(0.5,1,by=0.5/reps)
    }else{
        pi_item <- seq(0.5,1,by=0.5/reps)
        pi_feature <- rep(1,reps)
    }

    ############################
    ## calculate # of iterations needed for burn-in stage
    ############################
    if (100%%(pItem*100)==0){
        p_sample <-rep(pItem,1/pItem)
    }else{
        p_sample <-c(rep(pItem,floor(1/pItem)),1%%pItem)
    }

    if (100%%(pFeature*100)==0){
        p_sample_row <-rep(pFeature,1/pFeature)
    }else{
        p_sample_row <-c(rep(pFeature,floor(1/pFeature)),1%%pFeature)
    }
    lcm <-max(length(p_sample_row),length(p_sample))
    num_partition <-E*lcm ## number of partition needed for burn-in stage
    conf_record <-array(0,dim = c(reps+num_partition,n))

    ######################################
    ## BURN IN STAGE
    ##############################################
    message('Burn-in stage')
    sample_Burn <-sampleBurin(d,num_partition,p_sample,p_sample_row) ## 'sampleBurin' is in funmini.R
    n_burnin <-length(sample_Burn)
    subsample_o <-rep(0,ncol(d))
    subsample_f <-rep(0,nrow(d))
    if (adaptiveFeature==TRUE){
        feature_support  <- rep(0,nrow(d)) ## feature support
        names(feature_support) <- rownames(d)
    }

    CoAsso=mCount=ml <- matrix(c(0),ncol=n,nrow=n)

    for (i in seq_len(length(sample_Burn))){
        if(verbose){
            message(paste("  i =",i))
        }

        sample_x  <- sample_Burn[[i]]
        this_assignment  <- cluster_algo(submat = sample_x$submat,h = h,distance=distance,innerLinkage=innerLinkage)
        ## update records on minipatch sampling
        subsample_o  <- subsample_o+colnames(d)%in%colnames(d)[sample_x$subcols]
        subsample_f  <- subsample_f+rownames(d)%in%rownames(d)[sample_x$subrows]
        if (adaptiveFeature==TRUE){
            ##########################################
            ######## for feature importance
            ###########################################
            pvalu  <- get_pv(as.matrix(sample_x$submat),as.factor(this_assignment))
            pv  <- quantile(na.omit(pvalu),pp)
            pvalue  <- pvalu<=pv
            feature_support[rownames(sample_x$submat)[which(pvalue)]] <-1+ feature_support[rownames(sample_x$submat)[which(pvalue)]]
            feature_score <-feature_support/(subsample_f)
            feature_score[subsample_f==0] <-0
        }else{
            feature_support  <- feature_score<-NULL
        }
        ###################################
        mCount <- connectivityMatrix(rep(1,length(sample_x$subcols)),
                                     mCount,
                                     sample_x$subcols)

        ml <- connectivityMatrix(this_assignment,
                                 ml,
                                 sample_x$subcols)
        CoAsso  <- ml / mCount
        CoAsso[mCount==0] <-0

    }

    wi <-rep(1/n,n)
    wi_p <-rep(1/nrow(d),nrow(d))

    ############################################
    ## adaptive stage
    #####################################################
    ##message("adaptive stage")
    message('Adaptive stage')

    if (early_stop == FALSE){
        message('No early stopping')
        for (i in seq(1,reps)){
            update_IMPACC()
        }
    }else{
        conf_q <-NULL
        continue <-TRUE
        i <-1
        while (continue==TRUE & i<reps){
            update_IMPACC()
            conf_q <-c(conf_q,quantile(rowMeans(CoAsso*(1-CoAsso)),0.9)) ## take 90% quantile of the confusion
            if (i>num_unchange){
                ## find change of quantile
                cm  <- vapply(c(i:(i-(num_unchange-1))), function(x) abs(conf_q[x]-conf_q[x-1]),numeric(1))
                if (max(cm)<eps){
                    continue <-FALSE
                    message(paste0('Stop at iteration ',i+num_partition))

                }else{
                    continue <-TRUE
                }
            }
            i <-i+1
        }
    }

    #heatmap(CoAsso)
    labels  <- IMPACC_cluster(CoAsso,K=K,finalAlgorithm=finalAlgorithm,finalLinkage=finalLinkage)
    message('Done')

    return(list(ConsensusMatrix = CoAsso,
                labels=as.character(labels),
                feature_importance=feature_score,
                nIter = i))
}


IMPACC_cluster <-function(css=NULL,K=NULL,finalAlgorithm='hclust',finalLinkage='ward.D'){
    return(tryCatch({
        if(finalAlgorithm=='hclust'){
            hc <-hclust(as.dist(1-css),method=finalLinkage)
            ct  <- cutree(hc,K)
        }else if (finalAlgorithm=='kmedoid'){
            ct  <- cluster::pam(as.dist(1-css), K)$clustering
        }else if (finalAlgorithm=='spectral'){
            ct <- SNFtool::spectralClustering(css, K)
        }
        return(ct)
    }, error=function(e) NA))
}


######### HELPER FUNCTIONS
sampleBurin <- function(d,num_partition,
                        p_sample,
                        p_sample_row){
    list_col  <- lapply(seq_len(ceiling(num_partition/length(p_sample))), function(a) {
        g  <- sample(cut(seq(ncol(d)),ncol(d)*cumsum(c(0,p_sample))))
        partCol  <- split(seq_len(ncol(d)), g)
        names(partCol) <-NULL
        partCol
    })
    list_col  <- unlist(list_col, recursive = FALSE)
    list_row  <- lapply(seq_len(ceiling(num_partition/length(p_sample_row))), function(a) {
        g  <- sample(cut(seq(nrow(d)),nrow(d)*cumsum(c(0,p_sample_row))))
        partRow  <- split(seq_len(nrow(d)), g)
        names(partRow) <-NULL
        partRow
    })
    list_row  <- unlist(list_row, recursive = FALSE)
    res  <- lapply(seq_len(num_partition), function(i){
        list(submat=d[list_row[[i]],list_col[[i]]],
             subrows=list_row[[i]],
             subcols=list_col[[i]])
    })
    return(res)
}

########### update matrix that record joint clustering info (M^{h}),(N × N) connectivity matrix corresponding to dataset D(h)
connectivityMatrix <- function(clusterAssignments, m, sampleKey){
    ##input: named vector of cluster assignments, matrix to add connectivities
    ##output: connectivity matrix
    names(clusterAssignments) <- sampleKey
    #list samples by clusterId
    cls <- lapply(unique(clusterAssignments), function(i) as.numeric(names(clusterAssignments[clusterAssignments%in%i])))
    for ( i in seq_len(length(cls))){
        nelts <- seq_len(ncol(m))
        cl <- as.numeric( nelts %in% cls[[i]] ) ## produces a binary vector
        ## cl = 1 if this obs is in cluster i
        updt <- outer( cl, cl )
        #product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster
        ## (i,j) = 1 if (i,j) are in same cluster
        m <- m + updt
    }
    return(m)
}

### SPECIFY CLUSTERING MODELS
cluster_algo  <- function(submat,h,distance,innerLinkage){
    if(distance=="pearson" | distance=="spearman"){
        dis <- as.dist( 1-cor(submat,method=distance ))
    }else if( is(try(get(distance),silent=TRUE))=="function"){
        dis <- get(distance)(t(submat))
    }else{
        dis  <- dist(t(submat),method = distance)
    }

    this_cluster  <- hclust(dis, method=innerLinkage)
    this_cluster$height <- round(this_cluster$height, 6)
    this_assignment  <- cutree(this_cluster,h=quantile(this_cluster$height,h))
    return(this_assignment=this_assignment)
}




rowsums  <- function(xx){
    if (is.null(dim(xx))){
        xx
    }else{
        rowSums(xx)
    }
}
rowmeans  <- function(xx){
    if (is.null(dim(xx))){
        xx
    }else{
        rowMeans(xx)
    }
}
pv_anova <-function(X,y){
    n_y <-levels(y)
    ss_resi  <- rowsums(vapply(n_y, function(w) rowsums((X[,y==w]-rowmeans(X[,y==w]))^2),numeric(nrow(X))))
    ss_explained  <-  rowSums((X-rowMeans(X))^2)-ss_resi
    df1 <-length(n_y)-1
    df2  <- ncol(X)-length(n_y)
    FF  <- (ss_explained/df1)/(ss_resi/df2 )
    return(pf(FF, df1, df2, lower.tail = FALSE))
}


pv_anova_rank = function(X,y){
    X=as.matrix(X)
    n=length(y)
    n_y=levels(y)
    R = t(apply(X, 1, rank))
    TIES = apply(X, 1, table)
    TIEE = vapply(TIES, function(x) sum(x^3-x),numeric(1))
    STATISTIC = rowsums(vapply(n_y, function(w) rowsums(R[,y==w])^2/sum(y==w), numeric(nrow(X))))
    STATISTICS <- ((12 * STATISTIC / (n * (n + 1)) - 3 * (n + 1)) /
                       (1 - TIEE/ (n^3 - n)))
    PARAMETER <-  nlevels(y) - 1L
    PVAL <- vapply(STATISTICS, function(s) pchisq(s, PARAMETER, lower.tail = FALSE), numeric(1))
    return(PVAL)
}
pv_multinom = function(X,y){
    multinom_score = function(x,y){
        fit=multinom(y ~ x)
        pred = predict(fit, x, "probs")
        if (nlevels(y)==2){
            pred = cbind(1-pred,pred)
            predd = pred[cbind(seq_len(nrow(pred)), y)]
        }else{
            predd = pred[cbind(seq_len(nrow(pred)), y)]
        }
        return(mean(predd))
    }
    return(vapply(seq_len(nrow(X)), function(i) multinom_score(as.matrix(X)[i,],y), numeric(1)))
}


#################################
## SAMPLE MINIPATCHES
##################################

sample_MiniPatch <- function(d,
                             pSamp=NULL,
                             pRow=NULL,
                             weightsItem=NULL, ## vector
                             weightsFeature=NULL,
                             pi_item=1,
                             pi_feature=1,
                             qI=0.95,
                             qF=0.95 ){



    ## returns a list with the sample columns, as well as the sub-matrix & sample features (if necessary)
    space <-  ncol(d)
    sampleN <- floor(space*pSamp)

    if (pi_item <1){ ##adaptive subsampling
        upper  <- which(weightsItem>=quantile(weightsItem,qI))
        sampleN1  <- ceiling(min(sampleN,pi_item*length(upper)))
        if (length(upper)==1){
            sampCols1=upper
        }else{
            sampCols1 <- sort(sample(upper, sampleN1, replace = FALSE,prob = weightsItem[upper]))
        }
        sampCols2 <- sort(sample(seq_len(space)[-upper], sampleN-sampleN1, replace = FALSE))
        sampCols  <- sort(c(sampCols1,sampCols2))
    }else{ ## random sampling by probability
        sampCols  <- sort(sample(seq_len(space),sampleN,replace = FALSE,prob = weightsItem))
    }

    this_sample <- sampRows <- NA
    ## sample rows
    space  <- nrow(d)
    sampleN  <- floor(space*pRow)

    if (pi_feature <1){ ##adaptive subsampling
        if (qF<1&qF>0){
            upper  <- which(weightsFeature>=quantile(weightsFeature,qF))
        }else{
            upper  <- which(weightsFeature>=mean(weightsFeature)+qF*sd(weightsFeature))
        }


        sampleN1  <- ceiling(min(sampleN,pi_feature*length(upper)))
        if (sampleN1>1){
            sampRows1 <- sort(sample(upper, sampleN1, replace = FALSE,prob = weightsFeature[upper]))
        }else{
            sampRows1=upper
        }
        sampRows2 <- sort(sample(seq_len(space)[-upper], sampleN-sampleN1, replace = FALSE))
        sampRows  <- sort(c(sampRows1,sampRows2))
    }else{ ## random sampling
        sampRows <- sort(sample(seq_len(space), sampleN, replace = FALSE,prob = weightsFeature))
    }
    this_sample <- d[sampRows,sampCols]
    # dimnames(this_sample) <- NULL
    return( list(submat=this_sample,
                 subrows=sampRows,
                 subcols=sampCols ))
}

scalematrix <- function(data) {
    cm <- rowMeans(data)
    csd <- matrixStats::rowSds(data, center = cm)
    csd[which(csd==0)]=0.001
    (data - cm) / csd
}

