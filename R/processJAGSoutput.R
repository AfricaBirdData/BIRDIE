# This function and its helpers below are taken directly from the jagsUI package
# https://github.com/kenkellner/jagsUI/blob/master/R/processoutput.R

processJAGSoutput <- function(fit, DIC, params.omit, verbose = TRUE) {

    if(verbose){
        message('Calculating statistics.......')
    }

    # Get parameter names
    params <- colnames(fit[[1]])

    #Get number of chains
    m <- length(fit)

    #Collapse mcmc.lists into matrix
    mat = do.call(rbind,fit)

    #Get # of iterations / chain
    n <- dim(mat)[1] / m

    #Get parameter dimensions
    dim <- get.dim(params)

    #Create new parameter name vectors to handle non-scalar params
    expand <- sapply(strsplit(params, "\\["), "[", 1)
    params.simple <- unique(sapply(strsplit(params, "\\["), "[", 1))

    #Functions for statistics
    qs <- function(x,y){as.numeric(quantile(x,y, na.rm=TRUE))}
    #Overlap 0 function
    ov <- function(x){findInterval(0,sort(c(qs(x,0.025),qs(x,0.975))))==1}
    #f function (proportion of posterior with same sign as mean)
    gf <- function(x){if(mean(x, na.rm=TRUE)>=0){mean(x>=0, na.rm=TRUE)}else{mean(x<0, na.rm=TRUE)}}
    #n.eff function
    calcneff <- function(x,n,m){
        xp <- matrix(x,nrow=n,ncol=m)
        xdot <- apply(xp,2,mean, na.rm=TRUE)
        s2 <- apply(xp,2,var, na.rm=TRUE)
        W <- mean(s2)

        #Non-degenerate case
        if (is.na(W)){
            n.eff <- NA
        } else if ((W > 1.e-8) && (m > 1)) {
            B <- n*var(xdot)
            sig2hat <- ((n-1)*W + B)/n
            n.eff <- round(m*n*min(sig2hat/B,1),0)
            #Degenerate case
        } else {
            n.eff <- 1
        }
        n.eff
    }

    #Gelman diag function
    gd <- function(i,hold){
        r <- try(coda::gelman.diag(hold[,i], autoburnin=FALSE)$psrf[1], silent=TRUE)
        if(inherits(r, "try-error") || !is.finite(r)) {
            r <- NA
        }
        return(r)
    }

    #Make blank lists
    sims.list <- means <- rhat <- n.eff <- se <- as.list(rep(NA,length(params.simple)))
    q2.5 <- q25 <- q50 <- q75 <- q97.5 <- overlap0 <- f <- as.list(rep(NA,length(params.simple)))
    names(sims.list) <- names(means) <- names(rhat) <- names(n.eff) <- params.simple
    names(se) <- names(q2.5) <- names(q25) <- names(q50) <- names(q75) <- names(q97.5) <- params.simple
    names(overlap0) <- names(f) <- params.simple

    #This function modifies objects in global environment (output is discarded)
    #Calculates statistics for each parameter
    calc.stats <- function(prm){

        #If parameter is not a scalar (e.g. vector/array)
        if(!is.na(dim[prm][1])){

            #Get all samples
            sims.list[[prm]] <<- mat[,expand==prm,drop=FALSE]

            #if every iteration is NA, don't do anything else
            if(all(is.na(sims.list[[prm]]))){return(NA)}

            #If more than 1 chain, calculate rhat
            #Done separately for each element of non-scalar parameter to avoid errors
            if(m > 1 && (!prm %in% params.omit)){
                hold <- fit[, expand == prm, drop = FALSE]
                nelements <- sum(expand==prm)
                rhat.vals <- sapply(1:nelements,gd,hold=hold)
                names(rhat.vals) <- colnames(hold[[1]])
                rhat[[prm]] <<- populate(rhat.vals,dim[[prm]])
            } else if (m == 1){
                hold <- fit[,expand==prm]
                rhat[[prm]] <<- array(NA,dim=dim[[prm]])
            }

            #Calculate other statistics
            ld <- length(dim(sims.list[[prm]]))
            means[[prm]] <<- populate(colMeans(sims.list[[prm]], na.rm=TRUE),dim[[prm]])
            if(!prm%in%params.omit){
                se[[prm]] <<- populate(apply(sims.list[[prm]],c(2:ld),sd),dim=dim[[prm]])
                q2.5[[prm]] <<- populate(apply(sims.list[[prm]],c(2:ld),qs,0.025),dim=dim[[prm]])
                q25[[prm]] <<- populate(apply(sims.list[[prm]],c(2:ld),qs,0.25),dim=dim[[prm]])
                q50[[prm]] <<- populate(apply(sims.list[[prm]],c(2:ld),qs,0.5),dim=dim[[prm]])
                q75[[prm]] <<- populate(apply(sims.list[[prm]],c(2:ld),qs,0.75),dim=dim[[prm]])
                q97.5[[prm]] <<- populate(apply(sims.list[[prm]],c(2:ld),qs,0.975),dim=dim[[prm]])
                overlap0[[prm]] <<- populate(apply(sims.list[[prm]],c(2:ld),ov),dim=dim[[prm]])
                f[[prm]] <<- populate(apply(sims.list[[prm]],c(2:ld),gf),dim=dim[[prm]])
                n.eff[[prm]] <<- populate(apply(sims.list[[prm]],c(2:ld),calcneff,n,m),dim=dim[[prm]])
            }

            sims.list[[prm]] <<- populate(sims.list[[prm]],dim=dim[[prm]],simslist=T,samples=dim(mat)[1])

            #If parameter is a scalar
        } else {

            if(m > 1 && (!prm%in%params.omit)){rhat[[prm]] <<- coda::gelman.diag(fit[,prm],autoburnin=FALSE)$psrf[1]}

            sims.list[[prm]] <<- mat[,prm]

            if(all(is.na(sims.list[[prm]]))){return(NA)}

            means[[prm]] <<- mean(sims.list[[prm]], na.rm=TRUE)
            if(!prm%in%params.omit){
                se[[prm]] <<- sd(sims.list[[prm]], na.rm=TRUE)
                q2.5[[prm]] <<- qs(sims.list[[prm]],0.025)
                q25[[prm]] <<- qs(sims.list[[prm]],0.25)
                q50[[prm]] <<- qs(sims.list[[prm]],0.5)
                q75[[prm]] <<- qs(sims.list[[prm]],0.75)
                q97.5[[prm]] <<- qs(sims.list[[prm]],0.975)
                overlap0[[prm]] <<- ov(sims.list[[prm]])
                f[[prm]] <<- gf(sims.list[[prm]])
                n.eff[[prm]] <<- calcneff(sims.list[[prm]],n,m)}
        }

    }

    #Actually run function(nullout not used for anything)
    nullout <- sapply(params.simple, calc.stats)

    #Warn user if at least one Rhat value was NA
    rhat.sub <- unlist(rhat)[!is.na(unlist(means))]
    if(NA%in%rhat.sub&&verbose){
        options(warn=1)
        warning('At least one Rhat value could not be calculated.')
        options(warn=0,error=NULL)
    }

    #Do DIC/pD calculations if requested by user
    if(DIC & 'deviance' %in% params){
        dev <- matrix(data=mat[,'deviance'],ncol=m,nrow=n)
        pd <- numeric(m)
        dic <- numeric(m)
        for (i in 1:m){
            pd[i] <- var(dev[,i])/2
            dic[i] <- mean(dev[,i]) + pd[i]
        }
        pd <- mean(pd)
        dic <- mean(dic)

        #Return this list if DIC/pD requested
        return(list(sims.list=sims.list,mean=means,sd=se,q2.5=q2.5,q25=q25,q50=q50,q75=q75,q97.5=q97.5,overlap0=overlap0,
                    f=f,Rhat=rhat,n.eff=n.eff,pD=pd,DIC=dic))
    } else {
        #Otherwise return list without pD/DIC
        return(list(sims.list=sims.list,mean=means,sd=se,q2.5=q2.5,q25=q25,q50=q50,q75=q75,q97.5=q97.5,overlap0=overlap0,
                    f=f,Rhat=rhat,n.eff=n.eff))
    }

}


#This function gets the dimensions of non-scalar parameters
#for which the user has requested posterior distributions.

get.dim <- function(params){

    #Get all unique parameters (i.e., collapse indexed non-scalars)
    ps <- unique(sapply(strsplit(params, "\\["), "[", 1))
    #Slice indexes from non-scalar parameter entries
    test <- sapply(strsplit(params, "\\["), "[", 1)

    #Calculate dimension for each parameter i
    dim <- lapply(ps, function(i){

        #Extract indices from each element j of parameter i
        w <- params[test==i]
        getinds <- lapply(w,FUN=function(j){

            w2 <- strsplit(j,'\\[')[[1]][2]
            w3 <- strsplit(w2,"\\]")[[1]]
            w4 <- as.numeric(unlist(strsplit(w3,",")))
            return(w4)

        })

        #Get max value from each dimension of i
        collapsedinds <- do.call(rbind,getinds)
        apply(collapsedinds,2,max)

    })

    names(dim) = ps
    dim

}

populate <- function(input,dim,simslist=FALSE,samples=NULL){

    if(!simslist){

        charinds <- sub(".*\\[(.*)\\].*", "\\1", names(input), perl=TRUE)

        fill <- array(NA,dim=dim)

        for (i in 1:length(input)){

            ind <- lapply(strsplit(charinds[i], ','), as.integer)[[1]]
            fill[matrix(ind,1)] <- input[i]

        }
    } else {

        charinds <- sub(".*\\[(.*)\\].*", "\\1", colnames(input), perl=TRUE)

        fill <- array(NA,dim=c(samples,dim))

        for (i in 1:length(charinds)){

            #ind <- lapply(strsplit(charinds[i], ','), as.integer)[[1]]

            eval(parse(text=paste('fill[','1:',samples,',',charinds[i],']','<- input[,i]',sep="")))

        }
    }

    return(fill)

}
