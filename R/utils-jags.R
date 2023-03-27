#' Plot state series of a JAGS state-space model with two seasons
#'
#' @param fit A JAGS state-space model fitted to CWAC data
#' @param ssm_counts A data frame with the count data use to fit the
#' state-space model
#' @param linear If TRUE (default) abundance estimates and data are
#' transformed back to its original scale.
#' @param plot_options A list with two elements: colors - the colours of the
#' points that will appear in the plot (two values), and pers_theme  - A
#' personalized ggplot theme.
#'
#' @return A list with two elements: i) plot: a plot with summer and winter fitted states, as well as the long-term trend,
#' ii) data: the data used to create the individual plots. This is useful for extracting the data used by ggplot to render the plots
#' (e.g. for exporting to the dashboard)
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' counts <- barberspan
#' ssmcounts <- prepSsmData(counts, species = NULL)
#' fit <- fitCwacSsm(ssmcounts, mod_file = "mymodel.jags",
#' param = c("beta", "lambda", "sig.zeta",
#' "sig.w", "sig.eps", "sig.alpha", "sig.e", "mu_t", "mu_wt"))
#' plotSsm2ss(fit = fit, ssm_counts = ssmcounts)
#' }
plotJagsSsm2ss <- function(fit, ssm_counts, linear = TRUE,
                           plot_options = list(colors = NULL, pers_theme = NULL)){

    if(is.null(plot_options$colors)){
        plot_options$colors <- c("#71BD5E", "#B590C7")
    }

    # Data summary
    dat_summ <- ssm_counts %>%
        dplyr::mutate(season_est = as.integer(fit$mean$summer < 0.5) + 1L) %>%
        dplyr::group_by(site_id, year, season_est) %>%
        dplyr::summarise(count = log(mean(count+1, na.rm = TRUE))) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider(names_from = season_est, values_from = count) %>%
        dplyr::rename(count_s = "1",
                      count_w = "2")


    # Create a data frame with the posterior state
    post_stt <- data.frame(year = dat_summ$year,
                           site_id = dat_summ$site_id,
                           stt_s_est = as.vector(t(fit$mean$stt_s)),
                           stt_s_lb = as.vector(t(fit$q2.5$stt_s)),
                           stt_s_ub = as.vector(t(fit$q97.5$stt_s)),
                           stt_w_est = as.vector(t(fit$mean$stt_w)),
                           stt_w_lb = as.vector(t(fit$q2.5$stt_w)),
                           stt_w_ub = as.vector(t(fit$q97.5$stt_w)))

    # Pivot seasons
    dat_summ <- dat_summ %>%
        tidyr::pivot_longer(cols = c(count_s, count_w),
                            names_to = "season", values_to = "count") %>%
        dplyr::mutate(season = ifelse(grepl("_s", season), "summer", "winter"))

    post_stt <- post_stt %>%
        tidyr::pivot_longer(cols = -c(year, site_id),
                            names_to = "quantile", values_to = "value") %>%
        dplyr::mutate(season = ifelse(grepl("_s", quantile), "summer", "winter")) %>%
        dplyr::mutate(quantile = dplyr::case_when(grepl("_est", quantile) ~ "est",
                                                  grepl("_ub", quantile) ~ "ub",
                                                  grepl("_lb", quantile) ~ "lb"))

    # Add counts and location code
    post_stt <- post_stt %>%
        dplyr::left_join(dplyr::distinct(ssm_counts[,c("site_id", "LocationCode")]),
                         by = "site_id") %>%
        dplyr::rename(loc_code = LocationCode) %>%
        dplyr::left_join(dat_summ, by = c("site_id", "year", "season"))



    if(linear){
        post_stt <- post_stt %>%
            dplyr::mutate(value = exp(value) - 1,
                          count = exp(count) - 1) %>%
            dplyr::mutate(value = ifelse(value < 0, 0, value))

        abund_label <- "Abundance"

        # Cut axis when values are larger than 10 times the max count
        # ylims <- c(0, max(post_stt$count, na.rm = T) * 10)

    } else {

        abund_label <- "log abundance"
        ylims <- c(NA, NA)

    }


    site_plot_data <- vector("list", dplyr::n_distinct(post_stt$loc_code))


    for(i in seq_along(site_plot_data)){

        loc_code_sel <- post_stt %>%
            dplyr::filter(site_id == i) %>%
            dplyr::distinct(loc_code) %>%
            dplyr::pull(loc_code)

        # Abundance by season -----------------------------------------------------

        stt_plot <- post_stt %>%
            dplyr::filter(site_id == i) %>%
            ggplot() +
            geom_path(aes(x = year, y = value, linetype = quantile)) +
            geom_point(aes(x = year, y = count, col = season), show.legend = FALSE) +
            scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
            scale_colour_manual(values = plot_options$colors) +
            # coord_cartesian(ylim = ylims) +
            facet_wrap("season", ncol = 2) +
            xlab("Year") + ylab(abund_label) +
            plot_options$pers_theme


        # Trend plot --------------------------------------------------------------

        # Create a data frame with the posterior trend
        post_trd <- data.frame(beta_est = c(fit$mean$beta[i,], NA),
                               beta_lb = c(fit$q2.5$beta[i,], NA),
                               beta_ub = c(fit$q97.5$beta[i,], NA),
                               prop_est = fit$mean$lambda[i,],
                               prop_lb = fit$q2.5$lambda[i,],
                               prop_ub = fit$q97.5$lambda[i,],
                               year = unique(ssm_counts$year))

        if(linear){
            post_trd <- post_trd %>%
                dplyr::mutate(dplyr::across(.cols = -year,
                                            .fns = ~exp(.x)))
            rate_label <- "Growth rate"
            ratio_label <- "W/S ratio"

        } else {

            rate_label <- "log growth rate"
            ratio_label <- "log W/S ratio"

        }

        # Plot trend
        trd_plot <- post_trd %>%
            dplyr::select(-dplyr::starts_with("prop")) %>%
            tidyr::pivot_longer(cols = -year,
                                names_to = "quantile") %>%
            ggplot() +
            geom_path(aes(x = year, y = value, linetype = quantile)) +
            scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
            scale_colour_manual(values = plot_options$colors) +
            xlab("Year") + ylab(rate_label) +
            plot_options$pers_theme


        # Winter/summer ratio -----------------------------------------------------

        # Plot proportion summer/winter
        prop_plot <- post_trd %>%
            dplyr::select(-dplyr::starts_with("beta")) %>%
            tidyr::pivot_longer(cols = -year,
                                names_to = "quantile") %>%
            ggplot() +
            geom_path(aes(x = year, y = value, linetype = quantile)) +
            scale_linetype_manual(name = "", values = c(1, 2, 2), guide = NULL) +
            scale_colour_manual(values = plot_options$colors) +
            xlab("Year") + ylab(ratio_label) +
            plot_options$pers_theme


        # Save --------------------------------------------------------------------

        # Title
        plottitle <- paste(ifelse(unique(ssm_counts$spp) == "multi",
                                  "Multiple species",
                                  unique(ssm_counts$spp)),
                           "at site",
                           loc_code_sel)

        # To prevent opening multiple devices
        pfile <- tempfile()
        grDevices::png(pfile)
        p <- gridExtra::grid.arrange(stt_plot, trd_plot, prop_plot,
                                     layout_matrix = matrix(c(1,2,1,3), nrow = 2),
                                     top = plottitle)
        grDevices::dev.off()
        unlink(pfile)

        site_plot_data[[i]] <- list(plot = p,
                                    data = list(post_stt %>%
                                                    dplyr::filter(site_id == i),
                                                post_trd))

        names(site_plot_data)[[i]] <- loc_code_sel

    }

    return(site_plot_data)


}







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
