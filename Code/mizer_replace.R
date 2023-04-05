# Altered predation rate function
#	What is altered?
#		1. Two dimnames(pred_rate) lines to accommodate new list form
getPredRate <- function (params, n = initialN(params), n_pp = initialNResource(params), 
    n_other = initialNOther(params), t = 0, ...) 
{
    params <- validParams(params)
    no_sp <- dim(params@interaction)[1]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    f <- get(params@rates_funcs$PredRate)
    pred_rate <- f(params, n = n, n_pp = n_pp, n_other = n_other, 
        t = t, feeding_level = getFeedingLevel(params, n = n, 
            n_pp = n_pp, n_other = n_other, time_range = t))
    dimnames(pred_rate$benthic) <- list(sp = dimnames(params@initial_n)$sp, 
        w_prey = as.character(signif(params@w_full, 3)))
    dimnames(pred_rate$pelagic) <- list(sp = dimnames(params@initial_n)$sp, 
        w_prey = as.character(signif(params@w_full, 3)))
    pred_rate
}

# Altered predation mortality function
#	What is altered?
#		1. Nothing?? It just won't work unless I replicate the code
getPredMort <- function (object, n, n_pp, n_other, time_range, drop = TRUE, 
    ...) 
{
    if (is(object, "MizerParams")) {
        params <- validParams(object)
        if (missing(n)) 
            n <- params@initial_n
        if (missing(n_pp)) 
            n_pp <- params@initial_n_pp
        if (missing(n_other)) 
            n_other <- params@initial_n_other
        if (missing(time_range)) 
            time_range <- 0
        t <- min(time_range)
        f <- get(params@rates_funcs$PredMort)
        pred_mort <- f(params, n = n, n_pp = n_pp, n_other = n_other, 
            t = t, pred_rate = getPredRate(params, n = n, n_pp = n_pp, 
                n_other = n_other, t = t))
        dimnames(pred_mort) <- list(prey = dimnames(params@initial_n)$sp, 
            w_prey = dimnames(params@initial_n)$w)
        pred_mort
    } else {
        sim <- object
        if (missing(time_range)) {
            time_range <- dimnames(sim@n)$time
        }
        time_elements <- get_time_elements(sim, time_range)
        pred_mort_time <- plyr::aaply(which(time_elements), 1, 
            function(x) {
                n <- array(sim@n[x, , ], dim = dim(sim@n)[2:3])
                dimnames(n) <- dimnames(sim@n)[2:3]
                n_other <- sim@n_other[x, ]
                names(n_other) <- dimnames(sim@n_other)$component
                t <- as.numeric(dimnames(sim@n)$time[[x]])
                n_pp <- sim@n_pp[x, ]
                return(getPredMort(sim@params, n = n, n_pp = n_pp, 
                  n_other = n_other, time_range = t))
            }, .drop = FALSE)
        names(dimnames(pred_mort_time))[[1]] <- "time"
        pred_mort_time <- pred_mort_time[, , , drop = drop]
        return(pred_mort_time)
    }
}

# Altered plot predation mortality function
#	What is altered?
#		1. Also nothing, it just won't work unless I replicate the code
plotPredMort <- function (object, species = NULL, time_range, all.sizes = FALSE, 
    highlight = NULL, return_data = FALSE, ...) 
{
    assert_that(is.flag(all.sizes), is.flag(return_data))
    if (is(object, "MizerSim")) {
        if (missing(time_range)) {
            time_range <- max(as.numeric(dimnames(object@n)$time))
        }
        params <- object@params
    }
    else {
        params <- validParams(object)
    }
    pred_mort <- getPredMort(object, time_range = time_range, 
        drop = FALSE)
    if (length(dim(pred_mort)) == 3) {
        pred_mort <- apply(pred_mort, c(2, 3), mean)
    }
    species <- valid_species_arg(params, species)
    pred_mort <- pred_mort[as.character(dimnames(pred_mort)[[1]]) %in% 
        species, , drop = FALSE]
    plot_dat <- data.frame(w = rep(params@w, each = length(species)), 
        value = c(pred_mort), Species = species)
    if (!all.sizes) {
        for (sp in species) {
            plot_dat$value[plot_dat$Species == sp & (plot_dat$w < 
                params@species_params[sp, "w_min"] | plot_dat$w > 
                params@species_params[sp, "w_inf"])] <- NA
        }
        plot_dat <- plot_dat[complete.cases(plot_dat), ]
    }
    if (return_data) 
        return(plot_dat)
    p <- plotDataFrame(plot_dat, params, xlab = "Size [g]", 
        xtrans = "log10", highlight = highlight)
    suppressMessages(p <- p + scale_y_continuous(labels = prettyNum, 
        name = "Predation mortality [1/year]", limits = c(0, 
            max(plot_dat$value))))
    p
}

# Altered plot predation mortality function
#	What is altered?
#		1. Also nothing, it just won't work unless I replicate the code
setBevertonHolt <- function(params, R_factor = deprecated(), erepro, R_max, reproduction_level) {
    assert_that(is(params, "MizerParams"))
    no_sp <- nrow(params@species_params)
    num_args <- hasArg("erepro") + hasArg("R_max") + 
        hasArg("reproduction_level") + hasArg("R_factor")
    if (num_args > 1) {
        stop("You should only provide `params` and one other argument.")
    }
    if (num_args == 0) {
        erepro <- species_params(params)$erepro
    }
    if (!missing("erepro")) 
        values <- erepro
    if (hasArg("R_max")) 
        values <- R_max
    if (hasArg("reproduction_level")) 
        values <- reproduction_level
    if (hasArg("R_factor")) 
        values <- R_factor
    if (length(values) == 1 && is.null(names(values))) {
        values <- rep(values, no_sp)
    }
    if (length(values) != no_sp && is.null(names(values))) {
        stop("You need to supply a vector of length ", 
            no_sp, " or a single number or a named vector.")
    }
    if (is.null(names(values))) {
        names(values) <- params@species_params$species
    }
    values <- values[!is.na(values)]
    if (length(values) == 0) 
        return(params)
    if (!all(is.numeric(values))) {
        stop("You provided invalid non-numeric values.")
    }
    species <- valid_species_arg(params, names(values))
    sp_idx <- match(species, params@species_params$species)
    rdd <- getRDD(params)[species]
    rdi <- getRDI(params)[species]
    if (any(rdi == 0)) {
        stop("Some species have no reproduction.")
    }
    params@rates_funcs$RDD <- "BevertonHoltRDD"
    rdd_new <- getRequiredRDD(params)[species]
    if (!missing(erepro)) {
        erepro_new <- values
        erepro_old <- params@species_params$erepro[sp_idx]
        rdi_new <- rdi * erepro_new/erepro_old
        wrong <- rdi_new < rdd_new
        if (any(wrong)) {
            rdi_new[wrong] <- rdd_new[wrong]
            erepro_new[wrong] <- (erepro_old * rdi_new/rdi)[wrong]
            warning("For the following species `erepro` ", 
                "has been increased to the smallest ", 
                "possible value: ", paste0("erepro[", 
                  species[wrong], "] = ", signif(erepro_new[wrong], 
                    3), collapse = "; "))
        }
        r_max_new <- rdi_new * rdd_new/(rdi_new - rdd_new)
        r_max_new[is.nan(r_max_new)] <- Inf
        params@species_params$erepro[sp_idx] <- erepro_new
        params@species_params$R_max[sp_idx] <- r_max_new
        params@time_modified <- lubridate::now()
        return(params)
    }
    if (!missing(reproduction_level)) {
        if (!all(values >= 0 & values < 1)) {
            stop("The reproduction level must be smaller than 1 and non-negative.")
        }
        r_max_new <- rdd_new/values
    }
    if (!missing(R_factor)) {
        if (!all(values > 1)) {
            stop("The R_factor must be greater than 1.")
        }
        r_max_new <- rdd_new * values
    }
    if (!missing(R_max)) {
        wrong <- values < rdd_new
        if (any(wrong)) {
            warning("For the following species the requested `R_max` ", 
                "was too small and has been increased to give a ", 
                "reproduction level of 0.99: ", paste(species[wrong], 
                  collapse = ", "))
            values[wrong] <- rdd_new[wrong]/0.99
        }
        r_max_new <- values
    }
    rdi_new <- rdd_new/(1 - rdd_new/r_max_new)
    params@species_params$R_max[sp_idx] <- r_max_new
    params@species_params$erepro[sp_idx] <- params@species_params$erepro[sp_idx] * 
        rdi_new/rdi
    wrong <- params@species_params$erepro[sp_idx] > 1
    if (any(wrong)) {
        warning("The following species require an unrealistic reproductive ", 
            "efficiency greater than 1: ", paste(species[wrong], 
                collapse = ", "))
    }
    params@time_modified <- lubridate::now()
    return(params)
}

# Altered plot predation mortality function
#	What is altered?
#		1. Also nothing, it just won't work unless I replicate the code
getMort <- function (params, n = initialN(params), n_pp = initialNResource(params), 
    n_other = initialNOther(params), effort = getInitialEffort(params), 
    t = 0, ...) 
{
    params <- validParams(params)
    f <- get(params@rates_funcs$Mort)
    z <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t, 
        f_mort = getFMort(params, effort), pred_mort = getPredMort(params, 
            n = n, n_pp = n_pp, n_other = n_other, time_range = t))
    dimnames(z) <- list(prey = dimnames(params@initial_n)$sp, 
        w_prey = dimnames(params@initial_n)$w)
    return(z)
}

# Altered plot predation mortality function
#	What is altered?
#		1. Also nothing, it just won't work unless I replicate the code
getRequiredRDD <- function(params) {
    # Calculate required rdd
    mumu <- getMort(params)
    gg <- getEGrowth(params)
    rdd_new <- getRDD(params) # to get the right structure
    for (i in seq_len(nrow(params@species_params))) {
        gg0 <- gg[i, params@w_min_idx[i]]
        if (!(gg0 > 0)) {
            stop("Eggs of species ", params@species_params$species[i],
                 " have zero growth rate.")
        }
        mumu0 <- mumu[i, params@w_min_idx[i]]
        DW <- params@dw[params@w_min_idx[i]]
        n0 <- params@initial_n[i, params@w_min_idx[i]]
        if (!(n0 > 0)) {
            stop("Species ", params@species_params$species[i],
                 "appears to have no eggs.")
        }
        rdd_new[i] <- n0 * (gg0 + DW * mumu0)
    }
    rdd_new
}