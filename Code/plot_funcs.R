# Miscellaneous functions for plotting during calibration process

# Load packages
library(assertthat)

# Sets initial value of params to the no_t^th time step in sim
#	params: mizer params object
#	sim: mizer sim object
#	no_t: the time step of sim
#	returns modified params object with initial values to to no_t^th value in sim
myInitialValues = function (params, sim, no_t) {
    assert_that(is(params, "MizerParams"), is(sim, "MizerSim"))
    if (!identical(dim(sim@n)[2:3], dim(params@initial_n))) {
        stop("The consumer size spectrum of the simulation in `sim` has a ", 
            "different size from that in `params`.")
    }
    if (!identical(length(sim@n_pp[no_t, ]), length(params@initial_n_pp))) {
        stop("The resource size spectrum of the simulation in `sim` has a ", 
            "different size from that in `params`.")
    }
    if (!identical(length(sim@n_other[no_t, ]), length(params@initial_n_other))) {
        stop("The number of other components in the simulation in `sim` is ", 
            "different from that in `params`.")
    }
    if (!identical(length(sim@effort[no_t, ]), length(params@initial_effort))) {
        stop("The number of gears in the simulation in `sim` is ", 
            "different from that in `params`.")
    }
    if (!identical(dimnames(sim@effort)[[2]], names(params@initial_effort))) {
        stop("The gears in the simulation in `sim` have different names ", 
            "from those in `params`.")
    }
    params@initial_n[] <- sim@n[no_t, , ]
    params@initial_n_pp[] <- sim@n_pp[no_t, ]
    params@initial_n_other[] <- sim@n_other[no_t, ]
    params@initial_effort[] <- sim@effort[no_t, ]
    params@time_modified <- lubridate::now()
    params
}

# Gets growth curves for the no_t^th time step
#	object: mizer sim or params object
#	species: species for which to get growth curves
#	max_age: maximum age for growth
#	no_t: if object is sim, the time step for which to get growth
#	percentage: if TRUE, plots growth as a percentage of maximum
#	returns size at age
myGrowthCurves = function (object, species = NULL, max_age = 20, no_t, percentage = FALSE) {
    if (is(object, "MizerSim")) {
        params <- object@params
        params <- myInitialValues(params, object, no_t)
    }
    else if (is(object, "MizerParams")) {
        params <- validParams(object)
    }
    else {
        stop("The first argument to `getGrowthCurves()` must be a ", 
            "MizerParams or a MizerSim object.")
    }
    species <- valid_species_arg(params, species)
    idx <- which(params@species_params$species %in% species)
    species <- params@species_params$species[idx]
    age <- seq(0, max_age, length.out = 50)
    ws <- array(dim = c(length(species), length(age)), dimnames = list(Species = species, 
        Age = age))
    g <- getEGrowth(params)
    for (j in seq_along(species)) {
        i <- idx[j]
        g_fn <- stats::approxfun(c(params@w, params@species_params$w_inf[[i]]), 
            c(g[i, ], 0))
        myodefun <- function(t, state, parameters) {
            return(list(g_fn(state)))
        }
        ws[j, ] <- deSolve::ode(y = params@w[params@w_min_idx[i]], 
            times = age, func = myodefun)[, 2]
        if (percentage) {
            ws[j, ] <- ws[j, ]/params@species_params$w_inf[i] * 
                100
        }
    }
    return(ws)
}