# Universal feeding level scaled by hypoxia exposure
#	What is scaled?
#		1. Maximum intake level scaled by physiological hypoxia response
physoxiaFeedingLevel <- function(params, n, n_pp, n_other, t, encounter, ...) {
	return(encounter / (encounter + (params@intake_max * params@other_params$rs_maxin[t+1,,])))
}