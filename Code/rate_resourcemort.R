# Resource mortality
#	What is scaled?
#		1. Predation rate is applied to only the portion of the benthic resource exposed to predation due to preference
benthicResourceMort = function(params, n, n_pp, n_other, t, pred_rate, ...) {
	as.vector(params@species_params$interaction_resource %*% pred_rate$benthic)
}