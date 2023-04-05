# Predation mortality as affected by spatial occupancy
#	What is scaled?
#		1. Interaction matrix scaled by benthic occupancy
hypoxiaPredMort = function(params, n, n_pp, n_other, t, pred_rate, ...)  {
	
	# Species indices
	idx_sp = (length(params@w_full) - length(params@w) + 1):length(params@w_full)
	
	# Get benthic and pelagic mortality
	benthic_pm = (base::t(params@interaction) %*% pred_rate$benthic[,idx_sp,drop = FALSE])
	pelagic_pm = (base::t(params@interaction) %*% pred_rate$pelagic[,idx_sp,drop = FALSE])
	
	# Scale and return
	pm = params@other_params$benthic_occupancy[t+1,,] * benthic_pm + (1 - params@other_params$benthic_occupancy[t+1,,]) * pelagic_pm
	return(pm)
}