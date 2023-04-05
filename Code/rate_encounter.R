# Standard mizer encounter function with benthic resource only
#	What is scaled?
#		1. Benthic resource encounters scaled by occupancy
#		2. Prey species encounter scaled by occupancy
#		3. Overall encounter scaled by changes in search volume due to physiological hypoxia response
physoxiaEncounter = function(params, n, n_pp, n_other, t, ...) {
	
	# Calculate benthic and pelagic prey
	benthic_prey = mizerEncounter(params, n * params@other_params$benthic_occupancy[t+1,,], n_pp, n_other, t = t)
	pelagic_prey = mizerEncounter(params, n * (1 - params@other_params$benthic_occupancy[t+1,,]), n_other$n_pp_pelagic, n_other, t = t)
	
	# Scale and return total encounter
	encounter = params@other_params$benthic_occupancy[t+1,,] * benthic_prey + (1 - params@other_params$benthic_occupancy[t+1,,]) * pelagic_prey
	encounter = params@other_params$rs_search[t+1,,] * encounter
	return(encounter)
}