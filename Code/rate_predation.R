# Universal predation rate scaled by hypoxia exposure
#	What is scaled?
#		1. Feeding level has already been scaled, so only scale search volume to match with encounter
physoxiaPredRate <- function(params, n, n_pp, n_other, t, feeding_level, ...) {

	# Matrix dimensions
	no_sp = dim(params@interaction)[1]
	no_w = length(params@w)
	no_w_full = length(params@w_full)
    
	# If the feeding kernel does not have a fixed predator/prey mass ratio, use mizer version 0.3
	if (!is.null(comment(params@pred_kernel))) {
       	
		# Return list
		pred_rate = vector(mode = "list")
		
		# Total size bins
		n_total_in_size_bins_benthic = sweep(n*params@other_params$benthic_occupancy[t+1,,], 2, params@dw, '*', check.margin = FALSE)
		n_total_in_size_bins_pelagic = sweep(n*(1-params@other_params$benthic_occupancy[t+1,,]), 2, params@dw, '*', check.margin = FALSE)        

		# Predation rate scaled by hypoxia exposure
		pred_rate$benthic = sweep(params@pred_kernel, c(1,2),
					(1 - feeding_level) * 
					(params@search_vol * params@other_params$rs_search[t+1,,]) * 
					n_total_in_size_bins_benthic,
					"*", check.margin = FALSE)
		pred_rate$pelagic = sweep(params@pred_kernel, c(1,2),
					(1 - feeding_level) * 
					(params@search_vol * params@other_params$rs_search[t+1,,]) * 
					n_total_in_size_bins_pelagic,
					"*", check.margin = FALSE)
        
		# Integrate over all predator sizes
		pred_rate$benthic = colSums(aperm(pred_rate$benthic, c(2, 1, 3)), dims = 1)
		pred_rate$pelagic = colSums(aperm(pred_rate$pelagic, c(2, 1, 3)), dims = 1)
		return(pred_rate)
	} else {
    		
		# Return list
		pred_rate = vector(mode = "list")
		
		# Get indices of w_full that give w
		idx_sp = (no_w_full - no_w + 1):no_w_full
    
		# Express the result as a a convolution of Q and ft_pred_kernel_p
		Q_benthic = matrix(0, nrow = no_sp, ncol = no_w_full)
		Q_pelagic = Q_benthic
    
		# Fill the end of each row of Q with the proper values
		Q_benthic[, idx_sp] = sweep( (1 - feeding_level) * (params@search_vol * params@other_params$rs_search[t+1,,]) * (n * params@other_params$benthic_occupancy[t+1,,]), 2, params@dw, "*")
		Q_pelagic[, idx_sp] = sweep( (1 - feeding_level) * (params@search_vol * params@other_params$rs_search[t+1,,]) * (n * (1-params@other_params$benthic_occupancy[t+1,,])), 2, params@dw, "*")

		# Spectral integration in parallel over the different species
		pred_rate$benthic = Re(t(mvfft(t(params@ft_pred_kernel_p) * mvfft(t(Q_benthic)), inverse = TRUE))) / no_w_full
		pred_rate$pelagic = Re(t(mvfft(t(params@ft_pred_kernel_p) * mvfft(t(Q_pelagic)), inverse = TRUE))) / no_w_full
    
		# Due to numerical errors we might get negative or very small entries that should be 0
		pred_rate$benthic[pred_rate$benthic < 1e-18] = 0
		pred_rate$pelagic[pred_rate$pelagic < 1e-18] = 0
		
		# Entries beyond maximum size are zero
		pred_rate$benthic = pred_rate$benthic * params@ft_mask
		pred_rate$pelagic = pred_rate$pelagic * params@ft_mask
    
		return(pred_rate)
    	}
}