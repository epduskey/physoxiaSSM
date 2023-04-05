Duskey, E. 2023. Metabolic prioritization in response to environmental hypoxia: individual-level data informs community-level analyses. J Exp Mar Biol Ecol.

Data column headers in README_AllData.csv:

	- study: shorthand for each study constructed from first author's last name and year of publication
	- oxygen: oxygen level as reported
	- oxygen_unit: unit of reported oxygen level; either "percent_sat" for % air saturation or "kpa" for kilopascals
	- oxygen_mLL: oxygen level converted to mL*L^-1 using the site "https://loligosystems.azurewebsites.net/convert-oxygen-units" assuming 1 atm air pressure, and temperature and salinity as reported
	- temperature_C: temperature in degrees C
	- salinity_ppt: salinity in parts per thousand i.e. promille
	- activity: reported index or value of fish activity
	- activity_unit: the unit of activity, either NA for unit-less indices, "cm_s" for swim speed in cm*s^-1, "U_crit" for critical swim speed, "BL_s" for body lengths*s^-1, or "distance_m" for swim distance in m
	- SE_activity: reported standard error for activity
	- CV_activity: calculated coefficient of variation for reported activity level
	- reponse_activity: activity within each experiment as a proportion of the maximum reported value
	- ingestion: reported food ingestion
	- ingestion_unit: the unit of ingestion, either "g_day_fish" for average g*day^-1*fish^-1, "total_consumption" for total food consumption in g
	- SE_ingestion: reported or estimated standard error
	- CV_ingestion: calculated coefficient of variation for reported ingestion
	- response_ingestion: ingestion within each experiment as a proportion of the maximum reported value
	- fce: reported food conversion efficiency (FCE)
	- fce_unit: the unit of FCE, always "biomassgain_unitweight" i.e. biomass gain per unit weight of feed
	- SE_fce: reported or estimated standard error for FCE
	- CV_fce: calculated coefficient of variation for reported FCE
	- reponse_fce: FCE within each experiment as a proportion of the maximum reported value
	- nTreatment: number of experimental units
	- species: species used in experiment
	- age_class: age class of fish used in experiment, either "adult" or "juvenile"

Individual study notes:

	Chabot1999:
		
		- Chabot, D. and Dutil, J.D., 1999. Reduced growth of Atlantic cod in non‐lethal hypoxic conditions. Journal of Fish Biology, 55(3), pp.472-491. https://doi.org/10.1111/j.1095-8649.1999.tb00693.x
		- Activity index recorded as the number of fish filmed crossing a line drawn from center tank to edge in a high-visibility region
		- Activity scores tended to decrease with time in some treatments, therefore I calculated pooled mean and SE for each treatment across time periods

	Dutil2007:
		
		- Dutil, J.D., Sylvestre, E.L., Gamache, L., Larocque, R. and Guderley, H., 2007. Burst and coast use, swimming performance and metabolism of Atlantic cod Gadus morhua in sub‐lethal hypoxic conditions. Journal of Fish Biology, 71(2), pp.363-375. https://doi.org/10.1111/j.1095-8649.2007.01487.x
		- No second order statistics reported for FCE and ingestion in this study, so I used the mean of the coefficient of variation (CV) of all other studies I included in the meta-analysis for FCE and ingestion

	Herbert2005:
		
		- Herbert, N.A. and Steffensen, J.F., 2005. The response of Atlantic cod, Gadus morhua, to progressive hypoxia: fish swimming speed and physiological stress. Marine Biology, 147(6), pp.1403-1412. https://doi.org/10.1007/s00227-005-0003-8
		- Swim speed reported during steady and unsteady oxygen values, I only used those reported during steady oxygen

	Herbert2011:
		
		- Herbert, N.A., Skjæraasen, J.E., Nilsen, T., Salvanes, A.G. and Steffensen, J.F., 2011. The hypoxia avoidance behaviour of juvenile Atlantic cod (Gadus morhua L.) depends on the provision and pressure level of an O 2 refuge. Marine Biology, 158, pp.737-746. https://doi.org/10.1007/s00227-010-1601-7
		- Fish had a choice between a low and incrementally increasing oxygen habitat; mean oxygen and % time spent in each habitat was reported, so I used the latter to calculate a weighted mean of oxygen exposure
		- Oxygen and activity were also reported for an initial phase prior to the experiment, but I excluded these values
	
	Johansen2006:
		
		- Johansen, J.L., Herbert, N.A. and Steffensen, J.F., 2006. The behavioural and physiological response of Atlantic cod Gadus morhua L. to short‐term acute hypoxia. Journal of fish biology, 68(6), pp.1918-1924. https://doi.org/10.1111/j.1095-8649.2006.01080.x
		- Swim speed reported during a steady, high oxygen phase, and an experimental phase during which oxygen declined rapidly and acutely; I only included values during the experimental phase
	
	Koedijk2012:
		
		- Koedijk, R., Imsland, A.K., Folkvord, A., Stefansson, S.O., Jonassen, T.M. and Foss, A., 2012. Larval rearing environment influences the physiological adaptation in juvenile Atlantic cod, Gadus morhua. Aquaculture International, 20, pp.467-479. https://doi.org/10.1007/s10499-011-9478-0
		- Fish were exposed to variable oxygen only, or variable oxygen and high ammonia concentrations; I excluded the values with high ammonia
	
	Petersen2010

		- Petersen, L.H. and Gamperl, A.K., 2010. Effect of acute and chronic hypoxia on the swimming performance, metabolic capacity and cardiac function of Atlantic cod (Gadus morhua). Journal of Experimental Biology, 213(5), pp.808-819. https://doi.org/10.1242/jeb.033746
		- Critical swimming speed (U_crit) reported for fish acclimated to hypoxia or normoxia; I included values from both groups as neither were significantly different
	
	Schurmann1994
		
		- Schurmann, H. and Steffensen J., 1994. SPONTANEOUS SWIMMING ACTIVITY OF ATLANTIC COD GADUS MORHUA EXPOSED TO GRADED HYPOXIA AT THREE TEMPERATURES. Journal of Experimental Biology, 197(1): 129–142. https://doi.org/10.1242/jeb.197.1.129
		- No second order statistics reported for swim distance, but minimum and maximum were provided; I used the formula (max-min)/4 to calculate standard deviation, and these values to calculate standard error and coefficient of variation
	
	Skjæraasen2008
		
		- Skjæraasen, J.E., Nilsen, T., Meager, J.J., Herbert, N.A., Moberg, O., Tronci, V., Johansen, T. and Salvanes, A.G.V., 2008. Hypoxic avoidance behaviour in cod (Gadus morhua L.): The effect of temperature and haemoglobin genotype. Journal of Experimental Marine Biology and Ecology, 358(1), pp.70-77. https://doi.org/10.1016/j.jembe.2008.01.010
		- Fish had a choice between a hypoxic and relatively normoxic habitat; mean oxygen and % time spent in each habitat was reported, so I used the latter to calculate a weighted mean of oxygen exposure
	
	Thorarensen2017
		
		- Thorarensen, H., Gústavsson, A., Gunnarsson, S., Árnason, J., Steinarsson, A., Björnsdóttir, R. and Imsland, A.K., 2017. The effect of oxygen saturation on the growth and feed conversion of juvenile Atlantic cod (Gadus morhua L.). Aquaculture, 475, pp.24-28. https://doi.org/10.1016/j.aquaculture.2017.04.002
		- Feed conversion ratio (unit weight of feed*biomass gain^-1) reported rather than FCE, so I inverted the reported values to convert them to the same units as the other values