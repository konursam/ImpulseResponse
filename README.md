# ImpulseResponse
Impulse Response Research - Causality Models
Readme.txt
________________________________________________________________________
/base contains all add on functions used by the files below:
------------------------------------------------------------------------
..each of the below files produces IR estimate plots per specified model
------------------------------------------------------------------------
-LFP.m         :  2 channel model CT
-SPIKE.m       :  2 channel model DT
-LFP_3node.m   :  3 channel model CT
-SPIKE_3node.m :  3 channel model DT
-LFP_5node.m   :  5 channel model CT
-SPIKE_5node.m :  5 channel model DT

Within each of the above m-files, desired AR models/coefficients can be specified
_____________________________________
Impulse Response Polarity Estimation:
-------------------------------------
-Polar_est_LFP.m    :  estimate IR polarity by use of the first derivative method (lfp)
-Polar_est_spike.m. :  estimate IR polarity by use of the first derivative method (spk)
______________________
Significance Analysis:
----------------------
-Permutation_plots_LFP.m    :  produces significance analysis per fmax
-Permutation_plots_Spike.m  :  produces significance analysis per spike rate

-Perm_plots_PERPOINT_lfp.m. :  produces probability densities per IR point generated (per fmax)
-Perm_plots_PERPOINT_spk.m. :  produces probability densities per IR point generated (per spike rate)
