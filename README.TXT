README.TXT

This GitHub repository containes the code used in analysis of all datasets used in: 
"No state change in pelagic fish production and biodiversity during the Eocene-Oligocene Transition"

by Elizabeth C. Sibert, Michelle E. Zill, Ella T. Frigyik, and Ricahrd Norris

Code written and maintained by Elizabeth Sibert
Contact: esibert@fas.harvard.edu

Last updated; 12/27/2019

Code Files included: 
    setup.R - calls in datasets and does data management
    Diversity.R - calculations for fish biodiversity metrics reported in the paper
    figures.R - replicate the plots for the figures in the manuscript 
    functions.R - all internal functions used in the analyses
    mark.R - mark-recapture analyses
    scale_factor_for886.R - calcuation of how to scale ODP Site 886's IAR from >106um to >38um, in comparison to all other sites in the study


Data files included:
    EOT_dataset_comp.csv - all proxy and IAR datasets used in the study
    DSDP_596_EO_fish.csv - All teeth, ages, and identifications for DSDP Site 596
    ODP_689_EO_fish_csv.csv - All teeth, ages, and identifications for ODP Site 689
    temp-cramer.csv - Cramer et al 2011 compilation of temperature used in Figure 2
    sealevel-cramer.csv - Cramer et al 2011 compilation of sea level used in Figure 2

MARK directory includes all MARK input and output files generated by RMark package. Code is in mark.R

