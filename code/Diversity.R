#################################################
#                                               #
#   Code for cacluating Diversity metrics       #
#         Created: 12/27/2019                   #
#   *This code was in the Mark file originally  #
#                                               #
#################################################



##### 5. Calculate diversity indices ##### 
library(vegan)
diversity(counts.596, index = 'shannon')
diversity(counts.596, index = 'simpson')
rarefy(counts.596, sample = min(apply(counts.596, 1, sum))) 
rarecurve(counts.596, main = "DSDP Site 596 (S. Pacific)") 

diversity(counts.689, index = 'shannon')
diversity(counts.689, index = 'simpson')
rarefy(counts.689, min(apply(counts.689, 1, sum)))
rarecurve(counts.689, main = "ODP Site 689 (Antarctic)")

par(mfrow = c(1,2))
rarecurve(counts.596, main = "a) DSDP Site 596 (S. Pacific)") 
rarecurve(counts.689, main = "b) ODP Site 689 (Antarctic)")
