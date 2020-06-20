This repo contains scripts and other files used in data analysis used in Wutkowska et al., (2020?) 'Can root-associated' fungi mediate the impact of abiotic conditions on the growth of a High Arctic herb?'

List of files:

* flow.png - an overview of the data analysis procedures, this should clarify the flow (can be found in supplementary material to the manuscript)

BIOINFORMATICS:
* dada2_fencesA.R
* dada2_fencesB.R
* dada2_largeA.R
* dada2_largeB.R
* dada2_tempA.R
* dada2_tempB.R

* seqtabs_fencesA.rds
* seqtabs_fencesB.rds
* seqtabs_largeA.rds
* seqtabs_largeB.rds
* seqtabs_tempA.rds
* seqtabs_tempB.rds

* mergerunsRemovechimerasAssignTaxonomy.R

* asvtable_rarefied_nozerocolumns.txt - ASV table
* functionalannotation_taxtable_rarefied.txt - all the functional annotations for fungal ASVs included in the table (asvtable_rarefied_nozerocolumns.txt)

STATISTCS:
* meteodata.txt - concatenated meteorological data for all the localities in the study
* env_bistorta_temp.txt - input data for bistorta_SEM.R, includes also chosen meteorological data
* bistorta_SEM.R - R script generated for structural equation modeling of abiotic (edaphic, meteorological) and biotic variables (fungal and plant measurements).

