# Downing-et-al-2020a
Code and data for: Downing PA, Griffin AS. &amp; Cornwallis CK. 2020a. Group formation and the evolutionary pathway to social complexity in birds. Nature Ecology &amp; Evolution, 4: 479-486.

Number of supplementary items: four
1. GF_R_Code.R
2. GF_Tables_S1-S10.xlsx
3. GF_Data_Extraction.txt
4. GF_Supp_Methods.pdf


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: GF_Tables_S1-S10.xlsx

This Excel document contains the following sheets:

- Supplementary Table 1 (breeding systems)
	+ breeding systems and polyandry estimates for as many bird species as possible
	+ export as a csv.doc and read into R using: poly <- read.csv("data/polyandry.csv")
	+ 4884 data rows = 4730 species (some species have multiple polyandry estimates)
	+ column descriptions:\
		A. animal = latin binomial (matches the Jetz et al. nomenclature)\
		B. commonName = English name of each species\
		C. breedingSystem = whether each species is cooperative or non-cooperative\
		D. groupFormation = whether cooperative groups form as families or non-families\
		E. Riehl 2013 = Riehl's breeding system classification for each cooperative species\
		F. perEPP = percentage of chicks that result from extra-pair paternity\
		G. EPPN = number of chicks studied\
		H. perEPPbr = percentage of broods with extra-pair chicks\
		I. EPPbrN = number of broods studied\
		J. perEGP = percentage of chicks that result from extra-group paternity\
		K. EGPN = number of chicks studied\
		L. perEGPbr = percentage of broods with extra-group chicks\
		M. EGPbrN = number of broods studied\
		N. source = where the data were extracted from in each study\
		O. Cornwallis et al. 2017 = comparison with polyandry values in Cornwallis et al.\
		P. Brouwer & Griffith 2019 = comparison with polyandry values in Brouwer & Griffith\
		Q. parentage reference &/or breeding system\
		R. further references


- Supplementary Table 2 (group size)
	+ group size information for as many cooperatively breeding birds as possible
	+ export as a csv.doc and read into R using: groupSize <- read.csv("data/groupSize.csv")
	+ 172 data rows = 172 species (note that there are missing data)
	+ column descriptions:\
		A. animal = latin binomial (matches the Jetz et al. nomenclature)\
		B. commonName = English name of each species\
		C. groupFormation = whether cooperative groups form as families or non-families\
		D. mean = mean group size\
		E. median = median group size\
		F. mode = modal group size\
		G. min = minimum group size\
		H. max = maximum group size\
		I. nGroupYears = number of group years on which group size data based\
		J. notes = where the data were extracted from in each study\
		K. reference = study from which data were extracted


- Supplementary Table 3 (division of reproduction)
	+ data on how reproduction is divided in groups for as many cooperatively breeding birds as possible
	+ export as a csv.doc and read into R using: skew <- read.csv("data/skew.csv")
	+ 172 data rows = 172 species (note that there are missing data
	+ column descriptions:\
		A. animal = latin binomial (matches the Jetz et al. nomenclature)\
		B. commonName = English name of each species\
		C. breedingGroup = whether each species is cooperative or non-cooperative\
		D. breedingSystem = type of cooperative group\
		E. groupFormation = whether cooperative groups form as families or non-families\
		F. nGrMixedF = number of broods in multi-female groups with mixed maternity\
		G. nGrMonodF = number of broods in multi-female groups laid by a single female\
		H. perMixedF = percentage of broods in multi-female groups with mixed maternity\
		I. notes = how the data were extracted\
		J. nGrMixedM = number of broods in multi-male groups with mixed paternity\
		K. nGrMonodM = number of broods in multi-male groups sired by a single male\
		L. perMixedM = percentage of broods in multi-male groups with mixed paternity\
		M. notes = how the data were extracted\
		N. source = where the data were extracted from in each study\
		O. reference = study from which data were extracted


- Supplementary Table 4 (specialization)
	+ data on the correlation between group size and: maternal care, maternal fecundity and maternal survival for as many cooperatively breeding birds as possible (see effectSizes.txt for further info.)
	+ export as a csv.doc and read into R using: effects <- read.csv("data/effectSizes.csv")
	+ 114 data rows = 60 species (some species have multiple effect size estimates)
	+ column descriptions:\
		A. commonName = English name of each species\
		B. animal = latin binomial (matches the Jetz et al. nomenclature)\
		C. measurement = aspect of maternal care measured (brooding / feeding / incubation)\
		D. units = the units in which each measurement was made\
		E. zRcare = z-transformed correlation between maternal care and group size\
		F. RmeN = number of groups studied\
		G. varRme = sampling variance for zRcare (1 / (RmeN - 3))\
		H. notes = where the data were extracted from in each study\
		I. measurement = aspect of maternal fecundity measured (clutch size / re-nesting / egg size)\
		J. units = the units in which each measurement was made\
		K. zRfecundity = z-transformed correlation between fecundity and group size\
		L. RmiN = number of groups studied\
		M. varRmi = sampling variance for zRfecundity (1 / (RmiN - 3))\
		N. notes = where the data were extracted from in each study\
		O. measurement = survival probability of breeding females\
		P. units = the units in which each measurement was made\
		Q. zRsurvival = z-transformed correlation between breeding female survival and group size\
		R. RsxN = number of groups studied\
		S. varRsx = sampling variance for zRsurvival (1 / (RsxN - 3))\
		T. notes = where the data were extracted from in each study\
		U. groupFormation = whether cooperative groups form as families or non-families\
		V. care refs = study from which data for zRcare were extracted\
		W. fecundity refs = study from which data for zRfecundity were extracted\
		X. survival refs = study from which data for zRsurvival were extracted


- Supplementary Table 5 (results)
	+ parameter estimates from statistical models (see DowningetalRcode.R for further info.)
	+ 81 rows = output from 9 statistical models
	+ variable names:\
		beta = fixed or random effect parameter estimate from the model\
		lwr CI = lower 95% credible interval from the posterior distribution of the model\
		upr CI = upper 95% credible interval from the posterior distribution of the model


- Supplementary Table 9 (causality)
	+ summary of evidence that helper effects on breeder care and fecundity are causal
	+ 10 rows = 10 different studies examining causality of helper effects
	+ column descriptions:\
		A. commonName = English name of each species\
		B. animal = latin binomial (matches the Jetz et al. nomenclature)\
		C. measure = effect size measured (Zr care / Zr fecundity)\
		D. study type = how the study addressed causality\
		E. description = details of each experiment\
		F. finding = results of each experiment\
		G. conclusion = summary of helper effects\
		H. reference = study in which experiment conducted


- Supplementary Table 10 (stitchbird)
	+ parameter estimates from ancestral polyandry model with stitchbird assigned as non-cooperative
	+ 12 rows = output from one statistical model
	+ variable names:\
		beta = fixed or random effect parameter estimate from the model\
		lwr CI = lower 95% credible interval from the posterior distribution of the model\
		upr CI = upper 95% credible interval from the posterior distribution of the model


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: GF_R_Code.R

This R script contains all R code needed to replicate the analyses (including packages and functions)
Numbering is consistent with the methods section of the manuscript

- PART 1 - MAIN ANALYSES
	+ Part 1.1: ancestral state estimation
	+ Part 1.2 group formation and reproductive division of labour
	+ Part 1.3 group formation and group size
	+ Part 1.4 group formation and reproductive specialization


- PART 2 - SENSITIVITY ANALYSES
	+ Part 2.1 breeding system classification
	+ Part 2.2 different measures of fecundity and maternal care
	+ Part 2.3 uncertainty in ancestral state estimation and different polyandry measures


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: GF_Data_Extraction.txt

This plain text document contains details of how each effect size was calculated, organised by species


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: GF_Supp_Methods.pdf

This pdf document provides details of sensitivity analyses and the causality of helper effects and contains the following supplementary tables:

- Supplementary Table 6 (transitions between Riehl  2013 classifications)
	+ results of the transition analysis between breeding systems using Riehl's classifications

- Supplementary Table 7 (division of reproduction, group size and specialization using Riehl 2013 classifications)
	+ results of analyses of division of reproduction, group size and specialization using Riehl's classifications


- Supplementary Table 8 (uncertainty in ancestral state estimation using EPPbr and EGPbr)
	+ results from analyses of polyandry estimates for transitions from non cooperative (nonCoop) to family groups, non-family groups and non cooperative breeding systems with different levels of uncertainty in the assignment of ancestral breeding systems and using different measures of polyandry
	+ variable names:\
		% EPP br. = percentage of broods with extra-pair chicks\
		% EGP br. = percentage of broods with extra-group chicks

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
