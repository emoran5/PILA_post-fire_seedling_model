# PILA_post-fire_seedling_model

This repository contains data and code associated with "Sugar pine association genetics and performance in a post-fire restoration planting", by Emily V. Moran, Rainbow DeSilva, Jessica Wright, and Courtney Canning.

The seedling data files include:
Peavine_PILA_dec2_trim.txt   (2019 planting year/site)
Stumpy2018_PILAdec2_trim.txt (2018 planting year/site)
Stumpy2017_PILAdec2_trim.txt (2017 planting year/site)

Planting site climate data files include:
CWD_91_21.txt (annual Climatic Water Deficit)          Sno_91_21.txt (April snowpack by year) 
Precip_91_21.txt (annual precipitation)                JMin_91_21.txt (January minimum temperature by year)
JMax_91_21.txt  (July maximum temperature by year)

Seedling genotype data files include:
favor_allele_Sno2.txt         
favor_allele_Tmax2.txt        
favor_allele_Tmn2.txt
favor_allele_Cwd2.txt         
favor_allele_Ppt2.txt
These files show the genotypes for each genotype at loci associated with April snowpack, maximum July temperature, minimum January temperature, CWD, and annual precipitation, respectively. In these files, 2 = "has two of the allele found more frequently in locations where this climate variable is similar to the planting site (within * SD), AKA the 'locally favored allele'", 1 = "has one copy of the 'locally favored allele'", 0 = "has no copies of 'locally favored allele'", and 9 = "ungenotyped".



The R code files 'PILA_SeedlingDataOrg.r' and 'PILA_SeedlingDataOrg_Geno.r' organize and clean this data, preparing it as input for the statistical models. 
The first prepares the data for all seedlings (but not the genotype data) for the all-seedling source-elevation-only growth and survival models. It also makes the graphs of survival and growth by elevation and site, as seen in the paper. The second prepares the data for the genotyped-seedling-only models that may include either source elevation or genotype or both.



Only some examples of the tested model files are included, but the others can be reconstructed based on these examples.
'S8_PILA.r' is the survival model for all seedlings that includes Intercept, Site, Year/age, Height, Source elevation x Site interactions, and a climate variable.
'G7_PILA.r' is the growth model for all seedlings that includes Intercept, Site, Age, Height, Source elevation, and a climate variable.

'S9_Geno_PILA.r' is the survival model for genotyped seedlings that includes Intercept, Site, Year/age, Height, Genotype index (FG1 or FG2), and a climate variable.
'S7_Geno_PILA.r' is the survival model for genotyped seedlings that includes Intercept, Site, Year/age, Height, Source elevation x Site interactions, and a climate variable.

'G9_Geno_PILA.r' is the growth model for genotyped seedlings that includes Intercept, Site, Year/age, Height, Genotype index (FG1 or FG2), and a climate variable.
'G7_Geno_PILA.r' is the growth model for genotyped seedlings that includes Intercept, Site, Year/age, Height, Seed Zone-source elevation combos, and a climate variable.
