# Long-term selection in dairy ovine (Founders/Burn-in)
# for one trait: Milk Yield
# Simona Antonios
# source('simulation_sheep.r',echo=T)


# Clean
rm(list = ls())

# Load packages
library(package = "tidyverse")
library(package = "AlphaSimR")
library(degreenet)  # for simulation of herd sizes
library(data.table) # for fast data operations (reading, writing, ...)

# Source required functions
sourceFunctions = function() {
  paste0(source('functions_Mix99.R',echo=T))
}
sourceFunctions() # to ensure we have the latest version


# ---- Global parameters ----

nFemalesInLactation = 66600                                                        # no. of Females that are in lactation in the population
nEliteEwes          = 6000
nDamsOfFemales      = 60600
nLact1              = round(0.23 * nEliteEwes) + round(0.3 * nDamsOfFemales)       # no. of Ewes in lactation 1
nLact2              = round(0.35 * nEliteEwes) + 0.28 * round(nDamsOfFemales)      # no. of Ewes in lactation 2
nLact3              = round(0.26 * nEliteEwes) + 0.24 * round(nDamsOfFemales)      # no. of Ewes in lactation 3
nLact4              = round(0.16 * nEliteEwes) + 0.18 * round(nDamsOfFemales)      # no. of Ewes in lactation 4
nLambs              = trunc(((0.9 * 1.4 * 0.75 + 0.6 * 1.6 * 0.75))*nEwes/2)       # no. of Lamb, (1.5 prolificity rate, fertility rate of 0.6, and survival rate is of 0.7))
nYngFemales         = 0.48 * nLambs                                                # no. of successful born Lambs
nYngRams            = 900                                                          # no. of young Rams for PT scheme
nWtRams1            = 150                                                          # no. of waiting Rams for PT scheme
nEliteSires1        = 10                                                           # no. of elite sires selected every year (for the first year of AI)
nEliteSires2        = 10                                                           # no. of elite sires selected every year (for the second year of AI)
nEliteSires3        = 10                                                           # no. of elite sires selected every year (for the third year of AI)
nSiresOfFemales1    = 55                                                           # no. of elite sires of dams selected every year
nSiresOfFemales2    = 55                                                           # no. of elite sires of dams selected every year
nSiresOfFemales3    = 55                                                           # no. of elite sires of dams selected every year
nNaturalMatingRams  = 1000                                                         # no. of NM sires selected every year
startYear           = 1980                                                         # start year of the simulations
nHerds              = 237                                                          # number of herds
meanHerdSize        = 350                                                          # average herd size
sdHerdSize          = 129                                                          # sd of herd size
BaseNe              = 150                                                          # effective population size
nChr                = 26                                                           # no. of chromosomes
ChrSize             = 95 * 10e6                                                    # chr. size
nQTL                = 4500                                                         # no. of QTLs 
nQTLPerChr          = round(nQTL / nChr)                                           # no. of QTLs per chromosome
nSNPPerChr          = 4000                                                         # no. of segregation sites (markers) in the genome
RecRate             = 1.5e-8                                                       # Recombination rate
MutRate             = 1.5e-8                                                       # Mutation rate
nPTyrs              = 20                                                           # no. of progeny test years
nTraits             = 1                                                            # no. of traits (Milk Yield)

# Trait parameters:     
Addh2               = 0.34                                                         # heritability 
meangv              = 0                                                            # Mean genetic variance
addVar              = 1200                                                         # additive genetic variance
permVar             = 500                                                          # permanent environmental variance
resVar              = 1500                                                         # residual variance
phenVar             = 7000                                                         # specifying phenotypic variance for the trait
yearVar             = 1000                                                         # year variance
herdVar             = 1500                                                         # herd variance
herdYearVar         = 1000                                                         # herd/year covariance !!!

meanLac1            = 209                                                          # Mean of milk yield in Lactation 1
meanLac2            = 213                                                          # Mean of milk yield in Lactation 2
meanLac3            = 207                                                          # Mean of milk yield in Lactation 3
meanLac4            = 203                                                          # Mean of milk yield in Lactation 4

nonAddVar           = phenVar - addVar                                             # reduce additive variance from Phenotypic variance for non additive effects
trtMean             = 198                                                          # trait mean value 
trtStDev            = 83                                                           # trait standard deviation value 
trtVar              = 7000                                                         # trait variance


#Trait distribution
# Min.   : 10.05 
# 1st Qu.:135.34 
# Median :187.53 
# Mean   :197.52 
# 3rd Qu.:251.18 
# Max.   :828.39 


AncientNe = tibble(GenerationsAgo = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000, 7000, 8000, 9000, 10000),
                   Ne             = c(150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5500, 6000, 6500, 7000, 7500))

pdf("generations_ago_with_log.pdf")
with(AncientNe, plot(Ne ~ GenerationsAgo, type = "b", main = "Generation ago without log"))
with(AncientNe, plot(Ne ~ GenerationsAgo, type = "b", log = "xy", main = "Generation ago with log"))
dev.off()

# -------------------------------------- Founder Population ---------------------------------------

if (simStart == "FS") {
  runSim = 1 # Force start simulations from scratch
  if (file.exists("image_Year*.RData")) {
    system("rm image_Year*.RData")
  }
} else {
  lastImage = system("ls -ltr image_Year*.RData | tail -n -1", intern = T)
  print(lastImage)
  lastImage = tail(strsplit(lastImage, split=" ")[[1]],1)
  yearImage = as.numeric(gsub(".*?([0-9]+).*", "\\1", lastImage))
  if (yearImage < nPTyrs) { # Continue in PT period
    load(lastImage)
    runSim = 2
    }
  } else {
    print("simulation years are finished. Please increase simulation years if you want to run more years or start from scratch with command FS as the 6th argument")
    stop()
}

if (runSim == 1) {
Tmp = AncientNe %>%
  filter(Ne > BaseNe)
MaCSeNFlags = paste("-eN", Tmp$GenerationsAgo / (4 * BaseNe), Tmp$Ne / BaseNe, collapse = " ")
founderPop = runMacs(nInd = 10 * BaseNe,
                     nChr = 26,
                     segSites = 4000,
                     manualCommand = paste(as.integer(ChrSize),
                                           "-t", MutRate * 4 * BaseNe,
                                           "-r", RecRate * 4 * BaseNe,
                                           MaCSeNFlags),
                     manualGenLen = RecRate * ChrSize)
save.image("founder.RData")

# -------------------------------- AlphaSimR simulation parameters (SP) --------------------------------

SP = SimParam$new(founderPop)
SP$setSexes("yes_rand")
SP$setTrackPed(isTrackPed = TRUE)
SP$addTraitA(nQtlPerChr = nQTLPerChr, mean = trtMean, var = addVar)
# SP$setVarE() we assume that only ewes will have phenotypes, hence do not set this!!! Use setPhenoEwe() instead!
# SP$setCorE() we assume that only ewes will have phenotypes, hence do not set this!!! Use setPhenoEwe() instead!

# ---- Herds and herd and year effects ----

herdSize = simcmp(n = (10 * nHerds),
                  v = c(meanHerdSize, sdHerdSize))
# hist(herdSize)
herdSize = herdSize[herdSize > 9]
herdSize = sample(x = herdSize, size = nHerds)
# hist(herdSize)
herdEffect = cbind(rnorm(n = nHerds, mean = 0, sd = sqrt(herdVar[1])))
herdYearEffect = sampleHerdYearEffect(n = nHerds)
herds = list(
  herd = 1:nHerds,
  herdSize = herdSize,
  herdEffect = herdEffect,
  herdYearEffect = herdYearEffect
)

yearEffect = sampleYearEffect()

# ------------------------------------- Fill-in --------------------------------------

year = 0

# Generate initial founders
basePop = newPop(founderPop)


# Rams

# Elite Rams
# ... elite Rams in the 3rd year in service
eliteSires3 = randCross(pop = basePop, nCrosses = 10 * nEliteSires3)
eliteSires3 = selectWithinFam(pop = eliteSires3, nInd = 1, use = "rand", famType = "M")
eliteSires3@sex[] = "M"
eliteSires3@father[] = "0"
eliteSires3@mother[] = "0"
eliteSires3 = fillInMisc(pop = eliteSires3, year = startYear - 4.5)
# ... elite Rams in the 2nd year in service
eliteSires2 = randCross(pop = basePop, nCrosses = 10 * nEliteSires2)
eliteSires2 = selectWithinFam(pop = eliteSires2, nInd = 1, use = "rand", famType = "M")
eliteSires2@sex[] = "M"
eliteSires2@father[] = "0"
eliteSires2@mother[] = "0"
eliteSires2 = fillInMisc(pop = eliteSires2, year = startYear - 3.5)
# ... elite Rams in the 1st year in service
eliteSires1 = randCross(pop = basePop, nCrosses = 10 * nEliteSires1)
eliteSires1 = selectWithinFam(pop = eliteSires1, nInd = 1, use = "rand", famType = "M")
eliteSires1@sex[] = "M"
eliteSires1@father[] = "0"
eliteSires1@mother[] = "0"
eliteSires1 = fillInMisc(pop = eliteSires1, year = startYear - 2.5)
eliteSires = c(eliteSires3, eliteSires2, eliteSires1) # this is correct, they are the same in fill-in

#Rams of dams (AI)
# ... elite Rams in the 3rd year in service
SiresOfFemales3 = randCross(pop = basePop, nCrosses = 10 * nSiresOfFemales3)
SiresOfFemales3 = selectWithinFam(pop = SiresOfFemales3, nInd = 1, use = "rand", famType = "M")
SiresOfFemales3@sex[] = "M"
SiresOfFemales3@father[] = "0"
SiresOfFemales3@mother[] = "0"
SiresOfFemales3 = fillInMisc(pop = SiresOfFemales3, year = startYear - 4.5)
# ... elite Rams in the 2nd year in service
SiresOfFemales2 = randCross(pop = basePop, nCrosses = 10 * nSiresOfFemales2)
SiresOfFemales2 = selectWithinFam(pop = SiresOfFemales2, nInd = 1, use = "rand", famType = "M")
SiresOfFemales2@sex[] = "M"
SiresOfFemales2@father[] = "0"
SiresOfFemales2@mother[] = "0"
SiresOfFemales2 = fillInMisc(pop = SiresOfFemales2, year = startYear - 3.5)
# ... elite Rams in the 1st year in service
SiresOfFemales1 = randCross(pop = basePop, nCrosses = 10 * nSiresOfFemales1)
SiresOfFemales1 = selectWithinFam(pop = SiresOfFemales1, nInd = 1, use = "rand", famType = "M")
SiresOfFemales1@sex[] = "M"
SiresOfFemales1@father[] = "0"
SiresOfFemales1@mother[] = "0"
SiresOfFemales1 = fillInMisc(pop = SiresOfFemales1, year = startYear - 2.5)
SiresOfFemales = c(SiresOfFemales3, SiresOfFemales2, SiresOfFemales1) # this is correct, they are the same in fill-in


# Waiting Rams
# # ... 1.5 years old 
wtRams1  = randCross(basePop, nCrosses = nWtRams1)
wtRams1@sex[] = "M"
wtRams1@father[] = "0"
wtRams1@mother[] = "0"
wtRams1 = fillInMisc(pop = wtRams1, year = startYear - 1.5)

# Natural Mating Rams, 1.5
ntlMatingRams  = randCross(basePop, nCrosses = nNaturalMatingRams)
ntlMatingRams@sex[] = "M"
ntlMatingRams@father[] = "0"
ntlMatingRams@mother[] = "0"
ntlMatingRams = fillInMisc(pop = ntlMatingRams, year = startYear - 1.5)

# ... 1/12 year old
yngRams  = randCross(basePop, nCrosses = nYngRams)
yngRams@sex[] = "M"
yngRams@father[] = "0"
yngRams@mother[] = "0"
yngRams = fillInMisc(pop = yngRams, year = startYear -1/12)



# ewes
# Elite
# ... 4th lactation
EliteEwesLact4 = randCross(basePop, nCrosses = round(0.16 * nEliteEwes))
EliteEwesLact4@sex[] = "F"
EliteEwesLact4@father[] = "0"
EliteEwesLact4@mother[] = "0"
EliteEwesLact4 = fillInMisc(pop = EliteEwesLact4, herds = herds, permEnvVar = permVar,
                           year = startYear - 4)
# ... 3rd lactation
EliteEwesLact3 = randCross(basePop, nCrosses =  round(0.26 * nEliteEwes))
EliteEwesLact3@sex[] = "F"
EliteEwesLact3@father[] = "0"
EliteEwesLact3@mother[] = "0"
EliteEwesLact3 = fillInMisc(pop = EliteEwesLact3, herds = herds, permEnvVar = permVar,
                           year = startYear - 3)
# ... 2nd lactation
EliteEwesLact2 = randCross(basePop, nCrosses =  round(0.35 * nEliteEwes))
EliteEwesLact2@sex[] = "F"
EliteEwesLact2@father[] = "0"
EliteEwesLact2@mother[] = "0"
EliteEwesLact2 = fillInMisc(pop = EliteEwesLact2, herds = herds, permEnvVar = permVar,
                           year = startYear - 2)
# ... 1st lactation
EliteEwesLact1 = randCross(basePop, nCrosses =  round(0.23 * nEliteEwes))
EliteEwesLact1@sex[] = "F"
EliteEwesLact1@father[] = "0"
EliteEwesLact1@mother[] = "0"
EliteEwesLact1 = fillInMisc(pop = EliteEwesLact1, herds = herds, permEnvVar = permVar,
                           year = startYear - 1)

# Dams of females
# ... 4th lactation
DamofFemalesLact4 = randCross(basePop, nCrosses = round(0.18 * nDamsOfFemales))
DamofFemalesLact4@sex[] = "F"
DamofFemalesLact4@father[] = "0"
DamofFemalesLact4@mother[] = "0"
DamofFemalesLact4 = fillInMisc(pop = DamofFemalesLact4, herds = herds, permEnvVar = permVar,
                       year = startYear - 4)
# ... 3rd lactation
DamofFemalesLact3 = randCross(basePop, nCrosses = round(0.24 * nDamsOfFemales))
DamofFemalesLact3@sex[] = "F"
DamofFemalesLact3@father[] = "0"
DamofFemalesLact3@mother[] = "0"
DamofFemalesLact3 = fillInMisc(pop = DamofFemalesLact3, herds = herds, permEnvVar = permVar,
                       year = startYear - 3)
# ... 2nd lactation
DamofFemalesLact2 = randCross(basePop, nCrosses =  round(0.28 * nDamsOfFemales))
DamofFemalesLact2@sex[] = "F"
DamofFemalesLact2@father[] = "0"
DamofFemalesLact2@mother[] = "0"
DamofFemalesLact2 = fillInMisc(pop = DamofFemalesLact2, herds = herds, permEnvVar = permVar,
                       year = startYear - 2)
# ... 1st lactation
DamofFemalesLact1 = randCross(basePop, nCrosses =  round(0.3 * nDamsOfFemales))
DamofFemalesLact1@sex[] = "F"
DamofFemalesLact1@father[] = "0"
DamofFemalesLact1@mother[] = "0"
DamofFemalesLact1 = fillInMisc(pop = DamofFemalesLact1, herds = herds, permEnvVar = permVar,
                       year = startYear - 1)
                       
ewesLact1 = c(DamofFemalesLact1, EliteEwesLact1)
ewesLact2 = c(DamofFemalesLact2, EliteEwesLact2)
ewesLact3 = c(DamofFemalesLact3, EliteEwesLact3)
ewesLact4 = c(DamofFemalesLact4, EliteEwesLact4)

# young ewes
yngFemales = randCross(basePop, nCrosses = nYngFemales)
yngFemales@sex[] = "F"
yngFemales@father[] = "0"
yngFemales@mother[] = "0"
yngFemales = fillInMisc(pop = yngFemales, herds = herds, permEnvVar = permVar,
                     year = startYear - 1/12)

# new born
## From AI
n1 = round(400 * (nEliteSires1+nEliteSires2+nEliteSires3) * 0.75 * 0.6 * 1.6)
n2 = round(150 * 140 * 0.75 * 0.6 * 1.6)
n3 = round((nEwes/2 -(400 * (nEliteSires1+nEliteSires2+nEliteSires3) + 150 * 140)) * 0.75 * 0.6 * 1.6)
## From Natural mating
n4= round(1000 * 40 * 0.9 * 0.75 * 1.4)   
ramsId = c(sample(eliteSires@id, size = n1, replace = TRUE),
            sample(wtRams1@id, size = n2, replace = TRUE),
            sample(SiresOfFemales@id,size = n3, replace = TRUE),
            sample(ntlMatingRams@id,size = n4, replace = TRUE))

matingPlan = cbind(c(ewesLact1@id, ewesLact2@id, ewesLact3@id, ewesLact4@id), ramsId)
lambs = makeCross2(females = c(ewesLact1, ewesLact2, ewesLact3, ewesLact4),
                    males = c(eliteSires, wtRams1, SiresOfFemales, ntlMatingRams),
                    crossPlan = matingPlan)
lambs = selectInd(pop = lambs, nInd = nLambs, use = "rand")
lambs = fillInMisc(pop = lambs, mothers = c(ewesLact1, ewesLact2, ewesLact3, ewesLact4),
                    permEnvVar = permVar, year = startYear)

# Heterozygosity
# 
# lambsHet = calcHeterozygosity(pop = selectInd(lambs, nInd = ceiling(0.1 * nEwes), use = "rand"))
# ramsHet = calcHeterozygosity(pop = eliteSires)
# heterozygosity =
#   data.frame(year = startYear,
#              lambsHet = lambsHet$Heterozygosity, lambsHom = lambsHet$Homozygosity,
#              ramsHet = ramsHet$Heterozygosity, ramsHom = ramsHet$Homozygosity)
# write.table(x = heterozygosity, file = "heterozygosity.txt",
#             col.names = TRUE, row.names = FALSE, quote = FALSE)
# 

database = recordData( pop = eliteSires, year = startYear)
database = recordData( pop = SiresOfFemales, year = startYear)
database = recordData(database, pop = wtRams1, year = startYear)
database = recordData(database, pop = ntlMatingRams, year = startYear)
database = recordData(database, pop = yngRams, year = startYear)
database = recordData(database, pop = yngFemales, year = startYear)
database = recordData( pop = ewesLact4, year = startYear, lactation = 4)
database = recordData( pop = ewesLact3, year = startYear, lactation = 3)
database = recordData( pop = ewesLact2, year = startYear, lactation = 2)
database = recordData( pop = ewesLact1, year = startYear, lactation = 1)
database = recordData(database, pop = lambs, year = startYear)

save(x = database, file = "database.RData")

# Save image

save.image(file = paste("image_Year0_FillIn_simona.RData", sep = ""))

}

if (runSim == 1 | runSim == 2) {
  if (runSim == 1) {
    yearToDo = 1
  }
  if (runSim == 2) {
    yearToDo = year + 1
  }
# ---- Progeny testing stage ----

for (year in yearToDo:nPTyrs) {
  yearFull = startYear + year
  print(paste("Working on Progeny testing stage year", year,
              Sys.time(), "...", sep = " "))
  
  yearEffect = sampleYearEffect(n = 1)
  herds$herdYearEffect = sampleHerdYearEffect(n = nHerds)
  
  # ---- Phenotyping ----
  
  
  ewesLact4 = setPhenoEwe(ewesLact4, varE = resVar,
                          yearEffect = yearEffect, herds = herds,
                          traitMask = createTraitMask(ewesLact4))
  
  ewesLact3 = setPhenoEwe(ewesLact3, varE = resVar,
                          yearEffect = yearEffect, herds = herds,
                          traitMask = createTraitMask(ewesLact3))
  
  ewesLact2 = setPhenoEwe(ewesLact2, varE = resVar,
                          yearEffect = yearEffect, herds = herds,
                          traitMask = createTraitMask(ewesLact2))
  
  ewesLact1 = setPhenoEwe(ewesLact1, varE = resVar,
                          yearEffect = yearEffect, herds = herds,
                          traitMask = createTraitMask(ewesLact1))
  
  database = setDatabasePheno(database, pop = ewesLact1, trait = 1)
  database = setDatabasePheno(database, pop = ewesLact2, trait = 1)
  database = setDatabasePheno(database, pop = ewesLact3, trait = 1)
  database = setDatabasePheno(database, pop = ewesLact4, trait = 1)
  
  # ---- SELECTION BY CATEGORIES ----
  
  # ---- Rams ----
  # ---- Elite sires ----
  eliteSires3 = eliteSires2 # eliteSires3 are 4.5 years old here
  eliteSires2 = eliteSires1 # eliteSires2 are 3.5 years old here
  if (year == 1) {
    use = "rand"
  } else {
    use= "ebv"
  }
  # TODO: this if (!all... is needed because we set fathers to 0 - related to
  if (all(wtRams1@father == "0")) {
    eliteSires1 = selectInd(pop = wtRams1, nInd = nEliteSires1, # eliteSires1 are 2.5 years old here
                            use = use, trait = 1)
  } else {
    eliteSires1 = selectWithinFam(pop = wtRams1, nInd = 1, # eliteSires1 are 2.5 years old here
                                  use = use, trait = 1,
                                  famType = "M")
    eliteSires1 = selectInd(pop = eliteSires1, nInd = nEliteSires1,
                            use = use, trait = 1)
  }
  
  # ---- Sires of Dams ----
  SiresOfFemales3 = SiresOfFemales2 # SiresOfFemales3 are 4.5 years old here
  SiresOfFemales2 = SiresOfFemales1 # SiresOfFemales2 are 3.5 years old here
  if (all(wtRams1@father == "0")) {
    SiresOfFemales1 = selectInd(pop = wtRams1, nInd = nSiresOfFemales1, # SiresOfFemales1 are 2.5 years old here
                                use = use, trait = 1)
  } else {
    SiresOfFemales1 = selectWithinFam(pop = wtRams1, nInd = 1, # SiresOfFemales1 are 2.5 years old here
                                      use = use, trait = 1,
                                      famType = "M")
    SiresOfFemales1 = selectInd(pop = SiresOfFemales1, nInd = nSiresOfFemales1,
                                use = use, trait = 1)
  }
  
  # We need this sel and n before we update eliteSires
  # Note that we select yngRams only from lambs from proven rams (not waiting rams)
  sel = lambs@father %in% eliteSires@id
  sel1 = lambs@father %in% c(eliteSires@id,SiresOfFemales@id) 
  n = ceiling(nYngRams / length(eliteSires@id))
  eliteSires = c(eliteSires3, eliteSires2, eliteSires1)
  wtRams1 = selectInd(pop = yngRams, nInd = nWtRams1 ,
                          use = use, trait = 1) # wtRams1 are 1.5 years old here
  yngRams = selectWithinFam(pop = lambs[sel], nInd = n, # yngRams are 0.12 year old here
                             use = use, trait = 1,
                             sex = "M", famType = "M")
  yngRams = selectInd(pop = yngRams, nInd = nYngRams,
                       use = use, trait = 1)
  
  ntlMatingRams =  selectInd(pop = lambs[sel1], nInd = nNaturalMatingRams ,
                             use = use, trait = 1, sex = "M", famType = "M")
  # ---- ewes ----
  # ---Elites---
  
  EliteEwesLact4 = selectInd(ewesLact3, nInd = round(0.16 * nEliteEwes), use = "ebv") # EliteEwesLact4 are 5 years old here
  EliteEwesLact3 = selectInd(ewesLact2, nInd = round(0.26 * nEliteEwes), use = "ebv") # EliteEwesLact3 are 4 years old here
  EliteEwesLact2 = selectInd(ewesLact1, nInd = round(0.35 * nEliteEwes), use = "ebv") # EliteEwesLact2 are 3 years old here
  EliteEwesLact1 = selectInd(lambs[sel1], nInd = round(0.23 * nEliteEwes), use = "rand",sex ="F", trait = 1, famType = "M")   # EliteEwesLact1 are 2 years old here

  # Set phenotypes of the last lactation from the last generation to missing.
  # These are only copied phenotypes from the latest lactation.
  # Correct phenotypes will be added at the beginning of the next year.
  # Phenotypes of the last lactation are not part of the evaluation because lambs are selected before the
  # phenotype is recorded.
  EliteEwesLact4@pheno[,c(1:nTraits)] = NA
  EliteEwesLact3@pheno[,c(1:nTraits)] = NA
  EliteEwesLact2@pheno[,c(1:nTraits)] = NA
  EliteEwesLact1@pheno[,c(1:nTraits)] = NA
  
  # ---Dams of Females---
  DamofFemalesLact4 = selectInd(DamofFemalesLact3, round(0.18 * nDamsOfFemales), use = "rand") # DamofFemalesLact4 are 5 years old here
  DamofFemalesLact3 = selectInd(DamofFemalesLact2, round(0.24 * nDamsOfFemales), use = "rand") # DamofFemalesLact3 are 4 years old here
  DamofFemalesLact2 = selectInd(DamofFemalesLact1, round(0.28 * nDamsOfFemales), use = "rand") # DamofFemalesLact2 are 3 years old here
  DamofFemalesLact1 = selectInd(yngFemales, nInd = round(0.3 * nDamsOfFemales), use = "rand") # DamofFemalesLact1 are 2 years old here
  
  yngFemales = selectInd(lambs, nInd = nYngFemales, use = "rand", sex = "F") 
  
  # Set phenotypes of the last lactation from the last generation to missing.
  # These are only copied phenotypes from the latest lactation.
  # Correct phenotypes will be added at the beginning of the next year.
  # Phenotypes of the last lactation are not part of the evaluation because calves are selected before the
  # phenotype is recorded.
  DamofFemalesLact4@pheno[,c(1:nTraits)] = NA
  DamofFemalesLact3@pheno[,c(1:nTraits)] = NA
  DamofFemalesLact2@pheno[,c(1:nTraits)] = NA
  DamofFemalesLact1@pheno[,c(1:nTraits)] = NA
  
  # ---- lambs ----
  
  # See above comments on stage of animals and their use in the fill-in stage
  n1 = round(400 * (nEliteSires1+nEliteSires2+nEliteSires3) * 0.75 * 0.6 * 1.6)
  n2 = round(150 * 140 * 0.75 * 0.6 * 1.6)
  n3 = round((nEwes/2 -(400 * (nEliteSires1+nEliteSires2+nEliteSires3) + 150 * 140)) * 0.75 * 0.6 * 1.6)
  ## From Natural mating
  n4= round(1000 * 40 * 0.9 * 0.75 * 1.4)
  ramsId = c(sample(eliteSires@id, size = n1, replace = TRUE),
             sample(wtRams1@id, size = n2, replace = TRUE),
             sample(SiresOfFemales@id,size = n3, replace = TRUE),
             sample(ntlMatingRams@id,size = n4, replace = TRUE))
  matingPlan = cbind(c(ewesLact1@id, ewesLact2@id, ewesLact3@id, ewesLact4@id), ramsId)
  lambs = makeCross2(females = c(ewesLact1, ewesLact2, ewesLact3, ewesLact4),
                     males = c(eliteSires, wtRams1, SiresOfFemales, ntlMatingRams),
                     crossPlan = matingPlan)
  lambs = selectInd(pop = lambs, nInd = nLambs, use = "rand") # lambs are 0 years old here
  lambs = fillInMisc(pop = lambs, mothers = c(ewesLact1, ewesLact2, ewesLact3, ewesLact4),
                     permEnvVar = permVar, year = yearFull)

  # ---- Heterozygosity ----
  # 
  # lambsHet = calcHeterozygosity(pop = selectInd(lambs, nInd = ceiling(0.1 * nEwes), use = "rand"))
  # ramsHet = calcHeterozygosity(pop = eliteSires)
  # heterozygosity = rbind(
  #   heterozygosity,
  #   data.frame(year = yearFull,
  #              lambsHet = lambsHet$Heterozygosity, lambsHom = lambsHet$Homozygosity,
  #              ramsHet = ramsHet$Heterozygosity, ramsHom = ramsHet$Homozygosity)
  # )
  # write.table(x = heterozygosity, file = "heterozygosity.txt",
  #             col.names = TRUE, row.names = FALSE, quote = FALSE)
  # 

  # Data recording
  database = recordData( pop = eliteSires, year = yearFull)
  database = recordData( pop = SiresOfFemales, year = yearFull)
  database = recordData(database, pop = wtRams1, year = yearFull)
  database = recordData(database, pop = ntlMatingRams, year = yearFull)
  database = recordData(database, pop = yngRams, year = yearFull)
  database = recordData(database, pop = yngFemales, year = yearFull)
  database = recordData(database, pop = ewesLact4, year = yearFull, lactation = 4)
  database = recordData(database, pop = ewesLact3, year = yearFull, lactation = 3)
  database = recordData(database, pop = ewesLact2, year = yearFull, lactation = 2)
  database = recordData(database, pop = ewesLact1, year = yearFull, lactation = 1)
  database = recordData(database, pop = lambs, year = yearFull)
  
  save(x = database, file = "database.RData")
}
}
