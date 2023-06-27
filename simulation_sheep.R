
# Long-term selection in dairy ovine (Founders/Burn-in)
# for one trait: Milk Yield
# Simona Antonios
# source('simulation_sheep.r', echo = TRUE)

# Clean
rm(list = ls())
simStart = "FS"
# Load packages
# install.packages(pkg = c("tidyverse", "AlphaSimR", "degreenet", "data.table"))
library(tidyverse)
library(AlphaSimR)
library(degreenet)  # for simulation of herd sizes
library(data.table) # for fast data operations (reading, writing, ...)

# Source required functions
sourceFunctions = function() {
  paste0(source('functions_burnin.R',echo=T))
}
sourceFunctions() # to ensure we have the latest version

mainDir = getwd()
args = commandArgs(trailingOnly = TRUE)
rep = args[1]
dire = paste("/home/santonios/work/santonios/part2/simulation/rep", rep, sep = "")
unlink(dire, recursive = TRUE)
dir.create(path = dire, recursive = TRUE, showWarnings = FALSE)
setwd(dir = dire)

# ---- Global parameters ----

nEwes               = 80000
nFemalesInLactation = 66600                                                        # no. of Females that are in lactation in the population
nEliteEwes          = 6000
nDamsOfFemales      = 60600
nLact1              = round(0.23 * nEliteEwes) + round(0.3 * nDamsOfFemales)       # no. of Ewes in lactation 1
nLact2              = round(0.35 * nEliteEwes) + 0.28 * round(nDamsOfFemales)      # no. of Ewes in lactation 2
nLact3              = round(0.26 * nEliteEwes) + 0.24 * round(nDamsOfFemales)      # no. of Ewes in lactation 3
nLact4              = round(0.16 * nEliteEwes) + 0.18 * round(nDamsOfFemales)      # no. of Ewes in lactation 4
survRate            = 0.75                                                         # survival rate
NM_Fertility        = 0.9                                                          # fertility in NM
NM_prolificacy      = 1.4                                                          # prolificacy in NM
AI_Fertility        = 0.6                                                          # fertility in AI
AI_prolificacy      = 1.6                                                          # prolificacy in NM
nLambs              = trunc(((0.9 * 1.4 * 0.75 + 0.6 * 1.6 * 0.75))*nEwes/2)       # no. of Lamb, (1.5 prolificity rate, fertility rate of 0.6, and survival rate is of 0.7))
nYngFemales         = 0.48 * nLambs                                                # no. of successful born Lambs
nYngRams            = 900                                                          # no. of young Rams for PT scheme
nWtRams1            = 150                                                          # no. of waiting Rams for PT scheme
nWtRams2            = 150                                                          # no. of waiting Rams for PT scheme
nEliteSires1        = 10                                                           # no. of elite sires selected every year (for the first year of AI)
nEliteSires2        = 10                                                           # no. of elite sires selected every year (for the second year of AI)
nEliteSires3        = 10                                                           # no. of elite sires selected every year (for the third year of AI)
nsiresOfFemales1    = 55                                                           # no. of elite sires of dams selected every year
nsiresOfFemales2    = 55                                                           # no. of elite sires of dams selected every year
nsiresOfFemales3    = 55                                                           # no. of elite sires of dams selected every year
nNaturalMatingRams  = 1000                                                         # no. of NM sires selected every year
nEliteSireDose      = 400                                                          # no. of AI doses per elite sire
nwtRamsAIDose       = 85                                                           # no. of AI doses per wating ram
nNtlMatingDose      = 40                                                           # no. of Natural matings per NM ram
nDamsOfFemalesLact1 = round(0.30 * nDamsOfFemales)                                 # no. of dams of females in Lac1
nDamsOfFemalesLact2 = round(0.28 * nDamsOfFemales)                                 # no. of dams of females in Lac2
nDamsOfFemalesLact3 = round(0.24 * nDamsOfFemales)                                 # no. of dams of females in Lac3
nDamsOfFemalesLact4 = round(0.18 * nDamsOfFemales)                                 # no. of dams of females in Lac4
nEliteEwesLact1     = round(0.23 * nEliteEwes)                                     # no. of elite ewes in Lac1
nEliteEwesLact2     = round(0.35 * nEliteEwes)                                     # no. of elite ewes in Lac2
nEliteEwesLact3     = round(0.26 * nEliteEwes)                                     # no. of elite ewes in Lac3
nEliteEwesLact4     = round(0.16 * nEliteEwes)                                     # no. of elite ewes in Lac4
startYear           = 1980                                                         # start year of the simulations
nHerds              = 237                                                          # number of herds
meanHerdSize        = 320                                                          # average herd size
sdHerdSize          = 129                                                          # sd of herd size
BaseNe              = 150                                                          # effective population size
nChr                = 26                                                           # no. of chromosomes
ChrSize             = 95 * 1e6                                                     # chr. size
nQTL                = 4500                                                         # no. of QTLs
nQTLPerChr          = round(nQTL / nChr)                                           # no. of QTLs per chromosome
nSNPPerChr          = 4000                                                         # no. of segregation sites (markers) in the genome
RecRate             = 1.3e-8                                                       # Recombination rate
MutRate             = 1.5e-8                                                       # Mutation rate
nPTyrs              = 20                                                           # no. of progeny test years
nTraits             = 1                                                            # no. of traits (Milk Yield)

# Trait parameters:
Addh2               = 0.34                                                         # heritability
addVar              = 1000                                                         # additive genetic variance
permVar             = 400       # 500 (Old value) - 100 = 400 (new value)          # permanent environmental variance (old is what we started with, 100 was targeted domVar)
# addVar            = 1200                                                         # additive genetic variance
# permVar           = 500                                                          # permanent environmental variance
resVar              = 1500                                                         # residual variance
phenVar             = 7000                                                         # specifying phenotypic variance for the trait
yearVar             = 1000                                                         # year variance
herdVar             = 1500                                                         # herd variance
herdYearVar         = 1000                                                         # herd/year covariance !!!
domVar              = addVar * 0.1
meanDD              = 0.08
varDD               = 0.30

meanLac1            = 209                                                          # Mean of milk yield in Lactation 1
meanLac2            = 213                                                          # Mean of milk yield in Lactation 2
meanLac3            = 207                                                          # Mean of milk yield in Lactation 3
meanLac4            = 180                                                          # Mean of milk yield in Lactation 4

nonAddVar           = phenVar - addVar                                             # reduce additive variance from Phenotypic variance for non additive effects
trtMean             = 200                                                          # trait mean value
trtStDev            = 83                                                           # trait standard deviation value
trtVar              = 7000                                                         # trait variance

#  for AI
n1 = round(nEliteSireDose * (nEliteSires1+nEliteSires2+nEliteSires3) * survRate * AI_Fertility * AI_prolificacy) # 8640 progeny
n2 = round((nEwes/2 -(nEliteSireDose * (nEliteSires1+nEliteSires2+nEliteSires3) + nWtRams1 * nwtRamsAIDose)) * survRate * AI_Fertility * AI_prolificacy) # 10980 progeny,
n3 = round(nWtRams1 * nwtRamsAIDose * survRate * AI_Fertility * AI_prolificacy) # 9180 progeny
# From Natural mating
n4= round(nNaturalMatingRams * nNtlMatingDose * survRate * NM_Fertility * NM_prolificacy) # 37800 progeny

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
pdf("herd_density_PT_20_years.pdf")
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
  # founderPop = runMacs(nInd = 10 * BaseNe,
  #                      nChr = 26,
  #                      segSites = 4000,
  #                      manualCommand = paste(as.integer(ChrSize),
  #                                            "-t", MutRate * 4 * BaseNe,
  #                                            "-r", RecRate * 4 * BaseNe,
  #                                            MaCSeNFlags),
  #                      manualGenLen = RecRate * ChrSize)
  # save.image("founder.RData")

  founderPop = runMacs2(nInd = 10 * BaseNe,
                        nChr = 26,
                        segSites = 4000,
                        Ne = 150,
                        bp = ChrSize,
                        genLen =  RecRate * ChrSize *26 ,
                        mutRate = MutRate,
                        histNe = c(200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5500, 6000, 6500, 7000, 7500),
                        histGen = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000, 7000, 8000, 9000, 10000),
                       # returnCommand = TRUE
  )
  save.image("founder_runMacs2.RData")
  load("founder_runMacs2.RData")

  # -------------------------------- AlphaSimR simulation parameters (SP) --------------------------------

  SP = SimParam$new(founderPop)
  SP$setSexes("yes_rand") #yes_sys half and half
  SP$setTrackPed(isTrackPed = TRUE)

    # Fully additive trait
  # SP$addTraitA(nQtlPerChr = nQTLPerChr, mean = trtMean, var = addVar)
  # Trait with additive and dominance effects - we got the meanDD and varDD from this function
  # altAddTraitAD(nQtlPerChr = nQtlPerChr,
  #               mean = trtMean,
  #               varA = addVar,
  #               varD = domVar,
  #               limMeanDD = c(-1, 2),
  #               limVarDD = c(0, 2),
  #               inbrDepr = 70)
  # New trait called Trait1 was added
  # Dominance variance is 100.0
  # Inbreeding depression is 70.0
  # Used meanDD equals 0.08
  # Used varDD equals 0.28

  SP$addTraitAD(nQtlPerChr = nQTLPerChr, mean = trtMean, var = addVar, meanDD = meanDD, varDD = varDD)
  # SP$setVarE() we assume that only ewes will have phenotypes, hence do not set this!!! Use setPhenoEwe() instead!
  # SP$setCorE() we assume that only ewes will have phenotypes, hence do not set this!!! Use setPhenoEwe() instead!

  # ---- Herds and herd and year effects ----

  herdSize = simcmp(n = (10 * nHerds),
                    v = c(meanHerdSize, sdHerdSize)) # TODO: https://github.com/SimonaAntonios/sheep_simulation_AlphaSimR/issues/1
  # hist(herdSize)
  herdSize = herdSize[herdSize > 100]
  herdSize = sample(x = herdSize, size = nHerds)
  herdEffect = cbind(rnorm(n = nHerds, mean = 0,
                           sd = sqrt(herdVar[1])))
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
  eliteSires3 = fillInMisc(pop = eliteSires3, year = startYear - 5.5)
  # ... elite Rams in the 2nd year in service
  eliteSires2 = randCross(pop = basePop, nCrosses = 10 * nEliteSires2)
  eliteSires2 = selectWithinFam(pop = eliteSires2, nInd = 1, use = "rand", famType = "M")
  eliteSires2@sex[] = "M"
  eliteSires2@father[] = "0"
  eliteSires2@mother[] = "0"
  eliteSires2 = fillInMisc(pop = eliteSires2, year = startYear - 4.5)
  # ... elite Rams in the 1st year in service
  eliteSires1 = randCross(pop = basePop, nCrosses = 10 * nEliteSires1)
  eliteSires1 = selectWithinFam(pop = eliteSires1, nInd = 1, use = "rand", famType = "M")
  eliteSires1@sex[] = "M"
  eliteSires1@father[] = "0"
  eliteSires1@mother[] = "0"
  eliteSires1 = fillInMisc(pop = eliteSires1, year = startYear - 3.5)
  eliteSires = c(eliteSires3, eliteSires2, eliteSires1) # this is correct, they are the same in fill-in

  # Rams of dams (AI)
  # ... elite Rams in the 3rd year in service
  siresOfFemales3 = randCross(pop = basePop, nCrosses = 10 * nsiresOfFemales3)
  siresOfFemales3 = selectWithinFam(pop = siresOfFemales3, nInd = 1, use = "rand", famType = "M")
  siresOfFemales3@sex[] = "M"
  siresOfFemales3@father[] = "0"
  siresOfFemales3@mother[] = "0"
  siresOfFemales3 = fillInMisc(pop = siresOfFemales3, year = startYear - 5.5)
  # ... elite Rams in the 2nd year in service
  siresOfFemales2 = randCross(pop = basePop, nCrosses = 10 * nsiresOfFemales2)
  siresOfFemales2 = selectWithinFam(pop = siresOfFemales2, nInd = 1, use = "rand", famType = "M")
  siresOfFemales2@sex[] = "M"
  siresOfFemales2@father[] = "0"
  siresOfFemales2@mother[] = "0"
  siresOfFemales2 = fillInMisc(pop = siresOfFemales2, year = startYear - 4.5)
  # ... elite Rams in the 1st year in service
  siresOfFemales1 = randCross(pop = basePop, nCrosses = 10 * nsiresOfFemales1)
  siresOfFemales1 = selectWithinFam(pop = siresOfFemales1, nInd = 1, use = "rand", famType = "M")
  siresOfFemales1@sex[] = "M"
  siresOfFemales1@father[] = "0"
  siresOfFemales1@mother[] = "0"
  siresOfFemales1 = fillInMisc(pop = siresOfFemales1, year = startYear - 3.5)
  siresOfFemales = c(siresOfFemales3, siresOfFemales2, siresOfFemales1) # this is correct, they are the same in fill-in

  # Waiting Rams
  # ... 0.5 years old
  wtRams1  = randCross(basePop, nCrosses = nWtRams1)
  wtRams1@sex[] = "M"
  wtRams1@father[] = "0"
  wtRams1@mother[] = "0"
  wtRams1 = fillInMisc(pop = wtRams1, year = startYear - 0.5)
  # ... 1.5 years old
  wtRams2  = randCross(basePop, nCrosses = nWtRams2)
  wtRams2@sex[] = "M"
  wtRams2@father[] = "0"
  wtRams2@mother[] = "0"
  wtRams2 = fillInMisc(pop = wtRams2, year = startYear - 1.5)

  # Natural Mating Rams, 1.5
  ntlMatingRams  = randCross(basePop, nCrosses = nNaturalMatingRams)
  ntlMatingRams@sex[] = "M"
  ntlMatingRams@father[] = "0"
  ntlMatingRams@mother[] = "0"
  ntlMatingRams = fillInMisc(pop = ntlMatingRams, year = startYear - 1.5)

  # # ... 1/12 year old
  # yngRams  = randCross(basePop, nCrosses = nYngRams)
  # yngRams@sex[] = "M"
  # yngRams@father[] = "0"
  # yngRams@mother[] = "0"
  # yngRams = fillInMisc(pop = yngRams, year = startYear -1/12)

  # ewes
  # Elite
  # ... 4th lactation
  eliteEwesLact4 = randCross(basePop, nCrosses = nEliteEwesLact4)
  eliteEwesLact4@sex[] = "F"
  eliteEwesLact4@father[] = "0"
  eliteEwesLact4@mother[] = "0"
  eliteEwesLact4 = fillInMisc(pop = eliteEwesLact4, herds = herds, permEnvVar = permVar,
                              year = startYear - 4)
  # ... 3rd lactation
  eliteEwesLact3 = randCross(basePop, nCrosses = nEliteEwesLact3)
  eliteEwesLact3@sex[] = "F"
  eliteEwesLact3@father[] = "0"
  eliteEwesLact3@mother[] = "0"
  eliteEwesLact3 = fillInMisc(pop = eliteEwesLact3, herds = herds, permEnvVar = permVar,
                              year = startYear - 3)
  # ... 2nd lactation
  eliteEwesLact2 = randCross(basePop, nCrosses = nEliteEwesLact2)
  eliteEwesLact2@sex[] = "F"
  eliteEwesLact2@father[] = "0"
  eliteEwesLact2@mother[] = "0"
  eliteEwesLact2 = fillInMisc(pop = eliteEwesLact2, herds = herds, permEnvVar = permVar,
                              year = startYear - 2)
  # ... 1st lactation
  eliteEwesLact1 = randCross(basePop, nCrosses =  nEliteEwesLact1)
  eliteEwesLact1@sex[] = "F"
  eliteEwesLact1@father[] = "0"
  eliteEwesLact1@mother[] = "0"
  eliteEwesLact1 = fillInMisc(pop = eliteEwesLact1, herds = herds, permEnvVar = permVar,
                              year = startYear - 1)

  # Dams of females
  # ... 4th lactation
  damOfFemalesLact4 = randCross(basePop, nCrosses = nDamsOfFemalesLact4)
  damOfFemalesLact4@sex[] = "F"
  damOfFemalesLact4@father[] = "0"
  damOfFemalesLact4@mother[] = "0"
  damOfFemalesLact4 = fillInMisc(pop = damOfFemalesLact4, herds = herds, permEnvVar = permVar,
                                 year = startYear - 4)
  # ... 3rd lactation
  damOfFemalesLact3 = randCross(basePop, nCrosses = nDamsOfFemalesLact3)
  damOfFemalesLact3@sex[] = "F"
  damOfFemalesLact3@father[] = "0"
  damOfFemalesLact3@mother[] = "0"
  damOfFemalesLact3 = fillInMisc(pop = damOfFemalesLact3, herds = herds, permEnvVar = permVar,
                                 year = startYear - 3)
  # ... 2nd lactation
  damOfFemalesLact2 = randCross(basePop, nCrosses =  nDamsOfFemalesLact2)
  damOfFemalesLact2@sex[] = "F"
  damOfFemalesLact2@father[] = "0"
  damOfFemalesLact2@mother[] = "0"
  damOfFemalesLact2 = fillInMisc(pop = damOfFemalesLact2, herds = herds, permEnvVar = permVar,
                                 year = startYear - 2)
  # ... 1st lactation
  damOfFemalesLact1 = randCross(basePop, nCrosses =  nDamsOfFemalesLact1)
  damOfFemalesLact1@sex[] = "F"
  damOfFemalesLact1@father[] = "0"
  damOfFemalesLact1@mother[] = "0"
  damOfFemalesLact1 = fillInMisc(pop = damOfFemalesLact1, herds = herds, permEnvVar = permVar,
                                 year = startYear - 1)

  ewesLact1 = c(damOfFemalesLact1, eliteEwesLact1)
  ewesLact2 = c(damOfFemalesLact2, eliteEwesLact2)
  ewesLact3 = c(damOfFemalesLact3, eliteEwesLact3)
  ewesLact4 = c(damOfFemalesLact4, eliteEwesLact4)

  # young ewes
  # yngFemales = randCross(basePop, nCrosses = nYngFemales)
  # yngFemales@sex[] = "F"
  # yngFemales@father[] = "0"
  # yngFemales@mother[] = "0"
  # yngFemales = fillInMisc(pop = yngFemales, herds = herds, permEnvVar = permVar,
  #                         year = startYear - 1/12)


  # lambs
  matingPlan1 = cbind(c(eliteEwesLact1@id, eliteEwesLact2@id, eliteEwesLact3@id, eliteEwesLact4@id),
                      sample(eliteSires@id, size = nEliteEwes, replace = TRUE))
  n = n1 - nEliteEwes
  damOfFemalesId = c(damOfFemalesLact1@id, damOfFemalesLact2@id, damOfFemalesLact3@id, damOfFemalesLact4@id)
  damOfFemalesIdForElite = sample(damOfFemalesId, size = n)

  damOfFemalesIdForRest = damOfFemalesId[!damOfFemalesId %in% damOfFemalesIdForElite]
  matingPlan2 = cbind(damOfFemalesIdForElite,
                      sample(eliteSires@id, size = n, replace = TRUE))

  n = nDamsOfFemales - (n1 - nEliteEwes)
  matingPlan3 = cbind(damOfFemalesIdForRest,
                      sample(c(siresOfFemales@id, wtRams1@id, ntlMatingRams@id), size = n, replace = TRUE))

  matingPlan = rbind(matingPlan1, matingPlan2, matingPlan3)
  lambs = makeCross2(females = c(eliteEwesLact1, eliteEwesLact2, eliteEwesLact3, eliteEwesLact4,
                                 damOfFemalesLact1, damOfFemalesLact2, damOfFemalesLact3, damOfFemalesLact4),
                     males = c(eliteSires, siresOfFemales, wtRams1, ntlMatingRams),
                     crossPlan = matingPlan)
  lambs = fillInMisc(pop = lambs,
                     mothers = c(eliteEwesLact1, eliteEwesLact2, eliteEwesLact3, eliteEwesLact4,
                                 damOfFemalesLact1, damOfFemalesLact2, damOfFemalesLact3, damOfFemalesLact4),
                     permEnvVar = permVar, year = startYear)

  database = recordData(pop = eliteSires, year = startYear)
  database = recordData(database, pop = siresOfFemales, year = startYear)
  database = recordData(database, pop = wtRams1, year = startYear)
  database = recordData(database, pop = wtRams2, year = startYear)
  database = recordData(database, pop = ntlMatingRams, year = startYear)
  database = recordData(database, pop = damOfFemalesLact4, year = startYear, lactation = 4)
  database = recordData(database, pop = damOfFemalesLact3, year = startYear, lactation = 3)
  database = recordData(database, pop = damOfFemalesLact2, year = startYear, lactation = 2)
  database = recordData(database, pop = damOfFemalesLact1, year = startYear, lactation = 1)
  database = recordData(database, pop = eliteEwesLact4, year = startYear, lactation = 4)
  database = recordData(database, pop = eliteEwesLact3, year = startYear, lactation = 3)
  database = recordData(database, pop = eliteEwesLact2, year = startYear, lactation = 2)
  database = recordData(database, pop = eliteEwesLact1, year = startYear, lactation = 1)
  # database = recordData(database, pop = lambs, year = startYear)
  database = recordData(database, pop = lambsFemales, year = startYear)
  database = recordData(database, pop = lambsMalesNatMat, year = startYear)
  database = recordData(database, pop = lambsMalesElite, year = startYear)
  # database = recordData(database, pop = yngRams, year = startYear)
  # database = recordData(database, pop = yngFemales, year = startYear)

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
  #plot the herd size
  for (year in yearToDo:nPTyrs) {
    yearFull = startYear + year
    print(paste("Working on Progeny testing stage year", year,
                Sys.time(), "...", sep = " "))
    #  TODO: not necessary to use the fullyear, better to use year of the simulation: yearfull = year
    yearEffect = sampleYearEffect(n = 1)
    herds$herdYearEffect = sampleHerdYearEffect(n = nHerds)

    # ---- Phenotyping ----

    damOfFemalesLact4 = setPhenoEwe(damOfFemalesLact4, varE = resVar,
                                    mean = meanLac4, yearEffect = yearEffect, herds = herds)

    damOfFemalesLact3 = setPhenoEwe(damOfFemalesLact3, varE = resVar,
                                    mean = meanLac3, yearEffect = yearEffect, herds = herds)

    damOfFemalesLact2 = setPhenoEwe(damOfFemalesLact2, varE = resVar,
                                    mean = meanLac2, yearEffect = yearEffect, herds = herds)

    damOfFemalesLact1 = setPhenoEwe(damOfFemalesLact1, varE = resVar,
                                    mean = meanLac1, yearEffect = yearEffect, herds = herds)

    eliteEwesLact4 = setPhenoEwe(eliteEwesLact4, varE = resVar, mean = meanLac4,
                                 yearEffect = yearEffect, herds = herds)

    eliteEwesLact3 = setPhenoEwe(eliteEwesLact3, varE = resVar, mean = meanLac3,
                                 yearEffect = yearEffect, herds = herds)

    eliteEwesLact2 = setPhenoEwe(eliteEwesLact2, varE = resVar, mean = meanLac2,
                                 yearEffect = yearEffect, herds = herds)

    eliteEwesLact1 = setPhenoEwe(eliteEwesLact1, varE = resVar, mean = meanLac1,
                                 yearEffect = yearEffect, herds = herds)

    database = setDatabasePheno(database, pop = damOfFemalesLact1, trait = 1)
    database = setDatabasePheno(database, pop = damOfFemalesLact2, trait = 1)
    database = setDatabasePheno(database, pop = damOfFemalesLact3, trait = 1)
    database = setDatabasePheno(database, pop = damOfFemalesLact4, trait = 1)
    database = setDatabasePheno(database, pop = eliteEwesLact1, trait = 1)
    database = setDatabasePheno(database, pop = eliteEwesLact2, trait = 1)
    database = setDatabasePheno(database, pop = eliteEwesLact3, trait = 1)
    database = setDatabasePheno(database, pop = eliteEwesLact4, trait = 1)

    # ---- SELECTION BY CATEGORIES ----
    print(paste("Working on Rams selection year", year,
                Sys.time(), "...", sep = " "))
    # ---- Rams ----
    # ---- Elite sires ----
    print(paste("elite sires part", year,
                Sys.time(), "...", sep = " "))
    eliteSires3 = eliteSires2 # eliteSires3 are 4.5 years old here
    eliteSires2 = eliteSires1 # eliteSires2 are 3.5 years old here
    if (year == 1) {
      use = "rand"
    } else {
      use= "ebv"
    }
    # TODO: this if (!all... is needed because we set fathers to 0 - related to
    if (all(wtRams2@father == "0")) {
      eliteSires1 = selectInd(pop = wtRams2, nInd = nEliteSires1, # eliteSires1 are 2.5 years old here
                              use = use, trait = 1)
    } else {
      eliteSires1 = selectWithinFam(pop = wtRams2, nInd = 1, # eliteSires1 are 2.5 years old here
                                    use = use, trait = 1,
                                    famType = "M")
      eliteSires1 = selectInd(pop = eliteSires1, nInd = nEliteSires1,
                              use = use, trait = 1)
    }
    print(paste("sire of dams part", year,
                Sys.time(), "...", sep = " "))
    # ---- Sires of Dams ----
    siresOfFemales3 = siresOfFemales2 # siresOfFemales3 are 4.5 years old here
    siresOfFemales2 = siresOfFemales1 # siresOfFemales2 are 3.5 years old here
    if (year == 1) {
      use = "rand"
    } else {
      use= "ebv"
    }
    if (all(wtRams2@father == "0")) {
      siresOfFemales1 = selectInd(pop = wtRams2[!wtRams2@iid %in% eliteSires1@iid], nInd = nsiresOfFemales1, # siresOfFemales1 are 2.5 years old here
                                  use = use, trait = 1)
    } else {
      siresOfFemales1 = selectWithinFam(pop = wtRams2[!wtRams2@iid %in% eliteSires1@iid], nInd = 2, # siresOfFemales1 are 2.5 years old here
                                        use = use, trait = 1,
                                        famType = "M")
      siresOfFemales1 = selectInd(pop = siresOfFemales1, nInd = nsiresOfFemales1,
                                  use = use, trait = 1)
    }

    # TODO: change the nind in selectWithinFam to 2
    # TODO: better to use selectTop (TURE or FALSE) to choose th animals in selectInd

    eliteEwes = c(eliteEwesLact1, eliteEwesLact2, eliteEwesLact3, eliteEwesLact4)

    # We need this sel and n before we update eliteSires
    # Note that we select yngRams only from lambs from proven rams (not waiting rams)
    sel = lambs@father %in% eliteSires@id & lambs@mother %in% eliteEwes@id
    sel1 = lambs@father %in% c(eliteSires@id, siresOfFemales@id)
    # sel2 = ewesLact1@father %in% c(eliteSires@id, siresOfFemales@id)
    # sel3 = ewesLact2@father %in% c(eliteSires@id, siresOfFemales@id)
    # sel4 = ewesLact3@father %in% c(eliteSires@id, siresOfFemales@id)
    n = ceiling(nWtRams1 / length(eliteSires@id))
    eliteSires = c(eliteSires3, eliteSires2, eliteSires1)
    siresOfFemales = c(siresOfFemales3, siresOfFemales2, siresOfFemales1)
    wtRams2 = wtRams1
    wtRams1 = selectWithinFam(pop = lambs[sel], nInd = n, # yngRams are 0.12 year old here
                              use = use, trait = 1,
                              sex = "M", famType = "M", simParam = SP)
    wtRams1 = selectInd(pop = wtRams1, nInd = nWtRams1,
                        use = use, trait = 1)
    # wtRams1 = lambsMalesElite

    ntlMatingRams =  selectInd(pop = lambs[sel1], nInd = nNaturalMatingRams ,
                               use = use, trait = 1, sex = "M", famType = "M")
    # TODO: i need to exclude the lambs in ntlMAtig (sel1) from those that are in sel
    # ntlMatingRams = lambsMalesNatMat

    # ---- ewes ----
    # --- Elites ---
    print(paste("Working on ewes selection year", year,
                Sys.time(), "...", sep = " "))
    ewesLact1 = c(damOfFemalesLact1, eliteEwesLact1)
    ewesLact2 = c(damOfFemalesLact2, eliteEwesLact2)
    ewesLact3 = c(damOfFemalesLact3, eliteEwesLact3)
    ewesLact4 = c(damOfFemalesLact4, eliteEwesLact4)

    eliteEwesLact4 = selectInd(ewesLact3, nInd = round(0.16 * nEliteEwes), use = use) # eliteEwesLact4 are 4 years old here
    eliteEwesLact3 = selectInd(ewesLact2, nInd = round(0.26 * nEliteEwes), use = use) # eliteEwesLact3 are 3 years old here
    eliteEwesLact2 = selectInd(ewesLact1, nInd = round(0.35 * nEliteEwes), use = use) # eliteEwesLact2 are 2 years old here
    eliteEwesLact1 = selectInd(lambs[sel1], nInd = round(0.23 * nEliteEwes), use = use, sex ="F", trait = 1, famType = "M")   # eliteEwesLact1 are 1 years old here

    # Set phenotypes of the last lactation from the last generation to missing.
    # These are only copied phenotypes from the latest lactation.
    # Correct phenotypes will be added at the beginning of the next year.
    # Phenotypes of the last lactation are not part of the evaluation because lambs are selected before the
    # phenotype is recorded.
    # eliteEwesLact4@pheno[,c(1:nTraits)] = NA
    # eliteEwesLact3@pheno[,c(1:nTraits)] = NA
    # eliteEwesLact2@pheno[,c(1:nTraits)] = NA
    # eliteEwesLact1@pheno[,c(1:nTraits)] = NA

    # ---Dams of Females---
    damOfFemalesLact4 = selectInd(ewesLact3[!ewesLact3@iid %in% eliteEwesLact4@iid], nInd = round(0.18 * nDamsOfFemales), use = "rand") # damOfFemalesLact4 are 4 years old here
    damOfFemalesLact3 = selectInd(ewesLact2[!ewesLact2@iid %in% eliteEwesLact3@iid], nInd = round(0.24 * nDamsOfFemales), use = "rand") # damOfFemalesLact3 are 3 years old here
    damOfFemalesLact2 = selectInd(ewesLact1[!ewesLact1@iid %in% eliteEwesLact2@iid], nInd = round(0.28 * nDamsOfFemales), use = "rand") # damOfFemalesLact2 are 2 years old here
    damOfFemalesLact1 = selectInd(pop=lambs[!lambs@iid %in% eliteEwesLact1@iid], nInd=round(0.3 * nDamsOfFemales),  use = "rand", sex = "F") # damOfFemalesLact1 are 1 years old here

    #  TODO: better to keep only this part of the code ewesLact3[!ewesLact3@iid %in% eliteEwesLact4@iid] as below
    # damOfFemalesLact4 = !ewesLact3@iid %in% eliteEwesLact4@iid # damOfFemalesLact4 are 4 years old here
    # damOfFemalesLact3 = !ewesLact2@iid %in% eliteEwesLact3@iid # damOfFemalesLact3 are 3 years old here
    # damOfFemalesLact2 = !ewesLact1@iid %in% eliteEwesLact2@iid # damOfFemalesLact2 are 2 years old here
    # damOfFemalesLact1 = !lambs@iid %in% eliteEwesLact1@iid # damOfFemalesLact1 are 1 years old here


    # yngFemales = selectInd(lambs, nInd = nYngFemales, use = "rand", sex = "F")

    # Set phenotypes of the last lactation from the last generation to missing.
    # These are only copied phenotypes from the latest lactation.
    # Correct phenotypes will be added at the beginning of the next year.
    # Phenotypes of the last lactation are not part of the evaluation because calves are selected before the
    # phenotype is recorded.
    # damOfFemalesLact4@pheno[,c(1:nTraits)] = NA
    # damOfFemalesLact3@pheno[,c(1:nTraits)] = NA
    # damOfFemalesLact2@pheno[,c(1:nTraits)] = NA
    # damOfFemalesLact1@pheno[,c(1:nTraits)] = NA


    # ---- lambs ----
    print(paste("Working on lambs year", year,
                Sys.time(), "...", sep = " "))
    # See above comments on stage of animals and their use in the fill-in stage
    ramsId = c(sample(eliteSires@id, size = n1, replace = TRUE),
               sample(siresOfFemales@id,size = n2, replace = TRUE),
               sample(wtRams1@id, size = n3, replace = TRUE),
               sample(ntlMatingRams@id,size = n4, replace = TRUE))
    matingPlan = cbind(c(eliteEwesLact1@id, eliteEwesLact2@id, eliteEwesLact3@id, eliteEwesLact4@id, damOfFemalesLact1@id, damOfFemalesLact2@id, damOfFemalesLact3@id, damOfFemalesLact4@id), ramsId)
    #
    # lambsMalesElite ONLY 150, they need to be the offsrping of elite sires
    # lambsMalesNatMat ONLY 1000 they can
    #
    # REDUCE THE MATING PLAN BY HALF to have the females
    #  the ting that I can do is that I can do a vecot of the males ID and then shuttle it and then choose form it randomly to mate them with the females depending on what we want
    # lambsFemales = makeCross2(females = c(eliteEwesLact1, eliteEwesLact2, eliteEwesLact3, eliteEwesLact4, damOfFemalesLact1, damOfFemalesLact2, damOfFemalesLact3, damOfFemalesLact4),
    #                           males = c(eliteSires, siresOfFemales, wtRams1, ntlMatingRams),
    #                           crossPlan = matingPlan)
    # lambsFemales@sex = "F"
    #

    lambs = makeCross2(females = c(eliteEwesLact1, eliteEwesLact2, eliteEwesLact3, eliteEwesLact4, damOfFemalesLact1, damOfFemalesLact2, damOfFemalesLact3, damOfFemalesLact4),
                       males = c(eliteSires, siresOfFemales, wtRams1, ntlMatingRams),
                       crossPlan = matingPlan)
    #  TODO: better to shuffle the id of the females in each category for better mating 'data.frame(a[sample(nrow(a)),])'
    lambs = fillInMisc(pop = lambs, mothers = c(eliteEwesLact1, eliteEwesLact2, eliteEwesLact3, eliteEwesLact4, damOfFemalesLact1, damOfFemalesLact2, damOfFemalesLact3, damOfFemalesLact4),
                       permEnvVar = permVar)#, year = yearFull)

    # can I do this?
    # lambsFemales = selectInd(lambs, sex="F", nInd = sum(lambs@sex=="F"), use="rand")
    # lambsMalesElite = selectWithinFam(pop = lambs[sel], nInd = n,
    #                                   use = use, trait = 1,
    #                                   sex = "M", famType = "M")
    # lambsMalesElite = selectInd(lambs, sex="M", nInd= nWtRams1, use= use)
    # lambsMalesNatMat = selectInd(pop = lambs[sel1], nInd = nNaturalMatingRams ,
    #                              use = use, trait = 1, sex = "M", famType = "M")

    # Data recording
    database = recordData(database, pop = eliteSires, year = yearFull)
    database = recordData(database, pop = siresOfFemales, year = yearFull)
    database = recordData(database, pop = wtRams1, year = yearFull)
    database = recordData(database, pop = wtRams2, year = yearFull)
    database = recordData(database, pop = ntlMatingRams, year = yearFull)
    database = recordData(database, pop = eliteEwesLact4, year = yearFull, lactation = 4)
    database = recordData(database, pop = eliteEwesLact3, year = yearFull, lactation = 3)
    database = recordData(database, pop = eliteEwesLact2, year = yearFull, lactation = 2)
    database = recordData(database, pop = eliteEwesLact1, year = yearFull, lactation = 1)
    database = recordData(database, pop = damOfFemalesLact4, year = yearFull, lactation = 4)
    database = recordData(database, pop = damOfFemalesLact3, year = yearFull, lactation = 3)
    database = recordData(database, pop = damOfFemalesLact2, year = yearFull, lactation = 2)
    database = recordData(database, pop = damOfFemalesLact1, year = yearFull, lactation = 1)
    database = recordData(database, pop = lambs, year = yearFull)
    # database = recordData(database, pop = lambsFemales, year = startYear)
    # database = recordData(database, pop = lambsMalesNatMat, year = startYear)
    # database = recordData(database, pop = lambsMalesElite, year = startYear)

    # database = recordData(database, pop = yngRams, year = yearFull)
    # database = recordData(database, pop = yngFemales, year = yearFull)
    # database = recordData(database, pop = lambsMaleElite_wtrams, year = yearfull)

    save(x = database, file = "database.RData")
    # EBV
    print(paste("Setting EBV", year,
                Sys.time(), "...", sep = " "))

    OldDir = getwd()
    Dir = paste("Blup", year, sep = "_")
    unlink(paste(OldDir,Dir, sep = "/"), recursive = TRUE)
    dir.create(path = Dir, showWarnings = FALSE)
    setwd(dir = Dir)

    pedEbv = estimateBreedingValues(pedigree = SP$pedigree,
                                    database = database,
                                    trait = 1,
                                    vars = list(varA  = addVar,
                                                varPE = permVar,
                                                varHY = herdYearVar,
                                                varE  = resVar))

    setwd(dir = OldDir)

    # Set EBVs for every population
    eliteSires        = setEbv(eliteSires, ebv = pedEbv)
    siresOfFemales    = setEbv(siresOfFemales, ebv = pedEbv)
    wtRams1           = setEbv(wtRams1, ebv = pedEbv)
    wtRams2           = setEbv(wtRams2, ebv = pedEbv)
    ntlMatingRams     = setEbv(ntlMatingRams, ebv = pedEbv)
    eliteEwesLact4    = setEbv(eliteEwesLact4, ebv = pedEbv)
    eliteEwesLact3    = setEbv(eliteEwesLact3, ebv = pedEbv)
    eliteEwesLact2    = setEbv(eliteEwesLact2, ebv = pedEbv)
    eliteEwesLact1    = setEbv(eliteEwesLact1, ebv = pedEbv)
    damOfFemalesLact4 = setEbv(damOfFemalesLact4, ebv = pedEbv)
    damOfFemalesLact3 = setEbv(damOfFemalesLact3, ebv = pedEbv)
    damOfFemalesLact2 = setEbv(damOfFemalesLact2, ebv = pedEbv)
    damOfFemalesLact1 = setEbv(damOfFemalesLact1, ebv = pedEbv)
    lambs             = setEbv(lambs, ebv = pedEbv)

    database = setDatabaseEbv(database, pop = eliteSires)
    database = setDatabaseEbv(database, pop = siresOfFemales)
    database = setDatabaseEbv(database, pop = wtRams1)
    database = setDatabaseEbv(database, pop = wtRams2)
    database = setDatabaseEbv(database, pop = ntlMatingRams)
    database = setDatabaseEbv(database, pop = eliteEwesLact4)
    database = setDatabaseEbv(database, pop = eliteEwesLact3)
    database = setDatabaseEbv(database, pop = eliteEwesLact2)
    database = setDatabaseEbv(database, pop = eliteEwesLact1)
    database = setDatabaseEbv(database, pop = damOfFemalesLact4)
    database = setDatabaseEbv(database, pop = damOfFemalesLact3)
    database = setDatabaseEbv(database, pop = damOfFemalesLact2)
    database = setDatabaseEbv(database, pop = damOfFemalesLact1)
    TODO database = setDatabaseEbv(database, pop = lambs)
    #  here i need to add all the categories of lambs (wtrams, ntlmating, females)

    # yngRams@ebv             =as.matrix(pedEbv[pedEbv$IId %in% yngRams@iid, "EBV"])
    # yngFemales@ebv          =as.matrix(pedEbv[pedEbv$IId %in% yngFemales@iid, "EBV"])

    save.image(file = paste("image_Year", year, "_PT.RData", sep = ""))
    eliteEwes = c(eliteEwesLact1, eliteEwesLact2, eliteEwesLact3, eliteEwesLact4)
    damOfFemales=c(damOfFemalesLact1, damOfFemalesLact2, damOfFemalesLact3, damOfFemalesLact4)
    correlation = data.frame( eliteSires = cor(eliteSires@ebv, eliteSires@gv),
                              siresOfFemales = cor(siresOfFemales@ebv, siresOfFemales@gv),
                              eliteEwes = cor(eliteEwes@ebv, eliteEwes@gv),
                              damOfFemales = cor(damOfFemales@ebv, damOfFemales@gv),
                              wtRams1 = cor(wtRams1@ebv, wtRams1@gv),
                              wtRams2 = cor(wtRams2@ebv, wtRams2@gv),
                              ntlMatingRams = cor(ntlMatingRams@ebv, ntlMatingRams@gv),
                              row.names=paste("year", year, sep="")
    )
    write.table(x = correlation, file = paste("correlation", year, sep=""))
  }
}

save.image(file = paste("burnin", year, "_PT.RData", sep = ""))

setwd(dir = mainDir)

# "ELite_sire" "Sire_of_females" "eliteEwes" "damsOfFemales" "wtRams1" "wtRams2" "ntlMatingSires"
# cat correlation1 > cor
# sed -e '1d' correlation2 >> cor
# sed -e '1d' correlation3 >> cor
# sed -e '1d' correlation4 >> cor
# sed -e '1d' correlation5 >> cor
# sed -e '1d' correlation6 >> cor
# sed -e '1d' correlation7 >> cor
# sed -e '1d' correlation8 >> cor
# sed -e '1d' correlation9 >> cor
# sed -e '1d' correlation10 >>cor
# sed -e '1d' correlation11 >>cor
# sed -e '1d' correlation12 >>cor
# sed -e '1d' correlation13 >>cor
# sed -e '1d' correlation14 >>cor
# sed -e '1d' correlation15 >>cor
# sed -e '1d' correlation16 >>cor
# sed -e '1d' correlation17 >>cor
# sed -e '1d' correlation18 >>cor
# sed -e '1d' correlation19 >>cor
# sed -e '1d' correlation20 >>cor
