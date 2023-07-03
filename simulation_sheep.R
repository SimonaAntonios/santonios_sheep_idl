# Long-term selection in dairy ovine (Founders/Burn-in)
# for one trait: Milk Yield
# Simona Antonios
# source("simulation_sheep.r", echo = TRUE)

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
  paste0(source("functions.R", echo = TRUE))
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

histNe = c(200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5500, 6000, 6500, 7000, 7500) # historical Ne
histGen = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 5000, 7000, 8000, 9000, 10000) # historical generations
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

  if (FALSE) {
    founderPop = runMacs2(nInd = 10 * BaseNe,
                          nChr = nChr,
                          segSites = nSNPPerChr,
                          Ne = BaseNe,
                          bp = ChrSize,
                          genLen =  RecRate * ChrSize *26 ,
                          mutRate = MutRate,
                          histNe = histNe,
                          histGen = histGen,
                          # returnCommand = TRUE
    )
    save(founderPop, file = "founder_runMacs2.RData")
  }
  load("founder_runMacs2.RData")
  # load(file = paste("burnin", year, "_PT.RData", sep = ""))
  # sourceFunctions()

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

  SP$addTraitAD(nQtlPerChr = nQTLPerChr, mean = 0, var = addVar, meanDD = meanDD, varDD = varDD)
  # we will add lactation means later in setPhenoEwe
  # SP$setVarE() we assume that only ewes will have phenotypes, hence do not set this!!! Use setPhenoEwe() instead!
  # SP$setCorE() we assume that only ewes will have phenotypes, hence do not set this!!! Use setPhenoEwe() instead!

  # ---- Herds and herd and year effects ----

  herdSize = simcmp(n = (10 * nHerds),
                    v = c(round(meanHerdSize / 10), round(sdHerdSize / 10)))
  # hist(herdSize)
  herdSize = herdSize * 10
  herdSize = herdSize[herdSize > 100]
  herdSize = sample(x = herdSize, size = nHerds)
  herdEffect = cbind(rnorm(n = nHerds, mean = 0,
                           sd = sqrt(herdVar[1])))
  herdYearEffect = sampleHerdYearEffect(n = nHerds)
  herds = list(
    herd = 1:nHerds,
    herdSize = herdSize,
    herdEffect = as.matrix(herdEffect),
    herdYearEffect = as.matrix(herdYearEffect)
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

  database = recordData(pop = eliteSires3, year = startYear)
  database = recordData(database, pop = eliteSires2, year = startYear)
  database = recordData(database, pop = eliteSires1, year = startYear)
  database = recordData(database, pop = eliteSires, year = startYear)
  database = recordData(database, pop = siresOfFemales3, year = startYear)
  database = recordData(database, pop = siresOfFemales2, year = startYear)
  database = recordData(database, pop = siresOfFemales1, year = startYear)
  database = recordData(database, pop = siresOfFemales, year = startYear)
  database = recordData(database, pop = wtRams2, year = startYear)
  database = recordData(database, pop = wtRams1, year = startYear)
  database = recordData(database, pop = ntlMatingRams, year = startYear)
  database = recordData(database, pop = eliteEwesLact4, year = startYear, lactation = 4)
  database = recordData(database, pop = eliteEwesLact3, year = startYear, lactation = 3)
  database = recordData(database, pop = eliteEwesLact2, year = startYear, lactation = 2)
  database = recordData(database, pop = eliteEwesLact1, year = startYear, lactation = 1)
  database = recordData(database, pop = damOfFemalesLact4, year = startYear, lactation = 4)
  database = recordData(database, pop = damOfFemalesLact3, year = startYear, lactation = 3)
  database = recordData(database, pop = damOfFemalesLact2, year = startYear, lactation = 2)
  database = recordData(database, pop = damOfFemalesLact1, year = startYear, lactation = 1)
  database = recordData(database, pop = lambs, year = startYear)

  # Save current state
  # save(x = database, file = "database.RData")
  save.image(file = "image_Year0_FillIn.RData")
  # load(file = "image_Year0_FillIn.RData")
  # sourceFunctions()

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
    # year = yearToDo
    yearFull = startYear + year
    print(paste("Working on Progeny testing stage year", year, Sys.time(), "...", sep = " "))
    yearEffect = sampleYearEffect(n = 1)
    herds$herdYearEffect = sampleHerdYearEffect(n = nHerds)

    if (year <= 1) {
      use = "rand"
    } else {
      use= "ebv"
    }

    # ---- Phenotyping ----

    damOfFemalesLact4 = setPhenoEwe(damOfFemalesLact4, varE = resVar,
                                    mean = meanLac4, yearEffect = yearEffect, herds = herds)

    damOfFemalesLact3 = setPhenoEwe(damOfFemalesLact3, varE = resVar,
                                    mean = meanLac3, yearEffect = yearEffect, herds = herds)

    damOfFemalesLact2 = setPhenoEwe(damOfFemalesLact2, varE = resVar,
                                    mean = meanLac2, yearEffect = yearEffect, herds = herds)

    damOfFemalesLact1 = setPhenoEwe(damOfFemalesLact1, varE = resVar,
                                    mean = meanLac1, yearEffect = yearEffect, herds = herds)

    eliteEwesLact4 = setPhenoEwe(eliteEwesLact4, varE = resVar,
                                 mean = meanLac4, yearEffect = yearEffect, herds = herds)

    eliteEwesLact3 = setPhenoEwe(eliteEwesLact3, varE = resVar,
                                 mean = meanLac3, yearEffect = yearEffect, herds = herds)

    eliteEwesLact2 = setPhenoEwe(eliteEwesLact2, varE = resVar,
                                 mean = meanLac2, yearEffect = yearEffect, herds = herds)

    eliteEwesLact1 = setPhenoEwe(eliteEwesLact1, varE = resVar,
                                 mean = meanLac1, yearEffect = yearEffect, herds = herds)

    database = setDatabasePheno(database, pop = damOfFemalesLact4)
    database = setDatabasePheno(database, pop = damOfFemalesLact3)
    database = setDatabasePheno(database, pop = damOfFemalesLact2)
    database = setDatabasePheno(database, pop = damOfFemalesLact1)
    database = setDatabasePheno(database, pop = eliteEwesLact4)
    database = setDatabasePheno(database, pop = eliteEwesLact3)
    database = setDatabasePheno(database, pop = eliteEwesLact2)
    database = setDatabasePheno(database, pop = eliteEwesLact1)

    # ---- Rams ----
    print(paste("Working on Ram selection year", year, Sys.time(), "...", sep = " "))
    # ---- ... Elite sires ----
    eliteSires3 = eliteSires2 # eliteSires3 are 4.5 years old here
    eliteSires2 = eliteSires1 # eliteSires2 are 3.5 years old here

    if (all(wtRams2@father == "0")) {
      eliteSires1 = selectInd(pop = wtRams2, nInd = nEliteSires1, # eliteSires1 are 3.5 years old here
                              use = use)
    } else {
      eliteSires1 = selectWithinFam(pop = wtRams2, nInd = 1, # eliteSires1 are 3.5 years old here
                                    use = use,
                                    famType = "M")
      eliteSires1 = selectInd(pop = eliteSires1, nInd = nEliteSires1,
                              use = use)
    }

    # ---- ... Sires of Dams ----
    siresOfFemales3 = siresOfFemales2 # siresOfFemales3 are 5.5 years old here
    siresOfFemales2 = siresOfFemales1 # siresOfFemales2 are 4.5 years old here

    if (all(wtRams2@father == "0")) {
      siresOfFemales1 = selectInd(pop = wtRams2[!wtRams2@id %in% eliteSires1@id], nInd = nsiresOfFemales1, # siresOfFemales1 are 3.5 years old here
                                  use = use)
    } else {
      siresOfFemales1 = selectWithinFam(pop = wtRams2[!wtRams2@id %in% eliteSires1@id], nInd = 2, # siresOfFemales1 are 3.5 years old here
                                        use = use,
                                        famType = "M")
      siresOfFemales1 = selectInd(pop = siresOfFemales1, nInd = nsiresOfFemales1,
                                  use = use)
    }

    eliteEwes = c(eliteEwesLact1, eliteEwesLact2, eliteEwesLact3, eliteEwesLact4)

    # Note that we select yngRams only from lambs from proven rams (not waiting rams)
    selWtRams = lambs@father %in% eliteSires@id & lambs@mother %in% eliteEwes@id
    n = ceiling(nWtRams1 / length(eliteSires@id))
    wtRams2 = wtRams1
    wtRams1 = selectWithinFam(pop = lambs[selWtRams], nInd = n,
                              use = use,
                              sex = "M", famType = "M", simParam = SP)
    wtRams1 = selectInd(pop = wtRams1, nInd = nWtRams1,
                        use = use)

    selNtlRams = lambs@father %in% c(eliteSires@id, siresOfFemales@id) &
      !lambs@id %in% wtRams1@id
    ntlMatingRams =  selectInd(pop = lambs[selNtlRams], nInd = nNaturalMatingRams,
                               use = use, sex = "M", famType = "M")

    # Have to do it here since need to keep track of previous eliteSires in above steps
    eliteSires = c(eliteSires3, eliteSires2, eliteSires1)
    siresOfFemales = c(siresOfFemales3, siresOfFemales2, siresOfFemales1)

    # ---- Ewes ----
    print(paste("Working on Ewe selection year", year, Sys.time(), "...", sep = " "))
    # ---- ... Elite ewes ----
    ewesLact1 = c(damOfFemalesLact1, eliteEwesLact1)
    ewesLact2 = c(damOfFemalesLact2, eliteEwesLact2)
    ewesLact3 = c(damOfFemalesLact3, eliteEwesLact3)
    ewesLact4 = c(damOfFemalesLact4, eliteEwesLact4)

    # is it a good way to select the elite ewes, because they need to be the daughters of an AI male??
    eliteEwesLact4 = selectInd(ewesLact3, nInd = round(0.16 * nEliteEwes), use = use) # eliteEwesLact4 are 4 years old here
    eliteEwesLact3 = selectInd(ewesLact2, nInd = round(0.26 * nEliteEwes), use = use) # eliteEwesLact3 are 3 years old here
    eliteEwesLact2 = selectInd(ewesLact1, nInd = round(0.35 * nEliteEwes), use = use) # eliteEwesLact2 are 2 years old here
    eliteEwesLact1 = selectInd(lambs[selNtlRams], nInd = round(0.23 * nEliteEwes), use = use, sex ="F", famType = "M") # eliteEwesLact1 are 1 years old here

    # Set phenotypes to missing, because these are copied from the previous lactation.
    # We do this because of the recordData() call below - that would save wrong
    # phenotypes for these animals in this life stage! Correct phenotypes for these
    # animals in this life stage will be added at the beginning of the next year.
    eliteEwesLact4@pheno[,] = NA
    eliteEwesLact3@pheno[,] = NA
    eliteEwesLact2@pheno[,] = NA
    eliteEwesLact1@pheno[,] = NA

    # ---- ... Dams of Females ----

    damOfFemalesLact4 = selectInd(ewesLact3[!ewesLact3@id %in% eliteEwesLact4@id],
                                  nInd = round(0.18 * nDamsOfFemales), use = "rand") # damOfFemalesLact4 are 4 years old here
    damOfFemalesLact3 = selectInd(ewesLact2[!ewesLact2@id %in% eliteEwesLact3@id],
                                  nInd = round(0.24 * nDamsOfFemales), use = "rand") # damOfFemalesLact3 are 3 years old here
    damOfFemalesLact2 = selectInd(ewesLact1[!ewesLact1@id %in% eliteEwesLact2@id],
                                  nInd = round(0.28 * nDamsOfFemales), use = "rand") # damOfFemalesLact2 are 2 years old here
    damOfFemalesLact1 = selectInd(pop=lambs[!lambs@id %in% eliteEwesLact1@id],
                                  nInd = round(0.30 * nDamsOfFemales), use = "rand", sex = "F") # damOfFemalesLact1 are 1 years old here

    # Set phenotypes to missing, because these are copied from the previous lactation.
    # We do this because of the recordData() call below - that would save wrong
    # phenotypes for these animals in this life stage! Correct phenotypes for these
    # animals in this life stage will be added at the beginning of the next year.
    damOfFemalesLact4@pheno[,] = NA
    damOfFemalesLact3@pheno[,] = NA
    damOfFemalesLact2@pheno[,] = NA
    damOfFemalesLact1@pheno[,] = NA

    # ---- Lambs ----
    print(paste("Working on Lambs year", year, Sys.time(), "...", sep = " "))
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
                       permEnvVar = permVar, year = yearFull)

    # ---- Data recording & EBV ----

    database = recordData(database, pop = eliteSires3, year = yearFull)
    database = recordData(database, pop = eliteSires2, year = yearFull)
    database = recordData(database, pop = eliteSires1, year = yearFull)
    database = recordData(database, pop = eliteSires, year = yearFull)
    database = recordData(database, pop = siresOfFemales3, year = yearFull)
    database = recordData(database, pop = siresOfFemales2, year = yearFull)
    database = recordData(database, pop = siresOfFemales1, year = yearFull)
    database = recordData(database, pop = siresOfFemales, year = yearFull)
    database = recordData(database, pop = wtRams2, year = yearFull)
    database = recordData(database, pop = wtRams1, year = yearFull)
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
    # save(x = database, file = "database.RData")
    # load(file = "database.RData")
    # sourceFunctions()

    print(paste("Estimating breeding values", year, Sys.time(), "...", sep = " "))

    # We will remove male lambs from the evaluation - those that will never
    #   contribute to next generations
    sel = database$General$Pop == "lambs" & database$General$Sex == "M"
    maleLambs = database$General[sel, "IId"]
    sel = database$General$Pop %in% c("wtRams1", "ntlMatingRams")
    reproMales = database$General[sel, "IId"]
    culledMaleLambs = maleLambs$IId[!maleLambs$IId %in% reproMales$IId]
    removeCulledMaleLambs = row.names(SP$pedigree) %in% culledMaleLambs
    # sum(removeCulledMaleLambs); sum(!removeCulledMaleLambs)

    variances = list(varPE = permVar,
                     varA  = addVar,
                     varE  = resVar)
    # TODO
    # if (model == "X") {
    #   variances = list()
    # }
    pedEbv = estimateBreedingValues(pedigree = SP$pedigree,
                                    database = database,
                                    vars = variances,
                                    removeFromEvaluation = removeCulledMaleLambs)

    # Set EBVs for every population
    eliteSires3       = setEbv(eliteSires3, ebv = pedEbv)
    eliteSires2       = setEbv(eliteSires2, ebv = pedEbv)
    eliteSires1       = setEbv(eliteSires1, ebv = pedEbv)
    eliteSires        = setEbv(eliteSires, ebv = pedEbv)
    siresOfFemales3   = setEbv(siresOfFemales3, ebv = pedEbv)
    siresOfFemales2   = setEbv(siresOfFemales2, ebv = pedEbv)
    siresOfFemales1   = setEbv(siresOfFemales1, ebv = pedEbv)
    siresOfFemales    = setEbv(siresOfFemales, ebv = pedEbv)
    wtRams2           = setEbv(wtRams2, ebv = pedEbv)
    wtRams1           = setEbv(wtRams1, ebv = pedEbv)
    ntlMatingRams     = setEbv(ntlMatingRams, ebv = pedEbv)
    eliteEwesLact4    = setEbv(eliteEwesLact4, ebv = pedEbv)
    eliteEwesLact3    = setEbv(eliteEwesLact3, ebv = pedEbv)
    eliteEwesLact2    = setEbv(eliteEwesLact2, ebv = pedEbv)
    eliteEwesLact1    = setEbv(eliteEwesLact1, ebv = pedEbv)
    damOfFemalesLact4 = setEbv(damOfFemalesLact4, ebv = pedEbv)
    damOfFemalesLact3 = setEbv(damOfFemalesLact3, ebv = pedEbv)
    damOfFemalesLact2 = setEbv(damOfFemalesLact2, ebv = pedEbv)
    damOfFemalesLact1 = setEbv(damOfFemalesLact1, ebv = pedEbv)
    # ... we could use this, but we have excluded some culled lambs from the evaluation
    #     so we will just calculate parent average for all lambs - this works for
    #     pedigree BLUP, but not for genomic BLUP (if some lambs would be genotyped)
    # lambs             = setEbv(lambs, ebv = pedEbv)
    selM = match(x = lambs@mother, table = pedEbv$IId)
    selF = match(x = lambs@father, table = pedEbv$IId)
    lambs@ebv = as.matrix((pedEbv[selM, -1] + pedEbv[selF, -1]) / 2)

    # Data recording
    database = setDatabaseEbv(database, pop = eliteSires3)
    database = setDatabaseEbv(database, pop = eliteSires2)
    database = setDatabaseEbv(database, pop = eliteSires1)
    database = setDatabaseEbv(database, pop = eliteSires)
    database = setDatabaseEbv(database, pop = siresOfFemales3)
    database = setDatabaseEbv(database, pop = siresOfFemales2)
    database = setDatabaseEbv(database, pop = siresOfFemales1)
    database = setDatabaseEbv(database, pop = siresOfFemales)
    database = setDatabaseEbv(database, pop = wtRams2)
    database = setDatabaseEbv(database, pop = wtRams1)
    database = setDatabaseEbv(database, pop = ntlMatingRams)
    database = setDatabaseEbv(database, pop = eliteEwesLact4)
    database = setDatabaseEbv(database, pop = eliteEwesLact3)
    database = setDatabaseEbv(database, pop = eliteEwesLact2)
    database = setDatabaseEbv(database, pop = eliteEwesLact1)
    database = setDatabaseEbv(database, pop = damOfFemalesLact4)
    database = setDatabaseEbv(database, pop = damOfFemalesLact3)
    database = setDatabaseEbv(database, pop = damOfFemalesLact2)
    database = setDatabaseEbv(database, pop = damOfFemalesLact1)
    database = setDatabaseEbv(database, pop = lambs)

    eliteEwes = c(eliteEwesLact1, eliteEwesLact2, eliteEwesLact3, eliteEwesLact4)
    damOfFemales = c(damOfFemalesLact1, damOfFemalesLact2, damOfFemalesLact3, damOfFemalesLact4)

    # ---- Summarising EBV ----

    correlation = data.frame(year = year,
                             eliteSires3 = ebvAccuracy(eliteSires3),
                             eliteSires2 = ebvAccuracy(eliteSires2),
                             eliteSires1 = ebvAccuracy(eliteSires1),
                             eliteSires = ebvAccuracy(eliteSires),
                             siresOfFemales3 = ebvAccuracy(siresOfFemales3),
                             siresOfFemales2 = ebvAccuracy(siresOfFemales2),
                             siresOfFemales1 = ebvAccuracy(siresOfFemales1),
                             siresOfFemales = ebvAccuracy(siresOfFemales),
                             wtRams2 = ebvAccuracy(wtRams2),
                             wtRams1 = ebvAccuracy(wtRams1),
                             ntlMatingRams = ebvAccuracy(ntlMatingRams),
                             eliteEwesLact4 = ebvAccuracy(eliteEwesLact4),
                             eliteEwesLact3 = ebvAccuracy(eliteEwesLact3),
                             eliteEwesLact2 = ebvAccuracy(eliteEwesLact2),
                             eliteEwesLact1 = ebvAccuracy(eliteEwesLact1),
                             eliteEwes = ebvAccuracy(eliteEwes),
                             damOfFemalesLact4 = ebvAccuracy(damOfFemalesLact4),
                             damOfFemalesLact3 = ebvAccuracy(damOfFemalesLact3),
                             damOfFemalesLact2 = ebvAccuracy(damOfFemalesLact2),
                             damOfFemalesLact1 = ebvAccuracy(damOfFemalesLact1),
                             damOfFemales = ebvAccuracy(damOfFemales),
                             lambs = ebvAccuracy(lambs))
    if (year == 1) {
      add = FALSE
    } else {
      add = TRUE
    }
    write.table(x = correlation, file = "correlation.txt", append = add, col.names = !add)

    # Save image
    save.image(file = paste("image_Year", year, "_PT.RData", sep = ""))
    # load(file = paste("image_Year", year, "_PT.RData", sep = ""))
    # sourceFunctions()
  }
}

save.image(file = paste("burnin", year, "_PT.RData", sep = ""))
# load(file = paste("burnin", year, "_PT.RData", sep = ""))
# sourceFunctions()

setwd(dir = mainDir)
