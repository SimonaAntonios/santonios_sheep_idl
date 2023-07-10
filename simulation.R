#!/usr/bin/env Rscript
#
# Long-term selection in dairy sheep on milk yield with additive & dominance
#   effects to test Inbreeding Depression Load (IDL) model with or without
#   Optimal Contribution Selection (OCS) method.
#
# Simona Antonios, Jana Obsteter, Ivan Pocrnic, Gregor Gorjanc, Silvia Rodriguez-Ramilo, Zulma Vitezica (2023)
#
# This script runs one scenario at a time defined via arguments. Make sure this
#   script is executable:
#
# chmod a+x simulation.R
#
# Usage:
#
# ./simulation.R rep [burnin|scenario] [std|stdOCS|idl|idlOCS] [restart]
#
# where:
#
# std - standard pedigree-based BLUP model & truncation selection on BV
# idl - std extended with IDL & truncation selection on BV + IDL
# stdOCS - std model & OCS on BV
# idlOCS - idl model & OCS on BV + IDL
#
# examples:
#
# ./simulation.R 1 burnin
# ./simulation.R 1 scenario std
# ./simulation.R 1 scenario stdOCS

# ---- Setup -------------------------------------------------------------------

cat("Setup\n")

# Clean the working environment (mostly for interactive use)
rm(list = ls())

args = commandArgs(trailingOnly = TRUE)
# args = c("1", "burnin")
# args = c("1", "scenario", "std")
# args = c("1", "scenario", "stdOCS")

if (length(args) < 2) {
  stop("You must provide at least 2 arguments!")
}

# Argument 1 is replicate
rep = args[1]

# Argument 2: [burnin|scenario]
burninOrScenario = args[2]
if (burninOrScenario == "burnin") {
  burnin = TRUE
} else if (burninOrScenario == "scenario") {
  burnin = FALSE
} else {
  stop("Argument 2 must be either burnin or scenario!")
}
scenarios = !burnin

# Argument 3: [std|stdOCS|idl|idlOCS]
if (burnin) {
  scenario = NA
} else {
  if (length(args) < 3) {
    stop("You must provide at least 3 arguments when running scenarios!")
  }
  scenario = args[3]
  if (!scenario %in% c("std", "stdOCS", "idl", "idlOCS")) {
    stop("Argument 3 must be either std, stdOCS, idl, or idlOCS!")
  }
}

# Argument 3 or 4 (optional): [restart]
restart = FALSE
if (burnin) {
  if (length(args) > 2) {
    tmp = args[2]
    if (tmp == "restart") {
      restart = TRUE
    } else {
      stop("Argument 3 (when running burnin) must be restart or must be missing!")
    }
  }
} else {
  if (length(args) > 3) {
    tmp = args[3]
    if (tmp == "restart") {
      restart = TRUE
    } else {
      stop("Argument 4 (when running scenarios) must be restart or must be missing!")
    }
  }
}

cat("We are doing:\n")
cat("  * replicate:", rep, "\n")
cat("  *", burninOrScenario, "\n")
if (scenarios) {
  cat("  * scenario:", scenario, "\n")
}
cat("  * restarting:", restart, "\n")

# Load packages
# install.packages(pkg = c("tidyverse", "AlphaSimR", "degreenet", "data.table"))
library(tidyverse)
library(AlphaSimR)
library(degreenet)  # for simulation of herd sizes
library(data.table) # for fast data operations (reading, writing, ...)

if (interactive()) {
  # Simona's folder
  setwd(dir = "/home/santonios/work/santonios/part2/simulation")
  # Gregor's folder
  setwd(dir = "~/Storages/GitBox/sheep_simulation_AlphaSimR/")
}
# Source required functions
sourceFunctions = function(dir = NULL) {
  fileName = "functions.R"
  if (!is.null(dir)) {
    fileName = file.path(dir, fileName)
  }
  source(file = fileName)
}
sourceFunctions() # to ensure we have the latest version (countering save.image())

# ---- Global parameters -------------------------------------------------------

cat("Global parameters\n")

nEwes               = 80000
nFemalesInLactation = 66600                                                        # no. of females in lactation
nEliteEwes          = 6000
nDamsOfFemales      = 60600
survRate            = 0.75                                                         # survival rate
NM_Fertility        = 0.9                                                          # fertility in Natural Mating (NM)
NM_prolificacy      = 1.4                                                          # prolificacy in NM
AI_Fertility        = 0.6                                                          # fertility in Artificial Insemination (AI)
AI_prolificacy      = 1.6                                                          # prolificacy in AI
nLambs              = trunc(((NM_Fertility * NM_prolificacy * survRate +
                                AI_Fertility * AI_prolificacy * survRate)) * nEwes / 2) # no. of lambs (TODO: why divide by 2? because 50% the ewes are naturally mated and the other 50 are inseminated)

nWtRams1            = 150                                                          # no. of waiting rams for progeny testing
nWtRams2            = 150                                                          # no. of waiting rams for progeny testing
nEliteSires1        = 10                                                           # no. of elite sires selected every year (for the x-th year of AI)
nEliteSires2        = 10
nEliteSires3        = 10
nSiresOfFemales1    = 55                                                           # no. of elite sires of dams selected every year (for the x-th year)
nSiresOfFemales2    = 55
nSiresOfFemales3    = 55
nNaturalMatingRams  = 1000                                                         # no. of NM sires selected every year
nEliteSireDose      = 400                                                          # no. of AI doses per elite sire
nWtRamsAIDose       = 85                                                           # no. of AI doses per wating ram
nNtlMatingDose      = 40                                                           # no. of natural matings per NM ram

# Number of lambs from AI
n1 = round(nEliteSireDose * (nEliteSires1 + nEliteSires2 + nEliteSires3) * survRate * AI_Fertility * AI_prolificacy) # 8640 progeny
n2 = round((nEwes / 2 - (nEliteSireDose * (nEliteSires1 + nEliteSires2 + nEliteSires3) + nWtRams1 * nWtRamsAIDose)) * survRate * AI_Fertility * AI_prolificacy) # 10980 progeny
n3 = round(nWtRams1 * nWtRamsAIDose * survRate * AI_Fertility * AI_prolificacy) # 9180 progeny
# Number of lambs from NM
n4 = round(nNaturalMatingRams * nNtlMatingDose * survRate * NM_Fertility * NM_prolificacy) # 37800 progeny

pDamsOfFemalesLact1 = 0.30                                                         # prop. of dams of females in lactation X
pDamsOfFemalesLact2 = 0.28
pDamsOfFemalesLact3 = 0.24
pDamsOfFemalesLact4 = 0.18
nDamsOfFemalesLact1 = round(pDamsOfFemalesLact1 * nDamsOfFemales)                  # no. of dams of females in lactation X
nDamsOfFemalesLact2 = round(pDamsOfFemalesLact2 * nDamsOfFemales)
nDamsOfFemalesLact3 = round(pDamsOfFemalesLact3 * nDamsOfFemales)
nDamsOfFemalesLact4 = round(pDamsOfFemalesLact4 * nDamsOfFemales)

pEliteEwesLact1     = 0.23                                                         # prop. of elite ewes in lactation X
pEliteEwesLact2     = 0.35
pEliteEwesLact3     = 0.26
pEliteEwesLact4     = 0.16
nEliteEwesLact1     = round(pEliteEwesLact1 * nEliteEwes)                          # no. of elite ewes in lactation X
nEliteEwesLact2     = round(pEliteEwesLact2 * nEliteEwes)
nEliteEwesLact3     = round(pEliteEwesLact3 * nEliteEwes)
nEliteEwesLact4     = round(pEliteEwesLact4 * nEliteEwes)

startYear           = 1980                                                         # start year of the simulations

nHerds              = 237                                                          # number of herds
meanHerdSize        = 320                                                          # average herd size
sdHerdSize          = 129                                                          # sd of herd size

# Genome parameters:
BaseNe              = 150                                                          # effective population size
nChr                = 26                                                           # no. of chromosomes
ChrSize             = 95 * 1e6                                                     # chr. size
nQTL                = 4500                                                         # no. of QTLs
nQTLPerChr          = round(nQTL / nChr)                                           # no. of QTLs per chromosome
nSNPPerChr          = 4000                                                         # no. of markers per chromosome
RecRate             = 1.3e-8                                                       # recombination rate
MutRate             = 1.5e-8                                                       # mutation rate

# Historical Ne over generations
histNe  = c(200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5500, 6000, 6500, 7000,  7500)
histGen = c( 10,  20,  30,  40,  50,  60,  70,  80,  90, 100,  200,  300,  400,  500, 1000, 2000, 3000, 4000, 5000, 7000, 8000, 9000, 10000)
if (FALSE) {
  AncientNe = tibble(GenerationsAgo = c(1,      histGen),
                     Ne             = c(BaseNe, histNe))

  pdf("generations_ago_with_log.pdf")
  with(AncientNe, plot(Ne ~ GenerationsAgo, type = "b", main = "Generation ago with log", log = "xy"))
  dev.off()

  pdf("generations_ago_with.pdf")
  with(AncientNe, plot(Ne ~ GenerationsAgo, type = "b", main = "Generation ago without log"))
  dev.off()
}

# Trait parameters:
nTraits             = 1                                                            # no. of traits (milk yield)
meanLac1            = 209                                                          # mean of milk yield in lactation X
meanLac2            = 213
meanLac3            = 207
meanLac4            = 180

# addVar            = 1200                                                         # additive genetic variance
addVar              = 1000                                                         # additive genetic variance
# permVar           = 500                                                          # permanent environmental variance
permVar             = 400       # 500 (Old value) - 100 = 400 (new value)          # permanent environmental variance (we started with old, 100 was targeted domVar)
resVar              = 1500                                                         # residual variance
phenVar             = 7000                                                         # specifying phenotypic variance for the trait
yearVar             = 1000                                                         # year variance
herdVar             = 1500                                                         # herd variance
herdYearVar         = 1000                                                         # herd-year variance
domVar              = addVar * 0.1
meanDD              = 0.08
varDD               = 0.30
idlVar              = domVar # TODO: how do we get this? Estimate from the data?

nBurninYears        = 20                                                           # no. of years for burnin
nScenarioYears      = 20                                                           # no. of years for scenarios

# ---- Working directory -------------------------------------------------------

if (burnin) {
  dir = paste0("rep_", rep)
  if (!restart) {
    unlink(dir, recursive = TRUE)
    dir.create(path = dir)
  }
  setwd(dir)

  dir = "burnin"
  if (!restart) {
    unlink(dir, recursive = TRUE)
    dir.create(path = dir)
  }
  setwd(dir)
}

if (scenarios) {
  dir = paste0("rep_", rep)
  setwd(dir)
  if (!dir.exists("burnin")) {
    stop("Burnin folder missing for this rep!")
  }

  dir = paste0("sce_", scenario)
  if (!restart) {
    unlink(dir, recursive = TRUE)
    dir.create(path = dir)
  }
  setwd(dir)
}

cat(paste0("Working directory: ", getwd(), "\n"))

# ---- Restart -----------------------------------------------------------------

if (FALSE) {
  # Keeping this older code here should we want to enable restarts if some
  #   scenarios will be running very long or will be crashing ... This code
  #   will need adaptation to this script though!
  if (restart) {
    if (file.exists("image_Year*.RData")) {
      system("rm image_Year*.RData")
    }
  } else {
    lastImage = system("ls -ltr image_Year*.RData | tail -n -1", intern = TRUE)
    print(lastImage)
    lastImage = tail(strsplit(lastImage, split=" ")[[1]],1)
    yearImage = as.numeric(gsub(".*?([0-9]+).*", "\\1", lastImage))
    if (yearImage < nPTyrs) { # Continue in PT period
      load(lastImage)
      runSim = "continue"
    }
  }
}

yearToDo = 1

# ---- Burnin ------------------------------------------------------------------

if (burnin) {

  cat("Burnin\n")

  # ---- Founder genomes -------------------------------------------------------

  cat("Founder genomes", as.character(Sys.time()), "\n")

  if (FALSE) {
  founderPop = runMacs2(nInd = 10 * BaseNe,
                        nChr = nChr,
                        segSites = nSNPPerChr,
                        Ne = BaseNe,
                        bp = ChrSize,
                        genLen = RecRate * ChrSize,
                        mutRate = MutRate,
                        histNe = histNe,
                        histGen = histGen)
    # Save just the founderPop
  # save(founderPop, file = "founderPop.RData")
  }
  load("founderPop.RData")
  # load("../../founder_runMacs2.RData") # note that we are in rep_x/[burnin|scenario] folder at this stage
  sourceFunctions(dir = "../..") # to ensure we have the latest version (countering save.image(), if it was used)

  # ---- AlphaSimR simulation parameters ---------------------------------------

  cat("AlphaSimR simulation parameters", as.character(Sys.time()), "\n")

  SP = SimParam$new(founderPop)
  SP$setSexes(sexes = "yes_rand")
  SP$setTrackPed(isTrackPed = TRUE)

  # Fully additive trait
  # SP$addTraitA(nQtlPerChr = nQTLPerChr, mean = trtMean, var = addVar)

  # Trait with additive and dominance effects
  # ... we got the meanDD and varDD from this function and saved the values in
  #     global parameters
  if (FALSE) {
    altAddTraitAD(nQtlPerChr = nQtlPerChr,
                  mean = trtMean,
                  varA = addVar,
                  varD = domVar,
                  limMeanDD = c(-1, 2),
                  limVarDD = c(0, 2),
                  inbrDepr = 70)
    # New trait called Trait1 was added
    # Dominance variance is 100.0
    # Inbreeding depression is 70.0
    # Used meanDD equals 0.08
    # Used varDD equals 0.28
  }
  # ... and used these parameters here
  SP$addTraitAD(nQtlPerChr = nQTLPerChr, mean = 0, var = addVar, meanDD = meanDD, varDD = varDD)
  # ... we will add lactation means later in setPhenoEwe
  # SP$setVarE() we assume that only ewes will have phenotypes, hence do not set this!!! Use setPhenoEwe() instead!

  # ---- Herds & herd and year effects -----------------------------------------

  cat("Herds & herd and year effects", as.character(Sys.time()), "\n")

  # We are oversampling so we get a sufficient number of large enough herds.
  #   With these sheep parameters, the herds were not big enough, so we divided
  #   the inputs by 10 (the simcmp worked as expected with smaller inputs) and
  #   later rescaled herd sizes
  herdSize = simcmp(n = (10 * nHerds),
                    v = c(round(meanHerdSize / 10), round(sdHerdSize / 10)))
  # hist(herdSize); mean(herdSize); meanHerdSize - mean(herdSize); sd(herdSize); sdHerdSize - sd(herdSize)
  herdSize = herdSize * 10
  # hist(herdSize); mean(herdSize); meanHerdSize - mean(herdSize); sd(herdSize); sdHerdSize - sd(herdSize)
  # 314.7, 5.3, 128.8, 0.2
  herdSize = herdSize[herdSize > 100]
  # hist(herdSize); mean(herdSize); meanHerdSize - mean(herdSize); sd(herdSize); sdHerdSize - sd(herdSize)
  # 323.2, -3.2, 122.8, 6.2
  herdSize = sample(x = herdSize, size = nHerds)
  # hist(herdSize); mean(herdSize); meanHerdSize - mean(herdSize); sd(herdSize); sdHerdSize - sd(herdSize)
  # 327.2, -7.2, 115.8, 13.2
  herdEffect = cbind(rnorm(n = nHerds, mean = 0, sd = sqrt(herdVar)))
  herdYearEffect = sampleHerdYearEffect(n = nHerds)
  herds = list(
    herd = 1:nHerds,
    herdSize = herdSize,
    herdEffect = as.matrix(herdEffect),
    herdYearEffect = as.matrix(herdYearEffect)
  )
  yearEffect = sampleYearEffect()

  # ---- Fill-in ---------------------------------------------------------------

  cat("Fill-in", as.character(Sys.time()), "\n")

  # Generate initial founders
  basePop = newPop(founderPop)

  # nCrosses = 10 * nSomething below is to ensure we don't oversample founders

  # Elite sires (AI)
  # ... in the 3rd year in service
  eliteSires3 = randCross(pop = basePop, nCrosses = 10 * nEliteSires3)
  eliteSires3 = selectWithinFam(pop = eliteSires3, nInd = 1, use = "rand", famType = "M")
  eliteSires3@sex[] = "M"
  eliteSires3@father[] = "0"
  eliteSires3@mother[] = "0"
  eliteSires3 = fillInMisc(pop = eliteSires3, year = startYear - 5.5)
  # ... in the 2nd year in service
  eliteSires2 = randCross(pop = basePop, nCrosses = 10 * nEliteSires2)
  eliteSires2 = selectWithinFam(pop = eliteSires2, nInd = 1, use = "rand", famType = "M")
  eliteSires2@sex[] = "M"
  eliteSires2@father[] = "0"
  eliteSires2@mother[] = "0"
  eliteSires2 = fillInMisc(pop = eliteSires2, year = startYear - 4.5)
  # ... in the 1st year in service
  eliteSires1 = randCross(pop = basePop, nCrosses = 10 * nEliteSires1)
  eliteSires1 = selectWithinFam(pop = eliteSires1, nInd = 1, use = "rand", famType = "M")
  eliteSires1@sex[] = "M"
  eliteSires1@father[] = "0"
  eliteSires1@mother[] = "0"
  eliteSires1 = fillInMisc(pop = eliteSires1, year = startYear - 3.5)

  eliteSires = c(eliteSires3, eliteSires2, eliteSires1)

  # Sires of females (AI)
  # ... in the 3rd year in service
  siresOfFemales3 = randCross(pop = basePop, nCrosses = 10 * nSiresOfFemales3)
  siresOfFemales3 = selectWithinFam(pop = siresOfFemales3, nInd = 1, use = "rand", famType = "M")
  siresOfFemales3@sex[] = "M"
  siresOfFemales3@father[] = "0"
  siresOfFemales3@mother[] = "0"
  siresOfFemales3 = fillInMisc(pop = siresOfFemales3, year = startYear - 5.5)
  # ... in the 2nd year in service
  siresOfFemales2 = randCross(pop = basePop, nCrosses = 10 * nSiresOfFemales2)
  siresOfFemales2 = selectWithinFam(pop = siresOfFemales2, nInd = 1, use = "rand", famType = "M")
  siresOfFemales2@sex[] = "M"
  siresOfFemales2@father[] = "0"
  siresOfFemales2@mother[] = "0"
  siresOfFemales2 = fillInMisc(pop = siresOfFemales2, year = startYear - 4.5)
  # ... in the 1st year in service
  siresOfFemales1 = randCross(pop = basePop, nCrosses = 10 * nSiresOfFemales1)
  siresOfFemales1 = selectWithinFam(pop = siresOfFemales1, nInd = 1, use = "rand", famType = "M")
  siresOfFemales1@sex[] = "M"
  siresOfFemales1@father[] = "0"
  siresOfFemales1@mother[] = "0"
  siresOfFemales1 = fillInMisc(pop = siresOfFemales1, year = startYear - 3.5)

  siresOfFemales = c(siresOfFemales3, siresOfFemales2, siresOfFemales1)

  # Waiting rams
  # ... 1.5 years old
  wtRams2  = randCross(basePop, nCrosses = 10 * nWtRams2)
  wtRams2 = selectWithinFam(pop = wtRams2, nInd = 1, use = "rand", famType = "M")
  wtRams2@sex[] = "M"
  wtRams2@father[] = "0"
  wtRams2@mother[] = "0"
  wtRams2 = fillInMisc(pop = wtRams2, year = startYear - 1.5)
  # ... 0.5 years old
  wtRams1  = randCross(basePop, nCrosses = 10 * nWtRams1)
  wtRams1 = selectWithinFam(pop = wtRams1, nInd = 1, use = "rand", famType = "M")
  wtRams1@sex[] = "M"
  wtRams1@father[] = "0"
  wtRams1@mother[] = "0"
  wtRams1 = fillInMisc(pop = wtRams1, year = startYear - 0.5)

  # Natural mating rams (NM)
  ntlMatingRams  = randCross(basePop, nCrosses = nNaturalMatingRams)
  ntlMatingRams@sex[] = "M"
  ntlMatingRams@father[] = "0"
  ntlMatingRams@mother[] = "0"
  ntlMatingRams = fillInMisc(pop = ntlMatingRams, year = startYear - 1.5)

  # Elite ewes
  # ... in the 4th lactation
  eliteEwesLact4 = randCross(basePop, nCrosses = nEliteEwesLact4)
  eliteEwesLact4@sex[] = "F"
  eliteEwesLact4@father[] = "0"
  eliteEwesLact4@mother[] = "0"
  eliteEwesLact4 = fillInMisc(pop = eliteEwesLact4, year = startYear - 4,
                              herds = herds, permEnvVar = permVar)
  # ... in the 3rd lactation
  eliteEwesLact3 = randCross(basePop, nCrosses = nEliteEwesLact3)
  eliteEwesLact3@sex[] = "F"
  eliteEwesLact3@father[] = "0"
  eliteEwesLact3@mother[] = "0"
  eliteEwesLact3 = fillInMisc(pop = eliteEwesLact3, year = startYear - 3,
                              herds = herds, permEnvVar = permVar)
  # ... in the 2nd lactation
  eliteEwesLact2 = randCross(basePop, nCrosses = nEliteEwesLact2)
  eliteEwesLact2@sex[] = "F"
  eliteEwesLact2@father[] = "0"
  eliteEwesLact2@mother[] = "0"
  eliteEwesLact2 = fillInMisc(pop = eliteEwesLact2, year = startYear - 2,
                              herds = herds, permEnvVar = permVar)
  # ... in the 1st lactation
  eliteEwesLact1 = randCross(basePop, nCrosses =  nEliteEwesLact1)
  eliteEwesLact1@sex[] = "F"
  eliteEwesLact1@father[] = "0"
  eliteEwesLact1@mother[] = "0"
  eliteEwesLact1 = fillInMisc(pop = eliteEwesLact1, year = startYear - 1,
                              herds = herds, permEnvVar = permVar)

  # TODO: is this OK?
  eliteEwes = c(eliteEwesLact4, eliteEwesLact3, eliteEwesLact2, eliteEwesLact1)

  # Dams of females
  # ... in the 4th lactation
  damsOfFemalesLact4 = randCross(basePop, nCrosses = nDamsOfFemalesLact4)
  damsOfFemalesLact4@sex[] = "F"
  damsOfFemalesLact4@father[] = "0"
  damsOfFemalesLact4@mother[] = "0"
  damsOfFemalesLact4 = fillInMisc(pop = damsOfFemalesLact4, year = startYear - 4,
                                  herds = herds, permEnvVar = permVar)
  # ... in the 3rd lactation
  damsOfFemalesLact3 = randCross(basePop, nCrosses = nDamsOfFemalesLact3)
  damsOfFemalesLact3@sex[] = "F"
  damsOfFemalesLact3@father[] = "0"
  damsOfFemalesLact3@mother[] = "0"
  damsOfFemalesLact3 = fillInMisc(pop = damsOfFemalesLact3, year = startYear - 3,
                                  herds = herds, permEnvVar = permVar)
  # ... in the 2nd lactation
  damsOfFemalesLact2 = randCross(basePop, nCrosses =  nDamsOfFemalesLact2)
  damsOfFemalesLact2@sex[] = "F"
  damsOfFemalesLact2@father[] = "0"
  damsOfFemalesLact2@mother[] = "0"
  damsOfFemalesLact2 = fillInMisc(pop = damsOfFemalesLact2, year = startYear - 2,
                                  herds = herds, permEnvVar = permVar)
  # ... in the 1st lactation
  damsOfFemalesLact1 = randCross(basePop, nCrosses =  nDamsOfFemalesLact1)
  damsOfFemalesLact1@sex[] = "F"
  damsOfFemalesLact1@father[] = "0"
  damsOfFemalesLact1@mother[] = "0"
  damsOfFemalesLact1 = fillInMisc(pop = damsOfFemalesLact1, year = startYear - 1,
                                  herds = herds, permEnvVar = permVar)

  damsOfFemales = c(damsOfFemalesLact4, damsOfFemalesLact3, damsOfFemalesLact2, damsOfFemalesLact1)

  # Lambs
  matingPlan1 = cbind(eliteEwes@id,
                      sample(eliteSires@id, size = nEliteEwes, replace = TRUE))

  n = n1 - nEliteEwes
  damsOfFemalesId = damsOfFemales@id
  damsOfFemalesIdForElite = sample(damsOfFemalesId, size = n)
  matingPlan2 = cbind(damsOfFemalesIdForElite,
                      sample(eliteSires@id, size = n, replace = TRUE))

  n = nDamsOfFemales - (n1 - nEliteEwes)
  damsOfFemalesIdForRest = damsOfFemalesId[!damsOfFemalesId %in% damsOfFemalesIdForElite]
  # matingPlan3 = cbind(damsOfFemalesIdForRest,
  #                     sample(c(siresOfFemales@id, wtRams1@id, ntlMatingRams@id), size = n, replace = TRUE))
  # TODO: IS it good like this?
  matingPlan3 = cbind(damsOfFemalesIdForRest,
                      sample(c(siresOfFemales@id, size = n2, replace = TRUE),
                      sample(wtRams1@id, size = n3, replace = TRUE),
                      sample(ntlMatingRams@id, size = n4, replace = TRUE)))
  
  matingPlan = rbind(matingPlan1, matingPlan2, matingPlan3)
  lambs = makeCross2(females = c(eliteEwes, damsOfFemales),
                     males = c(eliteSires, siresOfFemales, wtRams1, ntlMatingRams),
                     crossPlan = matingPlan)
  lambs = fillInMisc(pop = lambs,
                     mothers = c(eliteEwes, damsOfFemales),
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
  # database = recordData(database, pop = eliteEwes, year = yearFull, lactation = ???) # can not save like this since lactations vary
  database = recordData(database, pop = damsOfFemalesLact4, year = startYear, lactation = 4)
  database = recordData(database, pop = damsOfFemalesLact3, year = startYear, lactation = 3)
  database = recordData(database, pop = damsOfFemalesLact2, year = startYear, lactation = 2)
  database = recordData(database, pop = damsOfFemalesLact1, year = startYear, lactation = 1)
  # database = recordData(database, pop = damsOfFemales, year = yearFull, lactation = ???) # can not save like this since lactations vary
  database = recordData(database, pop = lambs, year = startYear)

  # Save workspace image
  # save.image(file = "fillin.RData")
  # load(file = "fillin.RData")
  # sourceFunctions(dir = "../..") # to ensure we have the latest version (countering save.image())

  # ----- Loop breeding programme over years -----------------------------------

  cat("Loop breeding programme over years", as.character(Sys.time()), "\n")

  for (year in yearToDo:nBurninYears) {
    # year = yearToDo
    yearFull = startYear + year

    cat("Year", year, " (", yearFull, ") ", as.character(Sys.time()), "\n")

    yearEffect = sampleYearEffect()
    herds$herdYearEffect = sampleHerdYearEffect(n = nHerds)

    if (year <= 1) {
      use = "rand"
    } else {
      use= "ebv"
    }

    # ---- Phenotyping ----

    cat("Phenotyping", as.character(Sys.time()), "\n")

    eliteEwesLact4 = setPhenoEwe(eliteEwesLact4, varE = resVar,
                                 mean = meanLac4, yearEffect = yearEffect, herds = herds)

    eliteEwesLact3 = setPhenoEwe(eliteEwesLact3, varE = resVar,
                                 mean = meanLac3, yearEffect = yearEffect, herds = herds)

    eliteEwesLact2 = setPhenoEwe(eliteEwesLact2, varE = resVar,
                                 mean = meanLac2, yearEffect = yearEffect, herds = herds)

    eliteEwesLact1 = setPhenoEwe(eliteEwesLact1, varE = resVar,
                                 mean = meanLac1, yearEffect = yearEffect, herds = herds)

    damsOfFemalesLact4 = setPhenoEwe(damsOfFemalesLact4, varE = resVar,
                                     mean = meanLac4, yearEffect = yearEffect, herds = herds)

    damsOfFemalesLact3 = setPhenoEwe(damsOfFemalesLact3, varE = resVar,
                                     mean = meanLac3, yearEffect = yearEffect, herds = herds)

    damsOfFemalesLact2 = setPhenoEwe(damsOfFemalesLact2, varE = resVar,
                                     mean = meanLac2, yearEffect = yearEffect, herds = herds)

    damsOfFemalesLact1 = setPhenoEwe(damsOfFemalesLact1, varE = resVar,
                                     mean = meanLac1, yearEffect = yearEffect, herds = herds)

    database = setDatabasePheno(database, pop = eliteEwesLact4)
    database = setDatabasePheno(database, pop = eliteEwesLact3)
    database = setDatabasePheno(database, pop = eliteEwesLact2)
    database = setDatabasePheno(database, pop = eliteEwesLact1)
    database = setDatabasePheno(database, pop = damsOfFemalesLact4)
    database = setDatabasePheno(database, pop = damsOfFemalesLact3)
    database = setDatabasePheno(database, pop = damsOfFemalesLact2)
    database = setDatabasePheno(database, pop = damsOfFemalesLact1)

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
    
    pedEbv = estimateBreedingValues(pedigree = SP$pedigree,
                                    database = database,
                                    vars = variances,
                                    removeFromEvaluation = removeCulledMaleLambs)
    
    # Set EBVs for every population
    eliteSires3 = setEbv(eliteSires3, ebv = pedEbv)
    eliteSires2 = setEbv(eliteSires2, ebv = pedEbv)
    eliteSires1 = setEbv(eliteSires1, ebv = pedEbv)
    eliteSires  = setEbv(eliteSires, ebv = pedEbv)
    siresOfFemales3 = setEbv(siresOfFemales3, ebv = pedEbv)
    siresOfFemales2 = setEbv(siresOfFemales2, ebv = pedEbv)
    siresOfFemales1 = setEbv(siresOfFemales1, ebv = pedEbv)
    siresOfFemales = setEbv(siresOfFemales, ebv = pedEbv)
    wtRams2 = setEbv(wtRams2, ebv = pedEbv)
    wtRams1 = setEbv(wtRams1, ebv = pedEbv)
    ntlMatingRams = setEbv(ntlMatingRams, ebv = pedEbv)
    eliteEwesLact4 = setEbv(eliteEwesLact4, ebv = pedEbv)
    eliteEwesLact3 = setEbv(eliteEwesLact3, ebv = pedEbv)
    eliteEwesLact2 = setEbv(eliteEwesLact2, ebv = pedEbv)
    eliteEwesLact1 = setEbv(eliteEwesLact1, ebv = pedEbv)
    eliteEwes = setEbv(eliteEwes, ebv = pedEbv)
    damsOfFemalesLact4 = setEbv(damsOfFemalesLact4, ebv = pedEbv)
    damsOfFemalesLact3 = setEbv(damsOfFemalesLact3, ebv = pedEbv)
    damsOfFemalesLact2 = setEbv(damsOfFemalesLact2, ebv = pedEbv)
    damsOfFemalesLact1 = setEbv(damsOfFemalesLact1, ebv = pedEbv)
    damsOfFemales = setEbv(damsOfFemales, ebv = pedEbv)
    # ... we could use setEbv(), but we have excluded some culled lambs from the
    #     evaluation so we will just calculate parent average for all lambs - this
    #     works for pedigree BLUP, but not for genomic BLUP (if some lambs would
    #     have been genotyped)
    # lambs             = setEbv(lambs, ebv = pedEbv)
    selM = match(x = lambs@mother, table = pedEbv$IId)
    selF = match(x = lambs@father, table = pedEbv$IId)
    lambs@ebv = as.matrix((pedEbv[selM, -1] + pedEbv[selF, -1]) / 2)
    
    # Save current EBVs into the database
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
    # database = setDatabaseEbv(database, pop = eliteEwes) # we didn't save this pop in the database
    database = setDatabaseEbv(database, pop = damsOfFemalesLact4)
    database = setDatabaseEbv(database, pop = damsOfFemalesLact3)
    database = setDatabaseEbv(database, pop = damsOfFemalesLact2)
    database = setDatabaseEbv(database, pop = damsOfFemalesLact1)
    # database = setDatabaseEbv(database, pop = damsOfFemales) # we didn't save this pop in the database
    database = setDatabaseEbv(database, pop = lambs)
    
    # TODO: if we have done it above then we don't need to do it here again
    # eliteEwes = c(eliteEwesLact1, eliteEwesLact2, eliteEwesLact3, eliteEwesLact4)
    # damsOfFemales = c(damsOfFemalesLact1, damsOfFemalesLact2, damsOfFemalesLact3, damsOfFemalesLact4)
    
    # ---- Summarising results ----
    
    cat("... Summarising results", as.character(Sys.time()), "\n")
    
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
                             damsOfFemalesLact4 = ebvAccuracy(damsOfFemalesLact4),
                             damsOfFemalesLact3 = ebvAccuracy(damsOfFemalesLact3),
                             damsOfFemalesLact2 = ebvAccuracy(damsOfFemalesLact2),
                             damsOfFemalesLact1 = ebvAccuracy(damsOfFemalesLact1),
                             damsOfFemales = ebvAccuracy(damsOfFemales),
                             lambs = ebvAccuracy(lambs))
    add = ifelse(year == 1, FALSE, TRUE)
    write.table(x = correlation, file = "ebvAccuracy.txt", append = add, col.names = !add)
    
    # ---- Rams ----

    cat("Rams\n")

    # ---- ... Elite sires ----

    cat("... Elite sires", as.character(Sys.time()), "\n")

    eliteSires3 = eliteSires2 # eliteSires3 are 4.5 years old here
    eliteSires2 = eliteSires1 # eliteSires2 are 3.5 years old here

    if (all(wtRams2@father == "0")) {
      eliteSires1 = selectInd(pop = wtRams2, nInd = nEliteSires1, # eliteSires1 are 3.5 years old here
                              use = use)
    } else {
      eliteSires1 = selectWithinFam(pop = wtRams2, nInd = 1, # eliteSires1 are 3.5 years old here
                                    use = use, famType = "M")
      eliteSires1 = selectInd(pop = eliteSires1, nInd = nEliteSires1,
                              use = use)
    }

    # Don't create eliteSires here!

    # ---- ... Sires of females ----

    cat("... Sires of females", as.character(Sys.time()), "\n")

    siresOfFemales3 = siresOfFemales2 # siresOfFemales3 are 5.5 years old here
    siresOfFemales2 = siresOfFemales1 # siresOfFemales2 are 4.5 years old here

    if (all(wtRams2@father == "0")) {
      sel = !wtRams2@id %in% eliteSires1@id
      siresOfFemales1 = selectInd(pop = wtRams2[sel], nInd = nSiresOfFemales1, # siresOfFemales1 are 3.5 years old here
                                  use = use)
    } else {
      sel = !wtRams2@id %in% eliteSires1@id
      siresOfFemales1 = selectWithinFam(pop = wtRams2[sel], nInd = 2, # siresOfFemales1 are 3.5 years old here
                                        use = use, famType = "M")
      siresOfFemales1 = selectInd(pop = siresOfFemales1, nInd = nSiresOfFemales1,
                                  use = use)
    }

    # Note that we select wtRams only from lambs from elite rams (not waiting rams) and elite ewes
    selWtRams = lambs@father %in% eliteSires@id & lambs@mother %in% eliteEwes@id
    n = ceiling(nWtRams1 / length(eliteSires@id))
    wtRams2 = wtRams1
    wtRams1 = selectWithinFam(pop = lambs[selWtRams], nInd = n,
                              use = use, famType = "M", sex = "M")
    wtRams1 = selectInd(pop = wtRams1, nInd = nWtRams1,
                        use = use)

    selNtlRams = lambs@father %in% c(eliteSires@id, siresOfFemales@id) &
      !lambs@id %in% wtRams1@id
    ntlMatingRams =  selectInd(pop = lambs[selNtlRams], nInd = nNaturalMatingRams,
                               use = use, famType = "M", sex = "M")

    # Have to do it here since we need to keep track of previous eliteSires in
    #   the above steps
    eliteSires = c(eliteSires3, eliteSires2, eliteSires1)
    siresOfFemales = c(siresOfFemales3, siresOfFemales2, siresOfFemales1)

    # ---- Ewes ----

    cat("Ewes\n")

    # ---- ... Elite ewes ----

    cat("... Elite ewes", as.character(Sys.time()), "\n")

    ewesLact4 = c(damsOfFemalesLact4, eliteEwesLact4)
    ewesLact3 = c(damsOfFemalesLact3, eliteEwesLact3)
    ewesLact2 = c(damsOfFemalesLact2, eliteEwesLact2)
    ewesLact1 = c(damsOfFemalesLact1, eliteEwesLact1)

    # TODO: Selecting across all ewes, but should we focus on daughters of AI sires only?
    #       Hence selectInd(eliteEwesLact3, ...)?
    #  yes the elite ewes are the daughters of AI sires, but we were saying that since we are choosing them based 
    #  on their EBVs the majority have to be from elite sire no? or better to add this criterion (daughter of AI sires)? 
    eliteEwesLact4 = selectInd(ewesLact3, nInd = nEliteEwesLact4, use = use) # eliteEwesLact4 are 4 years old here
    eliteEwesLact3 = selectInd(ewesLact2, nInd = nEliteEwesLact3, use = use) # eliteEwesLact3 are 3 years old here
    eliteEwesLact2 = selectInd(ewesLact1, nInd = nEliteEwesLact2, use = use) # eliteEwesLact2 are 2 years old here
    eliteEwesLact1 = selectInd(lambs[selNtlRams], nInd = nEliteEwesLact1, use = use, sex = "F", famType = "M") # eliteEwesLact1 are 1 years old here

    # Set phenotypes to missing, because these are copied from the previous lactation.
    #   We do this because of the recordData() call below - that would save wrong
    #   phenotypes for these animals in this life stage! Correct phenotypes for these
    #   animals in this life stage will be added at the beginning of the next year.
    eliteEwesLact4@pheno[,] = NA
    eliteEwesLact3@pheno[,] = NA
    eliteEwesLact2@pheno[,] = NA
    eliteEwesLact1@pheno[,] = NA

    eliteEwes = c(eliteEwesLact4, eliteEwesLact3, eliteEwesLact2, eliteEwesLact1)

    # ---- ... Dams of females ----

    cat("... Dams of females", as.character(Sys.time()), "\n")

    damsOfFemalesLact4 = selectInd(ewesLact3[!ewesLact3@id %in% eliteEwesLact4@id],
                                  nInd = nDamsOfFemalesLact4, use = "rand") # damsOfFemalesLact4 are 4 years old here

    damsOfFemalesLact3 = selectInd(ewesLact2[!ewesLact2@id %in% eliteEwesLact3@id],
                                  nInd = nDamsOfFemalesLact3, use = "rand") # damsOfFemalesLact3 are 3 years old here

    damsOfFemalesLact2 = selectInd(ewesLact1[!ewesLact1@id %in% eliteEwesLact2@id],
                                  nInd = nDamsOfFemalesLact2, use = "rand") # damsOfFemalesLact2 are 2 years old here

    damsOfFemalesLact1 = selectInd(pop = lambs[!lambs@id %in% eliteEwesLact1@id],
                                  nInd = nDamsOfFemalesLact1, use = "rand", sex = "F") # damsOfFemalesLact1 are 1 years old here

    # TODO: we should revise these ages above - are damsOfFemalesLact1 really 1 year old here?

    # Set phenotypes to missing, because these are copied from the previous lactation.
    #   We do this because of the recordData() call below - that would save wrong
    #   phenotypes for these animals in this life stage! Correct phenotypes for these
    #   animals in this life stage will be added at the beginning of the next year.
    damsOfFemalesLact4@pheno[,] = NA
    damsOfFemalesLact3@pheno[,] = NA
    damsOfFemalesLact2@pheno[,] = NA
    damsOfFemalesLact1@pheno[,] = NA

    damsOfFemales = c(damsOfFemalesLact4, damsOfFemalesLact3, damsOfFemalesLact2, damsOfFemalesLact1)

    # ---- Lambs ----

    cat("Lambs", as.character(Sys.time()), "\n")

    matingPlan1 = cbind(eliteEwes@id,
                        sample(eliteSires@id, size = nEliteEwes, replace = TRUE))

    n = n1 - nEliteEwes
    damsOfFemalesId = damsOfFemales@id
    damsOfFemalesIdForElite = sample(damsOfFemalesId, size = n)
    matingPlan2 = cbind(damsOfFemalesIdForElite,
                        sample(eliteSires@id, size = n, replace = TRUE))

    n = nDamsOfFemales - (n1 - nEliteEwes)
    damsOfFemalesIdForRest = damsOfFemalesId[!damsOfFemalesId %in% damsOfFemalesIdForElite]
    # TODO: Is is good like this?
    matingPlan3 = cbind(damsOfFemalesIdForRest,
                        sample(c(siresOfFemales@id, size = n2, replace = TRUE),
                        sample(wtRams1@id, size = n3, replace = TRUE),
                        sample(ntlMatingRams@id, size = n4, replace = TRUE)))
    
    matingPlan = rbind(matingPlan1, matingPlan2, matingPlan3)
    lambs = makeCross2(females = c(eliteEwes, damsOfFemales),
                       males = c(eliteSires, siresOfFemales, wtRams1, ntlMatingRams),
                       crossPlan = matingPlan)
    lambs = fillInMisc(pop = lambs,
                       mothers = c(eliteEwes, damsOfFemales),
                       permEnvVar = permVar, year = yearFull)

    # ---- Data recording & estimating breeding values ----

    cat("... Data recording & estimating breeding values", as.character(Sys.time()), "\n")

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
    # database = recordData(database, pop = eliteEwes, year = yearFull, lactation = ???) # can not save like this since lactations vary
    database = recordData(database, pop = damsOfFemalesLact4, year = yearFull, lactation = 4)
    database = recordData(database, pop = damsOfFemalesLact3, year = yearFull, lactation = 3)
    database = recordData(database, pop = damsOfFemalesLact2, year = yearFull, lactation = 2)
    database = recordData(database, pop = damsOfFemalesLact1, year = yearFull, lactation = 1)
    # database = recordData(database, pop = damsOfFemales, year = yearFull, lactation = ???) # can not save like this since lactations vary
    database = recordData(database, pop = lambs, year = yearFull)

    # Save workspace image
    # fileName = paste("burnin_year", year, ".RData", sep = "")
    # save.image(file = fileName)
    # load(file = fileName)
    # sourceFunctions(dir = "../..") # to ensure we have the latest version (countering save.image())
  }

  # Save workspace image
  save.image(file = "burnin.RData")
  # load(file = "burnin.RData")
  # sourceFunctions(dir = "../..") # to ensure we have the latest version (countering save.image())
}

# ---- Scenarios ---------------------------------------------------------------

if (scenarios) {

  cat("Scenarios\n")

  fileName = "../burnin/burnin.RData"
  if (!file.exists(fileName)) {
    stop("Burnin data does not exist for this rep!")
  }
  load(file = fileName)
  sourceFunctions(dir = "../..") # to ensure we have the latest version (countering save.image())

  # ---- Loop breeding programme over years ------------------------------------

  cat("Loop breeding programme over years", as.character(Sys.time()), "\n")

  for (year in (nBurninYears + 1):nScenarioYears) {
    # year = yearToDo
    yearFull = startYear + year

    cat("Year", year, " (", yearFull, ") ", as.character(Sys.time()), "\n")

    yearEffect = sampleYearEffect()
    herds$herdYearEffect = sampleHerdYearEffect(n = nHerds)

    # ---- Phenotyping ----

    cat("Phenotyping", as.character(Sys.time()), "\n")

    eliteEwesLact4 = setPhenoEwe(eliteEwesLact4, varE = resVar,
                                 mean = meanLac4, yearEffect = yearEffect, herds = herds)

    eliteEwesLact3 = setPhenoEwe(eliteEwesLact3, varE = resVar,
                                 mean = meanLac3, yearEffect = yearEffect, herds = herds)

    eliteEwesLact2 = setPhenoEwe(eliteEwesLact2, varE = resVar,
                                 mean = meanLac2, yearEffect = yearEffect, herds = herds)

    eliteEwesLact1 = setPhenoEwe(eliteEwesLact1, varE = resVar,
                                 mean = meanLac1, yearEffect = yearEffect, herds = herds)

    damsOfFemalesLact4 = setPhenoEwe(damsOfFemalesLact4, varE = resVar,
                                     mean = meanLac4, yearEffect = yearEffect, herds = herds)

    damsOfFemalesLact3 = setPhenoEwe(damsOfFemalesLact3, varE = resVar,
                                     mean = meanLac3, yearEffect = yearEffect, herds = herds)

    damsOfFemalesLact2 = setPhenoEwe(damsOfFemalesLact2, varE = resVar,
                                     mean = meanLac2, yearEffect = yearEffect, herds = herds)

    damsOfFemalesLact1 = setPhenoEwe(damsOfFemalesLact1, varE = resVar,
                                     mean = meanLac1, yearEffect = yearEffect, herds = herds)

    database = setDatabasePheno(database, pop = eliteEwesLact4)
    database = setDatabasePheno(database, pop = eliteEwesLact3)
    database = setDatabasePheno(database, pop = eliteEwesLact2)
    database = setDatabasePheno(database, pop = eliteEwesLact1)
    database = setDatabasePheno(database, pop = damsOfFemalesLact4)
    database = setDatabasePheno(database, pop = damsOfFemalesLact3)
    database = setDatabasePheno(database, pop = damsOfFemalesLact2)
    database = setDatabasePheno(database, pop = damsOfFemalesLact1)
    
    # We will remove male lambs from the evaluation - those that will never
    #   contribute to next generations
    sel = database$General$Pop == "lambs" & database$General$Sex == "M"
    maleLambs = database$General[sel, "IId"]
    sel = database$General$Pop %in% c("wtRams1", "ntlMatingRams")
    reproMales = database$General[sel, "IId"]
    culledMaleLambs = maleLambs$IId[!maleLambs$IId %in% reproMales$IId]
    removeCulledMaleLambs = row.names(SP$pedigree) %in% culledMaleLambs
    # sum(removeCulledMaleLambs); sum(!removeCulledMaleLambs)
    
    if (scenario %in% c("std", "stdOCS")) {
      variances = list(varPE = permVar,
                       varA  = addVar,
                       varE  = resVar)
    } else if (scenario %in% c("idl", "idlOCS")) {
      variances = list(varPE = permVar,
                       varA  = addVar,
                       varIDL = idlVar,
                       varE  = resVar)
    }
    
    pedEbv = estimateBreedingValues(pedigree = SP$pedigree,
                                    database = database,
                                    vars = variances,
                                    removeFromEvaluation = removeCulledMaleLambs)
    
    # Set EBVs for every population
    eliteSires3 = setEbv(eliteSires3, ebv = pedEbv)
    eliteSires2 = setEbv(eliteSires2, ebv = pedEbv)
    eliteSires1 = setEbv(eliteSires1, ebv = pedEbv)
    eliteSires  = setEbv(eliteSires, ebv = pedEbv)
    siresOfFemales3 = setEbv(siresOfFemales3, ebv = pedEbv)
    siresOfFemales2 = setEbv(siresOfFemales2, ebv = pedEbv)
    siresOfFemales1 = setEbv(siresOfFemales1, ebv = pedEbv)
    siresOfFemales = setEbv(siresOfFemales, ebv = pedEbv)
    wtRams2 = setEbv(wtRams2, ebv = pedEbv)
    wtRams1 = setEbv(wtRams1, ebv = pedEbv)
    ntlMatingRams = setEbv(ntlMatingRams, ebv = pedEbv)
    eliteEwesLact4 = setEbv(eliteEwesLact4, ebv = pedEbv)
    eliteEwesLact3 = setEbv(eliteEwesLact3, ebv = pedEbv)
    eliteEwesLact2 = setEbv(eliteEwesLact2, ebv = pedEbv)
    eliteEwesLact1 = setEbv(eliteEwesLact1, ebv = pedEbv)
    eliteEwes = setEbv(eliteEwes, ebv = pedEbv)
    damsOfFemalesLact4 = setEbv(damsOfFemalesLact4, ebv = pedEbv)
    damsOfFemalesLact3 = setEbv(damsOfFemalesLact3, ebv = pedEbv)
    damsOfFemalesLact2 = setEbv(damsOfFemalesLact2, ebv = pedEbv)
    damsOfFemalesLact1 = setEbv(damsOfFemalesLact1, ebv = pedEbv)
    damsOfFemales = setEbv(damsOfFemales, ebv = pedEbv)
    # ... we could use setEbv(), but we have excluded some culled lambs from the
    #     evaluation so we will just calculate parent average for all lambs - this
    #     works for pedigree BLUP, but not for genomic BLUP (if some lambs would
    #     have been genotyped)
    # lambs             = setEbv(lambs, ebv = pedEbv)
    selM = match(x = lambs@mother, table = pedEbv$IId)
    selF = match(x = lambs@father, table = pedEbv$IId)
    lambs@ebv = as.matrix((pedEbv[selM, -1] + pedEbv[selF, -1]) / 2)
    
    # Save current EBVs into the database
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
    # database = setDatabaseEbv(database, pop = eliteEwes) # we didn't save this pop in the database
    database = setDatabaseEbv(database, pop = damsOfFemalesLact4)
    database = setDatabaseEbv(database, pop = damsOfFemalesLact3)
    database = setDatabaseEbv(database, pop = damsOfFemalesLact2)
    database = setDatabaseEbv(database, pop = damsOfFemalesLact1)
    # database = setDatabaseEbv(database, pop = damsOfFemales) # we didn't save this pop in the database
    database = setDatabaseEbv(database, pop = lambs)
    
    # TODO: if we have done it above then we don't need to do it here again
    # TODO: answer: we do not need it to update the values of the EBV for the accuracies below?
    # eliteEwes = c(eliteEwesLact1, eliteEwesLact2, eliteEwesLact3, eliteEwesLact4)
    # damsOfFemales = c(damsOfFemalesLact1, damsOfFemalesLact2, damsOfFemalesLact3, damsOfFemalesLact4)
    
    # ---- Summarising results ----
    
    cat("... Summarising results", as.character(Sys.time()), "\n")
    
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
                             damsOfFemalesLact4 = ebvAccuracy(damsOfFemalesLact4),
                             damsOfFemalesLact3 = ebvAccuracy(damsOfFemalesLact3),
                             damsOfFemalesLact2 = ebvAccuracy(damsOfFemalesLact2),
                             damsOfFemalesLact1 = ebvAccuracy(damsOfFemalesLact1),
                             damsOfFemales = ebvAccuracy(damsOfFemales),
                             lambs = ebvAccuracy(lambs))
    write.table(x = correlation, file = "ebvAccuracy.txt", append = TRUE, col.names = FALSE)
    
    # ---- Rams ----

    cat("Rams\n")

    # ---- ... Elite sires ----

    cat("... Elite sires", as.character(Sys.time()), "\n")

    eliteSires3 = eliteSires2 # eliteSires3 are 4.5 years old here
    eliteSires2 = eliteSires1 # eliteSires2 are 3.5 years old here

    if (all(wtRams2@father == "0")) {
      eliteSires1 = selectInd(pop = wtRams2, nInd = nEliteSires1, # eliteSires1 are 3.5 years old here
                              use = use)
    } else {
      eliteSires1 = selectWithinFam(pop = wtRams2, nInd = 1, # eliteSires1 are 3.5 years old here
                                    use = use, famType = "M")
      eliteSires1 = selectInd(pop = eliteSires1, nInd = nEliteSires1,
                              use = use)
    }

    # Don't create eliteSires here!
    # Why? 

    # ---- ... Sires of females ----

    cat("... Sires of females", as.character(Sys.time()), "\n")

    siresOfFemales3 = siresOfFemales2 # siresOfFemales3 are 5.5 years old here
    siresOfFemales2 = siresOfFemales1 # siresOfFemales2 are 4.5 years old here

    if (all(wtRams2@father == "0")) {
      sel = !wtRams2@id %in% eliteSires1@id
      siresOfFemales1 = selectInd(pop = wtRams2[sel], nInd = nSiresOfFemales1, # siresOfFemales1 are 3.5 years old here
                                  use = use)
    } else {
      sel = !wtRams2@id %in% eliteSires1@id
      siresOfFemales1 = selectWithinFam(pop = wtRams2[sel], nInd = 2, # siresOfFemales1 are 3.5 years old here
                                        use = use, famType = "M")
      siresOfFemales1 = selectInd(pop = siresOfFemales1, nInd = nSiresOfFemales1,
                                  use = use)
    }

    # Note that we select wtRams only from lambs from elite rams (not waiting rams) and elite ewes
    selWtRams = lambs@father %in% eliteSires@id & lambs@mother %in% eliteEwes@id
    n = ceiling(nWtRams1 / length(eliteSires@id))
    wtRams2 = wtRams1
    wtRams1 = selectWithinFam(pop = lambs[selWtRams], nInd = n,
                              use = use, famType = "M", sex = "M")
    wtRams1 = selectInd(pop = wtRams1, nInd = nWtRams1,
                        use = use)

    selNtlRams = lambs@father %in% c(eliteSires@id, siresOfFemales@id) &
      !lambs@id %in% wtRams1@id
    ntlMatingRams =  selectInd(pop = lambs[selNtlRams], nInd = nNaturalMatingRams,
                               use = use, famType = "M", sex = "M")

    # Have to do it here since we need to keep track of previous eliteSires in
    #   the above steps
    eliteSires = c(eliteSires3, eliteSires2, eliteSires1)
    siresOfFemales = c(siresOfFemales3, siresOfFemales2, siresOfFemales1)

    # ---- Ewes ----

    cat("Ewes\n")

    # ---- ... Elite ewes ----

    cat("... Elite ewes", as.character(Sys.time()), "\n")

    ewesLact4 = c(damsOfFemalesLact4, eliteEwesLact4)
    ewesLact3 = c(damsOfFemalesLact3, eliteEwesLact3)
    ewesLact2 = c(damsOfFemalesLact2, eliteEwesLact2)
    ewesLact1 = c(damsOfFemalesLact1, eliteEwesLact1)

    # TODO: Selecting across all ewes, but should we focus on daughters of AI sires only?
    #       Hence selectInd(eliteEwesLact3, ...)?
    #  yes the elite ewes are the daughters of AI sires, but we were saying that since we are choosing them based 
    #  on their EBVs the majority have to be from elite sire no? or better to add this criterion (daughter of AI sires)? 
    eliteEwesLact4 = selectInd(ewesLact3, nInd = nEliteEwesLact4, use = use) # eliteEwesLact4 are 4 years old here
    eliteEwesLact3 = selectInd(ewesLact2, nInd = nEliteEwesLact3, use = use) # eliteEwesLact3 are 3 years old here
    eliteEwesLact2 = selectInd(ewesLact1, nInd = nEliteEwesLact2, use = use) # eliteEwesLact2 are 2 years old here
    eliteEwesLact1 = selectInd(lambs[selNtlRams], nInd = nEliteEwesLact1, use = use, sex = "F", famType = "M") # eliteEwesLact1 are 1 years old here

    # Set phenotypes to missing, because these are copied from the previous lactation.
    #   We do this because of the recordData() call below - that would save wrong
    #   phenotypes for these animals in this life stage! Correct phenotypes for these
    #   animals in this life stage will be added at the beginning of the next year.
    eliteEwesLact4@pheno[,] = NA
    eliteEwesLact3@pheno[,] = NA
    eliteEwesLact2@pheno[,] = NA
    eliteEwesLact1@pheno[,] = NA

    eliteEwes = c(eliteEwesLact4, eliteEwesLact3, eliteEwesLact2, eliteEwesLact1)

    # ---- ... Dams of females ----

    cat("... Dams of females", as.character(Sys.time()), "\n")

    damsOfFemalesLact4 = selectInd(ewesLact3[!ewesLact3@id %in% eliteEwesLact4@id],
                                   nInd = nDamsOfFemalesLact4, use = "rand") # damsOfFemalesLact4 are 4 years old here

    damsOfFemalesLact3 = selectInd(ewesLact2[!ewesLact2@id %in% eliteEwesLact3@id],
                                   nInd = nDamsOfFemalesLact3, use = "rand") # damsOfFemalesLact3 are 3 years old here

    damsOfFemalesLact2 = selectInd(ewesLact1[!ewesLact1@id %in% eliteEwesLact2@id],
                                   nInd = nDamsOfFemalesLact2, use = "rand") # damsOfFemalesLact2 are 2 years old here

    damsOfFemalesLact1 = selectInd(pop = lambs[!lambs@id %in% eliteEwesLact1@id],
                                   nInd = nDamsOfFemalesLact1, use = "rand", sex = "F") # damsOfFemalesLact1 are 1 years old here

    # TODO: we should revise these ages above - are damsOfFemalesLact1 really 1 year old here?

    # Set phenotypes to missing, because these are copied from the previous lactation.
    #   We do this because of the recordData() call below - that would save wrong
    #   phenotypes for these animals in this life stage! Correct phenotypes for these
    #   animals in this life stage will be added at the beginning of the next year.
    damsOfFemalesLact4@pheno[,] = NA
    damsOfFemalesLact3@pheno[,] = NA
    damsOfFemalesLact2@pheno[,] = NA
    damsOfFemalesLact1@pheno[,] = NA

    damsOfFemales = c(damsOfFemalesLact4, damsOfFemalesLact3, damsOfFemalesLact2, damsOfFemalesLact1)

    # ---- Lambs ----

    cat("Lambs", as.character(Sys.time()), "\n")

    matingPlan1 = cbind(eliteEwes@id,
                        sample(eliteSires@id, size = nEliteEwes, replace = TRUE))

    n = n1 - nEliteEwes
    damsOfFemalesId = damsOfFemales@id
    damsOfFemalesIdForElite = sample(damsOfFemalesId, size = n)
    matingPlan2 = cbind(damsOfFemalesIdForElite,
                        sample(eliteSires@id, size = n, replace = TRUE))

    n = nDamsOfFemales - (n1 - nEliteEwes)
    damsOfFemalesIdForRest = damsOfFemalesId[!damsOfFemalesId %in% damsOfFemalesIdForElite]
    # matingPlan3 = cbind(damsOfFemalesIdForRest,
    #                     sample(c(siresOfFemales@id, wtRams1@id, ntlMatingRams@id), size = n, replace = TRUE))
    # TODO: here I added the number of contribution of each category is it correct like this?
    matingPlan3 = cbind(damsOfFemalesIdForRest,
                        sample(c(siresOfFemales@id, size = n2, replace = TRUE),
                        sample(wtRams1@id, size = n3, replace = TRUE),
                        sample(ntlMatingRams@id, size = n4, replace = TRUE)))

    matingPlan = rbind(matingPlan1, matingPlan2, matingPlan3)
    lambs = makeCross2(females = c(eliteEwes, damsOfFemales),
                       males = c(eliteSires, siresOfFemales, wtRams1, ntlMatingRams),
                       crossPlan = matingPlan)
    lambs = fillInMisc(pop = lambs,
                       mothers = c(eliteEwes, damsOfFemales),
                       permEnvVar = permVar, year = yearFull)

    # ---- Data recording & estimating breeding values ----

    cat("... Data recording & estimating breeding values", as.character(Sys.time()), "\n")

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
    # database = recordData(database, pop = eliteEwes, year = yearFull, lactation = ???) # can not save like this since lactations vary
    database = recordData(database, pop = damsOfFemalesLact4, year = yearFull, lactation = 4)
    database = recordData(database, pop = damsOfFemalesLact3, year = yearFull, lactation = 3)
    database = recordData(database, pop = damsOfFemalesLact2, year = yearFull, lactation = 2)
    database = recordData(database, pop = damsOfFemalesLact1, year = yearFull, lactation = 1)
    # database = recordData(database, pop = damsOfFemales, year = yearFull, lactation = ???) # can not save like this since lactations vary
    database = recordData(database, pop = lambs, year = yearFull)

 
    # Save workspace image
    # fileName = paste("scenario_year", year, ".RData", sep = "")
    # save.image(file = fileName)
    # load(file = fileName)
    # sourceFunctions(dir = "../..") # to ensure we have the latest version (countering save.image())
  }

  # Save workspace image
  save.image(file = "scenario.RData")
  # load(file = "scenario.RData")
  # sourceFunctions(dir = "../..") # to ensure we have the latest version (countering save.image())
}
