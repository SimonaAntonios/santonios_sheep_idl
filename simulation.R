#!/usr/bin/env Rscript
#
# This script runs one scenario at a time defined via arguments. Make sure this
#   script is executable:
#
# chmod a+x simulation.R
#
# Usage:
#
# ./simulation.R rep [burnin|scenario] [std|idl] [restart]
#
# where:
#
# std - standard pedigree-based BLUP model & truncation selection on BV
# idl - std extended with IDL & truncation selection on BV + IDL

#
# examples:
#
# ./simulation.R 1 burnin
# ./simulation.R 1 scenario std
# ./simulation.R 1 scenario XX (this one may counts for the mate allocation)

# ---- Setup -------------------------------------------------------------------

cat("Setup\n")

# Clean the working environment (mostly for interactive use)
rm(list = ls())

args = commandArgs(trailingOnly = TRUE)
# args = c("1", "burnin")
# args = c("1", "scenario", "std")
# args = c("1", "scenario", "XX")

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

# Argument 3: [std|idl|idl1|idl2|MA|EFI]
if (burnin) {
  scenario = NA
} else {
  if (length(args) < 3) {
    stop("You must provide at least 3 arguments when running scenarios!")
  }
  scenario = args[3]
  if (!scenario %in% c("std", "MA", "idl","idl1","idl2","EFI" )) {
    stop("Argument 3 must be either std, MA, or idl or ,idl1, idl2 !")
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
if (FALSE) {
  requiredPkgs = c("tidyverse", "AlphaSimR", "degreenet", "data.table", "pedigreemm")
  install.packages(pkg = requiredPkgs)
}
library(tidyverse)  # for data manipulation
library(degreenet)  # for simulation of herd sizes
library(data.table) # for fast data operations (reading, writing, ...)
library(pedigreemm) # for pedigree numerator relationship inverse matrix
library(AlphaSimR)  # for stochastic simulation of a breeding programme
library(pedigreeTools) # for new getPedNrmSubset function
library(optiSel) # for Optimal Contribution Selection strategy
library(lpSolve) # for the mate allocation 

if (interactive()) {
  # Simona's folder
  setwd(dir = "/home/santonios/work/santonios/part2/simulation")
  # Gregor's folder
  setwd(dir = "~/Storages/GitBox/HighlanderLab/santonios_sheep_idl/")
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

nEwes               = 80000/2                                                        # no. of all females
NM_Fertility        = 0.9                                                          # fertility in Natural Mating (NM)
NM_prolificacy      = 1.4                                                          # prolificacy in NM
AI_Fertility        = 0.6                                                          # fertility in Artificial Insemination (AI)
AI_prolificacy      = 1.6                                                          # prolificacy in AI
survRate            = 0.75                                                         # survival rate
NMLambRate = NM_Fertility * NM_prolificacy * survRate # 0.945
AILambRate = AI_Fertility * AI_prolificacy * survRate # 0.720
pEwesInAI = 0.5
nEwesInAI = pEwesInAI * nEwes
pEwesInNM = 1 - pEwesInAI
nEwesInNM = pEwesInNM * nEwes
# Based on these numbers we will get
if (FALSE) {
  NMLambRate * nEwesInNM + AILambRate* nEwesInAI # 66,600 lambs
}

# We are setting the no. of ewes from the above number of lambs (simulating only successful lactations!) and I am deviding only the females number by 2 to reduce the number of animals in the pedigree
nFemalesInLactation = 66600/2                                                        # no. of females in lactation
nDamsOfSires          = 6000/2
nDamsOfDams      = 60600/2

nWtRams1            = 150                                                          # no. of waiting rams for progeny testing
nWtRams2            = 150                                                          # no. of waiting rams for progeny testing
nAISiresOfSires1        = 10                                                           # no. of   AI Sires Of Sires selected every year (for the first year)
nAISiresOfSires2        = 10                                                           # no. of   AI Sires Of Sires selected every year (for the second year)
nAISiresOfSires3        = 10                                                           # no. of   AI Sires Of Sires selected every year (for the third year)
nAISiresOfSires         = nAISiresOfSires1 + nAISiresOfSires2 + nAISiresOfSires3
nAISiresOfDams1    = 55                                                           # no. of   AI Sires Of Siresof dams selected every year (for the first year)
nAISiresOfDams2    = 55                                                           # no. of   AI Sires Of Siresof dams selected every year (for the second year)
nAISiresOfDams3    = 55                                                           # no. of   AI Sires Of Siresof dams selected every year (for the third year)
nNMSires  = 1000                                                         # no. of NM sires selected every year
nAISiresOfSiresDose      = 400                                                          # no. of AI doses per AI Sire Of Sires
nWtRamsAIDose       = 85                                                           # no. of AI doses per wating ram
nNMDose      = 40                                                           # no. of natural matings per NM ram

# Number of lambs from different matings, I divided them by 2
nLambsFromAISiresOfSires = round(nAISiresOfSires * nAISiresOfSiresDose * AILambRate)/2 # 8,640/2
nLambsFromAIForPT = round(nWtRams1 * nWtRamsAIDose * AILambRate)/2 # 9,180/2
nLambsFromAIRest = round((nEwesInAI*2 - (nAISiresOfSires * nAISiresOfSiresDose + nWtRams1 * nWtRamsAIDose)) * AILambRate)/2 # 10,980/2
nLambsFromNM = round(nNMSires * nNMDose * NMLambRate)/2 # 37,800/2

pDamsOfDamsLact1 = 0.30                                                         # prop. of Dams Of Dams in lactation 1
pDamsOfDamsLact2 = 0.28                                                         # prop. of Dams Of Dams in lactation 2
pDamsOfDamsLact3 = 0.24                                                         # prop. of Dams Of Dams in lactation 3
pDamsOfDamsLact4 = 0.18                                                         # prop. of Dams Of Dams in lactation 4
nDamsOfDamsLact1 = round(pDamsOfDamsLact1 * nDamsOfDams)                  # no. of Dams Of Dams in lactation 1
nDamsOfDamsLact2 = round(pDamsOfDamsLact2 * nDamsOfDams)                  # no. of Dams Of Dams in lactation 2
nDamsOfDamsLact3 = round(pDamsOfDamsLact3 * nDamsOfDams)                  # no. of Dams Of Dams in lactation 3
nDamsOfDamsLact4 = round(pDamsOfDamsLact4 * nDamsOfDams)                  # no. of Dams Of Dams in lactation 4

pDamsOfSiresLact1     = 0.35                                                         # prop. of Dams Of Sires in lactation 2
pDamsOfSiresLact2     = 0.26                                                         # prop. of Dams Of Sires in lactation 3
pDamsOfSiresLact3     = 0.23                                                         # prop. of Dams Of Sires in lactation 1 
pDamsOfSiresLact4     = 0.16                                                         # prop. of Dams Of Sires in lactation 4
nDamsOfSiresLact1     = round(pDamsOfSiresLact1 * nDamsOfSires)                          # no. of Dams Of Sires in lactation 1
nDamsOfSiresLact2     = round(pDamsOfSiresLact2 * nDamsOfSires)                          # no. of Dams Of Sires in lactation 2
nDamsOfSiresLact3     = round(pDamsOfSiresLact3 * nDamsOfSires)                          # no. of Dams Of Sires in lactation 3
nDamsOfSiresLact4     = round(pDamsOfSiresLact4 * nDamsOfSires)                          # no. of Dams Of Sires in lactation 4

startYear           = 1980                                                         # start year of the simulations

nHerds              = 237                                                          # number of herds
meanHerdSize        = 320                                                          # average herd size
sdHerdSize          = 129                                                          # sd of herd size
genInt              = 4                                                            # generation interval

# Genome parameters:
BaseNe              = 150                                                          # effective population size
nChr                = 26                                                           # no. of chromosomes
ChrSize             = 95 * 1e6                                                     # chr. size
nQTL                = 239                                                          # no. of QTLs of MY,  4500 QTL for the 272 base traits
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
  with(AncientNe, plot(Ne ~ GenerationsAgo, type = "b", main = "Generation ago with log",  log = "xy"))
  dev.off()
  
  pdf("generations_ago_Ne.pdf")
  with(AncientNe, plot(Ne ~ GenerationsAgo, type = "b", main = "Generation ago without log", xlab = "Generations ago"))
  dev.off()
}

# Trait parameters:
meanLact1           = 209                                                          # mean of milk yield in lactation 1
meanLact2           = 213                                                          # mean of milk yield in lactation 2
meanLact3           = 207                                                          # mean of milk yield in lactation 3
meanLact4           = 180                                                          # mean of milk yield in lactation 4
meanLact = c(meanLact1, meanLact2, meanLact3, meanLact4)
nTmp1 = c(nDamsOfDamsLact1, nDamsOfDamsLact2, nDamsOfDamsLact3, nDamsOfDamsLact4)
nTmp2 = c(nDamsOfSiresLact1, nDamsOfSiresLact2, nDamsOfSiresLact3, nDamsOfSiresLact4)
meanLactAll = sum(c(nTmp1 * meanLact, nTmp2 * meanLact)) / nFemalesInLactation

varScale = 0.5                                                                     # downscale variances to get mostly positive phenotypes
yearVar             = 1000 * varScale                                              # year variance
herdVar             = 1200 * varScale                                              # herd variance
herdYearVar         = 1000 * varScale                                              # herd-year variance
addVar              = 1200 * varScale                                              # additive genetic variance
# addVar              = 1000 * varScale                                              # additive genetic variance
# permVar           = 500 * varScale                                               # permanent environmental variance
permVar             = 400 * varScale # 500 (Old value) - 100 = 400 (new value)     # permanent environmental variance (we started with old, 100 was targeted domVar)
resVar              = 1500 * varScale                                              # residual variance
domVar              = addVar * 0.1                                                 # dominance variance
fullInbreedDepress  = 70                                                           # depression with complete inbreeding
meanDD              = 0.48                                                         # dominance parameters - see use of altAddTraitAD() below
varDD               = 0.06                                                         # dominance parameters - see use of altAddTraitAD() below
idlVar              = 12923.0 * varScale
corrIDLAdd          = -0.1
covIDLAdd           = corrIDLAdd*(addVar*idlVar)^0.5
# TODO: how do we get this? Estimate from the data?
#       https://github.com/SimonaAntonios/santonios_sheep_idl/issues/17
# Based on these values we expect this phenotypic variance and ratios
if (FALSE) {
  # "Full" peno variance (still missing lactation means!)
  (phenVar = yearVar + herdVar + herdYearVar + permVar + addVar + domVar + resVar) # 3250
  permVar / phenVar # ~0.062
  addVar / phenVar # ~0.154
  domVar / phenVar # ~0.015
  
  # "Reduced" peno variance (still missing lactation means!)
  (phenVar = permVar + addVar + domVar + resVar) # 1500
  permVar / phenVar # ~0.133
  addVar / phenVar # ~0.333
  domVar / phenVar # ~0.033
}

nBurninYears        = 10                                                           # no. of years for burnin
nScenarioYears      = 15                                                           # no. of years for scenarios

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
  
  # if (FALSE) {
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
  save(founderPop, file = "founderPop.RData")
  # }
  # load("founderPop.RData")
  # load("../../founder_runMacs2.RData") # note that we are in rep_x/[burnin|scenario] folder at this stage
  # sourceFunctions(dir = "../..") # to ensure we have the latest version (countering save.image(), if it was used)
  
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
    source(file = "../../altAddTraitAD.R")
    # TODO: check if dominance params change with scaled down variances
    #       https://github.com/SimonaAntonios/santonios_sheep_idl/issues/18
    altAddTraitAD(nQtlPerChr = nQTLPerChr,
                  mean = meanLact1,
                  varA = addVar,
                  varD = domVar,
                  inbrDepr = fullInbreedDepress,
                  limMeanDD = c(-1, 2),
                  limVarDD = c(0, 2))
    # A new trait called Trait1 was added with  varScale = 0.5 and addVar = 1200 * varScale . 
    # varD = 59.31778 
    # inbrDepr = 70.01133 
    # meanDD = 0.4783672 
    # varDD = 0.05730109 
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
  herdEffect = sampleEffect(n = nHerds, var = herdVar)
  herdYearEffect = sampleEffect(n = nHerds, var = herdYearVar)
  herds = list(
    herd = 1:nHerds,
    herdSize = herdSize,
    herdEffect = herdEffect,
    herdYearEffect = herdYearEffect
  )
  yearEffect = sampleEffect(n = 1, var = yearVar)
  
  # ---- Fill-in ---------------------------------------------------------------
  
  cat("Fill-in", as.character(Sys.time()), "\n")
  
  # Generate initial founders
  basePop = newPop(founderPop)
  
  # nCrosses = 10 * nSomething below is to ensure we don't oversample founders
  
  #   AI Sires Of Sires(AI)
  # ... in the 3rd year in service
  AISiresOfSires3 = randCross(pop = basePop, nCrosses = 10 * nAISiresOfSires3)
  AISiresOfSires3 = selectWithinFam(pop = AISiresOfSires3, nInd = 1, use = "rand", famType = "M")
  AISiresOfSires3 = selectInd(pop = AISiresOfSires3, nInd = nAISiresOfSires3, use = "rand")
  AISiresOfSires3@sex[] = "M"
  AISiresOfSires3@father[] = "0"
  AISiresOfSires3@mother[] = "0"
  AISiresOfSires3 = fillInMisc(pop = AISiresOfSires3, year = startYear - 5.5)
  # ... in the 2nd year in service
  AISiresOfSires2 = randCross(pop = basePop, nCrosses = 10 * nAISiresOfSires2)
  AISiresOfSires2 = selectWithinFam(pop = AISiresOfSires2, nInd = 1, use = "rand", famType = "M")
  AISiresOfSires2 = selectInd(pop = AISiresOfSires2, nInd = nAISiresOfSires2, use = "rand")
  AISiresOfSires2@sex[] = "M"
  AISiresOfSires2@father[] = "0"
  AISiresOfSires2@mother[] = "0"
  AISiresOfSires2 = fillInMisc(pop = AISiresOfSires2, year = startYear - 4.5)
  # ... in the 1st year in service
  AISiresOfSires1 = randCross(pop = basePop, nCrosses = 10 * nAISiresOfSires1)
  AISiresOfSires1 = selectWithinFam(pop = AISiresOfSires1, nInd = 1, use = "rand", famType = "M")
  AISiresOfSires1 = selectInd(pop = AISiresOfSires1, nInd = nAISiresOfSires1, use = "rand")
  AISiresOfSires1@sex[] = "M"
  AISiresOfSires1@father[] = "0"
  AISiresOfSires1@mother[] = "0"
  AISiresOfSires1 = fillInMisc(pop = AISiresOfSires1, year = startYear - 3.5)
  
  AISiresOfSires = c(AISiresOfSires3, AISiresOfSires2, AISiresOfSires1)
  
  # AI Sire Of Dams (AI)
  # ... in the 3rd year in service
  AISiresOfDams3 = randCross(pop = basePop, nCrosses = 10 * nAISiresOfDams3)
  AISiresOfDams3 = selectWithinFam(pop = AISiresOfDams3, nInd = 1, use = "rand", famType = "M")
  AISiresOfDams3 = selectInd(pop = AISiresOfDams3, nInd = nAISiresOfDams3, use = "rand")
  AISiresOfDams3@sex[] = "M"
  AISiresOfDams3@father[] = "0"
  AISiresOfDams3@mother[] = "0"
  AISiresOfDams3 = fillInMisc(pop = AISiresOfDams3, year = startYear - 5.5)
  # ... in the 2nd year in service
  AISiresOfDams2 = randCross(pop = basePop, nCrosses = 10 * nAISiresOfDams2)
  AISiresOfDams2 = selectWithinFam(pop = AISiresOfDams2, nInd = 1, use = "rand", famType = "M")
  AISiresOfDams2 = selectInd(pop = AISiresOfDams2, nInd = nAISiresOfDams2, use = "rand")
  AISiresOfDams2@sex[] = "M"
  AISiresOfDams2@father[] = "0"
  AISiresOfDams2@mother[] = "0"
  AISiresOfDams2 = fillInMisc(pop = AISiresOfDams2, year = startYear - 4.5)
  # ... in the 1st year in service
  AISiresOfDams1 = randCross(pop = basePop, nCrosses = 10 * nAISiresOfDams1)
  AISiresOfDams1 = selectWithinFam(pop = AISiresOfDams1, nInd = 1, use = "rand", famType = "M")
  AISiresOfDams1 = selectInd(pop = AISiresOfDams1, nInd = nAISiresOfDams1, use = "rand")
  AISiresOfDams1@sex[] = "M"
  AISiresOfDams1@father[] = "0"
  AISiresOfDams1@mother[] = "0"
  AISiresOfDams1 = fillInMisc(pop = AISiresOfDams1, year = startYear - 3.5)
  
  AISiresOfDams = c(AISiresOfDams3, AISiresOfDams2, AISiresOfDams1)
  
  # Waiting rams
  # ... 1.5 years old
  wtRams2  = randCross(basePop, nCrosses = 10 * nWtRams2)
  wtRams2 = selectWithinFam(pop = wtRams2, nInd = 1, use = "rand", famType = "M")
  wtRams2 = selectInd(pop = wtRams2, nInd = nWtRams2, use = "rand")
  wtRams2@sex[] = "M"
  wtRams2@father[] = "0"
  wtRams2@mother[] = "0"
  wtRams2 = fillInMisc(pop = wtRams2, year = startYear - 1.5)
  # ... 0.5 years old
  wtRams1  = randCross(basePop, nCrosses = 10 * nWtRams1)
  wtRams1 = selectWithinFam(pop = wtRams1, nInd = 1, use = "rand", famType = "M")
  wtRams1 = selectInd(pop = wtRams1, nInd = nWtRams1, use = "rand")
  wtRams1@sex[] = "M"
  wtRams1@father[] = "0"
  wtRams1@mother[] = "0"
  wtRams1 = fillInMisc(pop = wtRams1, year = startYear - 0.5)
  
  # Natural mating rams (NM)
  NMSires  = randCross(basePop, nCrosses = nNMSires)
  NMSires@sex[] = "M"
  NMSires@father[] = "0"
  NMSires@mother[] = "0"
  NMSires = fillInMisc(pop = NMSires, year = startYear - 1.5)
  
  # Dams Of Sires
  # ... in the 4th lactation
  damsOfSiresLact4 = randCross(basePop, nCrosses = nDamsOfSiresLact4)
  damsOfSiresLact4@sex[] = "F"
  damsOfSiresLact4@father[] = "0"
  damsOfSiresLact4@mother[] = "0"
  damsOfSiresLact4 = fillInMisc(pop = damsOfSiresLact4, year = startYear - 4,
                                herds = herds, permEnvVar = permVar)
  # ... in the 3rd lactation
  damsOfSiresLact3 = randCross(basePop, nCrosses = nDamsOfSiresLact3)
  damsOfSiresLact3@sex[] = "F"
  damsOfSiresLact3@father[] = "0"
  damsOfSiresLact3@mother[] = "0"
  damsOfSiresLact3 = fillInMisc(pop = damsOfSiresLact3, year = startYear - 3,
                                herds = herds, permEnvVar = permVar)
  # ... in the 2nd lactation
  damsOfSiresLact2 = randCross(basePop, nCrosses = nDamsOfSiresLact2)
  damsOfSiresLact2@sex[] = "F"
  damsOfSiresLact2@father[] = "0"
  damsOfSiresLact2@mother[] = "0"
  damsOfSiresLact2 = fillInMisc(pop = damsOfSiresLact2, year = startYear - 2,
                                herds = herds, permEnvVar = permVar)
  # ... in the 1st lactation
  damsOfSiresLact1 = randCross(basePop, nCrosses =  nDamsOfSiresLact1)
  damsOfSiresLact1@sex[] = "F"
  damsOfSiresLact1@father[] = "0"
  damsOfSiresLact1@mother[] = "0"
  damsOfSiresLact1 = fillInMisc(pop = damsOfSiresLact1, year = startYear - 1,
                                herds = herds, permEnvVar = permVar)
  
  damsOfSires = c(damsOfSiresLact4, damsOfSiresLact3, damsOfSiresLact2, damsOfSiresLact1)
  
  # Dams Of Dams
  # ... in the 4th lactation
  damsOfDamsLact4 = randCross(basePop, nCrosses = nDamsOfDamsLact4)
  damsOfDamsLact4@sex[] = "F"
  damsOfDamsLact4@father[] = "0"
  damsOfDamsLact4@mother[] = "0"
  damsOfDamsLact4 = fillInMisc(pop = damsOfDamsLact4, year = startYear - 4,
                               herds = herds, permEnvVar = permVar)
  # ... in the 3rd lactation
  damsOfDamsLact3 = randCross(basePop, nCrosses = nDamsOfDamsLact3)
  damsOfDamsLact3@sex[] = "F"
  damsOfDamsLact3@father[] = "0"
  damsOfDamsLact3@mother[] = "0"
  damsOfDamsLact3 = fillInMisc(pop = damsOfDamsLact3, year = startYear - 3,
                               herds = herds, permEnvVar = permVar)
  # ... in the 2nd lactation
  damsOfDamsLact2 = randCross(basePop, nCrosses =  nDamsOfDamsLact2)
  damsOfDamsLact2@sex[] = "F"
  damsOfDamsLact2@father[] = "0"
  damsOfDamsLact2@mother[] = "0"
  damsOfDamsLact2 = fillInMisc(pop = damsOfDamsLact2, year = startYear - 2,
                               herds = herds, permEnvVar = permVar)
  # ... in the 1st lactation
  damsOfDamsLact1 = randCross(basePop, nCrosses =  nDamsOfDamsLact1)
  damsOfDamsLact1@sex[] = "F"
  damsOfDamsLact1@father[] = "0"
  damsOfDamsLact1@mother[] = "0"
  damsOfDamsLact1 = fillInMisc(pop = damsOfDamsLact1, year = startYear - 1,
                               herds = herds, permEnvVar = permVar)
  
  damsOfDams = c(damsOfDamsLact4, damsOfDamsLact3, damsOfDamsLact2, damsOfDamsLact1)
  females = c(damsOfSires, damsOfDams)
  # Lambs
  matingPlan1 = cbind(damsOfSires@id,
                      sample(AISiresOfSires@id, size = nDamsOfSires, replace = TRUE))
  
  n = nLambsFromAISiresOfSires - nDamsOfSires
  damsOfDamsId = damsOfDams@id
  damsOfDamsIdForElite = sample(damsOfDamsId, size = n)
  matingPlan2 = cbind(damsOfDamsIdForElite,
                      sample(AISiresOfSires@id, size = n, replace = TRUE))
  
  n = nDamsOfDams - (nLambsFromAISiresOfSires - nDamsOfSires)
  damsOfDamsIdForRest = damsOfDamsId[!damsOfDamsId %in% damsOfDamsIdForElite]
  matingPlan3 = cbind(damsOfDamsIdForRest,
                      sample(c(sample(AISiresOfDams@id, size = nLambsFromAIRest, replace = TRUE),
                               sample(wtRams1@id, size = nLambsFromAIForPT, replace = TRUE),
                               sample(NMSires@id, size = nLambsFromNM, replace = TRUE))))
  # Note that:
  # 1) sample(pop, size = n, replace = TRUE) gives us n contributions from pop
  # 2) sample(c(sample(), sample(), sample())) shuffles selected male contributions
  #    across female contributions (=random mating)
  
  matingPlan = rbind(matingPlan1, matingPlan2, matingPlan3)
  lambs = makeCross2(females = c(damsOfSires, damsOfDams),
                     males = c(AISiresOfSires, AISiresOfDams, wtRams1, NMSires),
                     crossPlan = matingPlan)
  lambs = fillInMisc(pop = lambs,
                     mothers = c(damsOfSires, damsOfDams),
                     permEnvVar = permVar, year = startYear)
  
  database = recordData(pop = AISiresOfSires3, year = startYear)
  database = recordData(database, pop = AISiresOfSires2, year = startYear)
  database = recordData(database, pop = AISiresOfSires1, year = startYear)
  database = recordData(database, pop = AISiresOfSires, year = startYear)
  database = recordData(database, pop = AISiresOfDams3, year = startYear)
  database = recordData(database, pop = AISiresOfDams2, year = startYear)
  database = recordData(database, pop = AISiresOfDams1, year = startYear)
  database = recordData(database, pop = AISiresOfDams, year = startYear)
  database = recordData(database, pop = wtRams2, year = startYear)
  database = recordData(database, pop = wtRams1, year = startYear)
  database = recordData(database, pop = NMSires, year = startYear)
  database = recordData(database, pop = damsOfSiresLact4, year = startYear, lactation = 4)
  database = recordData(database, pop = damsOfSiresLact3, year = startYear, lactation = 3)
  database = recordData(database, pop = damsOfSiresLact2, year = startYear, lactation = 2)
  database = recordData(database, pop = damsOfSiresLact1, year = startYear, lactation = 1)
  # database = recordData(database, pop = damsOfSires, year = yearFull, lactation = ???) # can not save like this since lactations vary
  database = recordData(database, pop = damsOfDamsLact4, year = startYear, lactation = 4)
  database = recordData(database, pop = damsOfDamsLact3, year = startYear, lactation = 3)
  database = recordData(database, pop = damsOfDamsLact2, year = startYear, lactation = 2)
  database = recordData(database, pop = damsOfDamsLact1, year = startYear, lactation = 1)
  # database = recordData(database, pop = damsOfDams, year = yearFull, lactation = ???) # can not save like this since lactations vary
  database = recordData(database, pop = lambs, year = startYear)
  
  # Save workspace image
  save.image(file = "fillin.RData")
  # load(file = "fillin.RData")
  # sourceFunctions(dir = "../..") # to ensure we have the latest version (countering save.image())
  
  # ----- Loop breeding programme over years -----------------------------------
  
  cat("Loop breeding programme over years", as.character(Sys.time()), "\n")
  
  for (year in yearToDo:nBurninYears) {
    # year = yearToDo
    yearFull = startYear + year
    
    cat("Year", year, " (", yearFull, ") ", as.character(Sys.time()), "\n")
    
    yearEffect = sampleEffect(n = 1, var = yearVar)
    herds$herdYearEffect = sampleEffect(n = nHerds, var = herdYearVar)
    
    if (year <= 1) {
      use = "rand"
    } else {
      use= "ebv"
    }
    
    # ---- Phenotyping ----
    
    cat("Phenotyping", as.character(Sys.time()), "\n")
    
    damsOfSiresLact4 = setPhenoEwe(damsOfSiresLact4, varE = resVar,
                                   mean = meanLact4, yearEffect = yearEffect, herds = herds)
    
    damsOfSiresLact3 = setPhenoEwe(damsOfSiresLact3, varE = resVar,
                                   mean = meanLact3, yearEffect = yearEffect, herds = herds)
    
    damsOfSiresLact2 = setPhenoEwe(damsOfSiresLact2, varE = resVar,
                                   mean = meanLact2, yearEffect = yearEffect, herds = herds)
    
    damsOfSiresLact1 = setPhenoEwe(damsOfSiresLact1, varE = resVar,
                                   mean = meanLact1, yearEffect = yearEffect, herds = herds)
    
    damsOfDamsLact4 = setPhenoEwe(damsOfDamsLact4, varE = resVar,
                                  mean = meanLact4, yearEffect = yearEffect, herds = herds)
    
    damsOfDamsLact3 = setPhenoEwe(damsOfDamsLact3, varE = resVar,
                                  mean = meanLact3, yearEffect = yearEffect, herds = herds)
    
    damsOfDamsLact2 = setPhenoEwe(damsOfDamsLact2, varE = resVar,
                                  mean = meanLact2, yearEffect = yearEffect, herds = herds)
    
    damsOfDamsLact1 = setPhenoEwe(damsOfDamsLact1, varE = resVar,
                                  mean = meanLact1, yearEffect = yearEffect, herds = herds)
    
    database = setDatabasePheno(database, pop = damsOfSiresLact4)
    database = setDatabasePheno(database, pop = damsOfSiresLact3)
    database = setDatabasePheno(database, pop = damsOfSiresLact2)
    database = setDatabasePheno(database, pop = damsOfSiresLact1)
    database = setDatabasePheno(database, pop = damsOfDamsLact4)
    database = setDatabasePheno(database, pop = damsOfDamsLact3)
    database = setDatabasePheno(database, pop = damsOfDamsLact2)
    database = setDatabasePheno(database, pop = damsOfDamsLact1)
    
    # ---- Genetic evaluation ----
    
    cat("Genetic evaluation", as.character(Sys.time()), "\n")
    
    # We will remove male lambs from the evaluation - those that will never
    #   contribute to next generations (we will still calculate their parent
    #   average, but later)
    # ... first we find all male lambs
    sel = database$General$Pop == "lambs" & database$General$Sex == "M"
    maleLambs = database$General[sel, "IId"]
    # ... then we find those that were selected for reproduction
    sel = database$General$Pop %in% c("wtRams1", "NMSires")
    reproMales = database$General[sel, "IId"]
    # ... the other were culled
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
    AISiresOfSires3 = setEbv(AISiresOfSires3, ebv = pedEbv)
    AISiresOfSires2 = setEbv(AISiresOfSires2, ebv = pedEbv)
    AISiresOfSires1 = setEbv(AISiresOfSires1, ebv = pedEbv)
    AISiresOfSires  = setEbv(AISiresOfSires, ebv = pedEbv)
    AISiresOfDams3 = setEbv(AISiresOfDams3, ebv = pedEbv)
    AISiresOfDams2 = setEbv(AISiresOfDams2, ebv = pedEbv)
    AISiresOfDams1 = setEbv(AISiresOfDams1, ebv = pedEbv)
    AISiresOfDams = setEbv(AISiresOfDams, ebv = pedEbv)
    wtRams2 = setEbv(wtRams2, ebv = pedEbv)
    wtRams1 = setEbv(wtRams1, ebv = pedEbv)
    NMSires = setEbv(NMSires, ebv = pedEbv)
    damsOfSiresLact4 = setEbv(damsOfSiresLact4, ebv = pedEbv)
    damsOfSiresLact3 = setEbv(damsOfSiresLact3, ebv = pedEbv)
    damsOfSiresLact2 = setEbv(damsOfSiresLact2, ebv = pedEbv)
    damsOfSiresLact1 = setEbv(damsOfSiresLact1, ebv = pedEbv)
    damsOfSires = setEbv(damsOfSires, ebv = pedEbv)
    damsOfDamsLact4 = setEbv(damsOfDamsLact4, ebv = pedEbv)
    damsOfDamsLact3 = setEbv(damsOfDamsLact3, ebv = pedEbv)
    damsOfDamsLact2 = setEbv(damsOfDamsLact2, ebv = pedEbv)
    damsOfDamsLact1 = setEbv(damsOfDamsLact1, ebv = pedEbv)
    damsOfDams = setEbv(damsOfDams, ebv = pedEbv)
    # ... we could use setEbv(lambs), but we have excluded culled lambs from the
    #     evaluation so we will just calculate parent average for all lambs - this
    #     works for pedigree BLUP, but not for genomic BLUP (if some lambs would
    #     have been genotyped)
    # lambs             = setEbv(lambs, ebv = pedEbv)
    selM = match(x = lambs@mother, table = pedEbv$IId)
    selF = match(x = lambs@father, table = pedEbv$IId)
    lambs@ebv = as.matrix((pedEbv[selM, -1] + pedEbv[selF, -1]) / 2)
    
    # Save current EBVs into the database
    database = setDatabaseEbv(database, pop = AISiresOfSires3)
    database = setDatabaseEbv(database, pop = AISiresOfSires2)
    database = setDatabaseEbv(database, pop = AISiresOfSires1)
    database = setDatabaseEbv(database, pop = AISiresOfSires)
    database = setDatabaseEbv(database, pop = AISiresOfDams3)
    database = setDatabaseEbv(database, pop = AISiresOfDams2)
    database = setDatabaseEbv(database, pop = AISiresOfDams1)
    database = setDatabaseEbv(database, pop = AISiresOfDams)
    database = setDatabaseEbv(database, pop = wtRams2)
    database = setDatabaseEbv(database, pop = wtRams1)
    database = setDatabaseEbv(database, pop = NMSires)
    database = setDatabaseEbv(database, pop = damsOfSiresLact4)
    database = setDatabaseEbv(database, pop = damsOfSiresLact3)
    database = setDatabaseEbv(database, pop = damsOfSiresLact2)
    database = setDatabaseEbv(database, pop = damsOfSiresLact1)
    # database = setDatabaseEbv(database, pop = damsOfSires) # we didn't save this pop in the database
    database = setDatabaseEbv(database, pop = damsOfDamsLact4)
    database = setDatabaseEbv(database, pop = damsOfDamsLact3)
    database = setDatabaseEbv(database, pop = damsOfDamsLact2)
    database = setDatabaseEbv(database, pop = damsOfDamsLact1)
    # database = setDatabaseEbv(database, pop = damsOfDams) # we didn't save this pop in the database
    database = setDatabaseEbv(database, pop = lambs)
    
    damsOfSires = c(damsOfSiresLact4, damsOfSiresLact3, damsOfSiresLact2, damsOfSiresLact1)
    damsOfDams = c(damsOfDamsLact4, damsOfDamsLact3, damsOfDamsLact2, damsOfDamsLact1)
    females = c(damsOfSires, damsOfDams)
    wtRams = c(wtRams1, wtRams2)
    males = c(AISiresOfSires, AISiresOfDams, wtRams, NMSires)
    correlationEbvTgv = data.frame(year = year,
                                   AISiresOfSires3 = calcAccuracyEbvVsTgv(AISiresOfSires3),
                                   AISiresOfSires2 = calcAccuracyEbvVsTgv(AISiresOfSires2),
                                   AISiresOfSires1 = calcAccuracyEbvVsTgv(AISiresOfSires1),
                                   AISiresOfSires = calcAccuracyEbvVsTgv(AISiresOfSires),
                                   AISiresOfDams3 = calcAccuracyEbvVsTgv(AISiresOfDams3),
                                   AISiresOfDams2 = calcAccuracyEbvVsTgv(AISiresOfDams2),
                                   AISiresOfDams1 = calcAccuracyEbvVsTgv(AISiresOfDams1),
                                   AISiresOfDams = calcAccuracyEbvVsTgv(AISiresOfDams),
                                   males = calcAccuracyEbvVsTgv(males),
                                   wtRams2 = calcAccuracyEbvVsTgv(wtRams2),
                                   wtRams1 = calcAccuracyEbvVsTgv(wtRams1),
                                   wtRams = calcAccuracyEbvVsTgv(wtRams),
                                   NMSires = calcAccuracyEbvVsTgv(NMSires),
                                   damsOfSiresLact4 = calcAccuracyEbvVsTgv(damsOfSiresLact4),
                                   damsOfSiresLact3 = calcAccuracyEbvVsTgv(damsOfSiresLact3),
                                   damsOfSiresLact2 = calcAccuracyEbvVsTgv(damsOfSiresLact2),
                                   damsOfSiresLact1 = calcAccuracyEbvVsTgv(damsOfSiresLact1),
                                   damsOfSires = calcAccuracyEbvVsTgv(damsOfSires),
                                   damsOfDamsLact4 = calcAccuracyEbvVsTgv(damsOfDamsLact4),
                                   damsOfDamsLact3 = calcAccuracyEbvVsTgv(damsOfDamsLact3),
                                   damsOfDamsLact2 = calcAccuracyEbvVsTgv(damsOfDamsLact2),
                                   damsOfDamsLact1 = calcAccuracyEbvVsTgv(damsOfDamsLact1),
                                   damsOfDams = calcAccuracyEbvVsTgv(damsOfDams),
                                   females = calcAccuracyEbvVsTgv(females),
                                   lambs = calcAccuracyEbvVsTgv(lambs))
    # add = ifelse(year == 1, FALSE, TRUE)
    # write.table(x = correlationEbvTgv, file = "calcAccuracyEbvVsTgv.txt", append = add, col.names = !add)
    if (!file.exists('calcAccuracyEbvVsTgv.txt')) {
      # If the file doesn't exist, write the column names along with the data
      write.table(correlationEbvTgv, file = "calcAccuracyEbvVsTgv.txt", col.names = TRUE, row.names = FALSE)
    } else {
      # If the file exists, append the data without writing column names
      write.table(correlationEbvTgv, file = "calcAccuracyEbvVsTgv.txt", col.names = FALSE, row.names = FALSE, append = TRUE)
    }
    
    
    correlationEbvTbv = data.frame(year = year,
                                   AISiresOfSires3 = calcAccuracyEbvVsTbv(AISiresOfSires3),
                                   AISiresOfSires2 = calcAccuracyEbvVsTbv(AISiresOfSires2),
                                   AISiresOfSires1 = calcAccuracyEbvVsTbv(AISiresOfSires1),
                                   AISiresOfSires = calcAccuracyEbvVsTbv(AISiresOfSires),
                                   AISiresOfDams3 = calcAccuracyEbvVsTbv(AISiresOfDams3),
                                   AISiresOfDams2 = calcAccuracyEbvVsTbv(AISiresOfDams2),
                                   AISiresOfDams1 = calcAccuracyEbvVsTbv(AISiresOfDams1),
                                   AISiresOfDams = calcAccuracyEbvVsTbv(AISiresOfDams),
                                   males = calcAccuracyEbvVsTbv(males),
                                   wtRams2 = calcAccuracyEbvVsTbv(wtRams2),
                                   wtRams1 = calcAccuracyEbvVsTbv(wtRams1),
                                   wtRams = calcAccuracyEbvVsTbv(wtRams),
                                   NMSires = calcAccuracyEbvVsTbv(NMSires),
                                   damsOfSiresLact4 = calcAccuracyEbvVsTbv(damsOfSiresLact4),
                                   damsOfSiresLact3 = calcAccuracyEbvVsTbv(damsOfSiresLact3),
                                   damsOfSiresLact2 = calcAccuracyEbvVsTbv(damsOfSiresLact2),
                                   damsOfSiresLact1 = calcAccuracyEbvVsTbv(damsOfSiresLact1),
                                   damsOfSires = calcAccuracyEbvVsTbv(damsOfSires),
                                   damsOfDamsLact4 = calcAccuracyEbvVsTbv(damsOfDamsLact4),
                                   damsOfDamsLact3 = calcAccuracyEbvVsTbv(damsOfDamsLact3),
                                   damsOfDamsLact2 = calcAccuracyEbvVsTbv(damsOfDamsLact2),
                                   damsOfDamsLact1 = calcAccuracyEbvVsTbv(damsOfDamsLact1),
                                   damsOfDams = calcAccuracyEbvVsTbv(damsOfDams),
                                   females = calcAccuracyEbvVsTbv(females),
                                   lambs = calcAccuracyEbvVsTbv(lambs))
    
    # add = ifelse(year == 1, FALSE, TRUE)
    # write.table(x = correlationEbvTbv, file = "calcAccuracyEbvVsTbv.txt", append = add, col.names = !add)
    if (!file.exists('calcAccuracyEbvVsTbv.txt')) {
      # If the file doesn't exist, write the column names along with the data
      write.table(correlationEbvTbv, file = "calcAccuracyEbvVsTbv.txt", col.names = TRUE, row.names = FALSE)
    } else {
      # If the file exists, append the data without writing column names
      write.table(correlationEbvTbv, file = "calcAccuracyEbvVsTbv.txt", col.names = FALSE, row.names = FALSE, append = TRUE)
    }
    
    inbreeding      = read.table('renf90.inb', header=F)
    inbreeding$year = year
    write.table(x = inbreeding, file = "inbreeding.txt", append = TRUE, col.names = FALSE)
    
    
    animal_data = data.frame(year        = year,
                             id_orig          = database$General$IId,
                             pop         = database$General$Pop,
                             YearOfBirth = database$General$YearOfBirth,
                             EBV_IDL         = database$Ebv,
                             GV          = database$Gv,
                             TBV         = database$General$Tbv.Trait1
    )
    
    animal_data = animal_data %>%
      filter(pop != "AISiresOfSires", pop != "AISiresOfDams" )
    
    # YearOfBirth = as.numeric(unlist(database$General$YearOfBirth))
    # animal_data$YearOfBirth =YearOfBirth
    # print(head(animal_data))
    #animal_data <- animal_data %>%
    #	select(id_orig, pop, YearOfBirth, EBV_IDL, GV, TBV)
    #animal_data <- animal_data[, c("year", "id_orig", "pop", "YearOfBirth", "EBV_IDL", "GV", "TBV")]
    
    write.table(x = animal_data, file = "animal_data.txt", append = TRUE, col.names = FALSE)
    
    animalData = data.frame(year        = year,
                            id_orig     = c(AISiresOfSires3@id,AISiresOfSires2@id,AISiresOfSires1@id,AISiresOfDams3@id,AISiresOfDams2@id,AISiresOfDams1@id,wtRams2@id,wtRams1@id,NMSires@id,damsOfSiresLact4@id, damsOfSiresLact3@id,damsOfSiresLact2@id,damsOfSiresLact1@id,damsOfDamsLact4@id, damsOfDamsLact3@id,damsOfDamsLact2@id,damsOfDamsLact1@id,lambs@id),
                            pop         = rep(c("AISiresOfSires3","AISiresOfSires2","AISiresOfSires1","AISiresOfDams3","AISiresOfDams2","AISiresOfDams1","wtRams2","wtRams1","NMSires","damsOfSiresLact4", "damsOfSiresLact3","damsOfSiresLact2","damsOfSiresLact1","damsOfDamsLact4", "damsOfDamsLact3","damsOfDamsLact2","damsOfDamsLact1","lambs"),
                                              c(AISiresOfSires3@nInd,AISiresOfSires2@nInd,AISiresOfSires1@nInd,AISiresOfDams3@nInd,AISiresOfDams2@nInd,AISiresOfDams1@nInd,wtRams2@nInd,wtRams1@nInd,NMSires@nInd,damsOfSiresLact4@nInd, damsOfSiresLact3@nInd,damsOfSiresLact2@nInd,damsOfSiresLact1@nInd,damsOfDamsLact4@nInd, damsOfDamsLact3@nInd,damsOfDamsLact2@nInd,damsOfDamsLact1@nInd,lambs@nInd)),
                            EBV_IDL     = c(AISiresOfSires3@ebv,AISiresOfSires2@ebv,AISiresOfSires1@ebv,AISiresOfDams3@ebv,AISiresOfDams2@ebv,AISiresOfDams1@ebv,wtRams2@ebv,wtRams1@ebv,NMSires@ebv,damsOfSiresLact4@ebv, damsOfSiresLact3@ebv,damsOfSiresLact2@ebv,damsOfSiresLact1@ebv,damsOfDamsLact4@ebv, damsOfDamsLact3@ebv,damsOfDamsLact2@ebv,damsOfDamsLact1@ebv,lambs@ebv),
                            GV          = c(AISiresOfSires3@gv,AISiresOfSires2@gv,AISiresOfSires1@gv,AISiresOfDams3@gv,AISiresOfDams2@gv,AISiresOfDams1@gv,wtRams2@gv,wtRams1@gv,NMSires@gv,damsOfSiresLact4@gv, damsOfSiresLact3@gv,damsOfSiresLact2@gv,damsOfSiresLact1@gv,damsOfDamsLact4@gv, damsOfDamsLact3@gv,damsOfDamsLact2@gv,damsOfDamsLact1@gv,lambs@gv),
                            TBV         = c(bv(AISiresOfSires3),bv(AISiresOfSires2),bv(AISiresOfSires1),bv(AISiresOfDams3),bv(AISiresOfDams2),bv(AISiresOfDams1),bv(wtRams2),bv(wtRams1),bv(NMSires),bv(damsOfSiresLact4), bv(damsOfSiresLact3),bv(damsOfSiresLact2),bv(damsOfSiresLact1),bv(damsOfDamsLact4), bv(damsOfDamsLact3),bv(damsOfDamsLact2),bv(damsOfDamsLact1),bv(lambs)))
    
    write.table(x = animalData, file = "animalData.txt", append = TRUE, col.names = FALSE)
    
    solution = read.table('solutions.orig', skip=1)
    solution$year = year
    if (!file.exists('solutions.orig.year')) {
      # If the file doesn't exist, write the column names along with the data
      write.table(solution, file = "solutions.orig.year", col.names = TRUE, row.names = FALSE)
    } else {
      # If the file exists, append the data without writing column names
      write.table(solution, file = "solutions.orig.year", col.names = FALSE, row.names = FALSE, append = TRUE)
    }
    
    # ---- Select rams ----
    # ---- ...   AI Sires Of Sires----
    
    cat("... AI Sire Of Siress", as.character(Sys.time()), "\n")
    
    AISiresOfSires3 = AISiresOfSires2 # AISiresOfSires3 are 4.5 years old here
    AISiresOfSires2 = AISiresOfSires1 # AISiresOfSires2 are 3.5 years old here
    
    if (all(wtRams2@father == "0")) {
      AISiresOfSires1 = selectInd(pop = wtRams2, nInd = nAISiresOfSires1, # AISiresOfSires1 are 3.5 years old here
                                  use = use)
    } else {
      AISiresOfSires1 = selectWithinFam(pop = wtRams2, nInd = 1, # AISiresOfSires1 are 3.5 years old here
                                        use = use, famType = "M")
      AISiresOfSires1 = selectInd(pop = AISiresOfSires1, nInd = nAISiresOfSires1,
                                  use = use)
    }
    
    # ---- ... AI Sire Of Dams ----
    
    cat("... AI Sire Of Dams", as.character(Sys.time()), "\n")
    
    AISiresOfDams3 = AISiresOfDams2 # AISiresOfDams3 are 5.5 years old here
    AISiresOfDams2 = AISiresOfDams1 # AISiresOfDams2 are 4.5 years old here
    
    if (all(wtRams2@father == "0")) {
      sel = !wtRams2@id %in% AISiresOfSires1@id
      AISiresOfDams1 = selectInd(pop = wtRams2[sel], nInd = nAISiresOfDams1, # AISiresOfDams1 are 3.5 years old here
                                 use = use)
    } else {
      sel = !wtRams2@id %in% AISiresOfSires1@id
      AISiresOfDams1 = selectWithinFam(pop = wtRams2[sel], nInd = 2, # AISiresOfDams1 are 3.5 years old here
                                       use = use, famType = "M")
      AISiresOfDams1 = selectInd(pop = AISiresOfDams1, nInd = nAISiresOfDams1,
                                 use = use)
    }
    
    cat("Select rams\n")
    
    # ---- ... Selecting males from lambs ----
    # Note that we select wtRams only from lambs from AI Sire Of Sires (not waiting rams) and Dams Of Sires
    selWtRams = lambs@father %in% AISiresOfSires@id & lambs@mother %in% damsOfSires@id
    n = ceiling(nWtRams1 / length(AISiresOfSires@id))
    wtRams2 = wtRams1
    wtRams1 = selectWithinFam(pop = lambs[selWtRams], nInd = n,
                              use = use, famType = "M", sex = "M")
    wtRams1 = selectInd(pop = wtRams1, nInd = nWtRams1,
                        use = use)
    
    selNtlLambs = lambs@father %in% c(AISiresOfSires@id, AISiresOfDams@id) &
      !lambs@id %in% wtRams1@id
    NMSires =  selectInd(pop = lambs[selNtlLambs], nInd = nNMSires,
                         use = use, famType = "M", sex = "M")
    
    
    AISiresOfSires = c(AISiresOfSires3, AISiresOfSires2, AISiresOfSires1)
    AISiresOfDams = c(AISiresOfDams3, AISiresOfDams2, AISiresOfDams1)
    males = c(AISiresOfSires, AISiresOfDams, wtRams1, wtRams2, NMSires)
    
    # ---- Select ewes ----
    
    cat("Select ewes\n")
    
    # ---- ... Dams Of Sires ----
    
    cat("... Dams Of Sires", as.character(Sys.time()), "\n")
    
    ewesLact4 = c(damsOfDamsLact4, damsOfSiresLact4)
    ewesLact3 = c(damsOfDamsLact3, damsOfSiresLact3)
    ewesLact2 = c(damsOfDamsLact2, damsOfSiresLact2)
    ewesLact1 = c(damsOfDamsLact1, damsOfSiresLact1)
   
    damsOfSiresLact4 = selectInd(damsOfSiresLact3, nInd = nDamsOfSiresLact4, use = use) # damsOfSiresLact4 are 4 years old here
    damsOfSiresLact3 = selectInd(damsOfSiresLact2, nInd = nDamsOfSiresLact3, use = use) # damsOfSiresLact3 are 3 years old here
    damsOfSiresLact2 = selectInd(damsOfSiresLact1, nInd = nDamsOfSiresLact2, use = use) # damsOfSiresLact2 are 2 years old here
    damsOfSiresLact1 = selectInd(lambs[selNtlLambs], nInd = nDamsOfSiresLact1, use = use, sex = "F", famType = "M") # damsOfSiresLact1 are 1 years old here
    
    
    # Note that damsOfSires are in reality selected only from AI sires only, hence
    #   we should use selectInd(damsOfSiresLact3, ...), but if we select on EBV we
    #   should grab the best females anyway!
    
    # Set phenotypes to missing, because these are copied from the previous lactation.
    #   We do this because of the recordData() call below - that would save wrong
    #   phenotypes for these animals in this life stage! Correct phenotypes for these
    #   animals in this life stage will be added at the beginning of the next year.
    damsOfSiresLact4@pheno[,] = NA
    damsOfSiresLact3@pheno[,] = NA
    damsOfSiresLact2@pheno[,] = NA
    damsOfSiresLact1@pheno[,] = NA
    
    damsOfSires = c(damsOfSiresLact4, damsOfSiresLact3, damsOfSiresLact2, damsOfSiresLact1)
    
    # ---- ... Dams Of Dams ----
    
    cat("... Dams Of Dams", as.character(Sys.time()), "\n")
    
    damsOfDamsLact4 = selectInd(ewesLact3[!ewesLact3@id %in% damsOfSiresLact4@id],
                                nInd = nDamsOfDamsLact4, use = 'rand') # damsOfDamsLact4 are 4 years old here
    
    damsOfDamsLact3 = selectInd(ewesLact2[!ewesLact2@id %in% damsOfSiresLact3@id],
                                nInd = nDamsOfDamsLact3, use = 'rand') # damsOfDamsLact3 are 3 years old here
    
    damsOfDamsLact2 = selectInd(ewesLact1[!ewesLact1@id %in% damsOfSiresLact2@id],
                                nInd = nDamsOfDamsLact2, use = 'rand') # damsOfDamsLact2 are 2 years old here
    
    damsOfDamsLact1 = selectInd(pop = lambs[!lambs@id %in% damsOfSiresLact1@id],
                                nInd = nDamsOfDamsLact1, use = 'rand', sex = "F") # damsOfDamsLact1 are 1 years old here
    # If you are confused about the ages note that:
    #   - these lambs were generated in previous year (year - 1)
    #   - now they are in 1 year old (year)
    #   - lactation phenotype will be generated next year (year + 2)
    #     (at the beginning of the year)
    
    # Set phenotypes to missing, because these are copied from the previous lactation.
    #   We do this because of the recordData() call below - that would save wrong
    #   phenotypes for these animals in this life stage! Correct phenotypes for these
    #   animals in this life stage will be added at the beginning of the next year.
    damsOfDamsLact4@pheno[,] = NA
    damsOfDamsLact3@pheno[,] = NA
    damsOfDamsLact2@pheno[,] = NA
    damsOfDamsLact1@pheno[,] = NA
    
    damsOfDams = c(damsOfDamsLact4, damsOfDamsLact3, damsOfDamsLact2, damsOfDamsLact1)
    females = c(damsOfSires, damsOfDams)
    # ---- Generate lambs ----
    # we are avoiding matings between siblings and between parent offspring
    cat("Generate lambs", as.character(Sys.time()), "\n")
    matingPlan1 <- matrix(NA, ncol = 2, nrow = nDamsOfSires)
    for (i in 1:nDamsOfSires) {
      repeat {
        sire <- sample(AISiresOfSires@id, size = 1)
        if (!areRelated(damsOfSires@id[i], sire, damsOfSires, AISiresOfSires)) {
          matingPlan1[i, ] <- c(damsOfSires@id[i], sire)
          break
        }
      }
    }
    
    n = nLambsFromAISiresOfSires - nDamsOfSires
    damsOfDamsId = damsOfDams@id
    damsOfDamsIdForElite = sample(damsOfDamsId, size = n)
    matingPlan2 <- matrix(NA, ncol = 2, nrow = n)
    for (i in 1:n) {
      repeat {
        sire <- sample(AISiresOfSires@id, size = 1)
        if (!areRelated(damsOfDamsIdForElite[i], sire, damsOfDams, AISiresOfSires)) {
          matingPlan2[i, ] <- c(damsOfDamsIdForElite[i], sire)
          break
        }
      }
    }
    
    n = nDamsOfDams - (nLambsFromAISiresOfSires - nDamsOfSires)
    matingPlan3 = matrix(NA, ncol = 2, nrow = n)
    damsOfDamsIdForRest = sample(damsOfDams@id[!damsOfDams@id %in% damsOfDamsIdForElite])
    siresAll= sample(c(sample(AISiresOfDams@id, size = nLambsFromAIRest, replace = TRUE),
                       sample(wtRams1@id, size = nLambsFromAIForPT, replace = TRUE),
                       sample(NMSires@id, size = nLambsFromNM, replace = TRUE)))
    siresPop = c(AISiresOfDams, wtRams1, NMSires)
    for (i in 1:n) {
      repeat {
        sire <- sample(siresAll, size = 1)
        if (!areRelated(damsOfDamsIdForRest[i], sire, damsOfDams, siresPop)) {
          matingPlan3[i, ] <- c(damsOfDamsIdForRest[i], sire)
          break
        }
      }
    }
    # Note that:
    # 1) sample(pop, size = n, replace = TRUE) gives us n contributions from pop
    # 2) sample(c(sample(), sample(), sample())) shuffles selected male contributions
    #    across female contributions (=random mating)
    
    matingPlan = rbind(matingPlan1, matingPlan2, matingPlan3)
    lambs = makeCross2(females = c(damsOfSires, damsOfDams),
                       males = c(AISiresOfSires, AISiresOfDams, wtRams1, NMSires),
                       crossPlan = matingPlan)
    lambs = fillInMisc(pop = lambs,
                       mothers = c(damsOfSires, damsOfDams),
                       permEnvVar = permVar, year = yearFull)
    
    # ---- Data recording ----
    
    cat("Data recording", as.character(Sys.time()), "\n")
    
    database = recordData(database, pop = AISiresOfSires3, year = yearFull)
    database = recordData(database, pop = AISiresOfSires2, year = yearFull)
    database = recordData(database, pop = AISiresOfSires1, year = yearFull)
    database = recordData(database, pop = AISiresOfSires, year = yearFull)
    database = recordData(database, pop = AISiresOfDams3, year = yearFull)
    database = recordData(database, pop = AISiresOfDams2, year = yearFull)
    database = recordData(database, pop = AISiresOfDams1, year = yearFull)
    database = recordData(database, pop = AISiresOfDams, year = yearFull)
    database = recordData(database, pop = wtRams2, year = yearFull)
    database = recordData(database, pop = wtRams1, year = yearFull)
    database = recordData(database, pop = NMSires, year = yearFull)
    database = recordData(database, pop = damsOfSiresLact4, year = yearFull, lactation = 4)
    database = recordData(database, pop = damsOfSiresLact3, year = yearFull, lactation = 3)
    database = recordData(database, pop = damsOfSiresLact2, year = yearFull, lactation = 2)
    database = recordData(database, pop = damsOfSiresLact1, year = yearFull, lactation = 1)
    # database = recordData(database, pop = damsOfSires, year = yearFull, lactation = ???) # can not save like this since lactations vary
    database = recordData(database, pop = damsOfDamsLact4, year = yearFull, lactation = 4)
    database = recordData(database, pop = damsOfDamsLact3, year = yearFull, lactation = 3)
    database = recordData(database, pop = damsOfDamsLact2, year = yearFull, lactation = 2)
    database = recordData(database, pop = damsOfDamsLact1, year = yearFull, lactation = 1)
    # database = recordData(database, pop = damsOfDams, year = yearFull, lactation = ???) # can not save like this since lactations vary
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
  
  for (year in (nBurninYears + 1):(nBurninYears + nScenarioYears)) {
    # year = yearToDo
    yearFull = startYear + year
    
    cat("Year", year, " (", yearFull, ") ", as.character(Sys.time()), "\n")
    
    yearEffect = sampleEffect(n = 1, var = yearVar)
    herds$herdYearEffect = sampleEffect(n = nHerds, var = herdYearVar)
    
    # ---- Phenotyping ----
    
    cat("Phenotyping", as.character(Sys.time()), "\n")
    
    damsOfSiresLact4 = setPhenoEwe(damsOfSiresLact4, varE = resVar,
                                   mean = meanLact4, yearEffect = yearEffect, herds = herds)
    
    damsOfSiresLact3 = setPhenoEwe(damsOfSiresLact3, varE = resVar,
                                   mean = meanLact3, yearEffect = yearEffect, herds = herds)
    
    damsOfSiresLact2 = setPhenoEwe(damsOfSiresLact2, varE = resVar,
                                   mean = meanLact2, yearEffect = yearEffect, herds = herds)
    
    damsOfSiresLact1 = setPhenoEwe(damsOfSiresLact1, varE = resVar,
                                   mean = meanLact1, yearEffect = yearEffect, herds = herds)
    
    damsOfDamsLact4 = setPhenoEwe(damsOfDamsLact4, varE = resVar,
                                  mean = meanLact4, yearEffect = yearEffect, herds = herds)
    
    damsOfDamsLact3 = setPhenoEwe(damsOfDamsLact3, varE = resVar,
                                  mean = meanLact3, yearEffect = yearEffect, herds = herds)
    
    damsOfDamsLact2 = setPhenoEwe(damsOfDamsLact2, varE = resVar,
                                  mean = meanLact2, yearEffect = yearEffect, herds = herds)
    
    damsOfDamsLact1 = setPhenoEwe(damsOfDamsLact1, varE = resVar,
                                  mean = meanLact1, yearEffect = yearEffect, herds = herds)
    
    database = setDatabasePheno(database, pop = damsOfSiresLact4)
    database = setDatabasePheno(database, pop = damsOfSiresLact3)
    database = setDatabasePheno(database, pop = damsOfSiresLact2)
    database = setDatabasePheno(database, pop = damsOfSiresLact1)
    database = setDatabasePheno(database, pop = damsOfDamsLact4)
    database = setDatabasePheno(database, pop = damsOfDamsLact3)
    database = setDatabasePheno(database, pop = damsOfDamsLact2)
    database = setDatabasePheno(database, pop = damsOfDamsLact1)
    
    
    # ---- Genetic evaluation ----
    
    cat("Genetic evaluation", as.character(Sys.time()), "\n")
    
    # We will remove male lambs from the evaluation - those that will never
    #   contribute to next generations (we will still calculate their parent
    #   average, but later)
    # ... first we find all male lambs
    sel = database$General$Pop == "lambs" & database$General$Sex == "M"
    maleLambs = database$General[sel, "IId"]
    # ... then we find those that were selected for reproduction
    sel = database$General$Pop %in% c("wtRams1", "NMSires")
    reproMales = database$General[sel, "IId"]
    # ... the other were culled
    culledMaleLambs = maleLambs$IId[!maleLambs$IId %in% reproMales$IId]
    removeCulledMaleLambs = row.names(SP$pedigree) %in% culledMaleLambs
    # sum(removeCulledMaleLambs); sum(!removeCulledMaleLambs)
    if ("scenario" %in% args) {
      scenario <- args[grep("^scenario$", args) + 1]
      
      if (scenario == "std"  || scenario == "EFI") {
        variances = list(varPE = permVar,
                         varA  = addVar,
                         varE  = resVar)
        
        pedEbv = estimateBreedingValues(pedigree = SP$pedigree,
                                        database = database,
                                        vars = variances,
                                        removeFromEvaluation = removeCulledMaleLambs)
        
      } else if (scenario == "idl"|| scenario == "idl1" || scenario == "idl2" || scenario == "MA") {
        variances = list(varPE = permVar,
                         varA  = addVar,
                         varIDL = idlVar,
                         varE  = resVar)
        # in the case of the scenarios idl, idl1 and idl2 the estimateBV_IDL will estimate 
        # the IDL, EBV+0.1IDL , EBV + 0.01 IDL, respectively so basically the EBV here is this combination
        pedEbv = estimateBV_IDL(pedigree = SP$pedigree,
                                database = database,
                                vars = variances,
                                removeFromEvaluation = removeCulledMaleLambs)
      } else {
        # Handle other scenarios if needed
        stop("Unknown scenario:", scenario)
      }
    }
    # Set EBVs for every population
    AISiresOfSires3 = setEbv(AISiresOfSires3, ebv = pedEbv)
    AISiresOfSires2 = setEbv(AISiresOfSires2, ebv = pedEbv)
    AISiresOfSires1 = setEbv(AISiresOfSires1, ebv = pedEbv)
    AISiresOfSires  = setEbv(AISiresOfSires, ebv = pedEbv)
    AISiresOfDams3 = setEbv(AISiresOfDams3, ebv = pedEbv)
    AISiresOfDams2 = setEbv(AISiresOfDams2, ebv = pedEbv)
    AISiresOfDams1 = setEbv(AISiresOfDams1, ebv = pedEbv)
    AISiresOfDams = setEbv(AISiresOfDams, ebv = pedEbv)
    wtRams2 = setEbv(wtRams2, ebv = pedEbv)
    wtRams1 = setEbv(wtRams1, ebv = pedEbv)
    NMSires = setEbv(NMSires, ebv = pedEbv)
    damsOfSiresLact4 = setEbv(damsOfSiresLact4, ebv = pedEbv)
    damsOfSiresLact3 = setEbv(damsOfSiresLact3, ebv = pedEbv)
    damsOfSiresLact2 = setEbv(damsOfSiresLact2, ebv = pedEbv)
    damsOfSiresLact1 = setEbv(damsOfSiresLact1, ebv = pedEbv)
    damsOfSires = setEbv(damsOfSires, ebv = pedEbv)
    damsOfDamsLact4 = setEbv(damsOfDamsLact4, ebv = pedEbv)
    damsOfDamsLact3 = setEbv(damsOfDamsLact3, ebv = pedEbv)
    damsOfDamsLact2 = setEbv(damsOfDamsLact2, ebv = pedEbv)
    damsOfDamsLact1 = setEbv(damsOfDamsLact1, ebv = pedEbv)
    damsOfDams = setEbv(damsOfDams, ebv = pedEbv)
    # ... we could use setEbv(lambs), but we have excluded culled lambs from the
    #     evaluation so we will just calculate parent average for all lambs - this
    #     works for pedigree BLUP, but not for genomic BLUP (if some lambs would
    #     have been genotyped)
    # lambs             = setEbv(lambs, ebv = pedEbv)
    selM = match(x = lambs@mother, table = pedEbv$IId)
    selF = match(x = lambs@father, table = pedEbv$IId)
    lambs@ebv = as.matrix((pedEbv[selM, -1] + pedEbv[selF, -1]) / 2)
    
    # Save current EBVs into the database
    database = setDatabaseEbv(database, pop = AISiresOfSires3)
    database = setDatabaseEbv(database, pop = AISiresOfSires2)
    database = setDatabaseEbv(database, pop = AISiresOfSires1)
    database = setDatabaseEbv(database, pop = AISiresOfSires)
    database = setDatabaseEbv(database, pop = AISiresOfDams3)
    database = setDatabaseEbv(database, pop = AISiresOfDams2)
    database = setDatabaseEbv(database, pop = AISiresOfDams1)
    database = setDatabaseEbv(database, pop = AISiresOfDams)
    database = setDatabaseEbv(database, pop = wtRams2)
    database = setDatabaseEbv(database, pop = wtRams1)
    database = setDatabaseEbv(database, pop = NMSires)
    database = setDatabaseEbv(database, pop = damsOfSiresLact4)
    database = setDatabaseEbv(database, pop = damsOfSiresLact3)
    database = setDatabaseEbv(database, pop = damsOfSiresLact2)
    database = setDatabaseEbv(database, pop = damsOfSiresLact1)
    # database = setDatabaseEbv(database, pop = damsOfSires) # we didn't save this pop in the database
    database = setDatabaseEbv(database, pop = damsOfDamsLact4)
    database = setDatabaseEbv(database, pop = damsOfDamsLact3)
    database = setDatabaseEbv(database, pop = damsOfDamsLact2)
    database = setDatabaseEbv(database, pop = damsOfDamsLact1)
    # database = setDatabaseEbv(database, pop = damsOfDams) # we didn't save this pop in the database
    database = setDatabaseEbv(database, pop = lambs)
    damsOfSires = c(damsOfSiresLact4, damsOfSiresLact3, damsOfSiresLact2, damsOfSiresLact1)
    damsOfDams = c(damsOfDamsLact4, damsOfDamsLact3, damsOfDamsLact2, damsOfDamsLact1)
    females = c(damsOfSires, damsOfDams)
    wtRams = c(wtRams1, wtRams2)
    
    # saving the inbreeding of each year
    inbreeding      = read.table('renf90.inb', header=F)
    inbreeding$year = year
    write.table(x = inbreeding, file = "inbreeding.txt", append = TRUE, col.names = FALSE)
    

    animal_data = data.frame(year        = year,
                             id_orig          = database$General$IId,
                             pop         = database$General$Pop,
                             YearOfBirth = database$General$YearOfBirth,
                             EBV_IDL         = database$Ebv,
                             GV          = database$Gv,
                             TBV         = database$General$Tbv.Trait1
    )
    animal_data = animal_data %>%
      filter(pop != "AISiresOfSires", pop != "AISiresOfDams" )
    write.table(x = animal_data, file = "animal_data.txt", append = TRUE, col.names = FALSE)

    animalData = data.frame(year        = year,
                            id_orig     = c(AISiresOfSires3@id,AISiresOfSires2@id,AISiresOfSires1@id,AISiresOfDams3@id,AISiresOfDams2@id,AISiresOfDams1@id,wtRams2@id,wtRams1@id,NMSires@id,damsOfSiresLact4@id, damsOfSiresLact3@id,damsOfSiresLact2@id,damsOfSiresLact1@id,damsOfDamsLact4@id, damsOfDamsLact3@id,damsOfDamsLact2@id,damsOfDamsLact1@id,lambs@id),
                            pop         = rep(c("AISiresOfSires3","AISiresOfSires2","AISiresOfSires1","AISiresOfDams3","AISiresOfDams2","AISiresOfDams1","wtRams2","wtRams1","NMSires","damsOfSiresLact4", "damsOfSiresLact3","damsOfSiresLact2","damsOfSiresLact1","damsOfDamsLact4", "damsOfDamsLact3","damsOfDamsLact2","damsOfDamsLact1","lambs"),
                                              c(AISiresOfSires3@nInd,AISiresOfSires2@nInd,AISiresOfSires1@nInd,AISiresOfDams3@nInd,AISiresOfDams2@nInd,AISiresOfDams1@nInd,wtRams2@nInd,wtRams1@nInd,NMSires@nInd,damsOfSiresLact4@nInd, damsOfSiresLact3@nInd,damsOfSiresLact2@nInd,damsOfSiresLact1@nInd,damsOfDamsLact4@nInd, damsOfDamsLact3@nInd,damsOfDamsLact2@nInd,damsOfDamsLact1@nInd,lambs@nInd)),
                            EBV_IDL     = c(AISiresOfSires3@ebv,AISiresOfSires2@ebv,AISiresOfSires1@ebv,AISiresOfDams3@ebv,AISiresOfDams2@ebv,AISiresOfDams1@ebv,wtRams2@ebv,wtRams1@ebv,NMSires@ebv,damsOfSiresLact4@ebv, damsOfSiresLact3@ebv,damsOfSiresLact2@ebv,damsOfSiresLact1@ebv,damsOfDamsLact4@ebv, damsOfDamsLact3@ebv,damsOfDamsLact2@ebv,damsOfDamsLact1@ebv,lambs@ebv),
                            GV          = c(AISiresOfSires3@gv,AISiresOfSires2@gv,AISiresOfSires1@gv,AISiresOfDams3@gv,AISiresOfDams2@gv,AISiresOfDams1@gv,wtRams2@gv,wtRams1@gv,NMSires@gv,damsOfSiresLact4@gv, damsOfSiresLact3@gv,damsOfSiresLact2@gv,damsOfSiresLact1@gv,damsOfDamsLact4@gv, damsOfDamsLact3@gv,damsOfDamsLact2@gv,damsOfDamsLact1@gv,lambs@gv),
                            TBV         = c(bv(AISiresOfSires3),bv(AISiresOfSires2),bv(AISiresOfSires1),bv(AISiresOfDams3),bv(AISiresOfDams2),bv(AISiresOfDams1),bv(wtRams2),bv(wtRams1),bv(NMSires),bv(damsOfSiresLact4), bv(damsOfSiresLact3),bv(damsOfSiresLact2),bv(damsOfSiresLact1),bv(damsOfDamsLact4), bv(damsOfDamsLact3),bv(damsOfDamsLact2),bv(damsOfDamsLact1),bv(lambs)))

    write.table(x = animalData, file = "animalData.txt", append = TRUE, col.names = FALSE)


    if ("scenario" %in% args) {
      scenario <- args[grep("^scenario$", args) + 1]
      # saving the solution of each year
      if (scenario == "std"  || scenario == "EFI") {
        solution = read.table('solutions.orig', skip=1)
        solution$year = year
        if (!file.exists('solutions.orig.year')) {
          # If the file doesn't exist, write the column names along with the data
          write.table(solution, file = "solutions.orig.year", col.names = TRUE, row.names = FALSE)
        } else {
          # If the file exists, append the data without writing column names
          write.table(solution, file = "solutions.orig.year", col.names = FALSE, row.names = FALSE, append = TRUE)
        }
      }
      else if (scenario == "idl"|| scenario == "idl1" || scenario == "idl2" || scenario == "MA") { 
        solution = read.table('solutions', skip=1)
        solution$year = year
        if (!file.exists('solutions.year')) {
          # If the file doesn't exist, write the column names along with the data
          write.table(solution, file = "solutions.year", col.names = TRUE, row.names = FALSE)
        } else {
          # If the file exists, append the data without writing column names
          write.table(solution, file = "solutions.year", col.names = FALSE, row.names = FALSE, append = TRUE)
        }
        
        # Calculate the correlation so I need to add the real ebv for the idl scenarios
        sol5 <-filter(solution, V2==5)
        colnames(sol5) = c("trait","effect","id_renum","EBV","SE","year")
        id <- read.table('id_renum_orig')
        colnames(id) = c('id_renum','id_orig')
        merged_data <- merge(sol5, id, by='id_renum')
        merged_data$id_orig=as.character(merged_data$id_orig)
        ad_id <- inner_join(animalData,merged_data, by ="id_orig")
        
        ### ----calculate of the accuracy---- ###
        cor_summary_ebv_gv <- ad_id %>%
          group_by(pop) %>%
          summarise(correlation = cor(EBV, GV, use = "complete.obs"))
        cor_summary_ebv_gv$year=year
        if (!file.exists('CorEBV_GV.txt')) {
          # If the file doesn't exist, write the column names along with the data
          write.table(cor_summary_ebv_gv, file = "CorEBV_GV.txt", col.names = TRUE, row.names = FALSE)
        } else {
          # If the file exists, append the data without writing column names
          write.table(cor_summary_ebv_gv, file = "CorEBV_GV.txt", col.names = FALSE, row.names = FALSE, append = TRUE)
        }
        
        cor_summary_ebv_tbv <- ad_id %>%
          group_by(pop) %>%
          summarise(correlation = cor(EBV, TBV, use = "complete.obs"))
        cor_summary_ebv_tbv$year=year
        if (!file.exists('CorEBV_TBV.txt')) {
          # If the file doesn't exist, write the column names along with the data
          write.table(cor_summary_ebv_tbv, file = "CorEBV_TBV.txt", col.names = TRUE, row.names = FALSE)
        } else {
          # If the file exists, append the data without writing column names
          write.table(cor_summary_ebv_tbv, file = "CorEBV_TBV.txt", col.names = FALSE, row.names = FALSE, append = TRUE)
        } 
      }
    }
    
     
    
  if (scenario =="std" || scenario =="EFI" ) {
      correlationEbvTgv = data.frame(year = year,
                                     AISiresOfSires3 = calcAccuracyEbvVsTgv(AISiresOfSires3),
                                     AISiresOfSires2 = calcAccuracyEbvVsTgv(AISiresOfSires2),
                                     AISiresOfSires1 = calcAccuracyEbvVsTgv(AISiresOfSires1),
                                     AISiresOfSires = calcAccuracyEbvVsTgv(AISiresOfSires),
                                     AISiresOfDams3 = calcAccuracyEbvVsTgv(AISiresOfDams3),
                                     AISiresOfDams2 = calcAccuracyEbvVsTgv(AISiresOfDams2),
                                     AISiresOfDams1 = calcAccuracyEbvVsTgv(AISiresOfDams1),
                                     AISiresOfDams = calcAccuracyEbvVsTgv(AISiresOfDams),
                                     wtRams2 = calcAccuracyEbvVsTgv(wtRams2),
                                     wtRams1 = calcAccuracyEbvVsTgv(wtRams1),
                                     wtRams = calcAccuracyEbvVsTgv(wtRams),
                                     NMSires = calcAccuracyEbvVsTgv(NMSires),
                                     damsOfSiresLact4 = calcAccuracyEbvVsTgv(damsOfSiresLact4),
                                     damsOfSiresLact3 = calcAccuracyEbvVsTgv(damsOfSiresLact3),
                                     damsOfSiresLact2 = calcAccuracyEbvVsTgv(damsOfSiresLact2),
                                     damsOfSiresLact1 = calcAccuracyEbvVsTgv(damsOfSiresLact1),
                                     damsOfSires = calcAccuracyEbvVsTgv(damsOfSires),
                                     damsOfDamsLact4 = calcAccuracyEbvVsTgv(damsOfDamsLact4),
                                     damsOfDamsLact3 = calcAccuracyEbvVsTgv(damsOfDamsLact3),
                                     damsOfDamsLact2 = calcAccuracyEbvVsTgv(damsOfDamsLact2),
                                     damsOfDamsLact1 = calcAccuracyEbvVsTgv(damsOfDamsLact1),
                                     damsOfDams = calcAccuracyEbvVsTgv(damsOfDams),
                                     lambs = calcAccuracyEbvVsTgv(lambs))
      write.table(x = correlationEbvTgv, file = "calcAccuracyEbvVsTgv.txt", append = TRUE, col.names = FALSE)
      # 
      correlationEbvTbv = data.frame(year = year,
                                     AISiresOfSires3 = calcAccuracyEbvVsTbv(AISiresOfSires3),
                                     AISiresOfSires2 = calcAccuracyEbvVsTbv(AISiresOfSires2),
                                     AISiresOfSires1 = calcAccuracyEbvVsTbv(AISiresOfSires1),
                                     AISiresOfSires = calcAccuracyEbvVsTbv(AISiresOfSires),
                                     AISiresOfDams3 = calcAccuracyEbvVsTbv(AISiresOfDams3),
                                     AISiresOfDams2 = calcAccuracyEbvVsTbv(AISiresOfDams2),
                                     AISiresOfDams1 = calcAccuracyEbvVsTbv(AISiresOfDams1),
                                     AISiresOfDams = calcAccuracyEbvVsTbv(AISiresOfDams),
                                     wtRams2 = calcAccuracyEbvVsTbv(wtRams2),
                                     wtRams1 = calcAccuracyEbvVsTbv(wtRams1),
                                     wtRams = calcAccuracyEbvVsTbv(wtRams),
                                     NMSires = calcAccuracyEbvVsTbv(NMSires),
                                     males = calcAccuracyEbvVsTbv(males),
                                     damsOfSiresLact4 = calcAccuracyEbvVsTbv(damsOfSiresLact4),
                                     damsOfSiresLact3 = calcAccuracyEbvVsTbv(damsOfSiresLact3),
                                     damsOfSiresLact2 = calcAccuracyEbvVsTbv(damsOfSiresLact2),
                                     damsOfSiresLact1 = calcAccuracyEbvVsTbv(damsOfSiresLact1),
                                     damsOfSires = calcAccuracyEbvVsTbv(damsOfSires),
                                     damsOfDamsLact4 = calcAccuracyEbvVsTbv(damsOfDamsLact4),
                                     damsOfDamsLact3 = calcAccuracyEbvVsTbv(damsOfDamsLact3),
                                     damsOfDamsLact2 = calcAccuracyEbvVsTbv(damsOfDamsLact2),
                                     damsOfDamsLact1 = calcAccuracyEbvVsTbv(damsOfDamsLact1),
                                     damsOfDams = calcAccuracyEbvVsTbv(damsOfDams),
                                     females = calcAccuracyEbvVsTbv(females),
                                     lambs = calcAccuracyEbvVsTbv(lambs))
      # add = ifelse(year == 1, FALSE, TRUE)
      write.table(x = correlationEbvTbv, file = "calcAccuracyEbvVsTbv.txt", append = TRUE , col.names = FALSE )
    }
    
    # ---- Select rams ----
    
    # ---- ...   AI Sires Of Sires----
    
    cat("... AI Sire Of Siress", as.character(Sys.time()), "\n")
    
    AISiresOfSires3 = AISiresOfSires2 # AISiresOfSires3 are 4.5 years old here
    AISiresOfSires2 = AISiresOfSires1 # AISiresOfSires2 are 3.5 years old here
    
    if ("scenario" %in% args) {
      scenario <- args[grep("^scenario$", args) + 1]
      # saving the solution of each year
      
      if (all(wtRams2@father == "0")) {
        AISiresOfSires1 = selectInd(pop = wtRams2, nInd = nAISiresOfSires1, # AISiresOfSires1 are 3.5 years old here
                                    use = use)
      } else  {
        if (scenario == "std" || scenario == "idl"|| scenario == "idl1" || scenario == "idl2" || scenario == "MA")  {
          
          
          AISiresOfSires1 = selectWithinFam(pop = wtRams2, nInd = 1, # AISiresOfSires1 are 3.5 years old here
                                            use = use, famType = "M")
          AISiresOfSires1 = selectInd(pop = AISiresOfSires1, nInd = nAISiresOfSires1,
                                      use = use)
        }
        
        else if ( scenario == "EFI"){
          # -------------------------------------
          # ---- AI sires of sires selection ----
          # -------------------------------------
          # collecting the id of the contemporary dams of sires and of the waiting rams 2
          dams_of_sire <- data.frame(damsOfSires@iid)
          wtrams2 <- data.frame(wtRams2@iid)
          fwrite(dams_of_sire, "dam", col.names = FALSE)
          fwrite(wtrams2, "sire", col.names = FALSE)
          # calculating the relationship coefficient between the waiting rams and all the contemporary females
          system(command = "ifort /save/user/santonios/bin/inbupgf90")
          system(command = "/save/user/santonios/bin/inbupgf90 --pedfile pedigree.txt --sire_file sire --dam_file dam --nrmcoeff")
          # reading the results of the inbupgf90, which represent all the possible combinations between the wtrams2 and the dams of sires 
          # and with their expected future relationship coefficient.
          nrmcoeff <- read.table('sire_dam.nrmcoeff')
          colnames(nrmcoeff) <- c('sire', 'dam', 'coan_coef')
          # summing the expected future relationship coefficient per sire
          nrmcoeff_sum <- nrmcoeff %>% group_by(sire) %>% summarize(mean_coan=mean(coan_coef))
          
          wtram2_ebv <- data.frame(sire = wtRams2@iid, ebv = wtRams2@ebv)
          # merging the ebv
          nrmcoeff_sum_ebv <- merge(nrmcoeff_sum,wtram2_ebv,by="sire")
          colnames(nrmcoeff_sum_ebv) <- c('sire', 'mean_coan', 'ebv')
          # getting the inbreeding depression param (b)
          sol <- read.table('solutions.orig', skip=1)
          b <- as.numeric(tail(sol$V5, n=1))
          # calculating EFI and the adjusted EBV
          nrmcoeff_ebv_EFI_ebva <- nrmcoeff_sum_ebv %>%
            dplyr::mutate(EFI = mean_coan / 2) %>%
            dplyr::mutate(ebva = ebv + (EFI / 2) * b) %>%
            dplyr::transmute(sire, mean_coan, ebv, EFI, ebva)  
          
          # to select the sires of sires 1 based on their adjusted EBV (ebva)
          top_10_iid <- nrmcoeff_ebv_EFI_ebva %>%
            arrange(desc(ebva)) %>%
            slice_head(n = 10)
          
          # Find the indices of the individuals that match top_10_iid
          indices <- which(wtRams2@iid %in% top_10_iid$sire)
          
          # Subset the individuals using the indices
          AISiresOfSires1 <- wtRams2[indices]
        }
      }
    }
    
    # ---- ... AI Sire Of Dams ----
    
    cat("... AI Sire Of Dams", as.character(Sys.time()), "\n")
    
    AISiresOfDams3 = AISiresOfDams2 # AISiresOfDams3 are 5.5 years old here
    AISiresOfDams2 = AISiresOfDams1 # AISiresOfDams2 are 4.5 years old here
    
    if (all(wtRams2@father == "0")) {
      sel = !wtRams2@id %in% AISiresOfSires1@id
      AISiresOfDams1 = selectInd(pop = wtRams2[sel], nInd = nAISiresOfDams1, # AISiresOfDams1 are 3.5 years old here
                                 use = use)
    } else {
      sel = !wtRams2@id %in% AISiresOfSires1@id
      AISiresOfDams1 = selectWithinFam(pop = wtRams2[sel], nInd = 2, # AISiresOfDams1 are 3.5 years old here
                                       use = use, famType = "M")
      AISiresOfDams1 = selectInd(pop = AISiresOfDams1, nInd = nAISiresOfDams1,
                                 use = use)
    }
    cat("Select rams\n")
    
    # ---- ... Selecting males from lambs ----
    # Note that we select wtRams only from lambs from AI Sire Of Sires (not waiting rams) and Dams Of Sires
    selWtRams = lambs@father %in% AISiresOfSires@id & lambs@mother %in% damsOfSires@id
    n = ceiling(nWtRams1 / length(AISiresOfSires@id))
    wtRams2 = wtRams1
    wtRams1 = selectWithinFam(pop = lambs[selWtRams], nInd = n,
                              use = use, famType = "M", sex = "M")
    wtRams1 = selectInd(pop = wtRams1, nInd = nWtRams1,
                        use = use)
    
    selAILambs = lambs@father %in% c(AISiresOfSires@id, AISiresOfDams@id) &
      !lambs@id %in% wtRams1@id
    NMSires =  selectInd(pop = lambs[selAILambs], nInd = nNMSires,
                         use = use, famType = "M", sex = "M")
    
    
    AISiresOfSires = c(AISiresOfSires3, AISiresOfSires2, AISiresOfSires1)
    AISiresOfDams = c(AISiresOfDams3, AISiresOfDams2, AISiresOfDams1)
    males = c(AISiresOfSires, AISiresOfDams, wtRams1, wtRams2, NMSires)
    
    # ---- Select ewes ----
    
    cat("Select ewes\n")
    
    # ---- ... Dams Of Sires ----
    
    cat("... Dams Of Sires", as.character(Sys.time()), "\n")
    
    ewesLact4 = c(damsOfDamsLact4, damsOfSiresLact4)
    ewesLact3 = c(damsOfDamsLact3, damsOfSiresLact3)
    ewesLact2 = c(damsOfDamsLact2, damsOfSiresLact2)
    ewesLact1 = c(damsOfDamsLact1, damsOfSiresLact1)
    
    damsOfSiresLact4 = selectInd(damsOfSiresLact3, nInd = nDamsOfSiresLact4, use = use) # damsOfSiresLact4 are 4 years old here
    damsOfSiresLact3 = selectInd(damsOfSiresLact2, nInd = nDamsOfSiresLact3, use = use) # damsOfSiresLact3 are 3 years old here
    damsOfSiresLact2 = selectInd(damsOfSiresLact1, nInd = nDamsOfSiresLact2, use = use) # damsOfSiresLact2 are 2 years old here
    damsOfSiresLact1 = selectInd(lambs[selNtlLambs], nInd = nDamsOfSiresLact1, use = use, sex = "F", famType = "M") # damsOfSiresLact1 are 1 years old here
    # Note that damsOfSires are in reality selected only from AI sires only, hence
    #   we should use selectInd(damsOfSiresLact3, ...), but if we select on EBV we
    #   should grab the best females anyway!
    
    # Set phenotypes to missing, because these are copied from the previous lactation.
    #   We do this because of the recordData() call below - that would save wrong
    #   phenotypes for these animals in this life stage! Correct phenotypes for these
    #   animals in this life stage will be added at the beginning of the next year.
    damsOfSiresLact4@pheno[,] = NA
    damsOfSiresLact3@pheno[,] = NA
    damsOfSiresLact2@pheno[,] = NA
    damsOfSiresLact1@pheno[,] = NA
    
    damsOfSires = c(damsOfSiresLact4, damsOfSiresLact3, damsOfSiresLact2, damsOfSiresLact1)
    
    # ---- ... Dams Of Dams ----
    
    cat("... Dams Of Dams", as.character(Sys.time()), "\n")
    
    damsOfDamsLact4 = selectInd(ewesLact3[!ewesLact3@id %in% damsOfSiresLact4@id],
                                nInd = nDamsOfDamsLact4, use = "rand") # damsOfDamsLact4 are 4 years old here
    
    damsOfDamsLact3 = selectInd(ewesLact2[!ewesLact2@id %in% damsOfSiresLact3@id],
                                nInd = nDamsOfDamsLact3, use = "rand") # damsOfDamsLact3 are 3 years old here
    
    damsOfDamsLact2 = selectInd(ewesLact1[!ewesLact1@id %in% damsOfSiresLact2@id],
                                nInd = nDamsOfDamsLact2, use = "rand") # damsOfDamsLact2 are 2 years old here
    
    damsOfDamsLact1 = selectInd(pop = lambs[!lambs@id %in% damsOfSiresLact1@id],
                                nInd = nDamsOfDamsLact1, use = "rand", sex = "F") # damsOfDamsLact1 are 1 years old here
    
    # Set phenotypes to missing, because these are copied from the previous lactation.
    #   We do this because of the recordData() call below - that would save wrong
    #   phenotypes for these animals in this life stage! Correct phenotypes for these
    #   animals in this life stage will be added at the beginning of the next year.
    damsOfDamsLact4@pheno[,] = NA
    damsOfDamsLact3@pheno[,] = NA
    damsOfDamsLact2@pheno[,] = NA
    damsOfDamsLact1@pheno[,] = NA
    
    damsOfDams = c(damsOfDamsLact4, damsOfDamsLact3, damsOfDamsLact2, damsOfDamsLact1)
    females = c(damsOfSires, damsOfDams)
    
    # ---- Generate lambs ----
    
    cat("Generate lambs", as.character(Sys.time()), "\n")
    matingPlan1 <- matrix(NA, ncol = 2, nrow = nDamsOfSires)
    if ("scenario" %in% args) {
      scenario <- args[grep("^scenario$", args) + 1]
      
      if (scenario == "idl" || scenario == "idl1" || scenario == "idl2" || scenario == "std"  || scenario == "EFI") {
        for (i in 1:nDamsOfSires) {
          repeat {
            sire <- sample(AISiresOfSires@id, size = 1)
            if (!areRelated(damsOfSires@id[i], sire, damsOfSires, AISiresOfSires)) {
              matingPlan1[i, ] <- c(damsOfSires@id[i], sire)
              break
            }
          }
        }
      }
      else if (scenario == "MA"){
    
datafile <- read.table(file="renadd05.ped", col.names = c("id_renum","sire","dam","v4","v5","v6","v7","v8","v9","id_orig"))
        # anim_id_pop = data.frame(id_orig     = c(AISiresOfSires3@id,AISiresOfSires2@id,AISiresOfSires1@id,AISiresOfDams3@id,AISiresOfDams2@id,AISiresOfDams1@id,wtRams2@id,wtRams1@id,NMSires@id,damsOfSiresLact4@id, damsOfSiresLact3@id,damsOfSiresLact2@id,damsOfSiresLact1@id,damsOfDamsLact4@id, damsOfDamsLact3@id,damsOfDamsLact2@id,damsOfDamsLact1@id,lambs@id),      
        #                          pop         = rep(c("AISiresOfSires3","AISiresOfSires2","AISiresOfSires1","AISiresOfDams3","AISiresOfDams2","AISiresOfDams1","wtRams2","wtRams1","NMSires","damsOfSiresLact4", "damsOfSiresLact3","damsOfSiresLact2","damsOfSiresLact1","damsOfDamsLact4", "damsOfDamsLact3","damsOfDamsLact2","damsOfDamsLact1","lambs"), 
        #                                            c(AISiresOfSires3@nInd,AISiresOfSires2@nInd,AISiresOfSires1@nInd,AISiresOfDams3@nInd,AISiresOfDams2@nInd,AISiresOfDams1@nInd,wtRams2@nInd,wtRams1@nInd,NMSires@nInd,damsOfSiresLact4@nInd, damsOfSiresLact3@nInd,damsOfSiresLact2@nInd,damsOfSiresLact1@nInd,damsOfDamsLact4@nInd, damsOfDamsLact3@nInd,damsOfDamsLact2@nInd,damsOfDamsLact1@nInd,lambs@nInd))
        #                          )
        anim_id_pop = data.frame(id_orig     = c(AISiresOfSires3@id,AISiresOfSires2@id,AISiresOfSires1@id,damsOfSiresLact4@id, damsOfSiresLact3@id,damsOfSiresLact2@id,damsOfSiresLact1@id), 
                                 pop         = rep(c("AISiresOfSires3","AISiresOfSires2","AISiresOfSires1","damsOfSiresLact4", "damsOfSiresLact3","damsOfSiresLact2","damsOfSiresLact1"),
                                                   c(AISiresOfSires3@nInd,AISiresOfSires2@nInd,AISiresOfSires1@nInd,damsOfSiresLact4@nInd, damsOfSiresLact3@nInd,damsOfSiresLact2@nInd,damsOfSiresLact1@nInd))
        )
        dams <- filter(anim_id_pop, pop=="damsOfSiresLact4" |  pop=="damsOfSiresLact3" |  pop=="damsOfSiresLact2" |  pop=="damsOfSiresLact1")
        ds <- damsOfSires@id
        write.table(dams, file = 'dams')
        write.table(ds,file='ds')

        datafile$id_renum = as.character(datafile$id_renum)
        datafile$id_orig  = as.character(datafile$id_orig)
        anim_id_pop$id_orig=as.character(anim_id_pop$id_orig)

        kfile <- read.table(file="K_id_renum.txt", col.names = c("id_renum","coef1","coef2","value1","anc1","value2","anc2","value3","anc3","value4","anc4","value5","anc5","value6","anc6","value7","anc7"))
        kfile$id_renum = as.character(kfile$id_renum)
        toto <- left_join(datafile,kfile, by = "id_renum") #%>%
        # left_join(., anim_id_pop, by = "id_orig")

        toto1 <- merge(toto, anim_id_pop, by = "id_orig")
        toto2 <- toto1 %>%
          select(c("id_renum","sire","dam","anc1","anc2","anc3","anc4","anc5","anc6","anc7","pop","id_orig"))

        males <- toto2 %>% filter(pop=="AISiresOfSires1"|pop=="AISiresOfSires2"|pop=="AISiresOfSires3")
        colnames(males) <- c("male_id","sire","dam","anc1","anc2","anc3","anc4","anc5","anc6","anc7","pop","male_id_orig")


        females <- toto2 %>% filter(pop=="damsOfSiresLact4"|pop=="damsOfSiresLact3"|pop=="damsOfSiresLact2"|pop=="damsOfSiresLact1")
        colnames(females) <- c("female_id","sire","dam","anc1","anc2","anc3","anc4","anc5","anc6","anc7","pop","female_id_orig")

        # Initialize an empty list to store results
        matings <- list()
        for (i in 1:nrow(males)) {
          for (j in 1:nrow(females)) {
            # Check if male and female are not related
            if (!(males[i, "dam"] == females[j, "female_id"] ||
                  males[i, "male_id"] == females[j, "sire"] ||
                  (males[i, "sire"] != 0 && females[j, "sire"] != 0 && males[i, "sire"] == females[j, "sire"] &&
                   males[i, "dam"] != 0 && females[j, "dam"] != 0 && males[i, "dam"] == females[j, "dam"]))) {

              # Extract ancestors starting from ancestor1
              male_ancestors <- unlist(replace(males[i, 4:10], is.na(males[i, 4:10]), 0))
              female_ancestors <- unlist(replace(females[j, 4:10], is.na(females[j, 4:10]), 0))

            } else if ((males[i, "dam"] == females[j, "female_id"] ||
                        males[i, "male_id"] == females[j, "sire"] ||
                        (males[i, "sire"] != 0 && females[j, "sire"] != 0 && males[i, "sire"] == females[j, "sire"] &&
                         males[i, "dam"] != 0 && females[j, "dam"] != 0 && males[i, "dam"] == females[j, "dam"]))) {

              # Set ancestors to 'A'
              male_ancestors <- rep('A', 7)
              female_ancestors <- rep('A', 7)
            }

            # Find common ancestors
            common <- intersect(male_ancestors, female_ancestors)

            # Create row data frame
            row <- data.frame(
              male_id = males[i, "male_id"],
              female_id = females[j, "female_id"],
              male_id_orig = males[i, "male_id_orig"],
              female_id_orig = females[j, "female_id_orig"],
              common_ancestor1 = ifelse(length(common) >= 1, common[1], "0"),
              common_ancestor2 = ifelse(length(common) >= 2, common[2], "0"),
              common_ancestor3 = ifelse(length(common) >= 3, common[3], "0"),
              common_ancestor4 = ifelse(length(common) >= 4, common[4], "0"),
              common_ancestor5 = ifelse(length(common) >= 5, common[5], "0"),
              common_ancestor6 = ifelse(length(common) >= 6, common[6], "0"),
              common_ancestor7 = ifelse(length(common) >= 7, common[7], "0"),
              stringsAsFactors = FALSE
            )

            # Append row to matings list
            matings <- c(matings, list(row))
          }
        }

        # Combine all rows into a single data frame
        matings <- do.call(rbind, matings)
        # Read the BLUP solutions file
        blup_sol <- read.table(file = "solutions", skip = 1)

        # Extract IDL information for animals (assuming effect number is 12)
        idl_ind <- blup_sol %>%
          filter(V2 == 12) %>%  # Change effect number to your animal effect number if different
          select(id_renum = V3, IDL = V4)
        # Initialize IDL columns in matings
        # if (!is.data.frame(matings)) {
        #   matings <- as.data.frame(matings)
        # }
        # Assign -100000 to the IDL columns
        matings$IDL1 <- -100000
        matings$IDL2 <- -100000
        matings$IDL3 <- -100000
        matings$IDL4 <- -100000
        matings$IDL5 <- -100000
        matings$IDL6 <- -100000
        matings$IDL7 <- -100000

        # Define idl_cols as a vector of IDL column names
        idl_cols <- paste0("IDL", 1:7)

        # Identify indices where common_ancestor1 is not 'A'
        non_A_indices <- matings$common_ancestor1 != "A"

        # Check if there are any such indices
        if (any(non_A_indices)) {
          # Merge IDL information into the matings data frame for all common ancestors
          for (i in 1:7) {
            matings[[idl_cols[i]]][non_A_indices] <- idl_ind$IDL[match(matings[[paste0("common_ancestor", i)]][non_A_indices], idl_ind$id_renum)]
          }
        }
        # Subset matings where common_ancestor1 is not 'A'
        non_A_indices <- matings$common_ancestor1 != "A"
        if (any(non_A_indices)) {
          # Merge IDL information into the matings data frame
          for (i in 1:7) {
            matings[[idl_cols[i]]][non_A_indices] <- idl_ind$IDL[match(matings[[paste0("common_ancestor", i)]][non_A_indices], idl_ind$id_renum)]
          }
        }


        # Replace NA with 0 for consistency
        matings[is.na(matings)] <- 0

        # Calculate total IDL
        matings_idl <- matings %>% transmute(male_id_orig,female_id_orig,total_idl = IDL1 + IDL2 + IDL3 + IDL4 + IDL5 + IDL6 + IDL7)

        matings_idl$total_idl[matings_idl$total_idl == 0] <- 100000
        write.table(matings_idl, 'matings_idl.txt')
        #matings_idl = matings %>% select('male_id','female_id','total_idl')
        obj.fun.expIDL <- matings_idl[,3]

        # save.image(file = "scenario1.RData")

        # nAISiresOfSires <- 30  ##number of males
        # nDamsOfSires <- 3000 ##optimise females, see parameters above

        nfpm <- nDamsOfSires/nAISiresOfSires ##(nfemales/nmales) #number of females per male (ratio)
        total_matings <- nDamsOfSires  #total number of matings equal to number of females

        constr.expIDL <- matrix(0, nDamsOfSires+nAISiresOfSires, nDamsOfSires*nAISiresOfSires)

        for(i in 1:nAISiresOfSires){
          for(j in 1:nDamsOfSires){
            constr.expIDL[i, nDamsOfSires*(i-1) + j] <- 1
            constr.expIDL[nAISiresOfSires+j, nDamsOfSires*(i-1) + j] <- 1
          }
        }

        ## adding constraints for number of matings equal to nfemales
        constr.expIDL <- rbind(constr.expIDL,rep(1,nDamsOfSires*nAISiresOfSires))

        ## adding direction at the end related to number of females to mate.
        constr.expIDL.dir <- c(rep("<=", nAISiresOfSires), rep("<=", nDamsOfSires),"=")

        ### adding constraint at the end related at the maximun number of females to mate.
        rhs.expIDL <- c(rep(nfpm, nAISiresOfSires), rep(1, nDamsOfSires), total_matings)

        #MINIMIZING EXPECTED IDL
        time <- proc.time()
        set.seed(123)  # For reproducibility
        options(scipen = 999)
        epsilon <- 1e-6  # Small value for perturbation
        random_factors <- 1 + runif(length(obj.fun.expIDL), -epsilon, epsilon)

        obj.fun.expIDL.perturbed <- obj.fun.expIDL * random_factors
        # save.image(file = "scenario2.RData")

        optim.expIDL <- lp("max", obj.fun.expIDL.perturbed,
                           constr.expIDL, constr.expIDL.dir, rhs.expIDL,
                           compute.sens = TRUE, all.bin=TRUE)
        print("optimization consumption time")

        proc.time() - time


        # save.image(file = "scenario3.RData")

        objval.expIDL <- optim.expIDL$objval #objective function value
        objval.expIDL

        sol.expIDL.matrix <- matrix(optim.expIDL$solution, nAISiresOfSires, nDamsOfSires, byrow=TRUE) #decision variables values, solution vector
        sire <- matings_idl$male_id_orig
        dam = matings_idl$female_id_orig[1:3000]
        v <- seq(1, (nDamsOfSires*nAISiresOfSires), by = nDamsOfSires)
        rownames(sol.expIDL.matrix) <- sire[v]
        colnames(sol.expIDL.matrix) <- dam[1:nDamsOfSires]
        positions <- which(sol.expIDL.matrix == 1, arr.ind = TRUE)
 # Create the resulting data frame
        matingPlan1_dataframe <- data.frame(
          Sire = rownames(sol.expIDL.matrix)[positions[, 1]],
          Dam = colnames(sol.expIDL.matrix)[positions[, 2]]
        )

        write.table(matingPlan1_dataframe, 'matingPlan1_dataframe')
        matingPlan1_dataframe <- matingPlan1_dataframe[, c("Dam", "Sire")]

        # Convert the reordered data frame to a matrix
        matingPlan1 <- as.matrix(matingPlan1_dataframe)

        # Remove row and column names
        rownames(matingPlan1) <- NULL
        colnames(matingPlan1) <- NULL

      }
    }


    n = nLambsFromAISiresOfSires - nDamsOfSires
    damsOfDamsId = damsOfDams@id
    damsOfDamsIdForElite = sample(damsOfDamsId, size = n)
    matingPlan2 <- matrix(NA, ncol = 2, nrow = n)
    for (i in 1:n) {
      repeat {
        sire <- sample(AISiresOfSires@id, size = 1)
        if (!areRelated(damsOfDamsIdForElite[i], sire, damsOfDams, AISiresOfSires)) {
          matingPlan2[i, ] <- c(damsOfDamsIdForElite[i], sire)
          break
        }
      }
    }
    
    n = nDamsOfDams - (nLambsFromAISiresOfSires - nDamsOfSires)
    matingPlan3 = matrix(NA, ncol = 2, nrow = n)
    damsOfDamsIdForRest = sample(damsOfDams@id[!damsOfDams@id %in% damsOfDamsIdForElite])
    siresAll= sample(c(sample(AISiresOfDams@id, size = nLambsFromAIRest, replace = TRUE),
                       sample(wtRams1@id, size = nLambsFromAIForPT, replace = TRUE),
                       sample(NMSires@id, size = nLambsFromNM, replace = TRUE)))
    siresPop = c(AISiresOfDams, wtRams1, NMSires)
    for (i in 1:n) {
      repeat {
        sire <- sample(siresAll, size = 1)
        if (!areRelated(damsOfDamsIdForRest[i], sire, damsOfDams, siresPop)) {
          matingPlan3[i, ] <- c(damsOfDamsIdForRest[i], sire)
          break
        }
      }
    }
    # Note that:
    # 1) sample(pop, size = n, replace = TRUE) gives us n contributions from pop
    # 2) sample(c(sample(), sample(), sample())) shuffles selected male contributions
    #    across female contributions (=random mating)
    
    matingPlan = rbind(matingPlan1, matingPlan2, matingPlan3)
    lambs = makeCross2(females = c(damsOfSires, damsOfDams),
                       males = c(AISiresOfSires, AISiresOfDams, wtRams1, NMSires),
                       crossPlan = matingPlan)
    lambs = fillInMisc(pop = lambs,
                       mothers = c(damsOfSires, damsOfDams),
                       permEnvVar = permVar, year = yearFull)
    
    # ---- Data recording ----
    
    cat("Data recording", as.character(Sys.time()), "\n")
    
    database = recordData(database, pop = AISiresOfSires3, year = yearFull)
    database = recordData(database, pop = AISiresOfSires2, year = yearFull)
    database = recordData(database, pop = AISiresOfSires1, year = yearFull)
    database = recordData(database, pop = AISiresOfSires, year = yearFull)
    database = recordData(database, pop = AISiresOfDams3, year = yearFull)
    database = recordData(database, pop = AISiresOfDams2, year = yearFull)
    database = recordData(database, pop = AISiresOfDams1, year = yearFull)
    database = recordData(database, pop = AISiresOfDams, year = yearFull)
    database = recordData(database, pop = wtRams2, year = yearFull)
    database = recordData(database, pop = wtRams1, year = yearFull)
    database = recordData(database, pop = NMSires, year = yearFull)
    database = recordData(database, pop = damsOfSiresLact4, year = yearFull, lactation = 4)
    database = recordData(database, pop = damsOfSiresLact3, year = yearFull, lactation = 3)
    database = recordData(database, pop = damsOfSiresLact2, year = yearFull, lactation = 2)
    database = recordData(database, pop = damsOfSiresLact1, year = yearFull, lactation = 1)
    # database = recordData(database, pop = damsOfSires, year = yearFull, lactation = ???) # can not save like this since lactations vary
    database = recordData(database, pop = damsOfDamsLact4, year = yearFull, lactation = 4)
    database = recordData(database, pop = damsOfDamsLact3, year = yearFull, lactation = 3)
    database = recordData(database, pop = damsOfDamsLact2, year = yearFull, lactation = 2)
    database = recordData(database, pop = damsOfDamsLact1, year = yearFull, lactation = 1)
    # database = recordData(database, pop = damsOfDams, year = yearFull, lactation = ???) # can not save like this since lactations vary
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
