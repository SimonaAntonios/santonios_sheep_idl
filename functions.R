# source(file = "functions.R")

sampleEffect = function(n, var) {
  as.matrix(rnorm(n = n, sd = sqrt(var)))
}

getHerd = function(pop) {
  # Getting herd information for individuals of a population
  # pop \code{\link{Pop-class}}
  ret = sapply(X = pop@misc, FUN = function(x) x$herd)
  return(ret)
}

getYob = function(pop) {
  # Getting year of birth information for individuals of a population
  # pop \code{\link{Pop-class}}
  ret = sapply(X = pop@misc, FUN = function(x) x$yearOfBirth)
  return(ret)
}

getPermEnvEffect = function(pop) {
  # Getting permanent environment effect for individuals of a population
  # pop \code{\link{Pop-class}}
  ret = sapply(X = pop@misc, FUN = function(x) x$permEnvEffect)
  return(ret)
}

getIDL = function(pop) {
  # Getting inbreeding depression load for individuals of a population
  # pop \code{\link{Pop-class}}
  ret = sapply(X = pop@misc, FUN = function(x) x$idl)
  return(ret)
}

getIIdPop = function(pop, popObject = NULL) {
  # pop \code{\link{Pop-class}}
  if (is.null(popObject)) {
    popObject = deparse(substitute(pop))
  }
  ret = paste(pop@id, popObject, sep = "-")
  return(ret)
}

areRelated <- function(id1, id2, pop1, pop2) {
  # this function is to retains individuals that are related like follow
  # if id1 that is from pop1 and id2 from pop2 are siblings or parent-offspring 
  mother1 <- pop1@mother[pop1@id == id1]
  father1 <- pop1@father[pop1@id == id1]
  
  mother2 <- pop2@mother[pop2@id == id2]
  father2 <- pop2@father[pop2@id == id2]
  
  # Check if individuals are siblings
  if (mother1 != 0 && mother2 != 0 && mother1 == mother2 &&
      father1 != 0 && father2 != 0 && father1 == father2) {
    # print('siblings') # Siblings
    return(TRUE)
  }
  # Check if individuals have a parent-offspring relationship
  if ((id1 %in% c(mother2, father2)) || (id2 %in% c(mother1, father1))) {
    # print('Parent-ofspring')
    # print(paste("id1", id1,"id2", id2))
    return(TRUE)  # Parent-ofspring
  }
  # print('non-related')# Not related
  return(FALSE)
}

fillInMisc = function(pop, mothers = NULL, herds = NULL, permEnvVar = NULL,
                      year = NA) {
  # Fill in the misc slot of a population
  # pop \code{\link{Pop-class}}
  # mothers population
  # herds list, with herd and herdSize nodes
  # startYear numeric, starting year
  # year numeric, current year
  # permEnvEffectVar numeric, permanent environment variance
  miscNULL = list(herd = NA, permEnvEffect = NA, yearOfBirth = NA)
  n = nInd(pop)
  if (!is.null(herds)) {
    herd = sample(x = herds$herd, size = n, replace = TRUE, prob = herds$herdSize)
  }
  if (!is.null(mothers)) {
    matchId = match(x = pop@mother, table = mothers@id)
    herd = getHerd(pop = mothers)[matchId]
  }
  if (!is.null(permEnvVar)) {
    permEnvEffect = rnorm(n = n, sd = sqrt(permEnvVar))
  }
  for (ind in 1:nInd(pop)) {
    pop@misc[[ind]] = miscNULL
    if (!is.null(herds) | !is.null(mothers)) {
      pop@misc[[ind]]$herd = herd[ind]
    }
    if (!is.null(permEnvVar)) {
      pop@misc[[ind]]$permEnvEffect = permEnvEffect[ind]
    }
    pop@misc[[ind]]$yearOfBirth = year
  }
  return(pop)
}

setPhenoEwe = function(pop, varE, mean, herds, yearEffect) {
  # Create a complex ewe phenotype
  # pop \code{\link{Pop-class}}
  # varE numeric, environmental variance
  # mean numeric, population mean (such as lactation mean)
  # herds list, holding herd, herdEffect (vector or matrix), and
  #   herd-year effect (vector or matrix) nodes
  # yearEffect numeric, year effect (scalar or vector)
  # traitMask matrix, specify which animals should have which pheno NA
  pop = setPheno(pop, varE = varE) # genetics + environment/residual
  herd = getHerd(pop)
  if (!is.matrix(mean)) {
    mean = as.matrix(mean)
  }
  if (!is.matrix(yearEffect)) {
    yearEffect = as.matrix(yearEffect)
  }
  pop@pheno = rep(1, times = nInd(pop)) %*% mean +
    pop@pheno +
    rep(1, times = nInd(pop)) %*% yearEffect +
    herds$herdEffect[herd, ] +
    herds$herdYearEffect[herd, ] +
    getPermEnvEffect(pop)
  return(pop)
}

recordData = function(database = NULL, pop = NULL, year, lactation = NA, label = NA) {
  # Record data from a population into a database
  # This database could be improved, but the logic here is that we are saving
  # records from an individual over the course of it's life so we will have multiple
  # records from different life stages. We will also be saving multiple estimates
  # of breeding values so we can see how these change over life stages!
  # pop \code{\link{Pop-class}}
  # TODO document other variables
  if (!is.null(pop)) {
    popObject = deparse(substitute(pop))
    yob = getYob(pop)
    herd = getHerd(pop)
    herdYear = paste(herd, year, sep = "_")
    tmp = list(General = data.table(IIdPop      = paste0(pop@id, "-", popObject),
                                    Pop         = popObject,
                                    Label       = label,
                                    IId         = pop@id,
                                    FId         = pop@father,
                                    MId         = pop@mother,
                                    Year        = year,
                                    YearOfBirth = yob,
                                    Sex         = pop@sex,
                                    Lactation   = lactation,
                                    Herd        = herd,
                                    HerdYear    = herdYear,
                                    Tbv   = bv(pop)),
               Pheno = pop@pheno,
               Gv    = pop@gv,
               Ebv   = pop@ebv)
    if (ncol(tmp$Ebv) == 0) {
      tmp$Ebv = matrix(data = as.numeric(NA),
                       nrow = nrow(pop@gv), ncol = ncol(pop@gv))
    }
    row.names(tmp$Pheno) = tmp$General$IIdPop
    row.names(tmp$Gv) = tmp$General$IIdPop
    row.names(tmp$Ebv) = tmp$General$IIdPop
    if (is.null(database)) {
      database = tmp
    } else {
      database = list(General = rbind(database$General, tmp$General),
                      Pheno = rbind(database$Pheno, tmp$Pheno),
                      Gv    = rbind(database$Gv, tmp$Gv),
                      Ebv   = rbind(database$Ebv, tmp$Ebv))
    }
    return(database)
  } else { # pop is NULL
    if (is.null(database)) {
      stop("Provide at least a non-NULL database or pop!")
    } else {
      return(database)
    }
  }
}

estimateBreedingValues = function(pedigree, database,
                                  trait = 1, na = -999, vars,
                                  removeFromEvaluation = NULL,
                                  inbLoadModel = FALSE) {
  # Estimate breeding values with external software - here we use blupf90.
  # pedigree SP$pedigree object from AlphaSimR
  # database list
  #   * General data.table with columns IId, Year, Herd, HerdYear
  #   * Pheno matrix with row.names equal to IId
  #   * Gv matrix with row.names equal to IId
  #   * Ebv matrix with row.names equal to IId
  # trait numeric, indicating which traits to analyse (one or more values,
  #   so, with two traits we have options: trait = 1, trait = 2, or trait = 1:2
  # na value used to denote missing value, say -999
  # vars list, variance components VarPE, VarA, VarIL, and VarE - vectors or matrices
  # removeFromEvaluation logical, vector of length equal to nrow(pedigree)
  #   indicating which animals should be removed from the evaluation (to speed
  #   things up)
  # inbLoadModel logical, are we running the standard model or model with
  #   inbreeding depression load (not yet used)

  # ---- Prepare pedigree file ----
  pedigree = data.table(IId = rownames(pedigree),
                        FId = pedigree[, "father"],
                        MId = pedigree[, "mother"])
  # Remove young males that will never have any progeny
  if (!is.null(removeFromEvaluation)) {
    pedigree = pedigree[!removeFromEvaluation, ]
  }
  fwrite(x = pedigree, file = "pedigree.txt", quote = FALSE,
    row.names = FALSE, col.names = FALSE, sep = " ")
  # rm(pedigree)

  # ---- Prepare phenotype file ----
  sel = rowSums(!is.na(database$Pheno[, trait, drop = FALSE])) > 0
  tmp = database$Pheno[sel, , drop = FALSE]
  phenotypes = data.table(IIdPop = row.names(tmp))
  phenotypes = cbind(phenotypes, tmp)
  colTrait = paste0("Pheno", trait)
  colnames(phenotypes) = c("IIdPop", colTrait)
  phenotypes = merge(x = phenotypes,
                     y = database$General[, c("IIdPop", "IId", "Year", "Herd", "Lactation")])
  sel = c("IId", "Year", "Herd", "Lactation", colTrait)
  fwrite(x = phenotypes[, ..sel], file = "phenotype.txt", quote = FALSE,
         row.names = FALSE, col.names = FALSE, sep = " ")
  rm(phenotypes)

  # ---- Prepare parameter file ----
  sink(file = "renum.par", type = "output")
  cat("# Renumf90 parameter file\n",
      "# herd-year\n",
      "COMBINE 6 2 3\n",
      "DATAFILE\n",
      "phenotype.txt\n",
      "TRAITS\n",
      "# milk yield\n",
      "5\n",
      "FIELDS_PASSED TO OUTPUT\n",
      "# id herd year lactation\n",
      "1 2 3 4\n",
      "WEIGHT(S)\n",
      "\n",
      "RESIDUAL_VARIANCE\n",
      vars$varE, "\n",
      "# herd\n",
      "EFFECT\n",
      "3 cross alpha\n",
      "# year\n",
      "EFFECT\n",
      "2 cross alpha\n",
      "# lactation\n",
      "EFFECT\n",
      "4 cross alpha\n",
      "# herd-year\n",
      "EFFECT\n",
      "6 cross alpha\n",
      "#animal\n",
      "EFFECT\n",
      "1 cross alpha\n",
      "RANDOM\n",
      "animal\n",
      "OPTIONAL\n",
      "pe\n",
      "FILE\n",
      "pedigree.txt\n",
      "FILE_POS\n",
      "1 2 3 0 0\n",
      "PED_DEPTH\n",
      "0\n",
      "INBREEDING\n",
      "pedigree.txt\n",
      "(CO)VARIANCES\n",
      vars$varA, "\n",
      "(CO)VARIANCES_PE\n",
      vars$varPE, "\n",
      "OPTION origID\n",
      "OPTION saveAinv\n",
      "OPTION saveAinvOrig\n",
    sep = "")
  sink()

  # ---- Run the command ----
  system(command = "echo renum.par | renumf90 | tee renum.log")
  system(command = "echo renf90.par | blupf90+ | tee blup.log")
  updateRenf90Par <- function(input_file, output_file) {
    # Read the existing parameter file
    parameters <- readLines(input_file)
    
    
    # Functions to update BLUPF90 parameter file with the new one including the IDL
    
    # this function is to replace certain values and lines
    ReplacingValuesRenf90 <- function(input_file, output_file, replacements) {
      # Read the existing parameter file
      parameters <- readLines(input_file)
      
      # Replace each old line with the corresponding new line
      for (replacement in replacements) {
        old_line <- replacement$old_line
        new_line <- replacement$new_line
        
        # Find the index of the line to be replaced
        line_to_replace_index <- which(parameters == old_line)
        
        # Replace the line with the new content
        parameters[line_to_replace_index] <- new_line
      }
      
      # Write the updated parameter file
      writeLines(parameters, output_file)
    }
    
    # Example usage
    input_file <- "renf90.par"
    output_file <- "updated_renf90.par"
    
    # Define replacements as a list of old and new lines
    replacements <- list(
      list(old_line = " renf90.dat", new_line = " renf90.dat.F"),
      list(old_line = "           6", new_line = "           10"),
      list(old_line = "     5", new_line = "     5  6"),
      list(old_line = "   500.00    ", new_line =paste("   500.00    ", round(covIDLAdd, digits = 2))),
      list(old_line = "     6", new_line = "     9")
      # Add more pairs of old and new lines as needed
    )
    
    ReplacingValuesRenf90(input_file, output_file, replacements)
    
    # this function is to add lines
    updateRenf90Par <- function(inputFile, outputFile) {
      # Read the existing parameter file
      parameters <- readLines(inputFile)
      
      # Find the position to add the NUMBER_OF_EFFECTS information
      numberOfEffectsPosition <- grep("NUMBER_OF_EFFECTS", parameters) + 1
      
      # Add the  information
      updatedParameters <- c(
        parameters[1:numberOfEffectsPosition],
        paste("10"),
        parameters[(numberOfEffectsPosition + 1):(numberOfEffectsPosition + 10)],
        paste("# cov1 nested in ancestor1 (in column 14)"),
        paste("13         0 cov 14"),
        paste("15         0 cov 16"),
        paste("17     ", nrow(pedigree), "cov 18", sep=""),
        paste("# permanent effect"),
        parameters[(numberOfEffectsPosition + 11):(numberOfEffectsPosition + 11)],
        paste("# total inbreeding"),
        paste("19         1 cov"),
        parameters[(numberOfEffectsPosition + 12):(numberOfEffectsPosition + 21)],
        paste("", round(covIDLAdd, digits = 2), idlVar, sep = "    "),
        parameters[(numberOfEffectsPosition + 22):length(parameters)]
      )
      # Write the updated parameter file
      writeLines(updatedParameters, outputFile)
    }
    
    
    # Example usage
    inputFile <- "updated_renf90.par"
    outputFile <- "renf90.par.idl"
    
    updateRenf90Par(inputFile, outputFile)
    
  
  # ---- Read in the solutions ----
  blup_sol = read_table(file = "solutions.orig",
    col_types = cols(
      .default = col_double(),
      trait = col_integer(),
      effect = col_integer(),
      level = col_integer(),
      original_id = col_character(),
      solution = col_double()))

  # ---- Extracting EBV from the file ----
  ebv_ind = blup_sol %>%
    filter(effect == 5) %>%            # Change effect number to your animal effect number
    select("original_id", "solution")
  colnames(ebv_ind) = c("IId", "EBV")
  return(ebv_ind)
}

setEbv = function(pop, ebv, trait = 1) {
  # Set EBV for a population
  # pop \code{\link{Pop-class}}
  # ebv data.table with columns IId and Ebv (one or more columns)
  # trait numeric, indicating which traits to set (one or more values
  #   so, with two traits options are: trait = 1, trait = 2, or trait = 1:2)
  if (ncol(pop@ebv) == 0) {
    pop@ebv = matrix(data = NA,
                     nrow = pop@nInd,
                     ncol = pop@nTraits)
  }
  matchId = match(x = pop@id, table = ebv$IId)
  traitDataTableCol = 1 + trait
  # pop@ebv[, trait] = as.matrix(ebv[matchId, ..traitDataTableCol]) # be careful ebv is data.table
  pop@ebv[, trait] = as.matrix(ebv[matchId, traitDataTableCol]) # be careful ebv is data.frame
  return(pop)
}

setDatabasePheno = function(database, pop, trait = 1) {
  # Takes phenotypes from the pop object and adds/updates them in database
  # database list
  # pop \code{\link{Pop-class}}
  # trait numeric, indicating which traits to set (one or more values,
  #   so, with two traits we have options: trait = 1, trait = 2, or trait = 1:2)
  # sel = getIIdPop(pop, popObject = "damOfFemalesLact4")
  sel = getIIdPop(pop, popObject = deparse(substitute(pop)))
  database$Pheno[sel, trait] = pop@pheno[, trait]
  return(database)
}

setDatabaseEbv = function(database, pop = NULL, trait = 1) {
  # Takes EBV from the pop object and adds/updates them in database
  # database list
  # pop \code{\link{Pop-class}}
  # trait numeric, indicating which traits to set (one or more values,
  #   so, with two traits we have options: trait = 1, trait = 2, or trait = 1:2)
  # sel = getIIdPop(pop, popObject = "eliteSires")
  sel = getIIdPop(pop, popObject = deparse(substitute(pop)))
  database$Ebv[sel, trait] = pop@ebv[, trait]
  return(database)
}

#' @rdname calcAccuracy
#' @title Calculate accuracy (=Pearson correlation) between two vectors
#'
#' @description Calculate accuracy (=Pearson correlation) between two vectors -
#'   the way breeder's do it
#'
#' @param x numeric
#' @param y numeric
#' @param digits numeric, round to digits for clearer presentation
#' @param use characther, see \code{\link{cor}}
#'
#' @return numeric
calcAccuracy = function(x, y, digits = 2, use = "complete.obs") {
  round(cor(x, y, use = use)[1, 1], digits = digits)
}

#' @rdname calcAccuracyEbvVsTgv
#' @title Calculate accuracy of Estimated Breeding Values versus True Genetic
#'   Values in a population
#'
#' @description Calculate accuracy of Estimated Breeding Values versus True
#'   Genetic Values in a population: cor(EBV, TGV)
#'
#' @param pop \code{\link{Pop-class}}
#' @param ... arguments passed to \code{calcAccuracy}
#'
#' @return numeric
calcAccuracyEbvVsTgv = function(pop, ...) {
  calcAccuracy(x = pop@ebv, y = gv(pop), ...)
}

#' @rdname calcAccuracyEbvVsTbv
#' @title Calculate accuracy of Estimated Breeding Values versus True Breeding
#'   Values in a population
#'
#' @description Calculate accuracy of Estimated Breeding Values versus True
#'   Breeding Values in a population: cor(EBV, TBV)
#'
#' @param pop \code{\link{Pop-class}}
#' @param basePop \code{\link{Pop-class}}, TODO
#' @param ... arguments passed to \code{calcAccuracy}
#'
#' @return numeric
calcAccuracyEbvVsTbv = function(pop, basePop = NULL, ...) {
  calcAccuracy(x = pop@ebv, y = bv(pop), ...)
}

#' @rdname calcAccuracyEbvEidlVsTgv
#' @title Calculate accuracy of Estimated Breeding Values plus Estimated
#'   Inbreeding Depression Load versus True Genetic Values in a population
#'
#' @description Calculate accuracy of Estimated Breeding Values plus Estimated
#'   Inbreeding Depression Load versus True Breeding Values in a population:
#'   cor(EBV + EIDL, TGV)
#'
#' @param pop \code{\link{Pop-class}}
#' @param ... arguments passed to \code{calcAccuracy}
#'
#' @return numeric
calcAccuracyEbvEidlVsTgv = function(pop, ...) {
  calcAccuracy(x = pop@ebv + getIDL(pop), y = gv(pop), ...)
}

#' @rdname calcAccuracyEbvEidlVsTbv
#' @title Calculate accuracy of Estimated Breeding Values plus Estimated
#'   Inbreeding Depression Load versus True Breeding Values in a population
#'
#' @description Calculate accuracy of Estimated Breeding Values plus Estimated
#'   Inbreeding Depression Load versus True Breeding Values in a population:
#'   cor(EBV + EIDL, TBV)
#'
#' @param pop \code{\link{Pop-class}}
#' @param basePop \code{\link{Pop-class}}, TODO
#' @param ... arguments passed to \code{calcAccuracy}
#'
#' @return numeric
calcAccuracyEbvEidlVsTbv = function(pop, basePop = NULL, ...) {
  calcAccuracy(x = pop@ebv + getIDL(pop), y = bv(pop), ...)
}

#' @rdname calcSummaryStat
#' @title Calculate summary statistic on a collection of populations
#'
#' @description Calculate summary statistic on a collection of populations
#'
#' @param pops list of \code{\link{Pop-class}}
#' @param FUN function, applied to each entry in \code{pops}
#' @param ... arguments passed to \code{FUN}
#'
#' @return data.frame with population name, statistic name, value
#'
#' @examples
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' pop1 = newPop(founderPop, simParam=SP)
#' pop2 = newPop(founderPop, simParam=SP)
#' calcAccuracyPvVsGv = function(pop) calcAccuracy(pheno(pop), gv(pop))
#' calcSummaryStat(pops = list(pop1 = pop1, pop2 = pop2), FUN = calcAccuracyPvVsGv)
calcSummaryStat = function(pops, FUN, statName = deparse(substitute(FUN)), ...) {
  if (!is.list(pops)) {
    pops = list(pops)
  }
  nPop = length(pops)
  ret = data.frame(population = names(pops),
                   statistic = statName,
                   value = NA)
  for (pop in 1:nPop) {
    ret[pop, "value"] = FUN(pops[[pop]], ...)
  }
  return(ret)
}

#' @rdname getPedNrmSubset
#' @title Get a subset of pedigree numerator relationship matrix efficiently
#'
#' @description Get a subset of pedigree numerator relationship matrix
#'   efficiently using the Colleau (2002) algorithm.
#'
#' @param pedNrmInv inverse pedigree numerator relationship matrix (ideally as a
#'   sparse matrix from a function like \code{\link[pedigreemm]{getAInv}})
#' @param ind character or numeric, names (if character) or position index (if
#'   numeric) of individuals in \code{pedNrmInv} for which we want the
#'   corresponding pedNrm
#' @param with character or numeric (see \code{ind}), the output will be
#'   \code{pedNrm[ind, with]}
#' @param sum character or numeric, individuals for which we want a sum
#'   (collapsed) parts of pedNrm instead of the individual rows and columns, the
#'   output will be |pedNrm[ind, ind] pedNrm[ind, s]|
#'                  |pedNrm[  s, ind] pedNrm[  s, s]|
#'
#' @return Pedigree numerator relationship matrix
#'
#' @reference Colleau, JJ. An indirect approach to the extensive calculation of
#'   relationship coefficients. Genet Sel Evol 34, 409 (2002).
#'   https://doi.org/10.1186/1297-9686-34-4-409
#'
#' @examples
#' # Example pedigree
#' #                         1   2  3   4  5  6  7  8
#' ped <- pedigree(sire = c(NA, NA, 1,  1, 4, 5, 5, 5),
#'                 dam  = c(NA, NA, 2, NA, 3, 2, 6, 6),
#'                 label= 1:8)
#'
#' # Pedigree numerator relationship matrix - covariance
#' (pedNrm <- getA(ped))
#'
#' # Pedigree numerator relationship matrix - precision (=inverse covariance)
#' (pedNrmInv <- getAInv(ped))
#'
#' # Subset of pedNrm[5:8, 5:8] - by subsetting pedNrm, which is inefficient
#' (pedNrmSubsetSlow <- pedNrm[5:8, 5:8])
#'
#' # Subset of pedNrm[5:8, 5:8] - via pedNrmInv, which is efficient
#' (pedNrmSubsetFast <- getPedNrmSubset(pedNrmInv, ind = 5:8))
#' pedNrmSubsetFast - pedNrmSubsetSlow
#' # above should be zero or near zero (~e-16)
#'
#' # Subset of pedNrm[5:8, 7:8]
#' (pedNrmSubsetWith <- getPedNrmSubset(pedNrmInv, ind = 5:8, with = 7:8))
#' pedNrmSubsetWith - pedNrm[7:8, 5:8]
#' # above should be zero or near zero (~e-16)
#'
#' # Collapsed/Sum subset of pedNrm
#' (pedNrmSubsetSum <- getPedNrmSubset(pedNrmInv, ind = 5:6, sum = 7:8))
#'
#' selInd <- as.character(5:6)
#' pedNrmSubsetFast[selInd, selInd] - pedNrmSubsetSum[selInd, selInd]
#' # above should be zero or near zero (~e-16)
#'
#' selSum <- as.character(7:8)
#' rowSums(pedNrmSubsetFast[selInd, selSum]) - pedNrmSubsetSum[selInd, "Sum"]
#' # above should be zero or near zero (~e-16)
#'
#' sum(pedNrmSubsetFast[selSum, selSum]) - pedNrmSubsetSum["Sum", "Sum"]
#' # above should be zero or near zero (~e-16)
#'
#' # If you need only diagonal (1 + inbreeding coefficient)
#' diag(pedNrmSubsetFast)
#'
#' # If you want values as a vector
#' c(as.matrix(pedNrmSubsetFast))
#'
#' @export
# getPedNrmSubset <- function(pedNrmInv, ind, with = ind, sum = NULL) {
#   indNames <- row.names(pedNrmInv)
#   if (is.character(ind)) {
#     ind <- match(x = ind, table = indNames, nomatch = 0)
#   }
#   if (is.character(with)) {
#     with <- match(x = with, table = indNames, nomatch = 0)
#     stop("WORK IN PROGRESS: Need to develop general approach for with!")
#   }
#   if (!is.null(sum)) {
#     if (is.character(sum)) {
#       sum <- match(x = sum, table = indNames, nomatch = 0)
#     }
#     test <- ind %in% sum
#     if (any(test)) {
#       warning("Some individuals in ind are present in sum\n",
#               "These individuals will be removed from ind, but kept in sum!")
#       ind <- ind[!test]
#     }
#     test <- with %in% sum
#     if (any(test)) {
#       warning("Some individuals in with are present in sum\n",
#               "These individuals will be removed from with, but kept in sum!")
#       with <- with[!test]
#     }
#   }
#   nIndAll <- nrow(pedNrmInv)
#   nInd <- length(ind)
#   # A x = y
#   # AInv A x = AInv y
#   # x = AInv y
#   # AInv y = x hence solving for y - matrix version of x and hence y
#   x <- sparseMatrix(i = ind, j = 1:nInd,
#                     dims = c(nIndAll, nInd))
#   x <- as(x, "dMatrix")
#   pedNrmSubset <- solve(pedNrmInv, x)[with, ]
#   dimnames(pedNrmSubset) <- list(as.character(with), as.character(ind))
#   if (!is.null(sum)) {
#     # As above A x = y etc. but here we are doing a vector of y
#     nSum <- length(sum)
#     x <- sparseMatrix(i = sum, j = rep(1, times = nSum),
#                       dims = c(nIndAll, 1))
#     x <- as(x, "dMatrix")
#     pedNrmSubsetSumTmp <- solve(pedNrmInv, x)
#     pedNrmSubsetSumVec <- pedNrmSubsetSumTmp[ind, , drop = FALSE]
#     dimnames(pedNrmSubsetSumVec) <- list(as.character(ind), "Sum")
# 
#     # We want x^TAx = y here, but we already have Ax = y so we just need x^Ty
#     pedNrmSubsetSumScalar <- crossprod(x, pedNrmSubsetSumTmp)
#     dimnames(pedNrmSubsetSumScalar) <- list("Sum", "Sum")
#     pedNrmSubset <- rbind(cbind(pedNrmSubset,          pedNrmSubsetSumVec),
#                           cbind(t(pedNrmSubsetSumVec), pedNrmSubsetSumScalar))
#   }
#   return(pedNrmSubset)
# }
# 
# 
# 

getPedNrmSubset <- function(pedNrmInv, ind, with = NULL, sum = NULL) {
  
  # ---- Setup ----
  
  indNames <- row.names(pedNrmInv)
  if (is.character(ind)) {
    ind <- match(x = ind, table = indNames, nomatch = 0)
  }
  if (is.null(with)) {
    with <- ind
  } else {
    stop("WORK IN PROGRESS: Need to develop general approach for with!")
    if (is.character(with)) {
      with <- match(x = with, table = indNames, nomatch = 0)
    }
  }
  nSumGroup <- 0
  if (!is.null(sum)) {
    if (!is.list(sum)) {
      sum <- list(sum)
    }
    nSumGroup <- length(sum)
    sumGroupName <- names(sum)
    if (is.null(sumGroupName)) {
      sumGroupName <- paste0("sumGroup", 1:nSumGroup)
    }
    for (sumGroup in 1:nSumGroup) {
      if (is.character(sum[[sumGroup]])) {
        sum[[sumGroup]] <- match(x = sum[[sumGroup]], table = indNames, nomatch = 0)
      }
      
      test <- ind %in% sum[[sumGroup]]
      if (any(test)) {
        warning("Some individuals in ind are present in",
                sumGroupName[sumGroup], "\n",
                "These individuals will be removed from ind, but kept in",
                sumGroupName[sumGroup], "!")
        ind <- ind[!test]
      }
      
      test <- with %in% sum[[sumGroup]]
      if (any(test)) {
        warning("Some individuals in with are present in",
                sumGroupName[sumGroup], "\n",
                "These individuals will be removed from with, but kept in",
                sumGroupName[sumGroup], "!")
        with <- with[!test]
      }
    }
  }
  
  # ---- Calculations ----
  
  nIndAll <- nrow(pedNrmInv)
  nInd <- length(ind)
  nIndPlusSumGroup <- nInd + nSumGroup
  pedNrmSubset <- sparseMatrix(x = 1, i = 1, j = 1,
                               dims = c(nIndPlusSumGroup, nIndPlusSumGroup),
                               symmetric = TRUE)
  # pedNrmSubset is dsCMatrix (symmetric sparse)
  if (is.null(sum)) {
    dimnames(pedNrmSubset) <- list(with, ind)
  } else {
    dimnames(pedNrmSubset) <- list(c(with, sumGroupName), c(ind, sumGroupName))
  }
  
  # First the standard A part for selected individuals
  # A x = y; if x is all 0s and a 1 in the k-th position then y is A[, k]
  # inv(A) A x = inv(A) y
  # inv(A) y = x; solve for y to get A[, k] - column
  # inv(A) Y = X; solve for Y to get A[, k] - matrix
  x <- sparseMatrix(x = 1, i = ind, j = 1:nInd,
                    dims = c(nIndAll, nInd))
  # x is dgCMatrix (general sparse)
  pedNrmSubset[1:nInd, 1:nInd] <- solve(pedNrmInv, x)[with, ]
  # solve() gives dgCMatrix (general sparse), but since pedNrmSubset is already
  # dsCMatrix (symmetric sparse) we retain this class;) Alternatively using
  #   pedNrmSubset[] <- or pedNrmSubset <- solve()
  # would give us dgCMatrix (general sparse), which we don't want and would have to cast
  #   pedNrmSubset <- as(pedNrmSubset, "symmetricMatrix") # dsCMatrix (symmetric sparse)
  
  # Now the sum groups
  if (!is.null(sum)) {
    # As above A x = y etc. but here we are using multiple 1s in y columns
    i <- unlist(sum)
    j <- rep(1:nSumGroup, times = sapply(X = sum, FUN = length))
    x <- sparseMatrix(x = 1, i = i, j = j,
                      dims = c(nIndAll, nSumGroup))
    # x is dgCMatrix (general sparse)
    pedNrmSubsetSumTmp <- solve(pedNrmInv, x)
    # pedNrmSubsetSumTmp is dgCMatrix (general sparse)
    
    # Creating temporary A_group matrix that will be added to the A matrix
    i <- rep(1:nInd, times = nSumGroup)
    pos <- (nInd + 1):nIndPlusSumGroup
    j <- rep(pos, each = nInd)
    tmp <- sparseMatrix(x = pedNrmSubsetSumTmp[ind, ], i = i, j = j,
                        dims = c(nIndPlusSumGroup, nIndPlusSumGroup),
                        symmetric = TRUE)
    # tmp is dsCMatrix (symmetric sparse)
    
    # We want x^TAx = y here, but we already have Ax = y so we just need x^Ty
    tmp[pos, pos] <- crossprod(x, pedNrmSubsetSumTmp)
    # crossprod() gives dgCMatrix (general sparse), but since tmp is already
    # dsCMatrix (symmetric sparse) we retain this class;)
    
    # Add the Ax = y and x^Ty to A
    pedNrmSubset <- pedNrmSubset + tmp
  }
  return(pedNrmSubset)
}