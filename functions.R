# source(file = "functions.R")

sampleEffect = function(n, var) {
  as.matrix(rnorm(n = n, sd = sqrt(var)))
}

getHerd = function(pop) {
  # Getting herd information for individuals of a population
  ret = sapply(X = pop@misc, FUN = function(x) x$herd)
  return(ret)
}

getYob = function(pop) {
  # Getting year of birth information for individuals of a population
  ret = sapply(X = pop@misc, FUN = function(x) x$yearOfBirth)
  return(ret)
}

getPermEnvEffect = function(pop) {
  # Getting permanent environment effect for individuals of a population
  ret = sapply(X = pop@misc, FUN = function(x) x$permEnvEffect)
  return(ret)
}

getIIdPop = function(pop, popObject = NULL) {
  if (is.null(popObject)) {
    popObject = deparse(substitute(pop))
  }
  ret = paste(pop@id, popObject, sep = "-")
  return(ret)
}

fillInMisc = function(pop, mothers = NULL, herds = NULL, permEnvVar = NULL,
                      year = NA) {
  # Fill in the misc slot of a population
  # pop population
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
  # pop population
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
                                    HerdYear    = herdYear),
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
  rm(pedigree)

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
  # pop Pop
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
  # pop Pop
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
  # pop Pop
  # trait numeric, indicating which traits to set (one or more values,
  #   so, with two traits we have options: trait = 1, trait = 2, or trait = 1:2)
  # sel = getIIdPop(pop, popObject = "eliteSires")
  sel = getIIdPop(pop, popObject = deparse(substitute(pop)))
  database$Ebv[sel, trait] = pop@ebv[, trait]
  return(database)
}

ebvAccuracy = function(pop, digits = 2 ) {
  round(cor(pop@ebv, pop@gv, use = "complete.obs")[1, 1], digits = digits)
}
