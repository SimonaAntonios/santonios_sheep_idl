#  source('functions_burnin.R',echo=T)

sampleHerdYearEffect = function(n) {
  cbind(rnorm(n = n, sd = sqrt(herdYearVar)))
}

sampleYearEffect = function(n = 1) {
  c(rnorm(n = n, sd = sqrt(yearVar)))
}

# sampleHerdYearEffect = function(n) {
#   cbind(rnorm(n = n, sd = sqrt(herdYearVar[1])),
#         rnorm(n = n, sd = sqrt(herdYearVar[2])))
# }
# 
# sampleYearEffect = function(n = 1) {
#   c(rnorm(n = n, sd = sqrt(yearVar[1])),
#     rnorm(n = n, sd = sqrt(yearVar[2])))
# }

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

createTraitMask = function(pop, herdsWithTrait2 = NULL) {
  # Create a logical trait mask - which animal will not have which trait
  # pop population
  # herdsWithTrait2 character, herds with that will record trait 2 (when NULL
  #   all herds/animals will have trait 2 mask set to TRUE - will not have
  #   trait 2 phenotype)
  traitMask = matrix(data = FALSE, nrow = nInd(pop), ncol = pop@nTraits)
  traitMask[, 1] = FALSE # don't mask the pheno for trait 1
  # traitMask[, 2] = TRUE # mask the pheno for trait 2
  # if (!is.null(herdsWithTrait2)) {
  #   herd = getHerd(pop)
  #   sel = herd %in% herdsWithTrait2
  #   traitMask[sel, 2] = FALSE # but not for these herds
  # }
  # return(traitMask)
}

setPhenoEwe = function(pop, varE, mean, herds, yearEffect){ #, traitMask) {
  # Create a complex cow phenotype
  # pop population
  # varE numeric, environmental variance
  # mean numeric, population mean (such as lactation mean)
  # herds list, holding herd, herdEffect (vector or matrix), and
  #   herd-year effect (vector or matrix) nodes
  # yearEffect numeric, year effect (scalar or vector)
  # traitMask matrix, specify which animals should have which pheno NA
  pop = setPheno(pop, varE = varE)
  herd = getHerd(pop)
  pop@pheno = mean + pop@pheno +
    yearEffect +
    herds$herdEffect[herd, ] +
    herds$herdYearEffect[herd, ] +
    getPermEnvEffect(pop)
  # if (any(traitMask)) {
  #   pop@pheno[traitMask] = NA
  # }
  return(pop)
}

setDatabasePheno = function(database, pop = NULL, trait = 1) {
  # Takes phenotypes and adds/updates them in database
  # database list
  # trait numeric, indicating which traits to set (one or more values,
  #   so, with two traits we have options: trait = 1, trait = 2, or trait = 1:2
  if(!is.null(pop)) {
    popName = deparse(substitute(pop))
    matchId <- match(x = paste(database$General$IId, database$General$Pop),
                     table = paste(pop@id, rep(popName, 1, length(pop@id))), nomatch = 0)
    database$Pheno[matchId != 0, trait] <- pop@pheno[matchId, trait]
  }
  return(database)
}

recordData = function(database = NULL, pop = NULL, year, lactation = NA, label = NA) {
  if (!is.null(pop)) {
    popObject = deparse(substitute(pop))
    yob = getYob(pop)
    herd = getHerd(pop)
    herdYear = paste(herd, year, sep = "_")
    tmp = list(General = data.table(Pop         = popObject,
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
    if (is.null(database)) {
      database = tmp
    } else {
      database = list(General = rbind(database$General,
                                      tmp$General),
                      Pheno = rbind(database$Pheno,
                                    tmp$Pheno),
                      Gv    = rbind(database$Gv,
                                    tmp$Gv),
                      Ebv   = rbind(database$Ebv,
                                    tmp$Ebv))
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


estimateBreedingValues = function(pedigree, database, genotypes = NULL,
                                  trait = 1, na = -999, vars, svd = FALSE,
                                  nCoreSvd = NULL, genVarPropSvd = NULL, ...) {
  
  # Estimate breeding values with other software - at the moment this is geared
  # towards Mix99, but we could have different code base for Mix99 and blupf90, say.
  # Pedigree SP$pedigree object from AlphaSimR
  # database list
  #   * General data.table with columns IId, Year, Herd, HerdYear
  #   * Pheno matrix
  #   * Gv matrix
  #   * Ebv matrix
  # genotypes character, string of genotype file name
  # trait numeric, indicating which traits to analyse (one or more values,
  #   so, with two traits we have options: trait = 1, trait = 2, or trait = 1:2
  # na value used to denote missing value, say -999
  # vars list, variance components VarA, VarPE, VarHY, and VarE - vectors or matrices
  # svd logical, sould we run SVD ssGBLUP
  
  # Prepare pedigree file
  pedigree = cbind(IId = rownames(pedigree),
                   FId = pedigree[, "father"],
                   MId = pedigree[, "mother"])
  write.table(
    x = pedigree,
    file = "pedigree",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = " "
  )
  rm(pedigree)
  
  # Prepare phenotype file
  nTrait = length(trait)
  multiTrait = nTrait > 1
  database$Pheno = database$Pheno[, trait, drop = FALSE]
  tmp = is.na(database$Pheno)
  database$Pheno[tmp] = na
  sel = rowSums(!tmp) > 0
  phenotypes = cbind(database$General[sel,],
                     database$Pheno[sel, trait])
  colNames = colnames(phenotypes)
  colTrait = paste("Pheno", trait, sep = "")
  n = length(colNames)
  n = seq(from = n - nTrait + 1,
          to = n,
          by = 1)
  colNames[n] = colTrait
  colnames(phenotypes) = colNames
  phenotypes$HerdYearId = as.numeric(factor(phenotypes$HerdYear))
  sel = c("IId", "Year", "Herd", "Lactation",  colTrait)
  fwrite(
    x = phenotypes[, ..sel],
    file = "performance",
    sep = " ",
    col.names = FALSE
  )
  rm(phenotypes)
  
  ## Count number of traits
  nTrait <- length(trait)
  
  # Create variance-covariance file
  if (nTrait == 1) {
    blupf90Var = paste(
      paste("1 1 1", vars$varHY),
      paste("2 1 1", vars$varPE),
      paste("3 1 1", vars$varA),
      paste("4 1 1", vars$varE),
      sep = "\n"
    )
    writeLines(text = blupf90Var, con = "blupf90.var", sep = "\n")
  } else if (nTrait == 2) {
    blupf90Var = paste(
      paste("1 1 1", vars$varHY[1]),
      paste("1 2 2", vars$varHY[2]),
      paste("1 2 1", (
        diag(sqrt(vars$varHY)) %*% vars$corHY %*% diag(sqrt(vars$varHY))
      )[1, 2]),
      paste("2 1 1", vars$varPE[1]),
      paste("2 2 2", vars$varPE[2]),
      paste("2 2 1", (
        diag(sqrt(vars$varPE)) %*% vars$corPE %*% diag(sqrt(vars$varPE))
      )[1, 2]),
      paste("3 1 1", vars$varA[1]),
      paste("3 2 2", vars$varA[2]),
      paste("3 2 1", (
        diag(sqrt(vars$varA)) %*% vars$corA %*% diag(sqrt(vars$varA))
      )[1, 2]),
      paste("4 1 1", vars$varE[1]),
      paste("4 2 2", vars$varE[2]),
      paste("4 2 1", (
        diag(sqrt(vars$varE)) %*% vars$corE %*% diag(sqrt(vars$varE))
      )[1, 2]),
      sep = "\n"
    )
    writeLines(text = blupf90Var, con = "blupf90.var", sep = "\n")
  } else if (nTrait > 2)
    (stop("You can not simulate more than two traits scenario."))
  
  ## Prepare parameter file
  prepare_par <- function() {
    sink("renum.par", type = "output")
    writeLines(
      "#renumf90 parametar file
            # herd
             COMBINE 6 3
            # year
             COMBINE 7  2
            # lactation
            COMBINE 8 4
            # herd year
            COMBINE 9 3 2
            DATAFILE
            performance
            TRAITS
            # milk yield
            5
            FIELDS_PASSED TO OUTPUT
            # official_animal_id cheptel campagne
            1 2 3
            WEIGHT(S)
            
            RESIDUAL_VARIANCE
            1500
            # cheptel
            EFFECT
            6 cross alpha
            # campagne
            EFFECT
            7 cross alpha
            # lactation
            EFFECT
            8 cross alpha
            # herd year
            EFFECT
            9 cross alpha
            #animal
            EFFECT
            1 cross alpha
            RANDOM
            animal
            OPTIONAL
            pe
            FILE
            pedigree
            FILE_POS
            1 2 3 0 0
            PED_DEPTH
            0
            INBREEDING
            pedigree
            (CO)VARIANCES
            1200
            (CO)VARIANCES_PE
            500
            OPTION origID
            OPTION saveAinv
            OPTION saveAinvOrig
                         "
    )
    sink()
  }
  
  prepare_par()
  
  system(command = "echo renum.par | /usr/local/bin/renumf90 | tee renum.log")
  system(command = "echo renf90.par | /usr/local/bin/blupf90+ | tee blup.log")
  
  blup_sol = read_table(
    "solutions.orig",
    col_names = FALSE,
    skip = 1,
    col_types = cols(
      .default = col_double(),
      X1 = col_double(),
      X2 = col_double(),
      X3 = col_double(),
      X4 = col_double(),
      X5 = col_double()
    )
  )
  colnames(blup_sol) = c("Trait", "Effect", "Level", "IId", "Solution")
  
  ## Extracting EBV from the file
  ebv_ind = blup_sol %>%
    filter(Trait == 1 &
             Effect == 5) %>%    # Change effect number to your animal effect number
    select("IId", "Solution")  # Instead of levels you will use whatever the column with orig ID is (same above)
  colnames(ebv_ind) = c("IId", "EBV")
  library(dplyr)
  # attach(ebv_ind)
  # ebv_ind <- ebv_ind[order(IId),]
  ebv_ind <- arrange(ebv_ind, IId)
  return(ebv_ind)
}

# renf90.dat
# milk_yield | effect1 | effect2 | effect3  | effect4 | new_id | old_id | year | herd

setEbv = function(pop, ebv, trait = 1) {
  # Set EBV for a population
  # pop population
  # ebv data.table with columns IId and Ebv (one or more columns)
  # trait numeric, indicating which traits to set (one or more values
  #       so, with two traits options are: trait = 1, trait = 2, or trait = 1:2
  if (!is.null(pop)) {
    matchId = match(x = pop@iid, table = ebv$IId)
    if (ncol(pop@ebv) == 0) {
      pop@ebv = matrix(data = NA,
                       nrow = pop@nInd,
                       ncol = pop@nTraits)
    }
    traitDataTableCol = 1 + trait
    # pop@ebv[, trait] = as.matrix(ebv[matchId, ..traitDataTableCol]) # be careful ebv is data.table
    pop@ebv[, trait] = as.matrix(ebv[matchId, traitDataTableCol]) # be careful ebv is data.frame
  }
  return(pop)
}

setDatabaseEbv = function(database, pop = NULL, trait = 1) {
  # Takes EBV from the pop object and adds/updates them in database
  # database list
  # trait numeric, indicating which traits to set (one or more values,
  #   so, with two traits we have options: trait = 1, trait = 2, or trait = 1:2
  if(!is.null(pop)) {
    popName = deparse(substitute(pop))
    matchId <- match(x = paste(database$General$IId, database$General$Pop),
                     table = paste(pop@id, rep(popName, 1, length(pop@id))), nomatch = 0)
    database$Ebv[matchId != 0, trait] <- pop@ebv[matchId, trait]
  }
  return(database)
}