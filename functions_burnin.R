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

setPhenoEwe = function(pop, varE, mean, herds, yearEffect, traitMask) {
  # Create a complex cow phenotype
  # pop population
  # varE numeric, environmental variance
  # herds list, holding herd, herdEffect (vector or matrix), and
  #   herd-year effect (vector or matrix) nodes
  # yearEffect numeric, year effect (scalar or vector)
  # traitMask matrix, specify which animals should have which pheno NA
  pop = setPheno(pop, varE = varE)
  herd = getHerd(pop)
  pop@pheno = pop@pheno +
    yearEffect +
    herds$herdEffect[herd, ] +
    herds$herdYearEffect[herd, ] +
    getPermEnvEffect(pop)
  if (any(traitMask)) {
    pop@pheno[traitMask] = NA
  }
  return(pop)
  mean = mean
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


}