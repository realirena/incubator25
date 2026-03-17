# functions for life expectancy decomposition
# adopted from: Aburto and Beltran-Sanchez (2019)
# Source: https://github.com/jmaburto/Violence-and-Lifespan-variation

#Some lifetable functions from LTuniform
AKm02a0 <- function(m0, sex = "male"){
  sex <- rep(sex, length(m0))
  ifelse(sex == "male", 
         ifelse(m0 < .0230, {0.14929 - 1.99545 * m0},
                ifelse(m0 < 0.08307, {0.02832 + 3.26201 * m0},.29915)),
         # f
         ifelse(m0 < 0.01724, {0.14903 - 2.05527 * m0},
                ifelse(m0 < 0.06891, {0.04667 + 3.88089 * m0}, 0.31411))
  )
}

# compute life expectancy
LifeExpectancy <- compiler::cmpfun(function(mx,sex = "female"){
  i.openage <- length(mx)
  OPENAGE   <- i.openage - 1
  RADIX     <- 1
  if (mx[i.openage] < 0.5 | is.na(mx[i.openage])) mx[i.openage] = mx[i.openage - 1]*1.1 
  ax        <- mx * 0 + .5
  # ax[1]     <- AKm02a0(m0 = mx[1], sex = sex)
  ax[i.openage] <- if (mx[i.openage] == 0) 0.5 else 1/mx[i.openage]
  
  qx        <- mx / (1 + (1 - ax) * mx)
  qx[i.openage]       <- ifelse(is.na(qx[i.openage]), NA, 1)
  
  
  px 				    <- 1 - qx
  px[is.nan(px)]      <- 0
  lx 			        <- c(RADIX, RADIX * cumprod(px[1:OPENAGE]))
  dx 				    <- lx * qx
  #dx[i.openage] <-0
  Lx 				    <- lx - (1 - ax) * dx
  Lx[i.openage ]	    <- dx[i.openage ] * ax[i.openage ]
  Lx[is.na(Lx)] <- 0
  Tx 				    <- c(rev(cumsum(rev(Lx[1:OPENAGE]))),0) + Lx[i.openage]
  ex 				    <- Tx / lx
  ex[is.na(ex)] <- 0
  ex[i.openage] <- if (ex[OPENAGE] == 0) 0 else ax[i.openage]
  ex[1]
})

e0frommxc <- function(mxcvec,sex, nages){
  dim(mxcvec) <- c(nages,length(mxcvec)/nages)
  mx          <- rowSums(mxcvec)
  LifeExpectancy(mx,sex)
}

#Decomposition function
mydecomp <- function (func, rates1, rates2, N, ...) {
  # y1 <- func(rates1, ...)
  # y2 <- func(rates2, ...)
  d <- rates2 - rates1
  n <- length(rates1)
  delta <- d/N
  x <- rates1 + d * matrix(rep(0.5:(N - 0.5)/N, length(rates1)), 
                           byrow = TRUE, ncol = N)
  cc <- matrix(0, nrow = n, ncol = N)
  for (j in 1:N) {
    for (i in 1:n) {
      z <- rep(0, n)
      z[i] <- delta[i]/2
      cc[i, j] <- func((x[, j] + z), ...) - func((x[, j] - 
                                                    z), ...)
    }
  }
  return(rowSums(cc))
}

# some helper function for the decomposition
# construct age x cause matrix per year
yr_mat <- function(df_yr, ages, causes){
  out <- matrix(0, length(ages), length(causes),
                dimnames = list(as.character(ages), causes))
  m <- as.matrix(df_yr[, ..causes])
  rownames(m) <- as.character(df_yr$age)
  out[rownames(m), colnames(m)] <- m
  out
}

# decomposition for state x sex x draw pair
decomp_grp <- function(df_grp, N=50, causes){
  df_grp <- copy(df_grp)
  setorder(df_grp, year, age)
  
  states <- unique(df_grp$state)
  sexes <- unique(df_grp$sex)
  draws <- unique(df_grp$draw)
  
  cat("state:", states, " sex:", sexes, " draw:", draws, "\n")
  
  ages <- sort(unique(df_grp$age))
  years <- sort(unique(df_grp$year))
  
  yr_mats <- lapply(years, function(y){
    yr_mat(df_grp[year==y], ages = ages, causes = causes)
  })
  names(yr_mats) <- years
  
  # list of decomposition results
  res_list <- lapply(seq_len(length(years)-1), function(k){
    y0 <- years[k]
    y1 <- years[k+1]
    
    mat0 <- yr_mats[[as.character(y0)]]
    mat1 <- yr_mats[[as.character(y1)]]
    
    contrib <- mydecomp(
      func = e0frommxc,
      rates1 = c(mat0), rates2 = c(mat1),
      N = N, sex = sexes, nages = length(ages)
    )
    dim(contrib) <- dim(mat0)
    dimnames(contrib) <- dimnames(mat0)
    
    tmp <- as.data.table(as.table(contrib))
    setnames(tmp,  names(tmp)[1:3], c("age", "cause", "contribution"))
    tmp[, `:=`(
      age = as.numeric(as.character(age)),
      state = states,
      sex = sexes,
      draw = draws,
      year_from = y0,
      year_to = y1
    )]
    tmp[, .(state, sex, draw, year_from, year_to, age, cause, contribution)]
  })
  rbindlist(res_list)
}
