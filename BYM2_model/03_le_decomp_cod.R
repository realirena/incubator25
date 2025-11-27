load("cod_data.RData")
load("posteriorBYM2_50draws.Rdata")
dat1 <- merge(draws_df, cod_data, 
              by = c("country", "state", "sex", "year", "age"),
              all.x = T)
library(data.table)
library(future)
library(future.apply)

dat1 <- as.data.table(dat1)

# sex code for LifeExpectancy
dat1[, sex_code := ifelse(sex == "female", "f", "m")]

# cause-specific mortality rates using "y_true"
dat1[, mx_hom   := y_true / pop]
dat1[, mx_other := other / pop]

dat1[is.na(mx_hom),   mx_hom := 0]
dat1[is.na(mx_other), mx_other := 0]

ages   <- sort(unique(dat1$age))
years  <- sort(unique(dat1$year))
causes <- c("mx_hom", "mx_other")
sexes     <- sort(unique(dat1$sex_code))
countries <- sort(unique(dat1$country))
states    <- sort(unique(dat1$state))
draws     <- sort(unique(dat1$draw))

Empty <- matrix(0, nrow = length(ages), ncol = length(causes),
                dimnames = list(ages, causes))



plan(multisession, workers = 10)   # parallel

Decomp.results <- list()

for (sx in sexes) {
  cat("Sex:", sx, "\n")
  sex_list <- list()
  sub_sex <- dat1[sex_code == sx]
  
  for (ctry in countries) {
    cat("  Country:", ctry, "\n")
    country_list <- list()
    sub_ctry <- sub_sex[country == ctry]
    
    for (st in states) {
      cat("    State:", st, "\n")
      sub_state <- sub_ctry[state == st]
      
      # subset draws for this state
      draws_state <- sort(unique(sub_state$draw))
      state_list <- list()
      
      for (dr in draws_state) {
        cat("      Draw:", dr, "\n")
        sub_draw <- sub_state[draw == dr]
        yrs <- sort(unique(sub_draw$year))
        
        # Need at least two years to decompose
        if (length(yrs) < 2) next
        
        # Build year matrices (age × cause)
        YearMatList <- lapply(yrs, function(yr) {
          tmp <- sub_draw[year == yr][order(age)]
          M <- as.matrix(tmp[, ..causes])
          rownames(M) <- tmp$age
          E <- Empty
          E[as.character(tmp$age), ] <- M
          E
        })
        names(YearMatList) <- yrs
        
        # Now decompose across time
        DecompYearList <- future_lapply(
          yrs[-length(yrs)],
          function(yr) {
            M1 <- YearMatList[[as.character(yr)]]
            M2 <- YearMatList[[as.character(yr + 1)]]
            
            contrib_vec <- mydecomp(
              func   = function(r) e0frommxc(r, sex_code = sx),
              rates1 = c(M1),
              rates2 = c(M2),
              N      = 50
            )
            
            contrib_mat <- matrix(
              contrib_vec,
              nrow = length(ages),
              ncol = length(causes),
              byrow = FALSE,
              dimnames = list(ages, causes)
            )
            
            contrib_mat
          }
        )
        names(DecompYearList) <- yrs[-length(yrs)]
        
        # Store per draw
        state_list[[as.character(dr)]] <- DecompYearList
      } # end draw loop
      
      # Store per state
      country_list[[st]] <- state_list
    } # end state loop
    
    # store per country
    sex_list[[ctry]] <- country_list
  } # end country loop
  
  # store per sex
  Decomp.results[[sx]] <- sex_list
} # end sex loop

plan(sequential)



Decomp.df <- extract_decomp_to_df(Decomp.results)
save(Decomp.df, file = "LEdecompResults.RData")
