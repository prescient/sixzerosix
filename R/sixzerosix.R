

#' Builds monthly survival curves for policy data using random forest survival regression
#' for observed survival periods and estimates the remaining LTV tail on a per policy basis
#' using an exponential linear regression.
#' @param id a vector to uniquely identify policies.
#' @param policy_start_date a vector of dates.
#' @param policy_end_date a vector of dates.
#' @param dead 0 indicates a live policy and 1 indicates a dead policy.
#' @param censor character string specifying the tyoe of censoring. Defaults to "left"
#' @param num_periods the number of periods to forecast.
#' @param x a data.frame with regressors for the random forest survival regression. All columns will be used in regression.
#' @param train_prop proportion of data set to use for training.
#' @param num_trees number of trees to use in random forest. Defaults to 500.
#' @param mtry number of variables to possibly split at each node.
#' @param min_node_size minimal node size. Defaults to 3.
#' @param importance importance type. Defaults to impurity. See ranger documentation for additional options.
#' @return a list with input data, random forest model, and survival probabilities

build_curves <- function(id,
                         policy_start_date,
                         policy_end_date,
                         dead,
                         censor = "left",
                         num_periods,
                         x,
                         train_prop = 0.5,
                         num_trees = 500,
                         mtry = NULL,
                         min_node_size = NULL,
                         max_depth = NULL,
                         importance = "impurity") {
  # sanity checks
  assertthat::assert_that(policy_start_date <= policy_end_date,
                          is.Date(policy_start_date),
                          is.Date(policy_end_date),
                          length(id) == length(policy_start_date),
                          length(policy_start_date) == length(policy_end_date),
                          length(dead) == length(policy_end_date),
                          nrow(x) == length(dead),
                          dead == 1 | dead == 0,
                          is.numeric(num_periods),
                          is.character(formula),
                          is.data.frame(x),
                          is.numeric(num_trees),
                          nrow(x[complete.cases(x), ]) == nrow(x))

  # prepare output list
  out                   <- list()
  out$id                <- id
  out$policy_start_date <- policy_start_date
  out$policy_end_date   <- policy_end_date
  out$dead              <- dead
  out$censor            <- censor
  out$num_periods       <- num_periods
  out$formula           <- formula
  out$x                 <- x
  out$train_prop        <- train_prop
  out$num_trees         <- num_trees
  out$mtry              <- mtry
  out$min_node_size     <- min_node_size
  out$max_depth         <- max_depth
  out$importance        <- importance
  out$num_obs           <- nrow(x)

  # determine train and test set
  train        <- sample(id, floor(length(id) * train_prop), replace = F)
  out$train_id <- ifelse(id %in% train, "train", "test")

  # calculate the number of months a policy has been alive. Adding 1 to the difference
  # provides a one month policy life for policies that start and stop in the same month
  # as it is assumed that the brokerage would receive a pro rata commission.
  months_alive <- (lubridate::year(policy_end_date) * 12 + lubridate::month(policy_end_date)) -
    (lubridate::year(policy_start_date) * 12 + lubridate::lubridate::month(policy_start_date)) + 1

  # months alive summary stats
  out$unique_months_alive     <- unique(months_alive)[order(unique(months_alive))]
  out$min_months_alive        <- min(months_alive)
  out$max_months_alive        <- max(months_alive)
  out$num_unique_months_alive <- length(unique(months_alive))
  out$span_months_alive       <- max(months_alive) - min(months_alive)
  out$surv_ob                 <- surv_ob
  out$fill_months             <- ifelse(length(unique(months_alive)) <
                                          max(months_alive) - min(months_alive) + 1, T, F)
  out$months_filled           <- 1:max(months_alive)[!(1:max_months_alive %in% out$unique_months_alive)]

  # create survival object for regression
  surv_ob <- survival::Surv(time   = months_alive[id %in% train],
                            event  = dead[id %in% train],
                            censor = censor)

  # build the model
  fm <- paste0("surv_ob ~ ", paste0(colnames(x), collapse = " + "))
  mdl <- ranger::ranger(formula       = as.formula(fm),
                        data          = x[id %in% train, ],
                        num.trees     = num_trees,
                        mtry          = mtry,
                        min.node.size = min_node_size,
                        max.depth     = max_depth,
                        importance    = importance)
  out$mdl <- mdl
  rf_survival_probs <- predict(mdl, x)$survival

  # survival model only predicts probabilities for months in which a death occurred
  # therefore we have to splice in missing values if they exist to model cash flows.
  # Spliced in values estimate the prior months survival probability.
  if(length(out$months_filled) > 0){
    months_filled <- out$months_filled
    unique_months <- out$unique_months_alive
    fill_list <- list()
    n <- 1
    for(i in months_filled){
      if(i == 1){
        fill_list[[n]] <- rep(1, nrow(x))
        n <- n + 1
      } else {
        fill_list[[n]] <- rf_survival_probs[,length(unique_months[unique_months < i])]
        n <- n + 1
      }
    }
    fill <- do.call(cbind, fill_list)
    all_probs <- cbind(rf_survival_probs, fill)
    current_order <- c(unique_months, months_filled)
    pos <- c(1:ncol(all_probs))
    all_probs <- all_probs[,pos[order(current_order)]]
    rf_survival_probs <- all_probs
    rm(months_filled, unique_months, fill_list, n, i, fill, current_order, pos, all_probs)
  }
  out$rf_survival_probs <- rf_survival_probs

  # need to model the rest of the tail. In future will model tail with a covariate for AEP
  if(ncol(out$rf_survival_probs) < num_periods){
    tail <- list()
    months <- 1:ncol(out$rf_survival_probs)
    tail_pred <- data.frame(months = (ncol(out$rf_survival_probs)+1):num_periods)
    for(i in 1:nrow(rf_survival_probs)){
      tail_df <- data.frame(y = out$rf_survival_probs[i,],
                            x = months)
      mdl <- lm(log(y) ~ x, data = tail_df)
      tail[[i]] <- exp(predict(mdl, tail_pred))
    }
    tail_mat <- matrix(unlist(tail), nrow = nrow(out$rf_survival_probs), byrow = T)
    rm(months, tail_pred, tail_df, mdl, tail)
  }




}
