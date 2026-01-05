#############################################X
#-------- Stopover Area Selection -----------X
#############################################X
# - fits an Integrated Point Process model  -X
#############################################X


# Load packages ----------------------

# required R packages
required_pkgs <- c("readr", 
                   "mgcv",
                   "dplyr", 
                   "ggplot2",
                   "rsample")

# load 
missing <- required_pkgs[!vapply(required_pkgs, 
                                 requireNamespace, 
                                 logical(1), 
                                 quietly = TRUE)]

# check if missing
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}



# Fit Integrated Point Process GAM ------------------


# load data
scale_df <- read_csv(("data/HUGO habitat/ipp_gam_scaled_df.csv"))


# fit
GAM_fit <- gam(test_col ~
                 is_biased + 
                 s(slope.mean, bs="ts", k=5) +
                 s(intOcc.sd, bs="ts", k=5) +
                 s(seasonal.LT, bs="ts", k=5) +
                 s(temp.LT, bs="ts", k=5) +
                 s(river.LT, bs="ts", k=5) +
                 s(cultivated.prp, bs="ts", k=5) +
                 s(lon, lat, bs="tp", k=30), 
               family = binomial, 
               method = "REML", 
               select = TRUE,
               data = scale_df
)

# review
summary(GAM_fit)

# inspect plots of variable effects
par(mfrow=c(2,3))
plot(GAM_fit, pages = 1, seWithMean = FALSE)
par(mfrow=c(1,1))


# Re-fit for cross-validation ------------------

# -----------------X
# create folds
# -----------------X

# split into five random folds
set.seed(429250)
n_folds <- 5
cv_df <- scale_df
folds <- vfold_cv(cv_df, v = 5, strata = test_col)

# create test table & vectors
test_table <- data.frame(test_col = cv_df$test_col)
test_table$preds <- numeric(nrow(test_table))

# -----------------X
# re-fit GAM
# -----------------X

for(k in seq_len(nrow(folds))) {
  
  # extract the training and testing indices 
  test_idx <- unlist(folds$splits[[k]]$in_id)  # testing set indices
  train_idx <- setdiff(seq_len(nrow(cv_df)), test_idx)
  
  # create train and test sets 
  testData <- cv_df[test_idx, ]
  trainData <- cv_df[train_idx, ]
  
  GAM_cv <- gam(test_col ~
                  is_biased + 
                  s(slope.mean, bs="ts", k=5) +
                  s(intOcc.sd, bs="ts", k=5) +
                  s(seasonal.LT, bs="ts", k=5) +
                  s(temp.LT, bs="ts", k=5) +
                  s(river.LT, bs="ts", k=5) +
                  s(cultivated.prp, bs="ts", k=5) +
                  s(lon, lat, bs="tp", k=30), 
                 family = binomial, 
                 method = "REML", 
                 select = TRUE,
                 data = cv_df
  )
  
  preds <- predict(GAM_cv, newdata = testData, type = "response")
  test_table$preds[test_idx] <- preds
  
  # clear
  rm(preds, GAM_cv, trainData, train_idx, testData, test_idx)
  gc()
  
}

# reformat table
test_table$pt_id <- scale_df$pt_id
test_table <- test_table %>% 
  mutate(test_col = as.factor(test_col)) %>% 
  dplyr::select(pt_id, test_col, preds)


# Evaluate performance ------------------

# via the Boyce Index
# approach adapted from:
# Liu et al. 2024 (https://doi.org/10.1111/ecog.07218)

# separate data
pred_all <- test_table$preds
pred_pres <- test_table$preds[test_table$test_col == 1]

# continuous boyce with spline fit
pred_all_scaled <- ((pred_all - min(pred_all))/(max(pred_all - min(pred_all))))
pred_pres_scaled <- ((pred_pres - min(pred_pres))/(max(pred_pres - min(pred_pres))))

# check range of the predicted values
range(pred_all_scaled)
range(pred_pres_scaled)
hist(pred_all_scaled, breaks = 50)

# function to compute boyce index
compute_boyce_spline <- function(fit, obs) { 
  
  # new bins
  brks <- unique(quantile(fit, probs = seq(0, 1, by = 0.05), na.rm = TRUE))
  nbin <- length(brks) - 1
  mids <- 0.5 * (brks[-1] + brks[-length(brks)])
  
  # new bin indices (keeps max inside last bin)
  i_fit <- findInterval(fit, brks, all.inside = TRUE)
  i_obs <- findInterval(obs, brks, all.inside = TRUE)
  
  # new expected and predicted proportions
  E_counts <- tabulate(i_fit, nbins = nbin); E_prop <- E_counts / sum(E_counts)
  P_counts <- tabulate(i_obs, nbins = nbin); P_prop <- P_counts / sum(P_counts)
  PE <- P_prop / E_prop
  PE[!is.finite(PE)] <- NA
  
  # new boyce calculation
  boyce_index <- suppressWarnings(cor(mids, PE, method = "spearman", use = "complete.obs"))
  
  # bootstrap CI across presences
  set.seed(606)
  B <- 400
  PE_boot <- replicate(B, {
    samp <- sample(obs, replace = TRUE)
    i_obs_b <- findInterval(samp, brks, all.inside = TRUE)
    P_b <- tabulate(i_obs_b, nbins = nbin); P_b <- P_b / sum(P_b)
    out <- P_b / E_prop
    out[!is.finite(out)] <- NA
    out
  })
  lo <- apply(PE_boot, 1, quantile, 0.025, na.rm = TRUE)
  hi <- apply(PE_boot, 1, quantile, 0.975, na.rm = TRUE)
  
  # new fit penalized spline
  df <- data.frame(mid = mids, pe = PE, lo = lo, hi = hi)
  keep <- is.finite(df$pe)
  gam_fit <- gam(pe ~ s(mid, bs = "ts", k = 20), data = df[keep, ], family = gaussian)
  df$pe_fit <- NA_real_
  df$pe_fit[keep] <- predict(gam_fit, newdata = df[keep, ])
  
  # new return
  list(data = df, gam = gam_fit, boyce_index = boyce_index)
}

# compute
result <- compute_boyce_spline(pred_all_scaled, pred_pres_scaled)
print(result)
hist(pred_all_scaled, breaks = 100)
quantile(pred_all_scaled, probs = seq(0,1,0.1), na.rm = TRUE)

# check P/E curve 
ggplot(result$data, aes(x = mid, y = pe)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(color = "darkgray") +
  geom_point(size = 1.1, color = "darkgray") +
  geom_line(aes(y = pe_fit), linewidth = 0.9) +
  scale_x_continuous(trans = "log10", breaks = c(0.0001, 0.001, 0.01, 0.1, 0.3, 0.5)) +
  annotate("text",
           x = min(result$data$mid, na.rm = TRUE),
           y = max(result$data$pe, na.rm = TRUE),
           label = paste0("Boyce Index = ", sprintf("%.2f", result$boyce_index)),
           hjust = -0.2, vjust = 0.01, size = 4.5) +
  labs(x = "Predicted Use",
       y = "Predicted-to-Expected (P/E) ratio") +
  theme_minimal(base_size = 12)

# check 99th quantile
quantile(result$data$mid, 0.99, na.rm = TRUE)

# find optimal threshold for binary high/low predictions
# via the reflection point
mid_seq <- seq(min(result$data$mid), max(result$data$mid), length.out = 1000)
pred_vals <- predict(result$gam, newdata = data.frame(mid = mid_seq))
slope <- diff(pred_vals) / diff(mid_seq)
flat_idx <- which(slope < 0.1 * max(slope))[1]  # where slope starts to flatten
reflection_point <- mid_seq[flat_idx + 1]
reflection_point # 0.0006215311
 
# check: which presences are rejected?
test_table %>% filter(test_col == 1 & preds < inflexion_point) %>% nrow() # 94
test_table %>% filter(test_col == 1 & preds < reflection_point) %>% nrow() # 4 presences occurred in lower-likelihood sites

# check: which absences are rejected?
test_table %>% filter(test_col == 0 & preds < inflexion_point) %>% nrow() # 49334
test_table %>% filter(test_col == 0 & preds < reflection_point) %>% nrow() # 10212 # 10212 absences occurred in higher-likelihood sites




