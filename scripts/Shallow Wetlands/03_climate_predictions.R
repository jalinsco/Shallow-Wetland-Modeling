#############################################X
#-------- Predict Shallow Wetlands ----------X
#############################################X
# -- a simple model w/ climatic variables ---X
#############################################X


# Load packages ----------------------

# required R packages
required_pkgs <- c("dplyr", 
                   "lme4",
                   "MuMIn",
                   "caret",
                   "tmap",
                   "sf",
                   "terra")


# load 
missing <- required_pkgs[!vapply(required_pkgs, 
                                 requireNamespace, 
                                 logical(1), 
                                 quietly = TRUE)]

# check if missing
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}



# Load data ----------------------

# reload 
t2 <- read_csv('data/Shallow Wetlands/gridded_climate_mod_t2.csv')
t3 <- read_csv('data/Shallow Wetlands/gridded_climate_mod_t3.csv')


# set ecoregion to factor
t2$ecoregion <- factor(t2$ecoregion)
lev <- levels(t2$ecoregion)
t3$ecoregion <- factor(t3$ecoregion, levels = lev)
 

# check for multicollinearity
num_vars <- t2 %>% dplyr::select(diff_dswe, diff_spei, diff_rangeT)
cor_mat <- cor(num_vars, use = "pairwise.complete.obs")
high_cor <- findCorrelation(cor_mat, cutoff = 0.5, verbose = TRUE, names = TRUE)


# Fit mixed model ----------------------

# fits mixed model w/ random effect for water/ecoregion
# for the first timestep (2021 to 2022)

mod <- lmer(diff_shallow ~ 
              diff_dswe + 
              diff_spei + 
              diff_rangeT + 
              (1 + diff_dswe | ecoregion), 
            data = t2,
            control = lme4::lmerControl(
              optimizer = "bobyqa",
              optCtrl = list(maxfun = 2e5)
            ))
summary(mod)

# check R squared
r.squaredGLMM(mod)


# Predict ----------------------

# check performance for the second time step (2022 to 2023)
t3$pred <- predict(mod, newdata = t3)
t3$resid <- t3$diff_shallow - t3$pred
caret::postResample(pred = t3$pred, obs = t3$diff_shallow)

ggplot(t3_unscaled, aes(x = diff_shallow, y = pred)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Observed", y = "Predicted")


# Residuals ----------------------

# rasterize residuals
t3_vect <- vect(st_as_sf(t3, coords = c("x", "y"), crs = "32614"))
t3_resid <- rasterize(t3_vect, r1, field = "resid", background = 0)

# set color palette
pcp_pal <- rev(diverging_hcl(100, palette = "Blue-Red"))
names(t3_resid) <- "Residual"

# plot residuals map 
tm_shape(t3_resid) +
  tm_raster(style = "cont",
            palette = pcp_pal, 
            midpoint = 0) +
  tm_shape(eco_excluded) +
  tm_polygons(col = "gray", border.col = "gray") +
  tm_layout(legend.outside = FALSE,
            legend.title.size = 0.8,
            legend.frame = FALSE,
            legend.position = c("left", "bottom"),
            legend.text.size = 0.6,
            frame = FALSE
  )
t