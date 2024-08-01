# Load packages ----
library(ggplot2)

# Formula 1-v-1 ----
## Build formula ----
formula <- "ENSG00000091513~ENSG00000197971"

## Find Bandwidth ----
bw <- gwr_bwSTE(gwr_method = "basic",
                formula = formula,
                m_sfe = msfe,
                sample_id = "JBO022",
                kernel = "bisquare",
                approach = "CV")

## Run GWR ----
gwr <- gwrSTE(gwr_method = "basic",
              formula = formula,
              m_sfe = msfe,
              sample_id = "JBO022",
              bw = bw,
              kernel = "bisquare")

## SF object ----
gwr_sf <- sf::st_as_sf(gwr$SDF)

## Summary statistics for Local_R2 ----
summary_stats <- summary(gwr_sf$Local_R2)
cat("Summary Statistics for Local R-squared Values:\n")
print(summary_stats)
cat("\n")

### Plot histogram of Local_R2 values ----
ggplot(gwr_sf, aes(x = Local_R2)) +
  geom_histogram(fill = "skyblue", color = "black") +
  labs(title = "Distribution of Local R-squared Values",
       x = "Local R-squared",
       y = "Frequency") +
  theme_classic()

### Create a spatial object ----
ggplot(gwr_sf) +
  geom_sf(aes(geometry = geometry, fill = Local_R2),
          colour = "grey5") +
  labs(fill = "Local R^2") +
  scale_fill_viridis_c() +
  theme_void()

# Interpreting a Local R-squared map
# A map of the local R-squared values visualizes the goodness of fit of the local
# regression model across different spatial locations in your tissue. Here's how
# you can interpret it:
#
#     High Local R-squared values: Areas with high local R-squared values indicate
# strong spatial relationships between the predictor variables (in your case,
# gene expression levels) and the response variable (dependent variable). This
# suggests that the local regression model provides a good fit to the data in
# these areas, meaning that the predictor variables explain a large portion of
# the variability in the response variable.
#
#     Low Local R-squared values: Conversely, areas with low local R-squared
# values indicate weak spatial relationships between the predictor variables
# and the response variable. This suggests that the local regression model may
# not accurately capture the variability in the data in these areas, and the
# predictor variables may not explain much of the variability in the response
# variable.
#
#     Spatial Patterns: You can also look for spatial patterns in the local
# R-squared values. Clusters of high or low values may indicate spatial
# heterogeneity in the relationship between the predictor variables and the
# response variable. For example, you may observe regions with consistently high
# R-squared values surrounded by regions with consistently low R-squared values,
# indicating spatial variation in the strength of the relationship.
#
#     Validation: It's important to validate the results by comparing them with
# known biological knowledge or conducting further analyses. For example, areas
# with high R-squared values may correspond to regions where the biological
# pathways represented by your gene sets are highly active or relevant, while
# areas with low R-squared values may indicate regions where other factors play
# a more significant role.
#
# Overall, interpreting a map of local R-squared values allows you to assess the
# spatial variation in the relationship between your predictor variables and the
# response variable, providing insights into the spatial structure of your data
# and potential spatial patterns of biological significance.

# What is the Intercept?
# The intercept in a regression model represents the value of the dependent
# variable when all independent variables are zero. In other words, it is the
# predicted value of the dependent variable when all predictor variables have
# a value of zero.
#
# In the context of spatial regression, such as GWR, the intercept represents
# the baseline value of the dependent variable at a specific location when all
# other predictors are absent or have no influence. It's the starting point of
# the regression line or plane.

# What is the R-squared?
# The R-squared value represents the proportion of variance in the dependent
# variable that is explained by the independent variables (gene expressions) in
# the local regression model for each location.

# What are the coefficients?
# In contrast, the coefficients (e.g., ENSG00000111796 and ENSG00000198851)
# represent the estimated effect of each independent variable on the dependent
# variable. These coefficients indicate the magnitude and direction of the
# relationship between each gene expression and the dependent variable.

## Coefficient map ----
k <- gwr$lm$rank





# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Formula 1-v-Many ----
## Build formula ----
formula_im <- "ENSG00000170476~ENSG00000111796+ENSG00000198851"

## Find Bandwidth ----
bw_im <- gwr_bwSTE(gwr_method = "basic",
                   formula = formula_im,
                   m_sfe = msfe,
                   sample_id = "JBO022",
                   kernel = "bisquare",
                   approach = "CV")

## Run GWR ----
gwr_im <- gwrSTE(gwr_method = "basic",
                 formula = formula_im,
                 m_sfe = msfe,
                 sample_id = "JBO022",
                 bw = bw_im,
                 kernel = "bisquare")

## SF object ----
gwr_im_sf <- sf::st_as_sf(gwr_im$SDF)

## Summary statistics for Local_R2 ----
summary_stats_im <- summary(gwr_im_sf$Local_R2)
cat("Summary Statistics for Local R-squared Values:\n")
print(summary_stats_im)
cat("\n")

### Plot histogram of Local_R2 values ----
ggplot(gwr_im_sf, aes(x = Local_R2)) +
  geom_histogram(fill = "skyblue", color = "black") +
  labs(title = "Distribution of Local R-squared Values",
       x = "Local R-squared",
       y = "Frequency") +
  theme_classic()

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Formula 1-v-Many test ----
## Test difference if remove dependent variables ----
## Build formula ----
formula_im2 <- "ENSG00000170476~ENSG00000111796"

## Find Bandwidth ----
bw_im2 <- gwr_bwSTE(gwr_method = "basic",
                   formula = formula_im2,
                   m_sfe = msfe,
                   sample_id = "JBO022",
                   kernel = "bisquare",
                   approach = "CV")

## Run GWR ----
gwr_im2 <- gwrSTE(gwr_method = "basic",
                 formula = formula_im2,
                 m_sfe = msfe,
                 sample_id = "JBO022",
                 bw = bw_im,
                 kernel = "bisquare")

## SF object ----
gwr_im2_sf <- sf::st_as_sf(gwr_im2$SDF)

## Summary statistics for Local_R2 ----
summary_stats_im2 <- summary(gwr_im2_sf$Local_R2)
cat("Summary Statistics for Local R-squared Values:\n")
print(summary_stats_im2)
cat("\n")

### Plot histogram of Local_R2 values ----
ggplot(gwr_im2_sf, aes(x = Local_R2)) +
  geom_histogram(fill = "skyblue", color = "black") +
  labs(title = "Distribution of Local R-squared Values",
       x = "Local R-squared",
       y = "Frequency") +
  theme_classic()


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Formula Many-v-Many ----
## Build formula ----
formula_mVm <- "ENSG00000163631+ENSG00000091513+ENSG00000160868+ENSG00000130649+ENSG00000141505+ENSG00000124253+ENSG00000257017+ENSG00000130707+ENSG00000130203~ENSG00000111796+ENSG00000198851"

## Find Bandwidth ----
bw_mVm <- gwr_bwSTE(gwr_method = "basic",
                   formula = formula_mVm,
                   m_sfe = msfe,
                   sample_id = "JBO022",
                   kernel = "bisquare",
                   approach = "CV")

## Run GWR ----
gwr_mVm <- gwrSTE(gwr_method = "basic",
                 formula = formula_mVm,
                 m_sfe = msfe,
                 sample_id = "JBO022",
                 bw = bw_mVm,
                 kernel = "bisquare")

## SF object ----
gwr_mVm_sf <- sf::st_as_sf(gwr_mVm$SDF)

## Summary statistics for Local_R2 ----
summary_stats_mVm <- summary(gwr_mVm_sf$Local_R2)
cat("Summary Statistics for Local R-squared Values:\n")
print(summary_stats_mVm)
cat("\n")

### Plot histogram of Local_R2 values ----
ggplot(gwr_mVm_sf, aes(x = Local_R2)) +
  geom_histogram(fill = "skyblue", color = "black") +
  labs(title = "Distribution of Local R-squared Values",
       x = "Local R-squared",
       y = "Frequency") +
  theme_classic()

