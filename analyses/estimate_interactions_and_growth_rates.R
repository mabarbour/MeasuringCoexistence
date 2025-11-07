
## Setup ----

# source data
source('code/manage_data.R')

# start a new session, see https://stackoverflow.com/a/66127391/2554330 for library(matlib)
options(rgl.useNULL = TRUE)
library(rgl)
library(matlib) # for inverse
library(visreg)
library(marginaleffects)
library(cowplot)
library(lmerTest)

# set plot theme
theme_set(theme_cowplot())

# organize data
RYG_data <- filter(aphid_data, Treatment %in% c("RYG")) # restrict to treatment with all 3 species
R_mono_data <- filter(aphid_data, Treatment %in% c("RP"))
Y_mono_data <- filter(aphid_data, Treatment %in% c("YP"))
G_mono_data <- filter(aphid_data, Treatment %in% c("GP"))

## Visualizations ----

ggplot(RYG_data, aes(x = week, group = ID)) +
  geom_line(aes(y = lnRt), color = "red") +
  geom_line(aes(y = lnYt), color = "yellow") +
  geom_line(aes(y = lnGt), color = "green") +
  facet_wrap(~ID)

# R monocultures
ggplot(R_mono_data, aes(x = R_cage_lag1, y = R_cage/R_cage_lag1)) +
  geom_point()

ggplot(R_mono_data, aes(x = R_cage_lag1, y = log(R_cage/R_cage_lag1))) +
  geom_point()

ggplot(R_mono_data, aes(x = log(R_cage_lag1), y = log(R_cage/R_cage_lag1))) +
  geom_point() +
  geom_smooth(method = "lm")

# Y monocultures
ggplot(Y_mono_data, aes(x = Y_cage_lag1, y = Y_cage/Y_cage_lag1)) +
  geom_point()

ggplot(Y_mono_data, aes(x = Y_cage_lag1, y = log(Y_cage/Y_cage_lag1))) +
  geom_point()

ggplot(Y_mono_data, aes(x = log(Y_cage_lag1), y = log(Y_cage/Y_cage_lag1))) +
  geom_point() +
  geom_smooth(method = "lm")

# G monocultures
ggplot(G_mono_data, aes(x = G_cage_lag1, y = G_cage/G_cage_lag1)) +
  geom_point()

ggplot(G_mono_data, aes(x = G_cage_lag1, y = log(G_cage/G_cage_lag1))) +
  geom_point()

ggplot(G_mono_data, aes(x = log(G_cage_lag1), y = log(G_cage/G_cage_lag1))) +
  geom_point() +
  geom_smooth(method = "lm")

## Estimate intrinsic growth rates and intraspecific interactions ----

# Red morph
lm_R <- lm(lnRt1 ~ lnRt, R_mono_data)
summary(lm_R)
visreg(lm_R)

# Yellow morph
lm_Y <- lm(lnYt1 ~ lnYt, Y_mono_data)
summary(lm_Y)
visreg(lm_Y)

# Green morph
lm_G <- lm(lnGt1 ~ lnGt, G_mono_data)
summary(lm_G)
visreg(lm_G)

## Estimate species interactions and intrinsic growth rates ----

## Simple approach: fit linear models for each morph

# Red morph
R_data <- filter(aphid_data, Treatment %in% c("RYG","RP")) # 
lm_R <- lm(lnRt1 ~ lnRt + lnYt + lnGt, R_data) #  
summary(lm_R)
#visreg(lm_R)
plot_predictions(lm_R, by = c("week","ID")) +
  geom_line(data = R_data, aes(x = week, y = lnRt1, color = ID), alpha = 0.5) +
  facet_wrap(~ID)

# Yellow morph
Y_data <- filter(aphid_data, Treatment %in% c("RYG","YP"))
lm_Y <- lm(lnYt1 ~ lnRt + lnYt + lnGt, Y_data)
summary(lm_Y)
#visreg(lm_Y)
plot_predictions(lm_Y, by = c("week","ID")) +
  geom_line(data = Y_data, aes(x = week, y = lnYt1, color = ID), alpha = 0.5) +
  facet_wrap(~ID)

# Green morph
G_data <- filter(aphid_data, Treatment %in% c("RYG","GP")) 
lm_G <- lm(lnGt1 ~ lnRt + lnYt + lnGt, G_data) 
summary(lm_G)
#visreg(lm_G)
plot_predictions(lm_G, by = c("week","ID")) +
  geom_line(data = G_data, aes(x = week, y = lnGt1, color = ID), alpha = 0.5) +
  facet_wrap(~ID)

# intrinsic growth rate vector
r_vector <- matrix(c(coef(lm_R)["(Intercept)"], 
                     coef(lm_Y)["(Intercept)"], 
                     coef(lm_G)["(Intercept)"]),
                   ncol = 1, nrow = 3,
                   dimnames = list(c("R","Y","G"),"r"))
r_vector

# interaction matrix
A_matrix <- matrix(c(coef(lm_R)["lnRt"], coef(lm_R)["lnYt"], coef(lm_R)["lnGt"],
                    coef(lm_Y)["lnRt"], coef(lm_Y)["lnYt"], coef(lm_Y)["lnGt"],
                    coef(lm_G)["lnRt"], coef(lm_G)["lnYt"], coef(lm_G)["lnGt"]),
                  ncol = 3, nrow = 3, byrow = T,
                  dimnames = list(c("R","Y","G"),c("R","Y","G")))
A_matrix

# identity matrix
I_matrix <- diag(1, nrow = 3, ncol = 3)

inv(I_matrix - A_matrix) %*% r_vector # equilibrium abundances on log scale
exp(inv(I_matrix - A_matrix) %*% r_vector) # raw scale

# for visualizing with StructuralCoexistence shiny app
# subtract the identity matrix to put intraspecific competition on same scale as interspecific interactions, then multiply by -1 to convert to formulation for Lotka-Volterra competition.
LV_competition_matrix <- round((A_matrix - I_matrix)*-1,2) 
LV_competition_matrix
round(r_vector,2)
