rm(list=ls())
library(R2jags)
library(runjags)
#library(mcmcplots)
library(tidyverse)
library(readxl)
library(EnvStats)
library(bayesplot)
library(truncnorm)
library(ggridges)
library(patchwork)
library(scales)
library(binom)
options(scipen=999)

#####generate observed sero estimates from citlali tables
##prop sero for all participants
sero_all <- data.frame(
  vacc_status = c("1 dose","2 doses","Overall","Unvaccinated"),
  positive = c(19, 70, 105,16),
  total = c(140, 181, 596,275)
)

sero_ciall <- binom.confint(sero_all$positive, sero_all$total, methods = "exact")
obs_all=cbind(sero_all, sero_ciall[, c("mean", "lower", "upper")])
obs_all=obs_all %>% filter(vacc_status=="Overall"|vacc_status=="Unvaccinated")
obs_all

sero_b472 <- data.frame(
  vacc_status = c("1 dose","2 doses","Overall","Unvaccinated"),
  positive = c(12, 37,62,13),
  total = c(35,57,181,89)
)

sero_cib472 <- binom.confint(sero_b472$positive, sero_b472$total, methods = "exact")
obs_b472=cbind(sero_b472, sero_cib472[, c("mean", "lower", "upper")])
obs_b472=obs_b472 %>% filter(vacc_status=="Overall"|vacc_status=="Unvaccinated")
obs_b472

sero_af72 <- data.frame(
  vacc_status = c("1 dose","2 doses","Overall","Unvaccinated"),
  positive = c(7,33,43,3),
  total = c(105,124,415,186)
)

sero_ciaf72 <- binom.confint(sero_af72$positive, sero_af72$total, methods = "exact")
obs_af72=cbind(sero_af72, sero_ciaf72[, c("mean", "lower", "upper")])
obs_af72=obs_af72 %>% filter(vacc_status=="Overall"|vacc_status=="Unvaccinated")
obs_af72

####extract the posterior cumulative incidence################
# Refit version: reads the extracted draws in outputs/summaries instead of the
# full .RData fits. Set to "refitted" or "published" and run once for each.
draws_prefix <- "refitted"

mcmc_matrixallWWw  <- readRDS(file.path("outputs", "summaries", paste0(draws_prefix, "_ww_key_draws.rds")))
mcmc_matrixallcass <- readRDS(file.path("outputs", "summaries", paste0(draws_prefix, "_case_key_draws.rds")))
mcmc_matrixallcomm <- readRDS(file.path("outputs", "summaries", paste0(draws_prefix, "_combined_key_draws.rds")))

# Extract raw posterior draws for cumulative incidence
Cuminc_ww_raw <- as.data.frame(mcmc_matrixallWWw[, grep("total_Cuminc", colnames(mcmc_matrixallWWw))])
Cuminc_cas_raw <- as.data.frame(mcmc_matrixallcass[, grep("total_Cuminc", colnames(mcmc_matrixallcass))])
Cuminc_com_raw <- as.data.frame(mcmc_matrixallcomm[, grep("total_Cuminc", colnames(mcmc_matrixallcomm))])

#######select end of outbreak
Cuminc_ww_finall <- Cuminc_ww_raw[, ncol(Cuminc_ww_raw)]
Cuminc_cas_finall <- Cuminc_cas_raw[, ncol(Cuminc_cas_raw)]
Cuminc_com_finall <- Cuminc_com_raw[, ncol(Cuminc_com_raw)]


#####convert cumulative incidence to seroprevalence
#pop_size <- 1278020 # catchment area population
pop_size <- 26100 # gbMSM population
pred_ww <- Cuminc_ww_finall / pop_size
pred_cas <- Cuminc_cas_finall / pop_size
pred_com <- Cuminc_com_finall / pop_size

posterior_violin <- bind_rows(
  data.frame(model = "WW-only", pred_sero = pred_ww),
  data.frame(model = "Case-only", pred_sero = pred_cas),
  data.frame(model = "Combined", pred_sero = pred_com)
)

vacc_status_levels <- c("1 dose", "2 doses", "Overall", "Unvaccinated")
vacc_status_levelsb <- c("Overall", "Unvaccinated")
vacc_status_levelsc <- "Unvaccinated"

posterior_replicated <- posterior_violin %>%
  tidyr::crossing(vacc_status = vacc_status_levelsc)

yy=posterior_replicated %>% filter(model=="WW-only")
quantile((yy$pred_sero)*100,prob=c(0.025,0.5,0.975))

# Expand observed data to match models (so point shows inside each violin)
obs_expandedall <- tidyr::crossing(
  model = unique(posterior_violin$model),
  obs_all %>% filter(vacc_status=="Unvaccinated")
)

# Combine observed seroprevalence into one dataset with a panel column
obs_all <- obs_expandedall %>% mutate(panel = "Overall population")
obs_af72 <- obs_af72 %>% 
  filter(vacc_status == "Unvaccinated") %>% 
  tidyr::crossing(model = unique(posterior_violin$model)) %>%
  mutate(panel = "Born after 1972")

obs_combined <- bind_rows(obs_all, obs_af72)

# Add panel column to posterior data so we can facet by it
# If you want model predictions same for both panels:
posterior_with_panel <- bind_rows(
  posterior_replicated %>% mutate(panel = "Overall population"),
  posterior_replicated %>% mutate(panel = "Born after 1972")
)

# # Single ggplot call with facet_wrap
# plot_combined <- ggplot(posterior_with_panel, aes(x = vacc_status, y = pred_sero*100, fill = model)) +
#   geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(width = 0.8)) +
#   
#   geom_point(data = obs_combined,
#              aes(x = vacc_status, y = mean * 100, group = model),
#              inherit.aes = FALSE,
#              position = position_dodge(width = 0.8),
#              size = 2.5, shape = 21, fill = "black") +
#   
#   geom_errorbar(data = obs_combined,
#                 aes(x = vacc_status,
#                     ymin = lower * 100,
#                     ymax = upper * 100,
#                     group = model),
#                 inherit.aes = FALSE,
#                 position = position_dodge(width = 0.8),
#                 width = 0.15) +
#   
#   scale_y_continuous("Seroprevalence (%)", limits = c(0, 10)) +
#   scale_fill_brewer(palette = "Set1", name = "Model") +
#   theme_minimal(base_size = 14) +
#   facet_wrap(~ panel) +
#   theme(
#     legend.position = "bottom",
#     axis.title.x = element_blank(),
#     panel.grid.major.x = element_blank()
#   )
# 
# plot_combined


off <- 0.12  # tweak 0.06–0.12

# Create a numeric x for BOTH posterior draws and observed summaries
posterior_plot <- posterior_with_panel %>%
  mutate(
    xnum = as.numeric(factor(vacc_status, levels = unique(vacc_status))),
    model_f = factor(model, levels = unique(model)),
    model_index = as.numeric(model_f) - (nlevels(model_f) + 1)/2,
    x = xnum + model_index * off
  )

obs_plot <- obs_combined %>%
  mutate(
    xnum = as.numeric(factor(vacc_status, levels = unique(posterior_with_panel$vacc_status))),
    model_f = factor(model, levels = unique(posterior_with_panel$model)),
    model_index = as.numeric(model_f) - (nlevels(model_f) + 1)/2,
    x = xnum + model_index * off
  )

# Store labels for x-axis
x_labs <- levels(factor(posterior_with_panel$vacc_status))

plot_combined <- ggplot(posterior_plot, aes(x = x, y = pred_sero * 100, fill = model)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  
  geom_errorbar(data = obs_plot,
                aes(x = x, ymin = lower * 100, ymax = upper * 100),
                inherit.aes = FALSE,
                width = 0.03) +
  
  geom_point(data = obs_plot,
             aes(x = x, y = mean * 100),
             inherit.aes = FALSE,
             size = 2.5, shape = 21, fill = "black") +
  
  scale_x_continuous(breaks = seq_along(x_labs), labels = x_labs) +
  scale_y_continuous("Seroprevalence (%)", limits = c(0, 10)) +
  scale_fill_brewer(palette = "Set1", name = "Model") +
  theme_minimal(base_size = 14) +
  # Refit version: both panels show the same all-gbMSM prediction, so the
  # post-1972 panel is labelled to make clear it is not a stratified prediction.
  facet_wrap(~ panel, labeller = labeller(
    panel = c(
      "Born after 1972" = "Model prediction (all gbMSM) vs post-1972 unvaccinated subgroup",
      "Overall population" = "Overall population"
    )
  )) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

plot_combined

ggsave(
  filename = paste0("outputs/figures/", draws_prefix, "_sero_mainfigure.pdf"),
  plot = plot_combined,
  width = 12,       # adjust width for your manuscript column/page size
  height = 7,      # adjust height accordingly
  units = "in",
  device = cairo_pdf # ensures high-quality text rendering
)

ggsave(
  filename = paste0("outputs/figures/", draws_prefix, "_sero_mainfigure.tiff"),
  plot = plot_combined,  # use your actual plot variable here
  width = 12,
  height = 7,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

##########combined model only##########

posterior_plotb=posterior_plot %>% filter(model=="Combined")
obs_plotb=obs_plot%>% filter(model=="Combined")

plot_combined <- ggplot(posterior_plotb, aes(x = x, y = pred_sero * 100, fill = model)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  
  geom_errorbar(data = obs_plotb,
                aes(x = x, ymin = lower * 100, ymax = upper * 100),
                inherit.aes = FALSE,
                width = 0.03) +
  
  geom_point(data = obs_plotb,
             aes(x = x, y = mean * 100),
             inherit.aes = FALSE,
             size = 2.5, shape = 21, fill = "black") +
  
  scale_x_continuous(breaks = seq_along(x_labs), labels = x_labs) +
  scale_y_continuous("Seroprevalence (%)", limits = c(0, 10)) +
  scale_fill_brewer(palette = "Set1", name = "Model") +
  theme_minimal(base_size = 14) +
  facet_wrap(~ panel) +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

plot_combined


plot_combined <- ggplot(posterior_plotb, 
                        aes(x = x, y = pred_sero * 100)) +
  
  # Model (violin)
  geom_violin(fill = "#D95F02", alpha = 0.5, color = "#D95F02", trim = FALSE) +
  
  # Observed CI
  geom_errorbar(data = obs_plotb,
                aes(x = x, ymin = lower * 100, ymax = upper * 100),
                inherit.aes = FALSE,
                width = 0.05,
                size = 0.8,
                color = "black") +
  
  # Observed point
  geom_point(data = obs_plotb,
             aes(x = x, y = mean * 100),
             inherit.aes = FALSE,
             size = 3,
             shape = 21,
             fill = "black") +
  
  # Arrow from observed → model (key addition)
  geom_segment(data = obs_plotb,
               aes(x = x, xend = x,
                   y = mean * 100,
                   yend = (mean * 100) + 1.5),
               arrow = arrow(length = unit(0.15, "cm")),
               inherit.aes = FALSE,
               color = "black") +
  
  # Label observed (optional but recommended)
  # geom_text(data = obs_plotb,
  #           aes(x = x, y = (mean * 100) - 1,
  #               label = "Observed"),
  #           size = 3.5,
  #           inherit.aes = FALSE) +
  geom_text(data = obs_plotb,
            aes(x = x, y = mean * 100,
                label = paste0(round(mean * 100, 1), "%")),
            nudge_x = 0.05,   # move slightly right (or -0.05 left)
            size = 3.5,
            inherit.aes = FALSE)+
  
  # Axis + labels
  scale_x_continuous(breaks = seq_along(x_labs), labels = x_labs) +
  scale_y_continuous("Seroprevalence (%)", limits = c(0, 10)) +
  
  # Clean facet labels
  facet_wrap(~ panel, labeller = labeller(
    panel = c(
      "Born after 1972" = "Born after 1972 (no prior immunity)",
      "Overall population" = "All participants"
    )
  )) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",   # remove legend (only one model)
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold")
  )

plot_combined




























############load the data with group estimates
# Refit version: same extracted draws as above, which include the group-level
# Cuminc columns.
mcmc_matrixallWW  <- mcmc_matrixallWWw
mcmc_matrixallcas <- mcmc_matrixallcass
mcmc_matrixallcom <- mcmc_matrixallcomm

# Select last column for each group (end of outbreak)
get_last_cuminc <- function(mcmc_matrix, n_groups = 3) {
  last_cols <- sapply(1:n_groups, function(g) {
    tail(grep(paste0("Cuminc\\[", g, ","), colnames(mcmc_matrix), value = TRUE), 1)
  })
  out <- mcmc_matrix[, last_cols]
  colnames(out) <- paste0("group", 1:n_groups)  # <— restore names
  out
}
Cuminc_cas_final <- get_last_cuminc(mcmc_matrixallcas, n_groups = 3)
Cuminc_ww_final  <- get_last_cuminc(mcmc_matrixallWW, n_groups = 3)
Cuminc_com_final <- get_last_cuminc(mcmc_matrixallcom, n_groups = 3)

# Example population sizes for 3 groups
group_pops <- c(15660,9135,1305)  # must sum to total ~26100

pred_ww_group <- as.data.frame(sweep(Cuminc_ww_final, 2, group_pops, "/"))
pred_cas_group <- as.data.frame(sweep(Cuminc_cas_final, 2, group_pops, "/"))
pred_com_group <- as.data.frame(sweep(Cuminc_com_final, 2, group_pops, "/"))


# Assuming pred_ww_group etc. are matrices with posterior draws per group
posterior_violin_group <- bind_rows(
  as.data.frame(pred_ww_group) %>% 
    mutate(draw = row_number()) %>%
    pivot_longer(cols = starts_with("group"), 
                 names_to = "group", 
                 values_to = "pred_sero") %>%
    mutate(model = "WW-only"),
  
  as.data.frame(pred_cas_group) %>%
    mutate(draw = row_number()) %>%
    pivot_longer(cols = starts_with("group"), 
                 names_to = "group", 
                 values_to = "pred_sero") %>%
    mutate(model = "Case-only"),
  
  as.data.frame(pred_com_group) %>%
    mutate(draw = row_number()) %>%
    pivot_longer(cols = starts_with("group"), 
                 names_to = "group", 
                 values_to = "pred_sero") %>%
    mutate(model = "Combined")
)


# Clean up group labels if necessary
posterior_violin_group <- posterior_violin_group %>%
  mutate(group = recode(group,
                        "group1" = "Group 1",
                        "group2" = "Group 2",
                        "group3" = "Group 3"))


yy=posterior_violin_group %>% filter(model=="Combined") %>% filter(group=="Group 3")
quantile((yy$pred_sero)*100,prob=c(0.025,0.5,0.975))





posterior_replicated_grp12 <- posterior_violin_group %>%
  filter(group %in% c("Group 1", "Group 2")) %>%
  tidyr::crossing(vacc_status = vacc_status_levelsc)

# Filter for group 3 separately
posterior_replicated_grp3 <- posterior_violin_group %>%
  filter(group == "Group 3") %>%
  tidyr::crossing(vacc_status = vacc_status_levelsc)

# Plot for groups 1 & 2 (shared y-axis)
plot_grp12 <- ggplot(posterior_replicated_grp12,
                     aes(x = vacc_status, y = pred_sero * 100, fill = model)) +
  geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(width = 0.8)) +
  geom_point(data = obs_combined,
             aes(x = vacc_status, y = mean * 100, group = model),
             inherit.aes = FALSE,
             position = position_dodge(width = 0.8),
             size = 2.5, shape = 21, fill = "black") +
  geom_errorbar(data = obs_combined,
                aes(x = vacc_status, ymin = lower * 100, ymax = upper * 100, group = model),
                inherit.aes = FALSE,
                position = position_dodge(width = 0.8),
                width = 0.15) +
  facet_grid(panel ~ group, scales = "fixed") +
  scale_y_continuous("Seroprevalence (%)", limits = c(0, 10)) +
  scale_fill_brewer(palette = "Set1", name = "Model") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

# Plot for group 3 (own scale)
plot_grp3 <- ggplot(posterior_replicated_grp3,
                    aes(x = vacc_status, y = pred_sero * 100, fill = model)) +
  geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(width = 0.8)) +
  geom_point(data = obs_combined,
             aes(x = vacc_status, y = mean * 100, group = model),
             inherit.aes = FALSE,
             position = position_dodge(width = 0.8),
             size = 2.5, shape = 21, fill = "black") +
  geom_errorbar(data = obs_combined,
                aes(x = vacc_status, ymin = lower * 100, ymax = upper * 100, group = model),
                inherit.aes = FALSE,
                position = position_dodge(width = 0.8),
                width = 0.15) +
  facet_grid(panel ~ group, scales = "fixed") +
  scale_y_continuous("Seroprevalence (%)", limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set1", name = "Model") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

# Combine plots with shared legend
fig_hybrid <- (plot_grp12+plot_grp3) +
  plot_layout(heights = c(1, 1))+
  #plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig_hybrid

ggsave(
  filename = paste0("outputs/figures/", draws_prefix, "_sero_suppfigs4.pdf"),
  plot = fig_hybrid,
  width = 12,       # adjust width for your manuscript column/page size
  height = 8,      # adjust height accordingly
  units = "in",
  device = cairo_pdf # ensures high-quality text rendering
)

ggsave(
  filename = paste0("outputs/figures/", draws_prefix, "_sero_suppfigs4.tiff"),
  plot = fig_hybrid,  # use your actual plot variable here
  width = 12,
  height = 8,
  dpi = 300,
  units = "in",
  device = "tiff",
  compression = "lzw"
)





































# Refit version: figures3 is not created anywhere in this script, so this save
# fails. Commented out.
# ggsave(
#   filename = "D:/Mpox25output/Figures/serob.tiff",
#   plot = figures3,  # use your actual plot variable here
#   width = 20,
#   height = 12,
#   dpi = 300,
#   units = "in",
#   device = "tiff",
#   compression = "lzw"
# )


# grp=ggplot(posterior_violin_group1, aes(x = group, y = pred_sero * 100, fill = model)) +
#   facet_wrap(~group)+
#   geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(width = 0.8)) +
#   scale_y_continuous("Predicted seroprevalence (%)", limits = c(0, NA)) +
#   scale_fill_brewer(palette = "Set1", name = "Model") +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "bottom",
#         panel.grid.major.x = element_blank())
# 
# 
# grp2=ggplot(posterior_violin_group2, aes(x = group, y = pred_sero * 100, fill = model)) +
#   facet_wrap(~group)+
#   geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(width = 0.8)) +
#   scale_y_continuous("Predicted seroprevalence (%)", limits = c(0, NA)) +
#   scale_fill_brewer(palette = "Set1", name = "Model") +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "bottom",
#         panel.grid.major.x = element_blank())
# 
# grp3=ggplot(posterior_violin_group3, aes(x = group, y = pred_sero * 100, fill = model)) +
#   facet_wrap(~group)+
#   geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(width = 0.8)) +
#   scale_y_continuous("Predicted seroprevalence (%)", limits = c(0, NA)) +
#   scale_fill_brewer(palette = "Set1", name = "Model") +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "bottom",
#         panel.grid.major.x = element_blank())

# grpall=ggplot(posterior_violin, aes(x = model, y = pred_sero * 100, fill = model)) +
#   geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(width = 0.8)) +
#   scale_y_continuous("Predicted seroprevalence (%)", limits = c(0, NA)) +
#   scale_fill_brewer(palette = "Set1", name = "Model") +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "bottom",
#         panel.grid.major.x = element_blank())

##############now add vaccine status to group 3 sero#####

















