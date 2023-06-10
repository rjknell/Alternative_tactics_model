#### Directional environmental change ----

#### Simulation runs for the alternative tactics model ----
# 
# 
# WARNING Running any substantial number of sims will take a long time (hours or days depending). The file "output_fixed_groups_100_reps_directional.csv" has the output from this

library(doParallel)
library(foreach)
library(dplyr)
library(ggplot2)
library(moments)
library(readr)
library(MASS)


# This will set the number of cores used. Pick a number suitable for your processor.
registerDoParallel(cores = 10)

# Set up grid of parameter values

K_cap <- c(250, 500, 1000)
group_size <- c(1,3,5,10)
max_offspring <- c(3,6,10)
threshold <- c(0, 0.25, 0.5, 1)
type <- c(TRUE, FALSE)
beta <- c(1,2,4)

# How many reps for each combination of parameter values?
rep <- 1:2

# Set up a data frame with all parameter combinations
params <-
  expand.grid(K_cap,
              group_size,
              max_offspring,
              threshold,
  						type,
  						beta,
              rep)



colnames(params) <-
  c("K_cap",
    "group_size",
    "max_offspring",
    "threshold",
  	"type",
  	"beta",
    "rep")

reps <- dim(params)[1]

### Load the model code

source("Alternative_tactics_model_fixed_groups.R")

### Run the sim ----

system.time(
	output1 <-
		foreach (
			i = 1:reps,
			.combine = "rbind",
			.packages = c("moments")
			) %dopar% adapt_sim_fixed_group(
				graphs = FALSE,
				e_var = "Directional",
				directional_rate = 0.005,
				time = 500,
				K_cap = params[i, 1],
				group_size = params[i, 2],
				max_offspring = params[i, 3],
				threshold = params[i, 4],
				type = params[i, 5],
				beta = params[i, 6]
			)
)


colnames(output1) <- c("extinct", "time_to_extinction")


sim_output <- data.frame(params, output1)

write_csv(sim_output, file = "choose_your_filename.csv")


#### Load the data and plot the heatmaps ----

#Load data

sim_output <- read.csv("output_fixed_groups_100_reps_directional.csv")

# Calculate proportion surviving for each set of parameter values
mean_ext <-
  sim_output %>% group_by(
    K_cap,
    group_size,
    max_offspring,
    threshold,
    type,
    beta
  ) %>% summarise(extinct_prop = mean(extinct))

#### Figure 2 Simultaneous ART (type = FALSE) ----

# Plot heatmap threshold vs group size, faceted by beta & K_cap. Separate plots for type = TRUE and FALSE. Only one value of max_offspring

#type = FALSE

# Set up labels

k_labs <-
  c("Carrying\ncapacity = 250",
    "Carrying\ncapacity = 500",
    "Carrying\ncapacity = 1000")

names(k_labs) <- c(250, 500, 1000)

b_labs <- c("\u03b2 = 1", "\u03b2 = 2", "\u03b2 = 4")
names(b_labs) <- c(1, 2, 4)

# Custom theme

theme_heatmap = function() {
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(
      color = "white", 
      fill = "white"
    ),
    axis.text = element_text(
      family = "sans",
      face = "bold",
      size = 10),
    strip.text = element_text(
      family = "sans",
      face = "bold",
      size = 12
    ),
    axis.title = element_text(
      family = "sans",
      face = "bold",
      size = 12
    ),
    legend.title = element_text(
      family = "sans",
      face = "bold",
      size = 10
    )
  )
}


# Plot heatmap

p1 <-
  ggplot(data = subset(mean_ext,
                       max_offspring == 6 & type == FALSE)) +
  aes(x = as.factor(group_size), y = as.factor(threshold)) +
  geom_tile(aes(fill = extinct_prop), colour = "grey90") +
  scale_fill_gradient2(
    low = "#ffeda0",
    high = "#f03b20",
    mid = "#feb24c",
    name = "Extinction\nprobability",
    limits = c(0, 1),
    guide = "colorbar" ,
    midpoint = 0.5
  ) +
  ylab("Mating tactic threshold") +
  xlab("Mating group size") +
  facet_grid(beta ~ K_cap,
             labeller = labeller(K_cap = k_labs, beta = b_labs)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_heatmap()


p1

# ggsave(
#   filename = "fig2.png",
#   plot = p1,
#   device = "png",
#   units = "cm",
#   width = 15,
#   height = 9
# )

#### Figure 4 Fixed ART (type = TRUE) ----
#type = TRUE

p2 <-
  ggplot(data = subset(mean_ext, max_offspring == 6 &
                         type == TRUE)) +
  aes(x = as.factor(group_size), y = as.factor(threshold)) +
  geom_tile(aes(fill = extinct_prop), colour = "grey90") +
  scale_fill_gradient2(
    low = "#ffeda0",
    high = "#f03b20",
    mid = "#feb24c",
    name = "Extinction\nprobability",
    limits = c(0, 1),
    guide = "colorbar" ,
    midpoint = 0.5
  ) +
  ylab("Mating tactic threshold") +
  xlab("Mating group size") +
  facet_grid(beta ~ K_cap,
             labeller = labeller(K_cap = k_labs, beta = b_labs)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_heatmap()


p2

# ggsave(
#   filename = "fig4.png",
#   plot = p2,
#   device = "png",
#   units = "cm",
#   width = 15,
#   height = 9
# )

#### Step change ----

#### Simulation runs for the alternative tactics model ----
# 
# 
# WARNING Running any substantial number of sims will take a long time (hours or days depending). The file "output_fixed_groups_100_reps_directional.csv" has the output from this

source("Alternative_tactics_model_fixed_groups.R")


# This will set the number of cores used. Pick a number suitable for your processor.
registerDoParallel(cores = 10)

# Set up grid of parameter values

K_cap <- c(500)
group_size <- c(1,3, 5, 10)
max_offspring <- c(3, 6, 10)
threshold <- c(0, 0.25, 0.5, 1)
type <- c(TRUE, FALSE)
beta <- c(1, 2, 4)

rep <- 1:100

params <-
  expand.grid(K_cap,
              group_size,
              max_offspring,
              threshold,
              type,
              beta,
              rep)



colnames(params) <-
  c("K_cap",
    "group_size",
    "max_offspring",
    "threshold",
    "type",
    "beta",
    "rep")

reps <- dim(params)[1]

system.time(
  return_time <-
    foreach (
      i = 1:reps,
      .combine = c,
      .packages = c("moments")
    ) %dopar% adapt_sim_fixed_group(
      e_var = "Step",
      graphs = FALSE,
      time = 600,
      K_cap = params[i, 1],
      group_size = params[i, 2],
      max_offspring = params[i, 3],
      threshold = params[i, 4],
      type = params[i, 5],
      beta = params[i, 6]
    )
)





sim_output2 <- data.frame(params, return_time)

write_csv(sim_output2, file = "choose_your_filename.csv")

#### Load data and plot graph ----

mean_return <- read.csv("output_step_return_100_reps.csv")

# Calculate proportion surviving for each set of parameter values
mean_return <-
  mean_return %>% group_by(
    K_cap,
    group_size,
    max_offspring,
    threshold,
    type,
    beta
  ) %>% summarise(
    mean_ret = mean(return_time, nrm = RUE),
    med_ret = median(return_time, na.rm = TRUE),
    Q1 = summary(return_time)[2],
    Q3 = summary(return_time)[5]
  )

min_return <- min(mean_return$mean_ret, na.rm = TRUE)
max_return <- max(mean_return$mean_ret, na.rm = TRUE)

mean_return$mean_return_stand <-
  (mean_return$mean_ret - min_return) / (max_return - min_return)

rm(min_return, max_return)

mean_return <- filter(mean_return, K_cap == 500)

# Set up labels and theme

mean_return$type <- ifelse(mean_return$type == TRUE, "Fixed", "Simultaneous")


offspring_labs <- c("Maximum\noffspring = 3", "Maximum\noffspring = 6", "Maximum\noffspring = 10")

names(offspring_labs) <- c(3,6,10)

b_labs <- c("\u03b2 = 1", "\u03b2 = 2", "\u03b2 = 4")
names(b_labs) <- c(1,2,4)

type_labs <- c("Simultaneous ART", "Fixed ART")
names(type_labs) <- c("Simultaneous", "Fixed")

# Custom theme

theme_linegraph = function() {
  theme(
    panel.background = element_rect(fill = "white", colour = NA), 
    panel.border = element_rect(fill = NA, colour = "grey20"),
    panel.grid = element_blank(),
    strip.background = element_rect(
      color = "white", 
      fill = "white"
    ),
    axis.text = element_text(
      family = "sans",
      face = "bold",
      size = 10),
    strip.text = element_text(
      family = "sans",
      face = "bold",
      size = 12
    ),
    axis.title = element_text(
      family = "sans",
      face = "bold",
      size = 12
    ),
    legend.title = element_text(
      family = "sans",
      face = "bold",
      size = 10
    )
  )
}

# Line graph both ARTs
p3 <- ggplot(data = subset(mean_return, max_offspring == 6)) +
  aes(x = group_size, y = med_ret, colour = as.factor(threshold), shape = as.factor(threshold)) +
  scale_x_continuous(breaks = c(1,3,5,10)) +
  scale_color_brewer(palette = "Set1") +
  # scale_colour_manual(values = pal2) +
  ylim(0, 350) +
  geom_line() +
  geom_point() +
  # theme_bw() +
  theme_linegraph() +
  labs(x = "Mating group size",
       y = "Return time (arbitrary units)",
       colour = "Mating tactic\nthreshold",
       shape = "Mating tactic\nthreshold") +
  facet_grid(beta ~ type,
             labeller = labeller(type = type_labs, beta = b_labs))

p3

# ggsave(
#   "fig5.png",
#   plot = p1,
#   device = "png",
#   width = 13,
#   height = 15,
#   units = "cm"
# )

#### CROEC sims ----

#### Simulation runs for the alternative tactics model ----
# 
# 
# WARNING Running any substantial number of sims will take a very long time (hours or days depending). The file "output_fixed_groups_100_reps_directional.csv" has the output from this

# Set up grid of parameter values

K_cap <- c(500)
group_size <- c(2,4,6,8,10)
max_offspring <- c(6)
threshold <- c(0, 0.25, 0.5)
type <- c(TRUE, FALSE)
beta <- c(1, 2, 4)
directional_rate <- c(0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008)

rep <- 1:2

params <-
  expand.grid(K_cap,
              group_size,
              max_offspring,
              threshold,
              type,
              beta,
              directional_rate,
              rep)



colnames(params) <-
  c("K_cap",
    "group_size",
    "max_offspring",
    "threshold",
    "type",
    "beta",
    "directional_rate",
    "rep")

reps <- dim(params)[1]

system.time(
  output3 <-
    foreach (
      i = 1:reps,
      .combine = "rbind",
      .packages = c("moments")
    ) %dopar% adapt_sim_fixed_group(
      e_var = "Directional",
      graphs = FALSE,
      time = 600,
      K_cap = params[i, 1],
      group_size = params[i, 2],
      max_offspring = params[i, 3],
      threshold = params[i, 4],
      type = params[i, 5],
      beta = params[i, 6],
      directional_rate = params[i,7]
    )
)

colnames(output3) <- c("extinct", "time_to_extinction")

sim_output3 <- data.frame(params, output3)

write_csv(sim_output3, file = "choose_your_filename.csv")

#Load data and calculate the CROEC values

sim_output4 <- read.csv("output_effect_size_100_reps.csv")

# Calculate LD50 rate for each combination of threshold and group size using dose.p from MASS

LD50_diffs_beta1 <- data.frame(
  group_size = c(2,4,6,8,10),
  LD50_diff_fixed_05_beta1 = numeric(5),
  LD50_diff_simultaneous_05_beta1 = numeric(5),
  LD50_diff_fixed_025_beta1 = numeric(5),
  LD50_diff_simultaneous_025_beta1 = numeric(5)
)

LD50_diffs_beta2 <- data.frame(
  group_size = c(2,4,6,8,10),	
  LD50_diff_fixed_05_beta2 = numeric(5),
  LD50_diff_simultaneous_05_beta2 = numeric(5),
  LD50_diff_fixed_025_beta2 = numeric(5),
  LD50_diff_simultaneous_025_beta2 = numeric(5)
)	

LD50_diffs_beta4 <- data.frame(
  group_size = c(2,4,6,8,10),
  LD50_diff_fixed_05_beta4 = numeric(5),
  LD50_diff_simultaneous_05_beta4 = numeric(5),
  LD50_diff_fixed_025_beta4 = numeric(5),
  LD50_diff_simultaneous_025_beta4 = numeric(5)
)

# Effect sizes for beta = 1

sim2 <- subset(sim_output4, beta == 1)

for(i in unique(sim2$group_size)) {
  mod1 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0",
      data = sim2
    )
  
  mod2 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0.25" &
        type == "FALSE",
      data = sim2
    )
  
  mod3 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0.25" &
        type == "TRUE",
      data = sim2
    )
  
  mod4 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0.5" &
        type == "FALSE",
      data = sim2
    )
  
  mod5 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0.5" &
        type == "TRUE",
      data = sim2
    )
  
  
  
  
  LD50_diffs_beta1$LD50_diff_fixed_05_beta1[i / 2] <-
    dose.p(mod5)[1] - dose.p(mod1)[1]
  
  LD50_diffs_beta1$LD50_diff_simultaneous_05_beta1[i / 2] <-
    dose.p(mod4)[1] - dose.p(mod1)[1]
  
  LD50_diffs_beta1$LD50_diff_fixed_025_beta1[i / 2] <-
    dose.p(mod3)[1] - dose.p(mod1)[1]
  
  LD50_diffs_beta1$LD50_diff_simultaneous_025_beta1[i / 2] <-
    dose.p(mod2)[1] - dose.p(mod1)[1]
}

# Effect sizes for beta = 2

sim2 <- subset(sim_output4, beta == 2)

for(i in unique(sim2$group_size)) {
  mod1 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0",
      data = sim2
    )
  
  mod2 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0.25" &
        type == "FALSE",
      data = sim2
    )
  
  mod3 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0.25" &
        type == "TRUE",
      data = sim2
    )
  
  mod4 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0.5" &
        type == "FALSE",
      data = sim2
    )
  
  mod5 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0.5" &
        type == "TRUE",
      data = sim2
    )
  
  
  
  
  LD50_diffs_beta2$LD50_diff_fixed_05_beta2[i / 2] <-
    dose.p(mod5)[1] - dose.p(mod1)[1]
  
  LD50_diffs_beta2$LD50_diff_simultaneous_05_beta2[i / 2] <-
    dose.p(mod4)[1] - dose.p(mod1)[1]
  
  LD50_diffs_beta2$LD50_diff_fixed_025_beta2[i / 2] <-
    dose.p(mod3)[1] - dose.p(mod1)[1]
  
  LD50_diffs_beta2$LD50_diff_simultaneous_025_beta2[i / 2] <-
    dose.p(mod2)[1] - dose.p(mod1)[1]
}



# Effect sizes for beta = 4

sim2 <- subset(sim_output4, beta == 4)

for(i in unique(sim2$group_size)) {
  mod1 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0",
      data = sim2
    )
  
  mod2 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0.25" &
        type == "FALSE",
      data = sim2
    )
  
  mod3 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0.25" &
        type == "TRUE",
      data = sim2
    )
  
  mod4 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0.5" &
        type == "FALSE",
      data = sim2
    )
  
  mod5 <-
    glm(
      Extinct ~ directional_rate,
      family = 'binomial',
      subset = sim2$group_size == i &
        sim2$threshold == "0.5" &
        type == "TRUE",
      data = sim2
    )
  
  
  
  
  LD50_diffs_beta4$LD50_diff_fixed_05_beta4[i / 2] <-
    dose.p(mod5)[1] - dose.p(mod1)[1]
  
  LD50_diffs_beta4$LD50_diff_simultaneous_05_beta4[i / 2] <-
    dose.p(mod4)[1] - dose.p(mod1)[1]
  
  LD50_diffs_beta4$LD50_diff_fixed_025_beta4[i / 2] <-
    dose.p(mod3)[1] - dose.p(mod1)[1]
  
  LD50_diffs_beta4$LD50_diff_simultaneous_025_beta4[i / 2] <-
    dose.p(mod2)[1] - dose.p(mod1)[1]
}



#For ggplot

effect_size <- c(
  LD50_diffs_beta1$LD50_diff_fixed_025_beta1,
  LD50_diffs_beta1$LD50_diff_fixed_05_beta1,
  LD50_diffs_beta1$LD50_diff_simultaneous_025_beta1,
  LD50_diffs_beta1$LD50_diff_simultaneous_05_beta1,
  LD50_diffs_beta2$LD50_diff_fixed_025_beta2,
  LD50_diffs_beta2$LD50_diff_fixed_05_beta2,
  LD50_diffs_beta2$LD50_diff_simultaneous_025_beta2,
  LD50_diffs_beta2$LD50_diff_simultaneous_05_beta2,
  LD50_diffs_beta4$LD50_diff_fixed_025_beta4,
  LD50_diffs_beta4$LD50_diff_fixed_05_beta4,
  LD50_diffs_beta4$LD50_diff_simultaneous_025_beta4,
  LD50_diffs_beta4$LD50_diff_simultaneous_05_beta4
)

effect_sizes <-
  data.frame(
    group_size = rep(c(2, 4, 6, 8, 10), times = 12),
    type = rep(rep(c(
      "Fixed", "Simultaneous"
    ), each = 10), times = 3),
    threshold = rep(rep(c(0.25, 0.5, 0.25, 0.5), each = 5), times = 3),
    beta = rep(c("Beta = 1", "Beta = 2", "Beta = 4"), each = 20),
    effect_size
  )

effect_sizes$type <- factor(effect_sizes$type, levels = c("Fixed", "Simultaneous"))
effect_sizes$beta <- as.factor(effect_sizes$beta)
effect_sizes$threshold <- as.factor(effect_sizes$threshold)

names(effect_sizes)[2] <- "type"

rm(effect_size)

# Labels
b_labs <- c("\u03b2 = 1", "\u03b2 = 2", "\u03b2 = 4")
names(b_labs) <- c("Beta = 1","Beta = 2","Beta = 4")

p1 <- ggplot(data = effect_sizes) +
  aes(x = group_size,
      y = effect_size,
      shape = type,
      col = type) +
  geom_point() +
  geom_line(aes(linetype = threshold)) +
  theme_bw() +
  # ylim(-0.0020, 0) +
  labs(
    x = "Mating group size",
    y = "Effect size",
    shape = "ART type",
    linetype = "Mating tactic\nthreshold",
    col = "ART type"
  )+
  scale_color_brewer(palette = "Set1") +
  geom_hline(yintercept = 0, alpha = 0.25) +
  facet_grid(. ~ beta, labeller = labeller(beta = b_labs)) +
  theme_linegraph()

p1

# ggsave(filename = "effect_sizes_figure.png", device = "png", units = "cm", width = 18, height = 14)