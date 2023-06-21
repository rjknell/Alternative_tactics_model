#### Directional environmental change ----

#### Simulation runs for the alternative tactics model ----
# 
# 
# WARNING Running any substantial number of sims will take a long time (hours or days depending). The file "output_fixed_groups_100_reps_directional.csv" has the output from this

# NB this script will be much easier to work thorugh if you 

library(doParallel)
library(foreach)
library(dplyr)
library(ggplot2)
library(moments)
library(readr)
library(MASS)


#### Custom themes ----

theme_linegraph = function() {
  theme(
    panel.background = element_rect(fill = "white", colour = NA), 
    panel.border = element_rect(fill = NA, colour = "grey20"),
    panel.grid = element_blank(),
    strip.background = element_rect(
      color = "white", 
      fill = "white"
    ),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
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

#### Load simulation function ----

source("Alternative_tactics_model.R")

# This will set the number of cores used. Pick a number suitable for your processor.
registerDoParallel(cores = 10)


#### Extinction under directional change ----

##### Simulation runs for directional change ----

###### Set up grid of parameter values ----

K_cap <- c(250, 500, 1000)
group_size <- c(1,3,5,10)
max_offspring <- c(3,6,10)
threshold <- c(0, 0.25, 0.5, 1)
type <- c("fixed", "simultaneous")
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


###### Run the sim ----

system.time(
	output1 <-
		foreach (
			i = 1:reps,
			.combine = "rbind",
			.packages = c("moments")
			) %dopar% ART_sim(
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


##### Load the data and plot the heatmaps ----

#Load data

sim_output1 <- read.csv("output_directional_100reps.csv")

# Calculate proportion surviving for each set of parameter values
mean_ext <-
  sim_output1 %>% group_by(
    K_cap,
    group_size,
    max_offspring,
    threshold,
    type,
    beta
  ) %>% summarise(extinct_prop = mean(extinct))

##### Figure 2 Simultaneous ART (type = "simultaneous") ----

# Plot heatmap threshold vs group size, faceted by beta & K_cap. Separate plots for type = "fixed" and "simultaneous". Only one value of max_offspring

#type = "simultaneous"

# Set up labels

k_labs <-
  c("Carrying\ncapacity = 250",
    "Carrying\ncapacity = 500",
    "Carrying\ncapacity = 1000")

names(k_labs) <- c(250, 500, 1000)

b_labs <- c("\u03b2 = 1", "\u03b2 = 2", "\u03b2 = 4")
names(b_labs) <- c(1, 2, 4)



# Plot heatmap

p1 <-
  ggplot(data = subset(mean_ext,
                       max_offspring == 6 & type == "simultaneous")) +
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

##### Figure 4 Fixed ART (type = "fixed") ----
#type = "fixed"

p2 <-
  ggplot(data = subset(mean_ext, max_offspring == 6 &
                         type == "fixed")) +
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


# 
# 
# WARNING Running any substantial number of sims will take a long time (hours or days depending). The file "output_fixed_groups_100_reps_directional.csv" has the output from this

source("Alternative_tactics_model.R")


# This will set the number of cores used. Pick a number suitable for your processor.
registerDoParallel(cores = 10)

##### Simulation runs for the alternative tactics model ----

###### Set up grid of parameter values ----

K_cap <- c(500)
group_size <- c(1,3, 5, 10)
max_offspring <- c(3, 6, 10)
threshold <- c(0, 0.25, 0.5, 1)
type <- c("fixed", "simultaneous")
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



###### Run the sim ----

system.time(
  return_time <-
    foreach (
      i = 1:reps,
      .combine = c,
      .packages = c("moments")
    ) %dopar% ART_sim(
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

##### Load data and plot graph ----

sim_output2 <- read.csv("output_step_return_100_reps.csv")

# Calculate proportion surviving for each set of parameter values
mean_return <-
  sim_output2 %>% group_by(
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

# Set up labels

mean_return$type <- ifelse(mean_return$type == "fixed", "Fixed", "Simultaneous")


offspring_labs <- c("Maximum\noffspring = 3", "Maximum\noffspring = 6", "Maximum\noffspring = 10")

names(offspring_labs) <- c(3,6,10)

b_labs <- c("\u03b2 = 1", "\u03b2 = 2", "\u03b2 = 4")
names(b_labs) <- c(1,2,4)

type_labs <- c("Simultaneous ART", "Fixed ART")
names(type_labs) <- c("Simultaneous", "Fixed")



# Line graph both ARTs
p3 <- ggplot(data = subset(mean_return, max_offspring == 6)) +
  aes(x = group_size, y = med_ret, colour = as.factor(threshold), shape = as.factor(threshold)) +
  scale_x_continuous(breaks = c(1,3,5,10)) +
  scale_color_brewer(palette = "Dark2") +
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

##### Simulation runs for the alternative tactics model ----
# 
# 
# WARNING Running any substantial number of sims will take a very long time (hours or days depending). The file "output_fixed_groups_100_reps_directional.csv" has the output from this

###### Set up grid of parameter values ----

K_cap <- c(500)
group_size <- c(2, 4, 6, 8, 10)
max_offspring <- c(6)
threshold <- c(0, 0.25, 0.5)
type <- c("fixed", "simultaneous")
beta <- c(1, 2, 4)
directional_rate <-
  c(
    0.002,
    0.003,
    0.004,
    0.005,
    0.006,
    0.007,
    0.008,
    0.00325,
    0.0035,
    0.00375,
    0.00425,
    0.0045,
    0.00475,
    0.00525,
    0.0055,
    0.00575
  )

rep <- 1:20

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

###### Run the sim ----

system.time(
  output3 <-
    foreach (
      i = 1:reps,
      .combine = "rbind",
      .packages = c("moments")
    ) %dopar% ART_sim(
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

sim_output3 <- read.csv("output_effect_size_100_reps.csv")

##### Calculate CROEC for each combination of parameter values using dose.p from MASS ----

# Set up data frame
crit_rates_ART <-
  expand.grid(
    type = c("fixed", "simultaneous"),
    group_size = c(2, 4, 6, 8, 10),
    threshold = c(0.25, 0.5),
    beta = c(1,2,4),
    CROEC = 0,
    effect_size = 0
  )

# Separate data frame for the CROEC rates when no ART because it makes life easier with GGplot
crit_rates_no_ART <-
  expand.grid(
    type = "None",
    group_size = c(2, 4, 6, 8, 10),
    threshold = 0,
    beta = c(1,2,4),
    CROEC= 0,
    effect_size = NA
  )

# Add empty variable for CROEC and effect size
# crit_rates_ART$CROEC <- numeric(length = nrow(crit_rates_ART))
# crit_rates_ART$effect_size <- numeric(length = nrow(crit_rates_ART))
# 
# crit_rates_no_ART$CROEC <- numeric(length = nrow(crit_rates_no_ART))
# crit_rates_no_ART$effect_size <- numeric(length = nrow(crit_rates_no_ART))

# Loop through each row of the output dataframe
for (i in 1:nrow(crit_rates_ART)) {
  # subset the main data frame on the basis of the values in the row of the output dataframe
  temp.data1 <-
    subset(
      sim_output3,
      type == crit_rates_ART$type[i] &
        group_size == crit_rates_ART$group_size[i] &
        threshold == crit_rates_ART$threshold[i] &
        beta == crit_rates_ART$beta[i]
        
    )
  
  # subset the main data frame to get data for threshold = 0 and the appropriate group size
  temp.data2 <-
    subset(
      sim_output3,
        group_size == crit_rates_ART$group_size[i] &
        threshold == 0 &
        beta == crit_rates_ART$beta[i]
    )
  
  # calculate glm for the present set of parameter values
  mod1 <-
    glm(extinct ~ directional_rate,
        data = temp.data1,
        family = "binomial")
  
  # calculate glm for threshold == 0 (so no ARTs)
  
  mod2 <-
    glm(extinct ~ directional_rate,
        data = temp.data2,
        family = "binomial")
  
  
  # calculate LD50 (CROEC)
  crit_rates_ART$CROEC[i] <- dose.p(mod1)[1]
  
  # calculate effect size
  if (crit_rates_ART$threshold[i] != 0) {
    crit_rates_ART$effect_size[i] <-  dose.p(mod1)[1] - dose.p(mod2)[1]
  }
  
  else crit_rates_ART$effect_size[i] <- NA
  
}

# Calculate CROEC rates for no ART

for (i in 1:nrow(crit_rates_no_ART)) {
  # subset the main data frame on the basis of the values in the row of the output dataframe
  temp.data1 <-
    subset(
      sim_output3,
        group_size == crit_rates_no_ART$group_size[i] &
        beta == crit_rates_no_ART$beta[i] &
        threshold == 0
    )
  
  # calculate glm for the present set of parameter values
  mod1 <-
    glm(extinct ~ directional_rate,
        data = temp.data1,
        family = "binomial")

 crit_rates_no_ART$CROEC[i] <- dose.p(mod1)[1] 
  
}  

# Change type so the legend has capitals
crit_rates_ART$type <- ifelse(crit_rates_ART$type == "fixed", "Fixed", "Simultaneous")

crit_rates <- rbind(crit_rates_ART, crit_rates_no_ART)


##### CROEC plot ----


# Labels
b_labs <- c("\u03b2 = 1", "\u03b2 = 2", "\u03b2 = 4")
names(b_labs) <- c(1,2,4)

# Palette
pal1 <- c( "black", "#1b9e77", "#d95f02")
pal2 <- c( "#1b9e77", "#d95f02", "#7570b3", "black")

p4 <- ggplot(data = crit_rates) +
  aes(x = group_size,
      y = CROEC,
      col = as.factor(threshold)) +
  geom_point() +
  geom_line(aes(linetype = type)) +
  geom_line(
    data = subset(crit_rates, type == "None"),
    aes(
      x = group_size,
      y = CROEC
    ),
    colour = "black",
    linetype = "dotted") +
  geom_point(data = subset(crit_rates, type == "None"),
             aes(
               x = group_size,
               y = CROEC
             ),
             colour = "black") +
  labs(x = "Mating group size", y = "Critical Rate of Environmental Change (arbitrary units)") +
  theme_linegraph() +
  scale_colour_manual(name = "Threshold",
                      breaks = c(0.25, 0.5),
                      values = pal2) +
  scale_linetype_manual(name = "Type of ART",
                        breaks = c("Fixed", "Simultaneous", "None"),
                        values = c("solid","longdash","dotted", "dotted")) +
  facet_grid(.~ beta, labeller = labeller(beta = b_labs)) 

p4

##### Effect size plot ----

p5 <- ggplot(data = subset(crit_rates, threshold !=0)) +
  aes(x = group_size,
      y = effect_size,
      colour = as.factor(threshold),
      linetype = type) +
  geom_point() +
  geom_line(aes(linetype = type)) +
  theme_linegraph() +
  # ylim(-0.0020, 0) +
  labs(
    x = "Mating group size",
    y = "Effect size",
    col = "Mating tactic\nthreshold",
    linetype = "ART type"
  )+
  scale_colour_manual(values = pal1[c(3,2)]) +
  geom_hline(yintercept = 0, alpha = 0.25, lty = 2) +
  facet_grid(. ~ beta, labeller = labeller(beta = b_labs))

p5

# ggsave(filename = "effect_sizes_figure.png", device = "png", units = "cm", width = 18, height = 14)


#### Increasing sneak factor sims ----

##### Simulation runs for changing the sneak factor ----

###### Set up grid of parameter values ----
K.cap <- c(500)
group.size <- c(2,4,6,8,10)
max.offspring <- c(6)
alpha <- c(1)
ART_success <- c(0.5, 0.7, 0.9)
init.threshold <- c(0, 0.25, 0.5)
type <- c("fixed", "simultaneous")
beta <- c(2)
c(
  0.002,
  0.003,
  0.004,
  0.005,
  0.006,
  0.007,
  0.008,
  0.00325,
  0.0035,
  0.00375,
  0.00425,
  0.0045,
  0.00475,
  0.00525,
  0.0055,
  0.00575
)

rep <- 1:100

params <-
  expand.grid(K.cap,
              group.size,
              max.offspring,
              alpha ,
              ART_success,
              init.threshold,
              obligate,
              beta,
              directional.rate,
              rep)



colnames(params) <-
  c("K.cap",
    "group.size",
    "max.offspring",
    "alpha" ,
    "ART_success",
    "init.threshold",
    "obligate",
    "beta",
    "directional.rate",
    "rep")

reps <- dim(params)[1]

###### Run the sim ----

system.time(
  output1 <-
    foreach (
      i = 1:reps,
      .combine = "rbind",
      .packages = c("moments")
    ) %dopar% ART_sim(
      mutation.thresh = 0,
      graphs = FALSE,
      e.var = "Directional",
      time = 500,
      K.cap = params[i, 1],
      group.size = params[i, 2],
      max.offspring = params[i, 3],
      alpha = params[i, 4],
      ART_success = params[i, 5],
      threshold = params[i, 6],
      obligate = params[i, 7],
      beta = params[i, 8],
      directional.rate = params[i, 9]
    )
)





sim_output5 <- data.frame(params, output1)

write_csv(sim_output5, file = "choose_your_filename.csv")

##### Load data ----

sim_output4 <- read.csv("output_ART_success_100reps.csv")


##### Calculate CROEC for each combination of parameter values using dose.p from MASS ----

# Set up data frame
crit_rates <-
  expand.grid(
    type = c("fixed", "simultaneous"),
    group_size = c(2, 4, 6, 8, 10),
    ART_success = c(0.5, 0.7, 0.9),
    threshold = c(0.25, 0.5),
    CROEC = 0,
    effect_size = 0
  )



# Calculate CROEC for no ART
crit_rates_no_ART <-
  expand.grid(
    type = "None",
    group_size = c(2,4,6,8,10),
    ART_success = 0,
    threshold = c(0.25, 0.5),
    CROEC = 0,
    effect_size = NA
  )

# Loop through each row of the output dataframe
for (i in 1:nrow(crit_rates)) {
  # subset the main data frame on the basis of the values in the row of the output dataframe
  temp.data1 <-
    subset(
      sim_output4,
      type == crit_rates$type[i] &
        group_size == crit_rates$group_size[i] &
        ART_success == crit_rates$ART_success[i] &
        threshold == crit_rates$threshold[i]
    )
  
  # subset the main data frame to get data for threshold = 0 and the appropriate group size
  temp.data2 <-
    subset(
      sim_output4,
        group_size == crit_rates$group_size[i] &
        threshold == 0
    )
  
  # calculate glm for the present set of parameter values
  mod1 <-
    glm(extinct ~ directional_rate,
        data = temp.data1,
        family = "binomial")
  
  # calculate glm for threshold == 0 (so no ARTs)
  
  mod2 <-
    glm(extinct ~ directional_rate,
        data = temp.data2,
        family = "binomial")
  
  
  # calculate LD50 (CROEC)
  crit_rates$CROEC[i] <- dose.p(mod1)[1]
  
  # CROEC for no ART
  crit_rates_no_ART$CROEC[which(crit_rates_no_ART$group_size == crit_rates$group_size[i])] <-
    dose.p(mod2)[1]
  
  # calculate effect size
  if (crit_rates$threshold[i] != 0) {
    crit_rates$effect_size[i] <-  dose.p(mod1)[1] - dose.p(mod2)[1]
  }
  
  else crit_rates$effect_size[i] <- NA
  
}

crit_rates$type <- ifelse(crit_rates$type == "fixed", "Fixed", "Simultaneous")

crit_rates <- rbind(crit_rates, crit_rates_no_ART)


##### Plot CROEC values ----
# 

pal1 <- c( "#1b9e77", "#d95f02", "#7570b3", "black")

thresh_labs <- c("Threshold = 0.25", "Threshold = 0.5")
names(thresh_labs) <- c(0.25, 0.5)

p5 <- ggplot(data = crit_rates) +
  aes(
    x = group_size,
    y = CROEC,
    colour = as.factor(ART_success)
  ) +
  geom_point() +
  geom_line(aes(x = group_size, y = CROEC, linetype = type)) +
  # geom_line(data = subset(crit_rates, type != "None"), aes(linetype = type)) +
  geom_line(
    data = subset(crit_rates, type == "None"),
    aes(x = group_size, y = CROEC),
    colour = "black",
    linetype = "dotted"
  ) +
  geom_point(
    data = subset(crit_rates, type == "None"),
    aes(x = group_size, y = CROEC),
    colour = "black"
  ) +
  theme_linegraph() +
  labs(x = "Mating group size", y = "Critical Rate of Environmental Change (Arbitrary Units)") +
  scale_colour_manual(name = "Minor male\nsuccess",
                      breaks = c(0.5, 0.7, 0.9),
                      values = pal1) +
  scale_linetype_manual(name = "Type of ART",
                        breaks = c("Fixed", "Simultaneous", "None"),
                        values = c("solid","longdash","dotted")) +
  # scale_linewidth_manual(name = "Type",
  #                        breaks = c("fixed", "simultaneous", "None"),
  #                        values = c(0.5,0.5,2)) +
  guides(shape = "none") +
  facet_grid(. ~ threshold, labeller = labeller(threshold = thresh_labs))

p5

##### Plot effect sizes ---- 

p6 <- ggplot(data = crit_rates[1:60, ]) +
  aes(x = group_size, y = effect_size, colour = as.factor(ART_success)) +
  geom_point() +
  labs(x = "Group size", y = "Effect size (arbitrary units)",
       col = "Minor male\nsuccess",
       linetype = "Type") +
  geom_line(aes(linetype = type)) +
  theme_linegraph() +
  scale_colour_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0, lty = 2, col = 'grey50') +
  facet_grid(.~ threshold, labeller = labeller(threshold = thresh_labs))

p6
