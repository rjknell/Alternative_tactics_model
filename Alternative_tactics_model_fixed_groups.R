

# Code for model of the effect of sexual selection on adaptation. This is based on the models detailed in Martinez-Ruiz, C. & Knell, R.J. (2017) Sexual selection can both increase and decrease extinction probability: reconciling demographic and evolutionary factors. The Journal of animal ecology, 86, 117-127 and Knell, R.J. & Martinez-Ruiz, C. (2017) Selective harvest focused on sexual signal traits can lead to extinction under directional environmental change. Proc. R. Soc. B, 284, 20171788. The main way that this model differs is that each individual is simply a row in the population data frame rather than an item in a list.


# The model is coded as a function called "adapt_sim_fixed_group". To run it, just paste the entire text file into the R console, press enter and then type adapt_sim_fixed_group(). The arguments to the function can be used to change the model parameters, so simulation_SS(K.cap = 500) will run the model for a carrying capacity of 500. All the model arguments are detailed below.

# Contact details: Rob Knell, School of Life Sciences, University of Hull, UK. email r.j.knell@hull.ac.uk

# This code is not currently published and is copyright Rob Knell 2020.

# Model parameters ----

#b_prob: probability of birth
#d_prob: probability of death
#K_cap: carrying capacity
#N0: starting number of individuals
#time: time of simulation
#env: environment
#max_offspring: Maximum offspring per female
#env_effect: relative importance of individual resources versus environment in determining condition
#I: Age of maturity
#group_size: size of groups that males are allocated to. Accounts for the strength of sexual selection (if 'group_size=1', no sexual selection)
#alpha: How strongly the display trait scales with condition in majors. Accounts for condition dependence in the trait
#beta: females preference strength
#display_inv: relative investment into display trait
#cost: cost of the sexual trait in males
#type: the ART is fixed if TRUE and do not compete for females in other ways, if FALSE it is simultaneous sensu Taborsky, M., Oliveira, R. F. & Brockmann, H. J. The evolution of alternative reproductive tactics: concepts and questions. Alternative reproductive tactics: an integrative approach 1, 21 (2008)
#mutation: "mutation rate" - sd of the random term added to the newly calculated e.genotype for offspring
#threshold: value for the major/minor condition threshold
#sneak_factor: determines the likelihood of minor males acquiring matings
#spermcomp: probability that an offspring of a female who mates with both a major and a minor male is fathered by the minor
#e_var: Defines different modes of environmental variation: Step, Random, Directional and None
#directional_rate: Strength of environmental change when it is directional - mean change per timestep
#random_sd1: Standard deviation of "ordinary" change when e_var is random
#random_sd2: Standard deviation of "rare event" change when e_var is random
#random_p_rare: Probability per timestep of a rare event when e_var is random
#step_size: size of the step change if the environmental change is set to step
#graphs: if TRUE then produce graphical output
#male_repro_data: if TRUE then return data on condition and reproductive success NB be careful this can generate very large data frames if K_cap and time are bigger numbers.

library(moments) #Needed to calculate skewness

adapt_sim_fixed_group <-
  function(b_prob = 1,
           d_prob = 0.08,
           K_cap = 500,
           N0 = 200,
           time = 500,
           env = 1,
           max_offspring = 3,
           env_effect = 1,
           I = 2,
           group_size = 6,
           alpha = 1,
           beta = 1,
           display_inv = 1,
           cost = 2,
  				 type  = TRUE,
           mutation = 0.05,
           threshold = 0.3,
           sneak_factor = 0.5,
           spermcomp = 0.5,
           e_var = "Random",
           directional_rate = 0.005,
           random_sd1 = 0.05,
           random_sd2 = 0.1,
           random_p_rare = 0.01,
  				 step_size = 0.25,
           graphs = TRUE,
  				 male_repro_data = FALSE) {
    
# Make sure all arguments are reasonable numbers ----

if (group_size < 1 | round(group_size) != group_size) {
  stop("group_size must be an integer equal or bigger than 1")
}

if (max_offspring < 1 | round(max_offspring) != max_offspring) {
  stop("max_offspring must be an integer equal or bigger than 1")
}

if (I < 1 | round(I) != I) {
  stop("I must be an integer equal or bigger than 1")
}

if (alpha < 0 | beta < 0 | cost < 0) {
  stop("alpha, beta and cost must be equal or bigger than 0")
}
    
# Set up initial population ----

    # All individuals alive to start
    # e_genotype = genotype determining fit to environment

population <- data.frame(alive = rep(1, length = N0),
                         sex = sample(c("F", "M"), size = N0, replace = TRUE),
                         age = round(runif(N0, min = 1, max = 10)),
                         e_genotype =rnorm(N0, mean = env, sd = 0.25),
                         resource = runif(n = N0))

#Condition for each individual ----
population$condition <- population$resource - env_effect * abs(population$e_genotype - env)


#If condition <= 0, they are dead and are removed
population <- population[population$condition > 0, ]

# Males allocated to a strategy ----
population$strategy <- ifelse(population$sex == "F", NA, "Major")

population$strategy[population$condition < threshold & population$sex  == "M"] <- "Minor"


# Display trait allocation ----
#Display trait (males) adjusted by condition dependence (alpha) and condition squared. 
  
population$display_trait <-
  ifelse(
    population$sex == "M" & population$age > I,
    display_inv * population$condition^alpha,
    0
  )

#Display trait for minors adjusted to zero
population$display_trait[population$strategy == "Minor"] <- 0


    
    
#Set up for population statistics: ----

lifespan <- NaN * (time/ 2)
mismatch_area <- NaN * (time/ 2)
strait <- NaN * (time/ 2)
e0 <- env
round <- 0
eff_pop <- 0

pop <- rep(0, times = time)
new_offspring <- rep(0, times = time)
males <- rep(0, times = time)
females <- rep(0, times = time)
env1 <- rep(0, times = time)
proportion_minors <-rep(0, times = time)
prop_offspring_minors <-rep(0, times = time)

condition_q1 <- rep(0, times = time)
condition_med <- rep(0, times = time)
condition_q3 <- rep(0, times = time)

e_genotype_q1 <- rep(0, times = time)
e_genotype_med <- rep(0, times = time)
e_genotype_q3 <- rep(0, times = time)

repro_skew <- rep(0, times = time)
  
    
    
    
    
# SIMULATION ----
    
for (i in seq(time)) {
      # Loop for each time increment


#Environmental change ----
  
  #Directional  change
  if (e_var == "Directional") {
    env<-
      ifelse(
        i < (time / 4),
        e0 + sample(c(-1, 1), 1) * runif(1, min = 0, max = 0.02),
        env<- env+ rnorm(1, mean = directional_rate, sd = 0.005)
      )
  }
  
  #Step change
  if (e_var == "Step") {
    e0 <- ifelse(i == 100, env + step_size, e0)
    env<- e0 + sample(c(-1, 1), 1) * runif(1, min = 0, max = 0.02)
  }
  
  #No change
  if (e_var == "None") {
    env<- e0 + sample(c(-1, 1), 1) * runif(1, min = 0, max = 0.02)
  }
  
  #Random change
  if (e_var == "Random") {
    env<-
      env+ sample(c(
        rnorm(1, mean = 0, sd = random_sd1),
        rnorm(1, mean = 0, sd = random_sd2)
      ),
      1,
      prob = c(1 - random_p_rare, random_p_rare))
  }
  
  if (e_var != "Random" &
      e_var != "None" & e_var != "Step" & e_var != "Directional") {
    stop("e_var must be either Directional, Step, Random or None")
  }
  
      
  # Ageing ----
  
  population$age <- population$age + 1 #Increment age by 1
  
  
  
  # Mature males
  is_male <- which(population$sex == "M" & population$age > I)
  
  # Mature females
  is_female <- which(population$sex == "F" & population$age > I)
  
  # Not males: females and immature males
  is_not_male <- which(population$sex == "F" | population$sex == "M" & population$age <= I)
  
  # Size of population
  
  pop_size <- nrow(population)
  
  # New condition, display etc. ----
  # 
  #Condition for each individual
  population$condition <- population$resource - env_effect * abs(population$e_genotype - env)
  
  #Display trait (males) adjusted by condition dependence (alpha) and condition squared.
  
  population$display_trait <-
    ifelse(
      population$sex == "M" & population$age > I,
      display_inv * population$condition^alpha,
      0
    )
  
  #Display trait for minors adjusted to zero
  population$display_trait[population$strategy == "Minor"] <- 0
  
  #Death: ----
  
  # pop_e <- length(ind) / K #Population effect on death
  pop_e <- pop_size / K_cap #Population effect on death
  
  #Set up vector for death probabilities
  prob_death <- numeric(length = pop_size)
  

  #Death probability. display_trait is zero in minors and in females
  prob_death <-
    d_prob * (pop_e + (abs(population$display_trait) * cost) + 
                (0.0154 * (population$age) ^ 2 - 
                   (0.169 * population$age) + 0.46))
                                     
  
  #Random numbers between 0 and 1
  random_death <- runif(n = pop_size)
  
  #If random number <= prob.death they die, otherwise they live
  population$alive <- ifelse(random_death <= prob_death, 0, 1)
  
  #If condition is negative they also die
  population$alive <- ifelse(population$condition <0, 0, population$alive)
  
  #Record ages at death
  age_at_death <- population$age[which(population$alive == 0)]
  
  #Remove dead individuals
  population <- population[-which(population$alive == 0), ]
  
  #Population size after deaths
  pop_size <- nrow(population)
  
  # Mature males after death applied
  is_male <- which(population$sex == "M" & population$age > I)
  
  # Mature females
  is_female <- which(population$sex == "F" & population$age > I)
  
  # Not males: females and immature males
  is_not_male <-
  	which(population$sex == "F" |
  					population$sex == "M" & population$age <= I)
  
  # Density dependence for new population size
  pop_e <- pop_size / K_cap

  # Extinction? ----
  
  #Break the main loop if all individuals are dead (length(ind)==0).
  
  extinct <- ifelse(pop_size == 0, 1, 0)
  
  extinction_t <- ifelse(extinct == 1, i, NA)
  
  if (extinct == 1) {
    break
  }
  

  
  # Mate choice and reproduction ----
  
  #This "if" prevents errors that arise when there are no mature males or no mature females left during a time step (length(is.male)==0)
  if (length(is_male) > 0 & length(is_female) > 0) {
    


# Reproduction ----
      
    # Generate new data frame with mature females only
    reproduction <- population[is_female, ]
    
    #Max number of offspring per female depending on condition
    #NB with b_prob set to 1 this is the final number of newborns
    reproduction$offspring <-
      round(max_offspring * reproduction$condition)
    
    #Offspring can't be negative
    reproduction$offspring <-
      ifelse(reproduction$offspring < 0, 0, reproduction$offspring)
    
    #Calculate the number of newborns
    reproduction$newborns <- numeric(length = nrow(reproduction))
    
    # Use a loop to make sure reproduction is working properly. So sue me.
    for (r in 1:nrow(reproduction)) {
      if (reproduction$offspring[r] == 0) {
        reproduction$newborns[r] <- 0
      }
      
      else {
        reproduction$newborns[r] <-
          sum(replicate(reproduction$offspring[r],
                        ifelse(runif(1) <= b_prob * (1 - pop_e),
                               1, 0)))
      }
    }
    
    #Which males are sexually mature majors or minors
    is_major <-
    	which(population$strategy == "Major" & population$age > I)
    is_minor <-
    	which(population$strategy == "Minor" & population$age > I)
    is_male <-
    	which(population$sex == "M" & population$age > I)
    is_female <- 
    	which(population$sex == "F" & population$age > I)
    
    # If >0 offspring and > 0 males----
    if(sum(reproduction$newborns) > 0 & length(is_male) > 0 & length(is_female) >0) {
    
    # Allocate major males ----
    

    
    #Modifies group_size if the number of available males is fewer than the group size.
    group_size2 <-
    	ifelse(group_size > length(is_major), length(is_major), group_size)
    
    
    #Allocate males if random mating or if no majors
    if (group_size == 1 | length(is_major) == 0) {
    	the_male <-
    		sample(is_male, size = length(is_female), replace = TRUE)
    }
    
    # If there are majors and mating is not random
    
    else {
    	
    	#Set up groups of males to sample from, one per female
    	
    	#If type == TRUE then only sample from majors
    	
    	if (type == TRUE) {

    		
    		# Sample majors for groups
    		major_sample <-
    			sample(is_major, size = length(is_major), replace = FALSE)
    		
    		# Remainder when #majors divided by group size
    		if (length(is_major) %% group_size2 > 0) {
    			# If groups doesn't divide into majors without a remainder, fill the last group with resampled males
    			major_sample <-
    				c(major_sample,
    					sample(
    						is_major,
    						size = group_size2 - (length(is_major) %% group_size2),
    						replace = FALSE
    					))
    			
    		}
    		
    		# Allocate majors to groups. byrow = TRUE reduces the probability of a male appearing twice in the same group 
    		male_groups <- matrix(data = major_sample, ncol = group_size2, byrow = TRUE)
    		
    	}

    	#If type == FALSE then sample from majors & minors
    	
    	if (type == FALSE) {
    		
    		# Sample males for groups
    		male_sample <-
    			sample(is_male, size = length(is_male), replace = FALSE)
    		
    		# Remainder when #majors divided by group size
    		if (length(is_male) %% group_size2 > 0) {
    			# If groups doesn't divide into majors without a remainder, fill the last group with resampled males
    			male_sample <-
    				c(
    					male_sample,
    					sample(
    						is_male,
    						size = group_size2 - (length(is_male) %% group_size2),
    						replace = FALSE
    					)
    				)
    			
    		}
    		
    		# Allocate males to groups
    		male_groups <-
    			matrix(data = male_sample, ncol = group_size2)
    		
    	}
        
    		
        #Get a value for display_trait for each male in the male.groups matrix.
        groups_display_traits <-
          matrix(data = population$display_trait[c(male_groups)], ncol = group_size2)
        
        #Order the males in each group by their display_trait values, highest in column 1
        
        for (v in 1:nrow(male_groups)) {
        	ord <- order(groups_display_traits[v, ], decreasing = TRUE)
        	male_groups[v, ] <- male_groups[v, ord]
        }
        
        # Allocate a male from each group to females
        
        the_male <- numeric(length = length(is_female))
        
        for(w in seq_along(is_female)) {

        	# Allocate each female to a group of males
        	allocated_group <- sample(nrow(male_groups), size = 1)
        	
        	# Calculate vector of probabilities for each male in the group based on rank
        	# Strength of female choice is included here: the higher beta, the more
        	# likely male 1 is to be chosen
        	male_probs <- (1/(1:group_size2)^beta)/(sum(1/(1:group_size2)^beta))
        	
        	# Choose male
        	chosen_one <- sample(1:group_size2, size = 1, prob = male_probs)
        	
        	# Sample males
        	the_male[w] <-
        		male_groups[allocated_group, chosen_one]

        }
        
        
        # #Generate matrix with first column as the row number and second as the column numbers
        # the_male_mat <- cbind(1:length(is_female), the_male_index)
        # 
        # #Use this matrix to extract the identity (row numbers from the population data frame) for each chosen male
        # the_male <- male_groups[the_male_mat]
        # 
      
    } #close else

      
# Allocate minor males ----

# Once matings with major males are established there is a probability that each female will also mate with a minor. This is determined by the sneak.factor parameter and by the proportion of minors

# Proportion of males that are minors
    prop_minors <-
      ifelse(
        length(is_male) > 0, 
             length(is_minor) / length(is_male), 
             NA
             )    
    
if(length(is_minor) > 0 & length(is_major) > 0) {      
      
no_females <- length(is_female)

# When a random number between 0 and 1 < sneak factor * prop_minors, a randomly selected  minor gets mated.    
succesful_minors <-
  ifelse(
    runif(no_females) < sneak_factor * prop_minors,
    sample(is_minor, size = no_females, replace = TRUE),
    0
  )

} # End if is_minor > 0
      
# Fill list of minors with 0s if no majors or no minors or random mating    
if (length(is_minor) == 0 | length(is_major) == 0) {
  
  succesful_minors <- rep(0, length = length(is_female))

}

# Offspring characteristics ----

      #For clarity put relevant data in a new data frame  
      new_features <-
          data.frame(reproduction$newborns, the_male, is_female, succesful_minors)

        names(new_features) <- c("newborns", "major", "female", "minor")
        
        #This "if" prevents errors that arise when there are no new individuals born in a time step
        # if (sum(new_features$newborns > 0)) {
          new_features <- new_features[which(new_features$newborns > 0), ]
          
          # Generates a vector with row numbers duplicated by the number of newborns
          dup <- rep(1:nrow(new_features), times = new_features$newborns)				
          
          # Use dup to generate a new matrix with a row for each new individual
          new_features <- new_features[dup,]
          
          # Number of newborns
          n_newborns <- nrow(new_features)
          
          # Select one male (major or minor) to be the father for each offspring. When a female has mated with a major and a minor there is a chance of either being the father as determined by the spermcomp parameter. 
          
          # If no majors present then the minors automatically win
          if(sum(new_features$major) == 0) {
            new_features$male <- new_features$minor
          }
          
          # If no minors present then the majors automatically win
          if(sum(new_features$minor) == 0) {
            new_features$male <- new_features$major
          }
          
          if(sum(new_features$minor) != 0 & sum(new_features$major) !=0) {
          
          prob_minor <- runif(n = n_newborns)            
            
          new_features$male <- ifelse(new_features$minor == 0 | prob_minor <= spermcomp, new_features$major, new_features$minor)
          
          rm(prob_minor)
          
          }
          
          # Record the proportion of offspring fathered by minors
          
          if(length(is_minor) > 0 & length(is_major) > 0) {
          	prop_offspring_minors[i] <-
          		sum(new_features$male %in% new_features$minor) /
          		(
          			sum(new_features$male %in% new_features$minor) + sum(new_features$male %in% new_features$major)
          		)
          }
          
          if (length(is_minor) == 0) {
          	prop_offspring_minors[i] <- 0
          }
          
          if (length(is_major) == 0 & length(is_minor) > 0) {
          	prop_offspring_minors[i] <- 1
          }
          
          # Generate new attributes for the newborns
          new_alive <- rep(1, length = n_newborns)
          
          new_sex <- sample(x = c("F", "M"), size = n_newborns, replace = TRUE)
          
          new_male <- which(new_sex == "M")
          new_female <- which(new_sex == "F")
          
          #Vectors with the genotype and sexual genotype of the reproductive couples.    
          mother_e_gen <- population$e_genotype[new_features$female]
          father_e_gen <- population$e_genotype[new_features$male]
          
          #The genotypes of each new individual are calculated as the 
          #mean of the genotypes of the parents + a random number 
          #drawn from a normal distribution (mean=0, sd =0.05)
          #
          new_e_genotype <- ((mother_e_gen + father_e_gen)/2) + rnorm(n_newborns, mean = 0, sd = mutation)

          new_resource = runif(n = n_newborns)

          #Condition for each individual
          new_condition <-new_resource - abs(new_e_genotype - env) * env_effect
          
          # Allocate new offspring to majors or minors ----
          
          new_strategy <- rep(NA, length = n_newborns)  
          
          new_strategy[new_sex == "M" & 
                         new_condition >= threshold] <- "Major"
          
          new_strategy[new_sex == "M" & 
                         new_condition < threshold] <- "Minor"
        
          #Display trait (males) adjusted by condition dependence (alpha).  

          new_display_trait <- numeric(length = n_newborns)
          
          new_majors <- which(new_strategy == "Major")
          
          if (length(new_majors) > 0) {
            
            new_display_trait[new_sex == "M" & new_strategy == "Major"] <-
              display_inv * new_condition[new_sex == "M" & new_strategy == "Major"] ^ alpha
          }
          
          new_age <- rep(0, length = n_newborns)
          
          # Put them all together in a new data frame
          new_population <- data.frame(new_alive, new_sex, new_age, new_e_genotype, new_resource, new_condition, new_strategy, new_display_trait)
          
          # Change the names
          names(new_population) <- c("alive", "sex", "age", "e_genotype", "resource", "condition", "strategy", "display_trait")
          
          #If condition <= 0, they are dead and are removed
          new_population <- new_population[new_population$condition > 0, ]
          
          #Add it on to the existing population data frame
          population <- rbind(population, new_population)
          
        }#close no newborns if
        
      }#Close no-males if
      

      
      #Population statistics ----
      
      # if (i >= round(t / 4) & i > round(3 / 4)) {
      #   round <- round + 1
      #   lifespan[round] <- mean(sapply(ind, function(x)
      #     x$age))
      #   mismatch.area[round] <-
      #     mean(sqrt(sapply(ind, function(x)
      #       x$mismatch)))
      #   strait[round] <-
      #     mean(sapply(ind, function(x)
      #       x$sex.trait), na.rm = T)
        # if(length(is.mature)>0){
        #Effective population size (Ne) calculated according to Wright's formulae: 4NmNf/(Nm+Nf)
        # eff_pop[round]<-(4*length(is.female)*length(is.male))/(length(is.female)+length(is.male))
        # }
      # }
      
      # Population stats for graphs ----
      
      pop[i] <- nrow(population)
      males[i] <- length(which(population$sex == "M" & population$age >I))
      females[i] <- length(which(population$sex == "F" & population$age >I))
      new_offspring[i] <- length(which(population$age == 0))                     
      env1[i] <- env
      
      quant_gene <-
        quantile(population$e_genotype) #quantiles of e_genotype
      e_genotype_q1[i] <- quant_gene[2]
      e_genotype_med[i] <- quant_gene[3]
      e_genotype_q3[i] <- quant_gene[4]
      
      quant_condition <-
        quantile(population$condition) #quantiles of condition
      condition_q1[i] <- quant_condition[2]
      condition_med[i] <- quant_condition[3]
      condition_q3[i] <- quant_condition[4]
      
      proportion_minors[i] <- prop_minors

      
      #Calculate the number of mature males which didn't reproduce
      no_repro <- length(which(!(is_male %in% new_features$male)))
      
      #Produce vector of numbers of offspring per male
      repro <- c(as.numeric(table(new_features$male)),
                 rep(0, length = no_repro)
      )
      
      #Calculate reproductive skew
      repro_skew[i] <- skewness(repro)
      
      if (male_repro_data == TRUE) {
      
      #Generate data frame of male id, condition and reproductive output
      
      repro_males <- as.numeric(names(table(new_features$male)))
      no_repro_males <- which(!(is_male %in% new_features$male))
      
      repro_males_id <- c(repro_males, no_repro_males)
      
      repro_males_cond <- population$condition[repro_males_id]
      
      
      male_repro_temp <- data.frame(time = rep(i, times = length(is_male)), male = repro_males_id, condition = repro_males_cond, offspring = repro)
      
      if(i == 1) {
      
      	male_repro <- male_repro_temp
      	
      }
      
      else {
      	
      	male_repro <- rbind(male_repro, male_repro_temp)
      	
      }
      
      }
      
}#close main loop
    
# return(repro_skew)



#Draw graphs ----
    if (graphs == TRUE) {
      
      repro_skew <- as.numeric(repro_skew)
      
      pop_stats <-
        data.frame(
          t1 = seq(time), 
          pop, 
          males, 
          females, 
          new_offspring,
          repro_skew)
      
      env_stats <-
        data.frame(
          t1 = seq(time),
          environment = env1,
          e_genotype_q1,
          e_genotype_med,
          e_genotype_q3
        )
      
      tactics_stats <-
        data.frame(
          t1 = seq(time),
          proportion_minors,
          prop_offspring_minors,
          condition_q1,
          condition_med,
          condition_q3
        )
      #return(pop_stats)
      
      par(mfrow = c(1, 3))
      
      plot(
        pop_stats$pop ~ pop_stats$t1,
        col = "#1b9e77",
        type = "l",
        xlab = "Time",
        ylab = "Population",
        lwd = 1,
        ylim = c(0, max(pop_stats$pop) * 1.1),
        xlim = c(0, time),
        font.lab = 2,
        cex.lab = 1.3,
        cex.axis = 1.2
      )
      points(
        pop_stats$males ~ pop_stats$t1,
        col = "#d95f02",
        type = "l",
        lwd = 0.5
      )
      points(
        pop_stats$females ~ pop_stats$t1,
        col = "#7570b3",
        type = "l",
        lwd = 0.5
      )
      
      text(x = 25, y = max(pop_stats$pop)*1.07, labels = "A", cex = 2.2, font = 2)
      
      # env_stats <- subset(env_stats, environment > 0)
      plot(
        environment[environment >0] ~ env_stats$t1[environment>0],
        col = "#d95f02",
        type = "l",
        lwd = 1,
        xlab = "Time",
        ylab = "Environment (orange) or genotype (green)",
        ylim = c(0.8, max(env_stats$environment) + 0.1),
        font.lab = 2,
        cex.lab = 1.3,
        cex.axis = 1.2,
        xlim = c(0, time),
        data = env_stats
      )
      
      points(
        e_genotype_q1[environment>0] ~ t1[environment>0],
        col = "#1b9e77",
        type = "l",
        lwd = 0.2,
        data = env_stats
      )
      
      points(
        e_genotype_med[environment>0] ~ t1[environment>0],
        col = "#1b9e77",
        type = "l",
        lwd = 1,
        data = env_stats
      )
      
      points(
        e_genotype_q3[environment>0] ~ t1[environment>0],
        col = "#1b9e77",
        type = "l",
        lwd = 0.2,
        data = env_stats
      )
      
      text(x = 25, y = max(env_stats$environment)*1.015, labels = "B", cex = 2.2, font = 2)
      
      # plot(
      #   tactics_stats$proportion_minors ~ tactics_stats$t1,
      #   col = "#d95f02",
      #   type = "l",
      #   xlab = "Time",
      #   ylab = "Proportion minors (#d95f02) or threshold",
      #   lwd = 1,
      #   ylim = c(0, 1)
      # )
      
      # points(
      # 	tactics_stats$prop_offspring_minors ~ tactics_stats$t1,
      # 	col = "darkgreen",
      # 	type = "l",
      # 	lwd = 0.2
      # 	)
      # 
      # 
      # points(
      #   tactics_stats$threshold_q1 ~ tactics_stats$t1,
      #   col = "#1b9e77",
      #   type = "l",
      #   lwd = 0.2
      # )
      # points(
      #   tactics_stats$threshold_med ~ tactics_stats$t1,
      #   col = "#1b9e77",
      #   type = "l",
      #   lwd = 1
      # )
      # points(
      #   tactics_stats$threshold_q3 ~ tactics_stats$t1,
      #   col = "#1b9e77",
      #   type = "l",
      #   lwd = 0.2
      # )
      
      # legend(
      #   "topleft",
      #   col = c("#d95f02", "#1b9e77"),
      #   lty = 1,
      #   legend = c("Proportion minors", "Threshold")
      # ) 
      
      plot(
        condition_q1[condition_med > 0] ~ t1[condition_med > 0],
        col = "#1b9e77",
        type = "l",
        xlab = "Time",
        ylab = "Condition (green) or proportion minors (orange)",
        lwd = 0.2,
        ylim = c(0,1),
        font.lab = 2,
        cex.lab = 1.3,
        cex.axis = 1.2,
        xlim = c(0, time),
        data = tactics_stats
      )
      
      points(
        condition_med[condition_med > 0] ~ t1[condition_med > 0],
        col = "#1b9e77",
        type = "l",
        lwd = 1,
        data = tactics_stats
      )
      
      points(
        condition_q3[condition_med > 0] ~ t1[condition_med > 0],
        col = "#1b9e77",
        type = "l",
        lwd = 0.2,
        data = tactics_stats
      )
      
      points(
      	proportion_minors[condition_med > 0] ~ t1[condition_med > 0],
      	col = "#d95f02",
      	type = "l",
      	data = tactics_stats
      )
      
      text(x = 25, y = 0.96, labels = "C", cex = 2.2, font = 2)
      
      # plot(pop_stats$repro_skew ~ pop_stats$t1,
      #      col = "darkgreen",
      #      type = "l",
      #      ylab = "Index of reproductive skew",
      #      # ylim = c(-1, 1)
      # )
      
      # plot(pop_stats$repro_skew ~ tactics_stats$proportion_minors,
      # 		 pch = 16,
      # 		 cex = 0.4,
      # 		 xlab = "Proportion minors",
      # 		 ylab = "Index of reproductive skew")
      
   par(mfrow = c(1,1))
      
    }
    
    # Statistics from run ----

    # total_mismatch <- sum(mismatch_area, na.rm = T)
    # median_lifespan <- median(lifespan, na.rm = T)
    # median_strait <- median(strait, na.rm = T)
    # # Ne<-mean(eff_pop, na.rm=T)
    # factors <- paste(K_cap, alpha, group_size, e_var, sep = "")
    # output <-
    #   c(total_mismatch,
    #     median_lifespan,
    #     median_strait,
    #     extinct,
    #     extinction_t,
    #     factors)
    #     

# For directional change return whether extinction occurred and if it did the time to extinction
    
if (e_var == "Directional" & male_repro_data == FALSE) {

		output <- c(extinct, extinction_t)
    return(output)

}

if (male_repro_data == TRUE) {

	return(male_repro)

}

# For step change give the return time

if(e_var == "Step") {
  output <- which(e_genotype_q3 > env1)[101:time][1]-100
  
  return(output)
}
    
    
    
}#Close function
