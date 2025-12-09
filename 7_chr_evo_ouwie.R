#Project: HISTORICAL BIOGEOGRAPHY AND NICHE EVOLUTION IN EUGENIINAE (MYRTACEAE)
#Author: Paulo Henrique Gaem
#University of Michigan, Ann Arbor

# 7. PERFORMING NICHE EVOLUTION ANALYSES WITH OUWIE

# Setting working directory
setwd("~/Documents/Projetos/3_PhD/1_Project/Eugenia Biogeography")

#loading packages
library(ape)
library(phytools)
library(OUwie)
library(dplyr)
library(openxlsx)

set.seed(747) # for reproducibility

#loading phylogenetic tree
eugtre <- read.tree('data/2_eugtree.tre')

#loading environmental variables
envirn1 <- read.csv('data/25_occs_envPCA.csv')

envirn <- envirn1 %>% #summarising environmental data (PCA eigenvectors) by mean and sd
  group_by(tip.labels) %>%
  summarise(
    sections = unique(sections),
    PC1.mean = mean(PC1, na.rm = TRUE),
    PC1.sd = sd(PC1, na.rm = TRUE),
    PC2.mean = mean(PC2, na.rm = TRUE),
    PC2.sd = sd(PC2, na.rm = TRUE),
    PC3.mean = mean(PC3, na.rm = TRUE),
    PC3.sd = sd(PC3, na.rm = TRUE)
  )

missing_labels <- setdiff(eugtre$tip.label, envirn$tip.labels)
print(missing_labels) # species for which no reliable occurrence points remained after cleaning

#prunning phylogeny in order to remove spp. without occurrence points
pruned_tree <- drop.tip(eugtre, setdiff(eugtre$tip.label, envirn$tip.labels))

#creating a vector with all extra-Neotropical tips
mrca_node <- getMRCA(pruned_tree, c("Eugenia.mespiloides_Wilson.10741", "Eugenia.tropophylla_Lowrey.2134")) #getting the Jossinia MRCA
joss.tree <- extract.clade(pruned_tree, mrca_node) #extracting the Jossinia clade
palaeo.tips <- unique(c(joss.tree$tip.label, 'Eugenia.cocosensis_VS.7857', 'Eugenia.pacifica_VS.7887',
                      'Eugenia.coronata_Flickinger.109', 'Eugenia.talbotii_Parmentier.5037'))

#creating OUwie data frame
traitData <- data.frame(tip.labels=envirn$tip.labels,
                        regimes=NA,
                        X.mean=NA,
                        X.sd=NA)

# Assigning regimes to trait data frame: setting 'Extra' for the tips in palaeo.tips, 'Neotr.' otherwise
traitData$regimes <- ifelse(traitData$tip.labels %in% palaeo.tips, "Extra", "Neotr")
traitData$regimes <- as.factor(traitData$regimes) #as.factor

# Performing ancestral state reconstruction in order to assign regimes to nodes in the phylogeny (required by OUwie)

#putting regimes in a vector
discrete_regimes <- as.character(traitData$regimes)
names(discrete_regimes) <- traitData$tip.labels

#running ancestral state reconstruction
ace_result <- ace(discrete_regimes, pruned_tree, model="ER", type="discrete")

#getting the most likely states for each node
most_likely_states <- colnames(ace_result$lik.anc)[apply(ace_result$lik.anc, 1, which.max)]

pruned_tree$node.label <- most_likely_states #putting ancestral states in tree

#visualising tree with ancestral state reconstructions as node labels

# Define colors for your node labels
colors_map <- setNames(c("red", "blue"), c("Extra", "Neotr"))

# Get colours for each tip based on its regime
tip_colors <- colors_map[as.character(traitData$regimes)]
names(tip_colors) <- traitData$tip.labels # Ensure names match for indexing

# Plot circular tree, add coloured node points
plot(pruned_tree, type = "fan", show.tip.label=F, edge.width=0.5, cex = 0.6)

nodelabels(pch = 21, bg = colors_map[pruned_tree$node.label], cex = 1.5)
tiplabels(pch = 21, bg = tip_colors[pruned_tree$tip.label], cex = 1.5) #phylo is ok

### APPROACH 1. Entire phylogeny, Neotropical vs. Palaeotropical regimes on OUwie
#running OUwie with all possible models for each PCA eingenvectors as traits

eigenvct <- c('PC1', 'PC2', 'PC3') #defining eigenvectors to loop through

summ_list <- list() #creating an empty list to save the three results

for (PC in eigenvct) {
  
  cat("RUNNING ANALYSES FOR:", PC, "\n")
  
  traitData$X.mean <- envirn[[paste0(PC, ".mean")]]
  traitData$X.sd <- envirn[[paste0(PC, ".sd")]]
  
  #distance of the minimum value to the origin (zero)
  displFactor <- abs(min(traitData$X.mean))
  
  # Shifting the mean value of the observed trait data in order to remove negative values (OUwie can't handle them)
  traitData[,'X.mean'] <- traitData[,'X.mean'] + displFactor
  
  #adding the average known tip fog value for species with unknown tip fog values (i.e., those which we have only one data point)
  traitData[,'X.sd'][is.na(traitData[,'X.sd'])] <- mean(traitData$X.sd, na.rm=T)
  
  #running OUwie
  
  # defining models to test (all possible models)
  ouwie_models <- c("BM1", "BMS", "OU1", "OUM", "OUMV") #, "OUMA", "OUMVA")
  
  # Create empty lists to store results
  results_list <- list()
  boot_list <- list()
  
  for (model in ouwie_models) {
    cat("Running model:", model, "\n")
    current_run <- OUwie(
      phy = pruned_tree, 
      data = traitData, 
      model = model, 
      simmap.tree = FALSE, 
      mserr = "known", 
      algorithm = "three.point",
      ub = c(90, 90, 90)
    )
    
    results_list[[model]] <- current_run
    
    if(model == "BM1" | model == "OU1"){
      n_regimes <- length(unique(pruned_tree$node.label))
    }else{
      n_regimes <- 1
    }
    
    boot.reps <- OUwie.boot(phy = pruned_tree,
                            data = traitData,
                            model= model,
                            nboot = 100,
                            alpha = rep(current_run$solution[1,], n_regimes), 
                            sigma.sq = rep(current_run$solution[2,], n_regimes),
                            theta = rep(current_run$theta[,1], n_regimes),
                            theta0 = current_run$theta[1,1],
                            simmap.tree = FALSE, 
                            algorithm = "three.point")
    
    boot.res <- as.data.frame(apply(boot.reps, 2, quantile, probs=c(0.025, 0.975), na.rm=TRUE))
    
    boot_list[[model]] <- boot.res
  }
  
  summ.ouwie <- data.frame(
    model = character(),
    alpha_1 = numeric(),
    alpha_1.LB  = numeric(),
    alpha_1.UB  = numeric(),
    alpha_2 = numeric(),
    alpha_2.LB = numeric(),
    alpha_2.UB = numeric(),
    sig.sq_1 = numeric(),
    sig.sq_1.LB = numeric(),
    sig.sq_1.UB = numeric(),
    sig.sq_2 = numeric(),
    sig.sq_2.LB = numeric(),
    sig.sq_2.UB = numeric(),
    optim_1 = numeric(),
    optim_1.LB = numeric(),
    optim_1.UB = numeric(),
    optim_2 = numeric(),
    optim_2.LB = numeric(),
    optim_2.UB = numeric(),
    lnL = numeric(),
    AIC = numeric(),
    AICc = numeric(),
    stringsAsFactors = F
  )
  
  #populating empty summary data.frame with the results
  
  #BM1 - single-rate (sigma-square) Brownian Motion
  summ.ouwie[1, "model"] <- results_list$BM1$model
  summ.ouwie[1, "sig.sq_1"] <- results_list$BM1$solution[2]
  summ.ouwie[1, "optim_1"] <- results_list$BM1$solution[3]-displFactor # -displFactor orrects for the shift of the mean value performed prior to running OUwie
  summ.ouwie[1, "lnL"] <- results_list$BM1$loglik
  summ.ouwie[1, "AIC"] <- results_list$BM1$AIC
  summ.ouwie[1, "AICc"] <- results_list$BM1$AICc
  
  summ.ouwie[1, "sig.sq_1.LB"] <- boot_list$BM1$sigma.sq_1[1]
  summ.ouwie[1, "sig.sq_1.UB"] <- boot_list$BM1$sigma.sq_1[2]
  summ.ouwie[1, "optim_1.LB"] <- boot_list$BM1$theta_root[1]
  summ.ouwie[1, "optim_1.UB"] <- boot_list$BM1$theta_root[2]
  
  #BMS - Brownian motion with a different sigma-square for each regime
  summ.ouwie[2, "model"] <- results_list$BMS$model
  summ.ouwie[2, "sig.sq_1"] <- results_list$BMS$solution[2]
  summ.ouwie[2, "sig.sq_2"] <- results_list$BMS$solution[5]
  summ.ouwie[2, "optim_1"] <- results_list$BMS$solution[3]-displFactor
  summ.ouwie[2, "lnL"] <- results_list$BMS$loglik
  summ.ouwie[2, "AIC"] <- results_list$BMS$AIC
  summ.ouwie[2, "AICc"] <- results_list$BMS$AICc
  
  summ.ouwie[2, "sig.sq_1.LB"] <- boot_list$BMS$sigma.sq_Extra[1]
  summ.ouwie[2, "sig.sq_1.UB"] <- boot_list$BMS$sigma.sq_Extra[2]
  summ.ouwie[2, "sig.sq_2.LB"] <- boot_list$BMS$sigma.sq_Neotr[1]
  summ.ouwie[2, "sig.sq_2.UB"] <- boot_list$BMS$sigma.sq_Neotr[2]
  summ.ouwie[2, "optim_1.LB"] <- boot_list$BMS$theta_root[1]-displFactor
  summ.ouwie[2, "optim_1.UB"] <- boot_list$BMS$theta_root[2]-displFactor
  
  #OU1 - Ornstein-Uhlenbeck model with a single alpha, sigma-square and optimum for both regimes
  summ.ouwie[3, "model"] <- results_list$OU1$model
  summ.ouwie[3, "alpha_1"] <- results_list$OU1$solution[1]
  summ.ouwie[3, "sig.sq_1"] <- results_list$OU1$solution[2]
  summ.ouwie[3, "optim_1"] <- results_list$OU1$solution[3]-displFactor
  summ.ouwie[3, "lnL"] <- results_list$OU1$loglik
  summ.ouwie[3, "AIC"] <- results_list$OU1$AIC
  summ.ouwie[3, "AICc"] <- results_list$OU1$AICc
  
  summ.ouwie[3, "alpha_1.LB"] <- boot_list$OU1$alpha_1[1]
  summ.ouwie[3, "alpha_1.UB"] <- boot_list$OU1$alpha_1[2]
  summ.ouwie[3, "sig.sq_1.LB"] <- boot_list$OU1$sigma.sq_1[1]
  summ.ouwie[3, "sig.sq_1.UB"] <- boot_list$OU1$sigma.sq_1[2]
  summ.ouwie[3, "optim_1.LB"] <- boot_list$OU1$theta_root[1]-displFactor
  summ.ouwie[3, "optim_1.UB"] <- boot_list$OU1$theta_root[2]-displFactor
  
  #OUM - Ornstein-Uhlenbeck model with a single alpha and sigma-square, but different optima for each regime
  summ.ouwie[4, "model"] <- results_list$OUM$model
  summ.ouwie[4, "alpha_1"] <- results_list$OUM$solution[1]
  summ.ouwie[4, "sig.sq_1"] <- results_list$OUM$solution[2]
  summ.ouwie[4, "optim_1"] <- results_list$OUM$solution[3]-displFactor
  summ.ouwie[4, "optim_2"] <- results_list$OUM$solution[6]-displFactor
  summ.ouwie[4, "lnL"] <- results_list$OUM$loglik
  summ.ouwie[4, "AIC"] <- results_list$OUM$AIC
  summ.ouwie[4, "AICc"] <- results_list$OUM$AICc
  
  summ.ouwie[4, "alpha_1.LB"] <- boot_list$OUM$alpha_Extra[1]
  summ.ouwie[4, "alpha_1.UB"] <- boot_list$OUM$alpha_Extra[2]
  summ.ouwie[4, "sig.sq_1.LB"] <- boot_list$OUM$sigma.sq_Extra[1]
  summ.ouwie[4, "sig.sq_1.UB"] <- boot_list$OUM$sigma.sq_Extra[2]
  summ.ouwie[4, "optim_1.LB"] <- boot_list$OUM$theta_Extra[1]-displFactor
  summ.ouwie[4, "optim_1.UB"] <- boot_list$OUM$theta_Extra[2]-displFactor
  summ.ouwie[4, "optim_2.LB"] <- boot_list$OUM$theta_Neotr[1]-displFactor
  summ.ouwie[4, "optim_2.UB"] <- boot_list$OUM$theta_Neotr[2]-displFactor
  
  
  #OUMV - Ornstein-Uhlenbeck model with a single alpha, but different sigma-squares and optima for each regime
  summ.ouwie[5, "model"] <- results_list$OUMV$model
  summ.ouwie[5, "alpha_1"] <- results_list$OUMV$solution[1]
  summ.ouwie[5, "sig.sq_1"] <- results_list$OUMV$solution[2]
  summ.ouwie[5, "sig.sq_2"] <- results_list$OUMV$solution[5]
  summ.ouwie[5, "optim_1"] <- results_list$OUMV$solution[3]-displFactor
  summ.ouwie[5, "optim_2"] <- results_list$OUMV$solution[6]-displFactor
  summ.ouwie[5, "lnL"] <- results_list$OUMV$loglik
  summ.ouwie[5, "AIC"] <- results_list$OUMV$AIC
  summ.ouwie[5, "AICc"] <- results_list$OUMV$AICc
  
  summ.ouwie[5, "alpha_1.LB"] <- boot_list$OUMV$alpha_Extra[1]
  summ.ouwie[5, "alpha_1.UB"] <- boot_list$OUMV$alpha_Extra[2]
  summ.ouwie[5, "sig.sq_1.LB"] <- boot_list$OUMV$sigma.sq_Extra[1]
  summ.ouwie[5, "sig.sq_1.UB"] <- boot_list$OUMV$sigma.sq_Extra[2]
  summ.ouwie[5, "sig.sq_2.LB"] <- boot_list$OUMV$sigma.sq_Neotr[1]
  summ.ouwie[5, "sig.sq_2.UB"] <- boot_list$OUMV$sigma.sq_Neotr[2]
  summ.ouwie[5, "optim_1.LB"] <- boot_list$OUMV$theta_Extra[1]-displFactor
  summ.ouwie[5, "optim_1.UB"] <- boot_list$OUMV$theta_Extra[2]-displFactor
  summ.ouwie[5, "optim_2.LB"] <- boot_list$OUMV$theta_Neotr[1]-displFactor
  summ.ouwie[5, "optim_2.UB"] <- boot_list$OUMV$theta_Neotr[2]-displFactor
  
  #OUMA - Ornstein-Uhlenbeck model with a single sigma-square, but different alpha and optima for each regime
  #summ.ouwie[6, "model"] <- results_list$OUMA$model
  #summ.ouwie[6, "alpha_1"] <- results_list$OUMA$solution[1]
  #summ.ouwie[6, "alpha_2"] <- results_list$OUMA$solution[4]
  #summ.ouwie[6, "sig.sq_1"] <- results_list$OUMA$solution[2]
  #summ.ouwie[6, "optim_1"] <- results_list$OUMA$solution[3]-displFactor
  #summ.ouwie[6, "optim_2"] <- results_list$OUMA$solution[6]-displFactor
  #summ.ouwie[6, "lnL"] <- results_list$OUMA$loglik
  #summ.ouwie[6, "AIC"] <- results_list$OUMA$AIC
  #summ.ouwie[6, "AICc"] <- results_list$OUMA$AICc
  
  #OUMVA - Ornstein-Uhlenbeck model with different alpha, sigma-square and optima for each regime
  #summ.ouwie[7, "model"] <- results_list$OUMVA$model
  #summ.ouwie[7, "alpha_1"] <- results_list$OUMVA$solution[1]
  #summ.ouwie[7, "alpha_2"] <- results_list$OUMVA$solution[4]
  #summ.ouwie[7, "sig.sq_1"] <- results_list$OUMVA$solution[2]
  #summ.ouwie[7, "sig.sq_2"] <- results_list$OUMVA$solution[5]
  #summ.ouwie[7, "optim_1"] <- results_list$OUMVA$solution[3]-displFactor
  #summ.ouwie[7, "optim_2"] <- results_list$OUMVA$solution[6]-displFactor
  #summ.ouwie[7, "lnL"] <- results_list$OUMVA$loglik
  #summ.ouwie[7, "AIC"] <- results_list$OUMVA$AIC
  #summ.ouwie[7, "AICc"] <- results_list$OUMVA$AICc
  
  #sort summary OUwie data table by AICc
  summ.ouwie <- summ.ouwie[order(summ.ouwie$AICc),]
  
  summ_list[[PC]] <- summ.ouwie #store for saving
}

#Save result list as an R object
save(results_list, file="data/26_ouwieResA.RData")

# Write to Excel workbook
wb <- createWorkbook()
for (PC in eigenvct) {
  addWorksheet(wb, PC)
  writeData(wb, PC, summ_list[[PC]])
}

saveWorkbook(wb, "results/3_nichevol_ouwieA.xlsx")

### APPROACH 2. A per-clade approach using BM1 models

#creating OUwie data frame
traitData2 <- data.frame(tip.labels=envirn$tip.labels,
                         sections=envirn$sections,
                         regimes=NA,
                         PC1.mean=envirn$PC1.mean,
                         PC1.sd=envirn$PC1.sd,
                         PC2.mean=envirn$PC2.mean,
                         PC2.sd=envirn$PC2.sd,
                         PC3.mean=envirn$PC3.mean,
                         PC3.sd=envirn$PC3.sd)

# Assigning regimes to trait data frame: setting 'Extra' for the tips in palaeo.tips, 'Neotr.' otherwise
traitData2$regimes <- ifelse(traitData2$tip.labels %in% palaeo.tips, "Extra", "Neotr")
traitData2$regimes <- as.factor(traitData2$regimes) #as.factor

# Create empty lists to store results

ouwie_list <- list()
boot2_list <- list()

eigenvct <- c('PC1', 'PC2', 'PC3') #defining eigenvectors to loop through

for (env_var in eigenvct){
  
  cat("\nPROCESSING EIGENVECTOR:", env_var, "\n")
  
  clades <- as.vector(unique(envirn$sections)) #vector containing names of clades

  for (clade in clades){
    
    cat("\nProcessing clade:", clade, "\n")
    
    clade_data <- traitData2[traitData2$sections == clade, ] #filtering only the clade being analysed
    
    cols_to_keep <- c("tip.labels", "sections", "regimes",
                      paste0(env_var, ".mean"), paste0(env_var, ".sd"))
    
    clade_data <- clade_data %>% select(all_of(cols_to_keep), -sections) #keeping only information of the eigenvector being analysed
    
    #OUwie cannot analyse traits that receive negative values. Let's shift all values so the minimum is zero.
    min_value <- min(clade_data[[paste0(env_var, ".mean")]], na.rm = TRUE)
    clade_data[[paste0(env_var, ".mean")]] <- clade_data[[paste0(env_var, ".mean")]] + abs(min_value)
    
    #adding the average known tip fog value for species with unknown tip fog values (i.e., those which we have only one data point)
    clade_data[[paste0(env_var, ".sd")]][is.na(clade_data[[paste0(env_var, ".sd")]])
    ] <- mean(clade_data[[paste0(env_var, ".sd")]], na.rm = TRUE)
    
    #prunning phylogeny to keep only tips that belong to a certain clade
    
    pruned_tree2 <- drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, clade_data$tip.labels))
    
    fit_BM <- OUwie(
      phy = pruned_tree2, 
      data = clade_data, 
      model = "BM1", 
      simmap.tree = FALSE, 
      mserr = "known", 
      algorithm = "three.point",
      ub = c(90, 90, 90)
    )
    
    ouwie_list[[paste(env_var, clade, sep = "_")]] <- fit_BM
    
    #  n_regimes <- length(unique(pruned_tree2$node.label))
    
    boot.reps <- OUwie.boot(phy = pruned_tree2,
                            data = clade_data,
                            model= "BM1",
                            nboot = 100,
                            alpha = rep(fit_BM$solution[1,], 2), 
                            sigma.sq = rep(max(fit_BM$solution[2,], 1e-8), 2),
                            theta = rep(fit_BM$theta[,1], 2),
                            theta0 = fit_BM$theta[1,1],
                            simmap.tree = FALSE, 
                            algorithm = "three.point")
    
    boot.res <- as.data.frame(apply(boot.reps, 2, quantile, probs=c(0.025, 0.975), na.rm=TRUE))
    
    boot2_list[[paste(env_var, clade, sep = "_")]] <- boot.res
  }
}

#Save result lists as an R object
save(ouwie_list, boot2_list, file="data/27_ouwieResB.RData")

# summarising results

summ.ouwie2 <- data.frame( #creating an empty data.frame 
  clade = character(),
  sig.sq = numeric(),
  sig.sq.LB = numeric(),
  sig.sq.UB = numeric(),
  lnL = numeric(),
  AIC = numeric(),
  AICc = numeric(),
  stringsAsFactors = F)

#populating empty summary data.frame with the OUwie and OUwie.boot results

res.rows <- names(ouwie_list) # creating a vector of objects inside ouwie_list (equal to boot2_list)

for (i in 1:length(ouwie_list)){ # iterating through objects within both lists
  summ.ouwie2[i, "clade"] <- names(ouwie_list)[i]
  summ.ouwie2[i, "sig.sq"] <- ouwie_list[[i]]$solution[2]
  summ.ouwie2[i, "lnL"] <- ouwie_list[[i]]$loglik
  summ.ouwie2[i, "AIC"] <- ouwie_list[[i]]$AIC
  summ.ouwie2[i, "AICc"] <- ouwie_list[[i]]$AICc
  
  summ.ouwie2[i, "sig.sq.LB"] <- boot2_list[[i]]$sigma.sq_1[1]
  summ.ouwie2[i, "sig.sq.UB"] <- boot2_list[[i]]$sigma.sq_1[2]
}

# sorting table by highest to lowers sig-sq value, within each eigenvector
PC_Group_Vector <- substr(summ.ouwie2$clade, start = 1, stop = 3)
summ.ouwie2 <- summ.ouwie2[order(PC_Group_Vector, -summ.ouwie2$sig.sq), ]

summ.ouwie3 <- summ.ouwie2 %>% mutate(
  sig.sq = as.numeric(format(sig.sq, scientific = FALSE)),
  sig.sq.LB = as.numeric(format(sig.sq.LB, scientific = FALSE)),
  sig.sq.UB = as.numeric(format(sig.sq.UB, scientific = FALSE)),
  lnL = as.numeric(format(lnL, scientific = FALSE)),
  AIC = as.numeric(format(AIC, scientific = FALSE)),
  AICc = as.numeric(format(AICc, scientific = FALSE))
)

write.xlsx(summ.ouwie3, file="results/4_nichevol_ouwieB.xlsx", rowNames=F)