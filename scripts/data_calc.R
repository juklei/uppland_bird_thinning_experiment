# dir("clean")

ldm <- read.csv("data/long_distance_migrants.csv")
bpo <- read.csv("clean/bpo_double.csv")
forest <- read.csv("clean/forest_experiment_data_JAGS.csv")
f_diff <- read.csv("clean/forest_experiment_data.csv")
# head(bpo)
# str(bpo)
# head(forest)

# ## Create matrix with seen at least once per year:
# T1 <- acast(bpo[, c("species", "n_obs", "obs_year")],
#             species ~ obs_year,
#             value.var = "n_obs",
#             fun.aggregate = function(x) sum(x) > 0)
# 
# ## Exclude species which where not seen at least once per year:
# red_names <- names(which(rowSums(T1) == 3))
# bpo <- droplevels(bpo[bpo$species %in% red_names, ])

# ## Only for one observer:
# bpo <- droplevels(bpo[bpo$observer == "jkn", ])
# forest <- droplevels(forest[forest$plot %in% bpo$plot, ])

## Create experiment variable where true control is not part of the experimental evaluation:
forest$e_within <- forest$experiment
forest$e_within[forest$treatment == "TC"] <- "before"

## Create treatment variable comparing true controls with intervention sites:
forest$TCvsNC <- ifelse(forest$treatment == "TC", "TC", "NC")

## Create post-treatments based on above or below mean reduction or increase of forest structural variables:

f_diff <- f_diff[f_diff$experiment == "difference" & f_diff$treatment %in% c("T", "URT"), ]

## For BA:
f_diff$BA_all <- ifelse(f_diff$BA > median(f_diff$BA), "above", "below")
f_diff$BA_dec <- ifelse(f_diff$BA_lov > median(f_diff$BA_lov), "above", "below")
f_diff$BA_spruce <- ifelse(f_diff$BA_gran > median(f_diff$BA_gran), "above", "below")
f_diff$BA_dw <- ifelse(f_diff$BA_dv > median(f_diff$BA_dv), "above", "below")

## For visibility:
f_diff$V <- ifelse(f_diff$laser_mean > median(f_diff$laser_mean), "above","below")

## Add to forest file and add TC and C to the respective treatments:
forest <- merge(forest, f_diff[, c(2,28:32)], by = "plot", all.x = TRUE)
forest[forest$treatment %in% c("C", "TC"), c(31:35)] <- 
  forest$treatment[forest$treatment %in% c("C", "TC")]

## Change characters to factors:
forest <- as.data.frame(unclass(forest))

## Make numeric levels for the treatment and experiment matrices:
forest$exp_num <- ifelse(forest$experiment == "before", 1, 2)
forest$e_within_num <- ifelse(forest$e_within == "before", 1, 2)
forest$treatment_num <- as.numeric(forest$treatment)
forest$TCvsNC_num <- as.numeric(forest$TCvsNC)
forest$BA_all_num <- as.numeric(forest$BA_all)
forest$BA_dec_num <- as.numeric(forest$BA_dec)
forest$BA_spruce_num <- as.numeric(forest$BA_spruce)
forest$BA_dw_num <- as.numeric(forest$BA_dw)
forest$V_num <- as.numeric(forest$V)

## Create data arrays for process model part:
exp_all <- acast(forest[, c("plot", "year", "exp_num")], year ~ plot)
e_within <- acast(forest[, c("plot", "year", "e_within_num")], year ~ plot)

## Create model data set:
data <- list(nobs = nrow(bpo),
             nspecies = max(as.numeric(bpo$species)),
             nyears = dim(exp_all)[1],
             nsites = dim(exp_all)[2],
             nblocks = nlevels(bpo$block),
             observed = bpo$n_obs,
             nvisits = bpo$n_visits,
             species = as.numeric(bpo$species),
             year = bpo$obs_year-(min(bpo$obs_year)-1),
             site = as.numeric(bpo$plot),
             block = as.numeric(bpo$block),
             observer = ifelse(bpo$observer == "jkn", 0, 1),
             exp = if(ref == "TC"){exp_all} else{e_within},
             treat = acast(unique(forest[, c("plot", paste0(TT, "_num"))]), . ~ plot)[1, ],
             eval = which(!levels(forest[, TT]) %in% c(ref, "TC")) ,
             ref = which(levels(forest[, TT]) == ref)) 

str(data)