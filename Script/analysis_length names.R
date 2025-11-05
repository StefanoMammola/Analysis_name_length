###########################################################################

# Species name length

###########################################################################

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.1.0) and R studio (v. 1.4.1103)
# Author: Stefano Mammola

# load R packages ---------------------------------------------------------

pacman::p_load("dplyr",
               "glmmTMB",
               "tidyverse",
               "ggplot2",
               "stringr",
               "glue",
               "MuMIn",
               "patchwork",
               "performance")

# Functions ------------------------------------------------------------

latin.readability <- function(name) {
  name <- tolower(name) # Convert to lowercase
  name_length <- nchar(gsub(" ", "", name))  # Count total characters (excluding spaces)
  consonant_clusters <- sum(gregexpr("[bcdfghjklmnpqrstvwxyz]{2,}", name)[[1]] > 0) # Count consonant clusters
  uncommon_letters <- sum(strsplit(name, NULL)[[1]] %in% c("j", "k", "x", "y", "z")) # Count uncommon letters
  
  readability <- 0.5 * name_length + 2 * consonant_clusters + uncommon_letters
  
  return(c(Length = name_length, Clusters = consonant_clusters, Uncommon = uncommon_letters, Readability = readability))
}


# Load data ---------------------------------------------------------------

db2  <- read.csv(file = "Data/full_data.csv", sep = ',', header = TRUE, as.is = FALSE)

str(db2)
head(db2,5)

# Data preparation --------------------------------------------------------

#Extracting description year
db2$Year_description <- sapply(stringr::str_extract_all(db2$author, "\\b\\d{4}\\b"), `[`, 1)

db2$Year_description <- as.numeric(db2$Year_description)

# Apply function for name length and readability and format results
db <- do.call(rbind, lapply(paste(db2$genus,db2$species), latin.readability))

db <- data.frame(kingdom = db2$kingdom, 
                 phylum = db2$phylum,
                 class =  db2$class,
                 order = db2$order,
                 family = db2$family,
                 genus = db2$genus,
                 species = db2$species, 
                 year = db2$Year_description,
                 db, 
                 citations = db2$Total_wos, 
                 wiki = db2$total_wiki_pgviews)

db$Readability <- db$Readability*10 #to get read of decimals

# Data exploration --------------------------------------------------------

GGally::ggpairs(data = results)

# Modelling ---------------------------------------------------------------

########################
# Citation model with length
########################

formula_m1 <- as.formula("citations ~ year + Length + (1 | phylum / class / order)")

# m1 <- glmmTMB(formula_m1, data = db, 
#               family = poisson,
#               control=glmmTMBControl(optimizer=optim,
#                                      optArgs=list(method="BFGS")))
# 
# performance::check_overdispersion(m1)

m1 <- glmmTMB(formula_m1, data = db, 
              family = nbinom2,
              control=glmmTMBControl(optimizer=optim,
                                     optArgs=list(method="BFGS")))

summary(m1)
# performance::check_model(m1)

# Model prediction
newdat <- data.frame(
  Length = seq(min(db$Length, na.rm = TRUE),
               max(db$Length, na.rm = TRUE),
               length.out = 100),
  year   = mean(db$year, na.rm = TRUE),   # year effects averaged
  phylum = NA, class = NA, order = NA     # random effects averaged
)

pred <- predict(m1, newdat, type = "response", se.fit = TRUE, re.form = NA, allow.new.levels=TRUE)
newdat$fit <- pred$fit
newdat$se  <- pred$se.fit
newdat$lcl <- newdat$fit - 1.96*newdat$se
newdat$ucl <- newdat$fit + 1.96*newdat$se

# Extract slope and p-value for Length
slope <- tidy(m1, effects = "fixed")$estimate[tidy_m1$term == "Length"]
pval  <- tidy(m1, effects = "fixed")$p.value[tidy_m1$term == "Length"]
rsq <- MuMIn::r.squaredGLMM(m1)[2]  # R² (conditional = fixed + random effects)

label_text <- glue::glue(
  "Slope = {round(slope, 3)}\n",
  "p = {format.pval(pval, digits = 3)}\n",
  "R² = {round(rsq, 3)}"
)

(plot_1 <- ggplot() +
  geom_point(data=db, aes(Length, citations), alpha=0.2) +
  geom_line(data=newdat, aes(Length, fit), color="blue", size=1.2) +
  geom_ribbon(data=newdat, aes(Length, ymin=lcl, ymax=ucl), alpha=0.2) +
  scale_y_continuous(trans = scales::pseudo_log_trans(), 
                     breaks = c(0, 1, 10, 100, 1000, 10000))+
    annotate("text",
             x = Inf, y = Inf,
             label = label_text,
             hjust = 1.1, vjust = 1.1,
             size = 4) +
  labs(x = "Species name length [number of characters]", 
       y = "Citations (log-scaled axis)")+
    theme_minimal(base_size = 12))

########################
# Citation model with Readability
########################

formula_m2 <- as.formula("citations ~ year + Readability + (1 | phylum / class / order)")

# m2 <- glmmTMB(formula_m2, data = db, 
#               family = poisson,
#               control=glmmTMBControl(optimizer=optim,
#                                      optArgs=list(method="BFGS")))
# 
# performance::check_overdispersion(m2)

m2 <- glmmTMB(formula_m2, data = db, 
              family = nbinom2,
              control=glmmTMBControl(optimizer=optim,
                                     optArgs=list(method="BFGS")))

summary(m2)
# performance::check_model(m2)

# Model prediction
newdat <- data.frame(
  Readability = seq(min(db$Readability, na.rm = TRUE),
               max(db$Readability, na.rm = TRUE),
               length.out = 100),
  year   = mean(db$year, na.rm = TRUE),   # year effects averaged
  class = NA, order = NA     # random effects averaged
)

pred <- predict(m2, newdat, type = "response", se.fit = TRUE, re.form = NA, allow.new.levels=TRUE)
newdat$fit <- pred$fit
newdat$se  <- pred$se.fit
newdat$lcl <- newdat$fit - 1.96*newdat$se
newdat$ucl <- newdat$fit + 1.96*newdat$se

# Extract slope and p-value for Length
slope <- tidy(m2, effects = "fixed")$estimate[tidy_m1$term == "Length"]
pval  <- tidy(m2, effects = "fixed")$p.value[tidy_m1$term == "Length"]
rsq <- MuMIn::r.squaredGLMM(m2)[2]  # R² (conditional = fixed + random effects)

label_text2 <- glue::glue(
  "Slope = {round(slope, 3)}\n",
  "p = {format.pval(pval, digits = 3)}\n",
  "R² = {round(rsq, 3)}"
)

(plot_2 <- ggplot() +
    geom_point(data=db, aes(Readability, citations), alpha=0.2) +
    geom_line(data=newdat, aes(Readability, fit), color="blue", size=1.2) +
    geom_ribbon(data=newdat, aes(Readability, ymin=lcl, ymax=ucl), alpha=0.2) +
    scale_y_continuous(trans = scales::pseudo_log_trans(), 
                       breaks = c(0, 1, 10, 100, 1000, 10000))+
    #scale_y_continuous(trans="log10") +
    annotate("text",
             x = Inf, y = Inf,
             label = label_text2,
             hjust = 1.1, vjust = 1.1,
             size = 4) +
    labs(x = "Species name readability", 
         y = "Citations (log-scaled axis)")+
    theme_minimal(base_size = 12))

########################
# Wiki model with Length
########################

formula_m3 <- as.formula("wiki ~ year + Length + (1 | phylum / class / order)")

# m3 <- glmmTMB(formula_m3, data = db, 
#               family = poisson,
#               control=glmmTMBControl(optimizer=optim,
#                                      optArgs=list(method="BFGS")))
# 
# performance::check_overdispersion(m3)

m3 <- glmmTMB(formula_m3, data = db, 
              family = nbinom2,
              control=glmmTMBControl(optimizer=optim,
                                     optArgs=list(method="BFGS")))

summary(m3)
# performance::check_model(m3)

# Model prediction
newdat <- data.frame(
  Length = seq(min(db$Length, na.rm = TRUE),
               max(db$Length, na.rm = TRUE),
               length.out = 100),
  year   = mean(db$year, na.rm = TRUE),   # year effects averaged
  phylum = NA, class = NA, order = NA     # random effects averaged
)

pred <- predict(m3, newdat, type = "response", se.fit = TRUE, re.form = NA, allow.new.levels=TRUE)
newdat$fit <- pred$fit
newdat$se  <- pred$se.fit
newdat$lcl <- newdat$fit - 1.96*newdat$se
newdat$ucl <- newdat$fit + 1.96*newdat$se

# Extract slope and p-value for Length
slope <- tidy(m3, effects = "fixed")$estimate[tidy_m1$term == "Length"]
pval  <- tidy(m3, effects = "fixed")$p.value[tidy_m1$term == "Length"]
rsq <- MuMIn::r.squaredGLMM(m3)[2]  # R² (conditional = fixed + random effects)

label_text3 <- glue::glue(
  "Slope = {round(slope, 3)}\n",
  "p = {format.pval(pval, digits = 3)}\n",
  "R² = {round(rsq, 3)}"
)

(plot_3 <- ggplot() +
    geom_point(data=db, aes(Length, wiki), alpha=0.2) +
    geom_line(data=newdat, aes(Length, fit), color="blue", size=1.2) +
    geom_ribbon(data=newdat, aes(Length, ymin=lcl, ymax=ucl), alpha=0.2) +
    scale_y_continuous(trans = scales::pseudo_log_trans(), 
                       breaks = c(0, 10, 100, 10000, 100000, 100000000))+
    annotate("text",
             x = Inf, y = Inf,
             label = label_text3,
             hjust = 1.1, vjust = 1.1,
             size = 4) +
    labs(x = "Species name length [number of characters]", 
         y = "Wikipedia views (log-scaled axis)")+
    theme_minimal(base_size = 12))

########################
# Wiki model with Readability
########################

formula_m4 <- as.formula("wiki ~ year + Readability + (1 | phylum / class / order)")

# m4 <- glmmTMB(formula_m4, data = db, 
#               family = poisson,
#               control=glmmTMBControl(optimizer=optim,
#                                      optArgs=list(method="BFGS")))
# 
# performance::check_overdispersion(m4)

m4 <- glmmTMB(formula_m4, data = db, 
              family = nbinom2,
              control=glmmTMBControl(optimizer=optim,
                                     optArgs=list(method="BFGS")))

summary(m4)
# performance::check_model(m4)

# Model prediction
newdat <- data.frame(
  Readability = seq(min(db$Readability, na.rm = TRUE),
               max(db$Readability, na.rm = TRUE),
               length.out = 100),
  year   = mean(db$year, na.rm = TRUE),   # year effects averaged
  phylum = NA, class = NA, order = NA     # random effects averaged
)

pred <- predict(m4, newdat, type = "response", se.fit = TRUE, re.form = NA, allow.new.levels=TRUE)
newdat$fit <- pred$fit
newdat$se  <- pred$se.fit
newdat$lcl <- newdat$fit - 1.96*newdat$se
newdat$ucl <- newdat$fit + 1.96*newdat$se

# Extract slope and p-value for Length
slope <- tidy(m4, effects = "fixed")$estimate[tidy_m1$term == "Length"]
pval  <- tidy(m4, effects = "fixed")$p.value[tidy_m1$term == "Length"]
rsq <- MuMIn::r.squaredGLMM(m4)[2]  # R² (conditional = fixed + random effects)

label_text4 <- glue::glue(
  "Slope = {round(slope, 3)}\n",
  "p = {format.pval(pval, digits = 3)}\n",
  "R² = {round(rsq, 3)}"
)

(plot_4 <- ggplot() +
    geom_point(data=db, aes(Readability, wiki), alpha=0.2) +
    geom_line(data=newdat, aes(Readability, fit), color="blue", size=1.2) +
    geom_ribbon(data=newdat, aes(Readability, ymin=lcl, ymax=ucl), alpha=0.2) +
    scale_y_continuous(trans = scales::pseudo_log_trans(), 
                       breaks = c(0, 10, 100, 10000, 100000, 100000000))+
    annotate("text",
             x = Inf, y = Inf,
             label = label_text4,
             hjust = 1.1, vjust = 1.1,
             size = 4) +
    labs(x = "Species name readability", 
         y = "Wikipedia views (log-scaled axis)")+
    theme_minimal(base_size = 12))

# Final plot --------------------------------------------------------------

pdf(file = "Figures/Figure_1.pdf", width = 8, height = 8)

ggpubr::ggarrange(plot_1,plot_2,plot_3,plot_4,
                 
                  common.legend = FALSE,
                  hjust = 0,
                  #align = "h",
                  labels = c("A", "B", "C", "D"),
                  ncol=2, nrow=2) 

dev.off()


