###########################################################################

# The power of naming: shorter and simpler species names draw more attention

###########################################################################

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.4.1) and R studio (v. 2023.12.0+369)
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
  
  readability <- 2 * consonant_clusters + uncommon_letters
  
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

# Data exploration --------------------------------------------------------

GGally::ggpairs(data = db, columns = c("Length", "Readability"))

cor.test(db$Readability,db$Length)

m0 <- glm(Readability ~ Length, family = "poisson", data = db)
performance::check_overdispersion(m0)
summary(m0)

db$resid_readability <- residuals(m0, type = "response")

#
db |>
  dplyr::arrange(desc(resid_readability)) |>  
  slice(1:3)

db |>
  dplyr::arrange(resid_readability) |>  
  slice(1:3)


# Figure S1 ---------------------------------------------------------------

(Fig_S1A <- ggplot(data=db, aes(Length)) +
   geom_histogram(fill = "grey5", bins = 25)+ 
  labs(x = "Species name length [number of characters]", 
       y = "Count")+
  theme_minimal(base_size = 12))

(Fig_S1B <- ggplot(data=db, aes(x = resid_readability)) +
    geom_histogram(fill = "grey5", bins = 25)+ 
    labs(x = "Reading difficulty [residuals]", 
         y = "")+
    theme_minimal(base_size = 12))


pdf(file = "Figures/Figure_S1.pdf", width = 8, height = 4)

ggpubr::ggarrange(Fig_S1A,Fig_S1B,
                  common.legend = FALSE,
                  hjust = 0,
                  #align = "h",
                  labels = c("A", "B"),
                  ncol=2, nrow=1) 

dev.off()

# Modelling ---------------------------------------------------------------

########################
# Citation model with length
########################

formula_m1 <- as.formula("citations ~ year + Length + resid_readability + (1 | phylum / class / order)")

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
performance::check_zeroinflation(m1)
performance::check_collinearity(m1)
# performance::check_model(m1)


# Model prediction
newdat <- data.frame(
  Length = seq(min(db$Length, na.rm = TRUE),
               max(db$Length, na.rm = TRUE),
               length.out = 100),
  resid_readability = mean(db$resid_readability, na.rm = TRUE), #readability effect is averaged
  year   = mean(db$year, na.rm = TRUE),   # year effects averaged
  phylum = NA, class = NA, order = NA     # random effects averaged
)

pred <- predict(m1, newdat, type = "response", se.fit = TRUE, re.form = NA, allow.new.levels=TRUE)
newdat$fit <- pred$fit
newdat$se  <- pred$se.fit
newdat$lcl <- newdat$fit - 1.96*newdat$se
newdat$ucl <- newdat$fit + 1.96*newdat$se

# Extract slope and p-value for Length
slope <- summary(m1)$coefficients$cond[3,1]
pval  <- summary(m1)$coefficients$cond[3,4]
rsq <- MuMIn::r.squaredGLMM(m1)[2]  # R² (conditional = fixed + random effects)

label_text <- glue::glue(
  "Slope = {round(slope, 3)}\n",
  "p = {format.pval(pval, digits = 3)}\n",
  "R² = {round(rsq, 3)}"
)

(plot_1 <- ggplot() +
  geom_point(data=db, aes(Length, citations), alpha=0.2) +
  geom_line(data=newdat, aes(Length, fit), color="blue", linewidth = 1.2) +
  geom_ribbon(data=newdat, aes(Length, ymin=lcl, ymax=ucl), fill = "blue", alpha=0.2) +
  scale_y_continuous(trans = scales::pseudo_log_trans(), 
                     breaks = c(0, 1, 10, 100, 1000, 10000))+
    annotate("text",
             x = Inf, y = Inf,
             label = label_text,
             hjust = 1.1, vjust = 1.1,
             size = 4) +
  labs(x = "", 
       y = "Citations (log-scaled axis)")+
    theme_minimal(base_size = 12))

# plot readability

# Model prediction
newdat <- data.frame(
  resid_readability = seq(min(db$resid_readability, na.rm = TRUE),
                          max(db$resid_readability, na.rm = TRUE),
                          length.out = 100),
  Length   = mean(db$Length, na.rm = TRUE),   # Length effects averaged
  year   = mean(db$year, na.rm = TRUE),   # year effects averaged
  phylum = NA, class = NA, order = NA     # random effects averaged
)

pred <- predict(m1, newdat, type = "response", se.fit = TRUE, re.form = NA, allow.new.levels=TRUE)
newdat$fit <- pred$fit
newdat$se  <- pred$se.fit
newdat$lcl <- newdat$fit - 1.96*newdat$se
newdat$ucl <- newdat$fit + 1.96*newdat$se

# Extract slope and p-value for Length
slope <- summary(m1)$coefficients$cond[4,1]
pval  <- summary(m1)$coefficients$cond[4,4]

label_text4 <- glue::glue(
  "Slope = {round(slope, 3)}\n",
  "p = {format.pval(pval, digits = 3)}\n"
)

(plot_2 <- ggplot() +
    geom_point(data=db, aes(resid_readability, citations), alpha=0.2) +
    geom_line(data=newdat, aes(resid_readability, fit), color="blue", linewidth=1.2) +
    geom_ribbon(data=newdat, aes(resid_readability, ymin=lcl, ymax=ucl),  fill = "blue", alpha=0.2) +
    scale_y_continuous(trans = scales::pseudo_log_trans(), 
                       breaks = c(0, 10, 100, 10000, 100000, 100000000))+
    annotate("text",
             x = Inf, y = Inf,
             label = label_text4,
             hjust = 1.1, vjust = 1.1,
             size = 4) +
    labs(x = "", 
         y = "")+
    theme_minimal(base_size = 12))

########################
# Wiki model with Length
########################

formula_m3 <- as.formula("wiki ~ year + Length + resid_readability + (1 | phylum / class / order)")

m3 <- glmmTMB(formula_m3, data = db, 
              family = nbinom2,
              control=glmmTMBControl(optimizer=optim,
                                     optArgs=list(method="BFGS")))

summary(m3)
performance::check_zeroinflation(m3)
performance::check_collinearity(m3)
# performance::check_model(m3)

# Model prediction
newdat <- data.frame(
  Length = seq(min(db$Length, na.rm = TRUE),
               max(db$Length, na.rm = TRUE),
               length.out = 100),
  resid_readability = mean(db$resid_readability, na.rm = TRUE), # readability effect averaged
  year   = mean(db$year, na.rm = TRUE),   # year effects averaged
  phylum = NA, class = NA, order = NA     # random effects averaged
)

pred <- predict(m3, newdat, type = "response", se.fit = TRUE, re.form = NA, allow.new.levels=TRUE)
newdat$fit <- pred$fit
newdat$se  <- pred$se.fit
newdat$lcl <- newdat$fit - 1.96*newdat$se
newdat$ucl <- newdat$fit + 1.96*newdat$se

# Extract slope and p-value for Length
slope <- summary(m3)$coefficients$cond[3,1]
pval  <- summary(m3)$coefficients$cond[3,4]
rsq <- MuMIn::r.squaredGLMM(m3)[2]  # R² (conditional = fixed + random effects)

label_text3 <- glue::glue(
  "Slope = {round(slope, 3)}\n",
  "p = {format.pval(pval, digits = 3)}\n",
  "R² = {round(rsq, 3)}"
)

(plot_3 <- ggplot() +
    geom_point(data=db, aes(Length, wiki), alpha=0.2) +
    geom_line(data=newdat, aes(Length, fit), color="blue", linewidth = 1.2) +
    geom_ribbon(data=newdat, aes(Length, ymin=lcl, ymax=ucl), fill = "blue", alpha=0.2) +
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

# plot with Readability

# Model prediction
newdat <- data.frame(
  resid_readability = seq(min(db$resid_readability, na.rm = TRUE),
                          max(db$resid_readability, na.rm = TRUE),
                          length.out = 100),
  Length   = mean(db$Length, na.rm = TRUE),   # Length effects averaged
  year   = mean(db$year, na.rm = TRUE),   # year effects averaged
  phylum = NA, class = NA, order = NA     # random effects averaged
)

pred <- predict(m3, newdat, type = "response", se.fit = TRUE, re.form = NA, allow.new.levels=TRUE)
newdat$fit <- pred$fit
newdat$se  <- pred$se.fit
newdat$lcl <- newdat$fit - 1.96*newdat$se
newdat$ucl <- newdat$fit + 1.96*newdat$se

# Extract slope and p-value for Length
slope <- summary(m3)$coefficients$cond[4,1]
pval  <- summary(m3)$coefficients$cond[4,4]

label_text4 <- glue::glue(
  "Slope = {round(slope, 3)}\n",
  "p = {format.pval(pval, digits = 3)}\n"
)

(plot_4 <- ggplot() +
    geom_point(data=db, aes(resid_readability, wiki), alpha=0.2) +
    geom_line(data=newdat, aes(resid_readability, fit), color="blue", linewidth=1.2) +
    geom_ribbon(data=newdat, aes(resid_readability, ymin=lcl, ymax=ucl),fill = "blue",  alpha=0.2) +
    scale_y_continuous(trans = scales::pseudo_log_trans(),
                       breaks = c(0, 10, 100, 10000, 100000, 100000000))+
    annotate("text",
             x = Inf, y = Inf,
             label = label_text4,
             hjust = 1.1, vjust = 1.1,
             size = 4) +
    labs(x = "Reading difficulty", 
         y = "")+
    theme_minimal(base_size = 12))

# Final plot --------------------------------------------------------------

pdf(file = "Figures/Figure_1.pdf", width = 9, height = 8)

ggpubr::ggarrange(plot_1,plot_2,plot_3,plot_4,
                  common.legend = FALSE,
                  hjust = 0,
                  #align = "h",
                  labels = c("A", "B", "C", "D"),
                  ncol=2, nrow=2) 

dev.off()


# db$wiki01 <- ifelse(db$wiki>0, 1, 0)
# table(db$wiki01)
# 
# formula_m5 <- as.formula("wiki01 ~ year + Length + (1 | phylum / class / order)")
# 
# m5 <- glmmTMB(formula_m5, data = db, 
#               family = binomial,
#               control=glmmTMBControl(optimizer=optim,
#                                      optArgs=list(method="BFGS")))
# summary(m5)

# formula_m6 <- as.formula("wiki01 ~ year + Readability + (1 | phylum / class / order)")
# 
# m6 <- glmmTMB(formula_m6, data = db,
#               family = binomial,
#               control=glmmTMBControl(optimizer=optim,
#                                      optArgs=list(method="BFGS")))
# summary(m6)

# Data transformation -----------------------------------------------------

# Balancing levels IUCN
db2$IUCN_rec <- db2$IUCN

levels(db2$IUCN_rec) <- c("Endangered", "Unknown", "Endangered", "Endangered", "Least concern","Unknown","Least concern","Endangered")
table(db2$IUCN_rec)

db2$IUCN_rec <- relevel(db2$IUCN_rec, "Unknown") #setting baseline

# Balancing levels Domain
db2$domain_rec <- db2$domain

levels(db2$domain_rec) <- c("freshwater","multiple","multiple","multiple",
                           "marine", "multiple", "terrestrial", "multiple", "terrestrial")

db2$domain_rec <- relevel(db2$domain_rec, "multiple") #setting baseline

table(db2$domain_rec)

# Homogenize distribution
db2 <- db2 %>% 
  dplyr::mutate(log_uniqueness_family = log(uniqueness_family+1),
                log_uniqueness_genus = log(uniqueness_genus+1),
                log_range_size = log(range_size+1),
                log_size_avg = log(size_avg+1),
                log_distance_hu = log(mean_divergence_time_Mya+1))

db2 <- db2 %>% 
  dplyr::mutate(scaled_uniqueness_family = scale(log_uniqueness_family, center = TRUE, scale = TRUE),
                scaled_log_distance_human = scale(log_distance_hu, center = TRUE, scale = TRUE),
                scaled_range_size = scale(log_range_size, center = TRUE, scale = TRUE),
                scaled_size = scale(log_size_avg, center = TRUE, scale = TRUE))

db3 <- data.frame(db,
            db2 |> dplyr::select(biogeography,
            scaled_size,
            colorful,
            color_blu,
            color_red,
            scaled_range_size,
            domain_rec,
            IUCN_rec,
            scaled_uniqueness_family,
            common_name,
            human_use,
            harmful_to_human,
            scaled_log_distance_human))  

random <- "(1 | phylum/class/order) + (1 | biogeography)"

colnames(db3)

#formula
model.formula.db3 <- as.formula(paste0("Length ~",
                                       paste(colnames(db3)[16:ncol(db3)], collapse = " + "),
                                       "+",
                                       random))

M1 <- glmmTMB::glmmTMB(model.formula.db3,
                       family = poisson, 
                       data = db3)

# Model validation
performance::check_overdispersion(M1) #Model is overdispersed
performance::check_model(M1)
summary(M1)


#formula
model.formula.db3 <- as.formula(paste0("resid_readability ~",
                                       paste(colnames(db3)[16:ncol(db3)], collapse = " + "),
                                       "+",
                                       random))

M1 <- lme4::lmer(model.formula.db3,
                       data = db3)

# Model validation
parameters::parameters(M1)
performance::check_model(M1)
