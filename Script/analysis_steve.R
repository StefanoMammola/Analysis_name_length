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
               "performance")

# Settings ----------------------------------------------------------------

set.seed(42)

theme_set(theme_bw())

theme_update(
  legend.position = "right", #No legend
  plot.background = element_blank(), #No background
  panel.grid = element_blank(), #No gridlines
  axis.text.x = element_text(size = 10, colour = "black"),
  axis.title = element_text(size = 14, colour = "black")#Font 10 for x axis
)

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
                 dbwiki = db2$total_wiki_pgviews)

# Data exploration --------------------------------------------------------

GGally::ggpairs(data = results)

# Modelling ---------------------------------------------------------------

formula_m1 <- as.formula("citations ~ year + Length + (1 | phylum / class / order)")

m1 <- glmmTMB(formula_m1, data = db, 
              family = poisson,
              control=glmmTMBControl(optimizer=optim,
                                     optArgs=list(method="BFGS")))
performance::check_overdispersion(m1)

m1 <- glmmTMB(formula_m1, data = db, 
              family = nbinom1,
              control=glmmTMBControl(optimizer=optim,
                                     optArgs=list(method="BFGS")))

summary(m1)


# Plot 1: Length vs log(citations + 1)
plot1 <- ggplot(db, aes(x = Length, y = citations)) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # Regression line
  labs(title = "Length vs log(Citations + 1)",
       x = "Length",
       y = "log(Citations + 1)")

# Plot 2: Length vs log(dbwiki + 1)
plot2 <- ggplot(results, aes(x = Length, y = log(dbwiki + 1))) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "red") +  # Regression line
  labs(title = "Length vs log(DBWiki + 1)",
       x = "Length",
       y = "log(DBWiki + 1)")

# Plot 3: Readability vs log(citations + 1)
plot3 <- ggplot(results, aes(x = Readability, y = log(citations + 1))) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "green") +  # Regression line
  labs(title = "Readability vs log(Citations + 1)",
       x = "Readability",
       y = "log(Citations + 1)")

# Plot 4: Readability vs log(dbwiki + 1)
plot4 <- ggplot(results, aes(x = Readability, y = log(dbwiki + 1))) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "purple") +  # Regression line
  labs(title = "Readability vs log(DBWiki + 1)",
       x = "Readability",
       y = "log(DBWiki + 1)")

library(patchwork)

# Arrange plots in a 2x2 grid
combined_plots <- (plot1 + plot2) / (plot3 + plot4)
print(combined_plots)


# Apply function and format results

results2 <- results[results$Length < 40,]

# Plot 1: Length vs log(citations + 1)
plot1 <- ggplot(results2, aes(x = Length, y = log(citations + 1))) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # Regression line
  labs(title = "Length vs log(Citations + 1)",
       x = "Length",
       y = "log(Citations + 1)")

# Plot 2: Length vs log(dbwiki + 1)
plot2 <- ggplot(results2, aes(x = Length, y = log(dbwiki + 1))) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "red") +  # Regression line
  labs(title = "Length vs log(DBWiki + 1)",
       x = "Length",
       y = "log(DBWiki + 1)")

# Plot 3: Readability vs log(citations + 1)
plot3 <- ggplot(results2, aes(x = Readability, y = log(citations + 1))) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "green") +  # Regression line
  labs(title = "Readability vs log(Citations + 1)",
       x = "Readability",
       y = "log(Citations + 1)")

# Plot 4: Readability vs log(dbwiki + 1)
plot4 <- ggplot(results2, aes(x = Readability, y = log(dbwiki + 1))) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = TRUE, color = "purple") +  # Regression line
  labs(title = "Readability vs log(DBWiki + 1)",
       x = "Readability",
       y = "log(DBWiki + 1)")

library(patchwork)

# Arrange plots in a 2x2 grid
combined_plots2 <- (plot1 + plot2) / (plot3 + plot4)
print(combined_plots2)
