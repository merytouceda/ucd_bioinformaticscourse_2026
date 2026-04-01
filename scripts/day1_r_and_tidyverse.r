
# installing and loading packages
install.packages("tidyverse") #you only have to do this once

library(tidyverse)

# loading prepackaged data
# note that after running the following line you can see a Data object (a dataframe specifically) on the rightside pannel --> 
# click on it to view it
data("starwars") 

# ------------------------------------------------------------------------------
# 1. Just Enough Base R to Survive
# ------------------------------------------------------------------------------
# Create an object and assign it a value
x <- 5

# better naming of variables
my_favorite_number <- 5
greeting <- "Hello world!"

# Examine variable content
# this will make the contents of the variables be printed to the "Console", which is below 
# Note that you can also look to the rightside pannel and see the objects under "Values" --> 
print(my_favorite_number)
print(greeting)


# There are several data types in R
# numeric: numbers (integers and decimals)
my_age <- 30
pi_approx <- 3.14159

# character: text, always wrapped in quotes
my_name <- "Obi-Wan Kenobi"
my_species <- "Human"

# logical: TRUE or FALSE (always uppercase in R)
is_jedi <- TRUE
is_sith <- FALSE

# you can check what type an object is with class()
class(my_age)        # "numeric"
class(my_name)       # "character"
class(is_jedi)       # "logical"

# why does this matter? some functions only work on certain types
# for example, you can do math on numeric but not character
my_age + 10          # works fine
# my_name + 10       # would give an error — uncomment to see it



# ------------------------------------------------------------------------------
# 2. Vectors: the fundamental unit of R
# ------------------------------------------------------------------------------
# a vector is a sequence of values of the same type
# think of it as a single column in a spreadsheet
heights <- c(172, 167, 96, 202, 150)
names_vec <- c("Luke", "Leia", "R2-D2", "Darth Vader", "Yoda")
is_human <- c(TRUE, TRUE, FALSE, TRUE, FALSE)

# vectors can be named
names(heights) <- names_vec
heights

# you can do math on numeric vectors — it applies to every element
heights / 100          # convert cm to meters
heights - mean(heights) # how far from average height?

# check the length of a vector
length(heights)



# ------------------------------------------------------------------------------
# 3. Objects vs Functions
# ------------------------------------------------------------------------------
# Everything that exists is an object. Everything that happens is a function.

# objects just sit there holding a value
my_number <- 42
my_vector <- c(1, 2, 3, 4, 5)

# functions DO something — they take an input and give something back
# the parentheses are the giveaway: if it has parentheses, it's a function
sqrt(my_number)         # give it a number, get a number back
sum(my_vector)          # give it a vector, get a number back
length(my_vector)       # give it a vector, get its length back
toupper(greeting)       # give it a character, get it in uppercase

# you can save what a function gives back into a new object
square_root <- sqrt(my_number)
square_root

# this is exactly what happens when you load data:
# read_csv() is a function — you hand it a file path, it gives back a data frame
# starwars is an object — a data frame sitting in your environment


# ------------------------------------------------------------------------------
# 4. A first look at a data frame
# ------------------------------------------------------------------------------
# a data frame is a table: rows are observations, columns are variables
# each column is a vector — that's why data types matter

# how many rows and columns?
nrow(starwars)
ncol(starwars)
dim(starwars)

# get a quick overview — note the data type listed under each column name
glimpse(starwars)

# summary statistics for every column
summary(starwars)

# look at the first few rows
head(starwars)

# access a single column with $ — this gives you a vector
starwars$name
starwars$height

# basic stats on a column
mean(starwars$height, na.rm = TRUE)   # na.rm = TRUE tells R to ignore missing values
max(starwars$height, na.rm = TRUE)
min(starwars$height, na.rm = TRUE)

# na.rm is worth explaining: some characters have no recorded height (feel free to explore the data with View)
# R is strict — if there are NAs and you don't tell it what to do, it returns NA
mean(starwars$height)           # returns NA
mean(starwars$height, na.rm = TRUE)  # returns the mean of non-missing values


# ------------------------------------------------------------------------------
# 5. The dplyr Core 5
# ------------------------------------------------------------------------------
# dplyr is part of the tidyverse and gives us a set of verbs to work with data frames
# each verb does one thing, takes a data frame, and gives a data frame back
# remember: a data frame goes in, a data frame comes out

# first, a note on the pipe operator %>%
# the pipe takes what is on the left and passes it to the function on the right
# think of it as "and then"
# instead of writing:
summary(starwars)
# you can write:
starwars %>% summary()
# they do the same thing — but the pipe becomes powerful when you chain steps together

# the keyboard shortcut for %>% is Ctrl + Shift + M (Windows) or Cmd + Shift + M (Mac)
# go to Tools -> Global Options -> Code and check "Use native pipe operator" 
# if you want |> instead, but we will use %>% throughout this course


# ----- filter(): keep rows that match a condition -----

# keep only human characters
humans <- starwars %>% 
  filter(species == "Human")

humans
nrow(humans)

# keep only characters taller than 180cm
tall_characters <- starwars %>% 
  filter(height > 180)

# filter with multiple conditions
# , between conditions means AND — both must be true
tall_humans <- starwars %>% 
  filter(species == "Human", height > 180)

# | means OR — at least one must be true
humans_or_droids <- starwars %>% 
  filter(species == "Human" | species == "Droid")

# useful shortcut for multiple values of the same variable
humans_or_droids <- starwars %>% 
  filter(species %in% c("Human", "Droid"))

# filter out missing values
# is.na() checks if a value is missing
# ! means NOT — so !is.na() means "is not missing"
no_missing_height <- starwars %>% 
  filter(!is.na(height))


# ----- select(): keep columns you want -----

# keep only name, height, and species
starwars %>% 
  select(name, height, species)

# remove a column with -
starwars %>% 
  select(-films, -vehicles, -starships)

# select a range of columns
starwars %>% 
  select(name:eye_color)

# select columns that start with a pattern
starwars %>% 
  select(starts_with("s"))  # skin_color, species, starships, sex

# combine filter and select with the pipe
# "give me the name and height of all humans"
starwars %>% 
  filter(species == "Human") %>% 
  select(name, height)


# ----- mutate(): create or modify a column -----

# convert height from cm to meters
starwars %>% 
  mutate(height_m = height / 100) %>% 
  select(name, height, height_m)

# create a BMI column
starwars %>% 
  mutate(
    height_m = height / 100,
    bmi = mass / height_m^2
  ) %>% 
  select(name, height_m, mass, bmi)

# mutate with a condition using if_else()
# if_else(condition, value if TRUE, value if FALSE)
starwars %>% 
  mutate(is_tall = if_else(height > 180, "tall", "not tall")) %>% 
  select(name, height, is_tall)


# if_else() works great for two outcomes (tall / not tall)
# but what if you have more than two categories?
# that's where case_when() comes in — it's like if_else() but for multiple conditions
# think of it as: "when THIS is true, give me THAT"
starwars %>% 
  mutate(height_category = case_when(
    height < 100              ~ "short",
    height >= 100 & height < 180 ~ "average",
    height >= 180             ~ "tall",
    is.na(height)             ~ "unknown"
  )) %>% 
  select(name, height, height_category)


# ----- arrange(): sort rows -----

# sort by height, shortest first (ascending by default)
starwars %>% 
  arrange(height) %>% 
  select(name, height)

# sort descending — tallest first
starwars %>% 
  arrange(desc(height)) %>% 
  select(name, height)

# sort by multiple columns
starwars %>% 
  arrange(species, desc(height)) %>% 
  select(name, species, height)


# ----- summarize(): collapse rows to a single summary -----

# average height across all characters
starwars %>% 
  summarize(mean_height = mean(height, na.rm = TRUE))

# multiple summaries at once
starwars %>% 
  summarize(
    mean_height = mean(height, na.rm = TRUE),
    median_height = median(height, na.rm = TRUE),
    min_height = min(height, na.rm = TRUE),
    max_height = max(height, na.rm = TRUE),
    n = n()   # n() counts the number of rows
  )


# ------------------------------------------------------------------------------
# 6. Group Operations
# ------------------------------------------------------------------------------
# group_by() + summarize() is where dplyr gets really powerful
# it splits the data frame into groups, applies summarize to each, and combines the results

# average height by species
starwars %>% 
  group_by(species) %>% 
  summarize(
    mean_height = mean(height, na.rm = TRUE),
    n = n()
  )


# filter to species with more than 2 characters first
# then summarize — chaining multiple steps
starwars %>% 
  filter(!is.na(species)) %>% 
  group_by(species) %>% 
  summarize(
    mean_height = mean(height, na.rm = TRUE),
    mean_mass = mean(mass, na.rm = TRUE),
    n = n()
  ) %>% 
  filter(n > 2) %>%          # only species with more than 2 characters
  arrange(desc(mean_height))


# count() is a shortcut for group_by() + summarize(n = n())
# how many characters per species?
starwars %>% 
  count(species, sort = TRUE)   # sort = TRUE orders by n descending

# count by two variables
starwars %>% 
  count(species, gender, sort = TRUE)




# ------------------------------------------------------------------------------
# 7. First ggplot
# ------------------------------------------------------------------------------
# ggplot2 builds plots in layers
# the grammar: ggplot(data, aes(x, y)) + geom_*()
# aes() = aesthetics: which columns map to which visual properties
# geom_*() = geometry: what type of plot

# ----- bar chart from a summarize result -----

# first build the summary
species_counts <- starwars %>% 
  count(species, sort = TRUE) %>% 
  filter(n > 1)                 # only species with more than 1 character

# then plot it
ggplot(species_counts, aes(x = species, y = n)) +
  geom_col()

# hard to read — flip the axes
ggplot(species_counts, aes(x = n, y = species)) +
  geom_col()

# reorder bars by count — much more readable
ggplot(species_counts, aes(x = n, y = reorder(species, n))) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Number of characters per species",
    x = "Count",
    y = "Species"
  )

# ----- scatter plot -----

# height vs mass
ggplot(starwars, aes(x = height, y = mass)) +
  geom_point()

# there is an outlier — who is it?
starwars %>% 
  filter(mass > 1000) %>% 
  select(name, height, mass, species)
# Jabba the Hutt — worth removing for a cleaner plot

starwars %>% 
  filter(mass < 1000) %>%        # pipe directly into ggplot
  ggplot(aes(x = height, y = mass)) +
  geom_point()

# map a third variable to color
starwars %>% 
  filter(mass < 1000, !is.na(gender)) %>% 
  ggplot(aes(x = height, y = mass, color = gender)) +
  geom_point(size = 3, alpha = 0.7) +   # size = point size, alpha = transparency
  labs(
    title = "Height vs mass in Star Wars characters",
    x = "Height (cm)",
    y = "Mass (kg)",
    color = "Gender"
  )

# add a trend line
starwars %>% 
  filter(mass < 1000, !is.na(gender)) %>% 
  ggplot(aes(x = height, y = mass, color = gender)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +  # lm = linear model, se = confidence band
  labs(
    title = "Height vs mass in Star Wars characters",
    x = "Height (cm)",
    y = "Mass (kg)",
    color = "Gender"
  )




# ----- a note on plot types and geom_() -----
# the geom you choose depends on what you want to show:
# geom_col() / geom_bar()  — counts or values per category
# geom_point()             — relationship between two numeric variables
# geom_line()              — trends over time or ordered sequences
# geom_boxplot()           — distribution of a numeric variable per group
# geom_histogram()         — distribution of a single numeric variable
# geom_violin()            — like boxplot but shows full distribution shape

# quick examples:

# boxplot: height distribution by gender
starwars %>% 
  filter(!is.na(gender)) %>% 
  ggplot(aes(x = gender, y = height, fill = gender)) +
  geom_boxplot() +
  labs(title = "Height distribution by gender")

# histogram: distribution of heights
starwars %>% 
  filter(!is.na(height)) %>% 
  ggplot(aes(x = height)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  labs(title = "Distribution of character heights")

# ggplot always wants a tidy data frame
# this is exactly what dplyr produces — that's why the two fit together perfectly
# a typical pattern you will use constantly:
starwars %>%                              # start with data
  filter(!is.na(height)) %>%             # clean it
  group_by(species) %>%                  # group
  summarize(mean_height = mean(height),  # summarize
            n = n()) %>% 
  filter(n > 2) %>%                      # filter the summary
  ggplot(aes(x = reorder(species, mean_height), # plot
             y = mean_height)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Mean height by species (species with >2 characters)",
    x = "Species",
    y = "Mean height (cm)"
  )

