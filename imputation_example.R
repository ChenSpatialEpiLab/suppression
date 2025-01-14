
# This is an example of how to do imputation for suppressed counts in the CDC Wonder data.
# The example data has mortality by state, year, race/ethnicity, and gender.
# We impute by drawing from integers between 1 and 9 with equal probability, and store the imputed datasets
# in a list object.
# We then run an analysis (in this example, computing age standardized rates), and
# combine the results using Rubin's rules.

# dependencies ------------------------------------------------------------
library(tidyverse)
library(mice) # Needed to combine results using Rubin's rules


# read in CDC Wonder data -------------------------------------------------

df_hispanic <- read.delim("cdcwonder_2015-2020_Hispanic.txt") |>
  janitor::clean_names() |>
  mutate(deaths = as.numeric(deaths),
         population = as.numeric(population),
         race = "Hispanic or Latino") |>
  filter(!state=="")


df_nonhispanic <- read.delim("cdcwonder_2015-2020_Non-Hispanic.txt") |>
  janitor::clean_names() |>
  mutate(deaths = as.numeric(deaths),
         population = as.numeric(population)) |>
  filter(!state=="")

df_combined <- bind_rows(df_nonhispanic, df_hispanic) |>
  select(-hispanic_origin, -hispanic_origin_code) |>
  rename(agecat = ten_year_age_groups) |>
  filter(!agecat == "Not Stated")


# create a list of imputed datasets  -------------------------------

# specify the number of imputations
n_imputations <- 10

# create a list object to hold the imputed datasets
df_imputed <- list()

# generate imputed datasets and store them in the list
for (iter in c(1:n_imputations)) {
  df_imputed[[iter]] <- df_combined |>
      mutate(deaths = replace(deaths, is.na(deaths), sample(1:9, sum(is.na(deaths)), replace = TRUE)))
}


# as an example: implement direct age standardization ---------------------

# In this example, we are using 11 age categories
these_agecat <- c('< 1 year',
                  paste0(c('1-4', '5-14','15-24','25-34',
                               '35-44','45-54','55-64','65-74','75-84','85+'), ' years'))

# read SEER standard population data. The 19 age categories need to be collapsed
# into 13 age categories to match what is in the mortality data.
seer_age19 <- read_fwf("stdpop.19ages.txt",
                       fwf_widths(c(3,3,8),
                                  c("standard","age","std_raw"))) %>%
  # In this example, we will use standard "201", which corresponds to
  # 2000 U.S. Std Million (19 age groups) according to the documentation
  filter(standard=="201") |>
  # In this example, we need to collapse some of the categories to match what we have available in the mortality data
  mutate(agecat=recode(age,
                       '000'="< 1 year",                   
                       '001'="1-4 years",
                       '002'="5-14 years",
                       '003'="5-14 years",
                       '004'="15-24 years",
                       '005'="15-24 years",
                       '006'="25-34 years",
                       '007'="25-34 years",
                       '008'="35-44 years",
                       '009'="35-44 years",
                       '010'="45-54 years",
                       '011'="45-54 years",
                       '012'="55-64 years",
                       '013'="55-64 years",
                       '014'="65-74 years",
                       '015'="65-74 years",
                       '016'="75-84 years",
                       '017'="75-84 years",
                       '018'="85+ years"),
         std.pop=as.numeric(std_raw)) |>
  group_by(agecat) |>
  summarise(std=sum(std.pop)) |>
  mutate(agecat = factor(agecat, levels = these_agecat))


# write a function to perform age standardization

f_age_standardize <- function(data, standard) {
  # We have to merge the dataset with the age standard, by age category
  df_std <- data |>
    # merge with the SEER age standard
    left_join(standard, by = "agecat") |>
    # calculate stratum-specific rates
    # and the variance of the stratum-specific rates
    mutate(rate = deaths / population,
           var_rate = deaths / (population^2)) |>
    # to aggregate over age strata, we need to group
    group_by(state, year, race, gender) |>
    # summarize by taking weighted sums, where the weights are the 
    # weights from the SEER standard
    # Don't forget to divide by the sum of the weights
    dplyr::summarize(age_std_rate = sum(std * rate)/sum(std),
                     # Note that the variance of the age-standardized rate is itself 
                     # a weighted sum, where the weights  are 
                     var_age_std_rate = sum(std^2 * var_rate)/(sum(std)^2)) |>
    ungroup() |>
    # create a row index variable which we will use later for matching to the pooled results
    mutate(row_index = row_number())
}



# use purrr::map to run the analysis in each element of the list ----------

list_agestd_rates <- purrr::map(
  .x = df_imputed,
  .f = ~ f_age_standardize(.x, seer_age19)
)


# use pool.scalar() from the mice package to combine results --------------

# Suppose each of the data frames in df_list
# has columns: estimate, variance
# They are all the same length (3050 rows).

pooled_results <- map_dfr(
  .x = seq_len(nrow(list_agestd_rates[[1]])),   # each element of list_agestd_rates should have the same number of rows
  .f = function(i) {
    # Get the vector of estimates across all list elements for row i
    Q_vec <- map_dbl(list_agestd_rates, ~ .x$age_std_rate[i])
    # Get the vector of variances across all list elements for row i
    U_vec <- map_dbl(list_agestd_rates, ~ .x$var_age_std_rate[i])
    
    # Apply Rubin's rules for univariate estimates
    # n can be set to the original sample size or NULL if not known
    pooled <- pool.scalar(Q = Q_vec, U = U_vec, n = NULL)
    
    # Create a single-row tibble for the results
    tibble(
      row_index  = i,
      Qbar       = pooled$qbar,  # pooled estimate
      T          = pooled$t      # pooled total variance
    )
  }
)

# put back the identifying variables for each row
pooled_results <- pooled_results |>
  left_join(list_agestd_rates[[1]] |> select(row_index, state, year, race, gender),
            by = "row_index")


# Note that Qbar is the point estimate and T is the pooled variance
# The analytic code will need to be amended for other kinds of analyses. This is just an illustration.

