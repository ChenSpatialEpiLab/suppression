# suppression
Code for dealing with data suppression in CDC WONDER data

CDC Wonder mortality data suppresses death counts between 1 and 9.

We want code that multiply imputes by drawing from a random uniform distribution between 1 and 9. We create
multiply imputed datasets that fill in the missing values, and store this in a list object.

We then use purrr::map to run an analysis in each imputed dataset, and store the results in another list object.

We illustrate using an example of computing age-standardized rates by state, year, race/ethnicity, and gender.

Then, we use purrr::map_dfr to use mice::pool.scalar to combine the estimates using Rubin's rules. 
