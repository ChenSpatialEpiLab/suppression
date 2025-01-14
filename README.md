# suppression
Code for dealing with data suppression in CDC WONDER data

CDC Wonder mortality data suppresses death counts between 1 and 5.

We want code that multiply imputes by drawing from a random uniform distribution between 1 and 5.
