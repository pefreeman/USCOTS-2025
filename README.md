# USCOTS-2025
Materials for the Multinomial Simulation Poster at USCOTS 2025

The `R Markdown` files contained in this repo include

1. `Analyze_Die_Tosses.Rmd`: given the results of a die-tossing experiment,
what are the $p$-values that arise when we perform the chi-square goodness
of fit test and conduct multinomial simulations, under the null hypothesis
that the die is fair?

2. `CI_Multinomial.Rmd`: given the results of a set of multinomial simulations
(an estimated $p$-value), can we construct a 95% confidence interval for the 
true $p$-value? (Short answer: yes.)

3. `P_Value_Uniformity.Rmd`: are the $p$-values that arise when we perform 
the chi-square goodness-of-fit test and conduct multinomial simulations, 
under the null hypothesis that the die is fair, uniformly distributed (as
we would expect them to be when the null hypothesis is correct)?

To be clear: none of these materials exist to "convince" you, the instructor,
that multinomial simulations are superior to the chi-square goodness-of-fit
test; rather, they are provided as starter materials for exercises that you
might want to have your students complete. Feel free to download and amend
these files!

Much of the material in these files are also packaged in the `R Shiny` app
`multinomial_app.R`.
