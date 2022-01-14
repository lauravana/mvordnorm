# mvordnorm

R package mvordnorm implements composite likelihood estimation in a multivariate regression model of continuous and ordinal variables.
For the ordinal responses, marginally ordinal probit regression model is employed. A regression model with normally distributed errors
is employed for the continuous responses. The dependence among the multiple responses is modeled by assuming that the normally distributed variables
and the latent variables underlying the ordinal responses come from a multivariate distribution with a general correlation structure.
In addition, regression coefficients (and threshold parameters for the ordinal responses) are varying across the multiple response dimensions.
