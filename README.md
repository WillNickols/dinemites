# dinemites

To compile package
```
roxygen2::roxygenize()
```

To install package:
```
conda env create -f dinemites.yml
```

```
install.packages('devtools')
devtools::install_github('WillNickols/dinemites', type = 'source')
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```
