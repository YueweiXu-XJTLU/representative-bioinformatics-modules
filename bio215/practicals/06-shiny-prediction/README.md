# Shiny prediction practical

[Notebook](Practical_6.Rmd) · [saved HTML](Practical_6.html)

`app.R` is the application also present in the original deployment folder. `app_by_hand.R` and `app_only_practical_6.R` are distinct development variants. The test input and `img_ai.png` are retained; the teaching dataset is shared in `../../data/`.

The fitted random forest and example CSV are identical to the assets in [m6APrediction](https://github.com/YueweiXu-XJTLU/m6APrediction/tree/23fa952f1f85604d6436532d4afe7c31095ac1bb/inst/extdata). From this directory, install the package and restore its assets before running the app:

```r
remotes::install_github("YueweiXu-XJTLU/m6APrediction@23fa952f1f85604d6436532d4afe7c31095ac1bb")
assets <- c("rf_fit.rds", "m6A_input_example.csv")
file.copy(file.path(system.file("extdata", package = "m6APrediction"), assets), ".")
shiny::runApp(".")
```

The restored assets are excluded from this repository by `.gitignore`. Application dependencies are named in the source files; the applications were not launched during this reorganization.

***This repository contains coursework for educational purposes only and should not be treated as a rigorously validated academic reference.***
