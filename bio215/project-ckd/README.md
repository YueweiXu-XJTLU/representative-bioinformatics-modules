# Chronic kidney disease prediction

[Submitted report](reports/report.pdf) · [figures](figures/) · [Shiny application](ShinyApp/)

The project combines preprocessing, GLM/random-forest/SVM comparisons, model evaluation, SHAP interpretation, and an interactive Shiny interface.

## Running the application

- The package is maintained separately in [ckdPrediction](https://github.com/YueweiXu-XJTLU/ckdPrediction). Five model/data files in the application were identical to its [package assets](https://github.com/YueweiXu-XJTLU/ckdPrediction/tree/6282df900bce0f7c2d3cf30e65f3d340db1b2d09/inst/extdata), so they are restored from that package. From `ShinyApp/`:

```r
remotes::install_github("YueweiXu-XJTLU/ckdPrediction@6282df900bce0f7c2d3cf30e65f3d340db1b2d09")
assets <- c("glm_fit.rds", "rf_fit.rds", "svm_fit.rds", "train_df.rds", "test_df.rds")
file.copy(file.path(system.file("extdata", package = "ckdPrediction"), assets), ".")
shiny::runApp(".")
```

- The Shiny web application's [source code](app://-/ShinyApp/) is included in this repository, and the [deployed application](https://yueweixu.shinyapps.io/ckdPrediction) is available online.

***This repository contains coursework for educational purposes only and should not be treated as a rigorously validated academic reference.***
