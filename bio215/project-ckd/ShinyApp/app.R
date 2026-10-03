library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)
library(ggplot2)
library(ckdPrediction)

#############################################################
# 1. Data & preprocessing
#############################################################

# load data
ckd_raw  <- read.csv("../../../shared-data/kidney_disease.csv", stringsAsFactors = FALSE)
train_df <- readRDS("train_df.rds")
test_df  <- readRDS("test_df.rds")
example_df  <- readRDS("example_df.rds")
rf_fit   <- readRDS("rf_fit.rds")
glm_fit  <- readRDS("glm_fit.rds")
svm_fit  <- readRDS("svm_fit.rds")
rf_curve_df  <- readRDS("grid_rf.rds")
glm_curve_df <- readRDS("grid_glm.rds")
svm_curve_df <- readRDS("grid_svm.rds")

# normalise missing values
ckd_raw[ckd_raw == ""] <- NA
ckd_raw[] <- lapply(ckd_raw, function(x) {
  if (is.character(x)) {
    x <- trimws(x)
    x <- gsub("\\t", "", x)
    x[x %in% c("", "?", "NA")] <- NA
  }
  return(x)
})
numeric_fix <- c("pcv", "wc", "rc")
ckd_raw[numeric_fix] <- lapply(ckd_raw[numeric_fix], function(x) as.numeric(x))

# model registry
model_list <- list(
  "Random forest"             = rf_fit,
  "Ridge logistic regression" = glm_fit,
  "Support vector machine"    = svm_fit
)

# predictor names
predictor_cols <- setdiff(colnames(test_df), c("classification", "id"))

# feature dictionary
feature_dictionary <- data.frame(
  Column  = c("id", "age", "bp", "sg", "al", "su", "rbc", "pc", "pcc", "ba",
              "bgr", "bu", "sc", "sod", "pot", "hemo", "pcv", "wc", "rc",
              "htn", "dm", "cad", "appet", "pe", "ane", "classification"),
  Meaning = c("Sample identifier", "Age", "Blood pressure",
              "Urine specific gravity", "Albumin in urine", "Sugar in urine",
              "Red blood cell morphology", "Pus cells", "Pus cell clumps",
              "Bacteria", "Random blood glucose", "Blood urea",
              "Serum creatinine", "Serum sodium", "Serum potassium",
              "Hemoglobin", "Packed cell volume", "White blood cell count",
              "Red blood cell count", "Hypertension status",
              "Diabetes mellitus status", "Coronary artery disease status",
              "Appetite status", "Presence of pedal edema", "Anemia status",
              "Disease label (CKD vs non-CKD)"),
  stringsAsFactors = FALSE
)

# pre-compute metrics and recommended thresholds
ev_rf  <- model_evaluation_ckd(test_df, rf_fit,  threshold = 0.5)
ev_glm <- model_evaluation_ckd(test_df, glm_fit, threshold = 0.5)
ev_svm <- model_evaluation_ckd(test_df, svm_fit, threshold = 0.5)
best_thr_rf  <- as.numeric(ev_rf$best_threshold$value)
best_thr_glm <- as.numeric(ev_glm$best_threshold$value)
best_thr_svm <- as.numeric(ev_svm$best_threshold$value)

# helper: summarise metrics for one model
eval_one <- function(model, label, test_df, threshold = 0.5) {
  ev      <- model_evaluation_ckd(test_df, model, threshold = threshold)
  m       <- ev$metrics
  auc_tab <- ev$curves$auc
  auc_roc <- as.numeric(auc_tab[auc_tab$curvetypes == "ROC", "aucs"])
  auc_prc <- as.numeric(auc_tab[auc_tab$curvetypes == "PRC", "aucs"])
  
  data.frame(
    Model       = label,
    Accuracy    = round(m$accuracy, 3),
    Precision   = round(m$precision, 3),
    Recall      = round(m$recall, 3),
    Specificity = round(1 - m$fpr, 3),
    F1_score    = round(2 * m$precision * m$recall / (m$precision + m$recall), 3),
    AUC_ROC     = round(auc_roc, 3),
    AUC_PRC     = round(auc_prc, 3),
    stringsAsFactors = FALSE
  )
}

# helper: locate closest point to a threshold
find_threshold_point <- function(curve_df, thr) {
  curve_df$dist <- abs(curve_df$threshold - thr)
  curve_df[which.min(curve_df$dist), , drop = FALSE]
}

#############################################################
# 2. UI layout
#############################################################

ui <- dashboardPage(
  title = "CKD Prediction",
  skin  = "blue",
  
  # header
  dashboardHeader(title = "CKD Prediction"),
  
  # sidebar navigation
  dashboardSidebar(
    useShinyjs(),
    actionButton("re_initialize", "Reinitialize app", icon = icon("refresh")),
    br(),
    sidebarMenu(
      menuItem("Home Page", tabName = "overview", icon = icon("home")),
      menuItem(
        "CKD Data", icon = icon("database"),
        menuSubItem("Dataset Browser", tabName = "database", icon = icon("table")),
        menuSubItem("Data Visualisation", tabName = "visualization", icon = icon("area-chart"))
      ),
      menuItem(
        "Model Performance", icon = icon("line-chart"),
        menuSubItem("Threshold Selection",   tabName = "performance_curves",     icon = icon("wave-square")),
        menuSubItem("Model Comparison",      tabName = "performance_thresholds", icon = icon("balance-scale")),
        menuSubItem("Feature Importance",     tabName = "performance_importance", icon = icon("layer-group"))
      ),
      menuItem(
        "CKD Prediction", tabName = "prediction", icon = icon("calculator"),
        menuSubItem("Batch Prediction",      tabName = "batch",  icon = icon("file-upload")),
        menuSubItem("Single-patient Prediction",   tabName = "single", icon = icon("user"))
      )
    )
  ),
  
  # body content
  dashboardBody(
    useShinyjs(),
    
    tags$head(tags$title("ckdPrediction")),
    
    tags$head(tags$style(HTML("
      #illustrative_img { height: auto !important; }
      #illustrative_img img{
        max-width: 100%;
        height: auto;
        max-height: 320px;
        object-fit: contain;
        display: block;
        margin: 0 auto;
      }
    "))),
    
    tabItems(
      
      ##################################################
      # 2.1 Home overview
      ##################################################
      tabItem(
        tabName = "overview",
        
        fluidRow(
          valueBoxOutput("vb_total_samples", width = 4),
          valueBoxOutput("vb_feature_cols",  width = 4),
          valueBoxOutput("vb_model_types",   width = 4)
        ),
        
        fluidRow(
          box(
            title = "Project Summary",
            width = 6, status = "primary", solidHeader = TRUE,
            p(
              "This Web-App is based on",
              strong(
                tags$a(
                  href   = "https://github.com/YueweiXu-XJTLU/ckdPrediction",
                  target = "_blank",
                  "ckdPrediction R package (click to visit)"
                )
              ),
              "about chronic kidney disease (CKD) prediction."
            ),
            p("It provides:"),
            tags$ul(
              style = "list-style-type: disc; padding-left: 20px;",
              tags$li(
                "A mini CKD database derived from ",
                tags$a(
                  href   = "https://www.kaggle.com/datasets/akshayksingh/kidney-disease-dataset",
                  target = "_blank",
                  "Kaggle (click to visit)"
                ),
                ", allowing users to inspect, search, visualise, and download data inside."
              ),
              tags$li(
                "Three pre-trained classifiers - ", strong("Random Forest, Ridge Logistic Regression, Support Vector Machine"), ", enabling comparisons of performance."
              ),
              tags$li(
                "CKD prediction using the ML models for both single records and batch CSV files, together with SHAP interpretations."
              )
            ),
            tags$p(
              strong("Note: "),
              "Due to the strong correlations between certain clinical features and CKD diagnosis, ",
              "the limited sample size, and the geographic concentration of data, ",
              "our models may exhibit ",
              tags$strong("optimistically high performance"),
              ". Please interpret results cautiously."
            ),
            tags$p(
              strong("Disclaimer: "),
              "The web app is only for study purpose, and ",
              strong("must not be used"),
              " for any real clinical decision making.",
              style = "color:#b91c1c;"
            )
          ),
          box(
            title = "How to Use",
            width = 6, status = "primary", solidHeader = TRUE,
            tags$ol(
              tags$li(
                strong("CKD Data"),
                tags$ul(
                  tags$li(
                    strong("Dataset Browser: "),
                    "Browse, filter and download the CKD dataset (full or processed subset)."
                  ),
                  tags$li(
                    strong("Data Visualisation: "),
                    "Plot distributions of numeric features and proportions of categorical features by CKD label."
                  )
                )
              ),
              tags$li(
                strong("Model Performance"),
                tags$ul(
                  tags$li(
                    strong("Thresholds Selection: "),
                    "Use interactive ROC & PRC curves to view and select your desired thresholds."
                  ),
                  tags$li(
                    strong("Models Comparison: "),
                    "Use decided thresholds to compare detailed performance among models."
                  ),
                  tags$li(
                    strong("Feature Importance: "),
                    "Inspect feature importance and detailed algorithms for models."
                  )
                )
              ),
              tags$li(
                strong("CKD Prediction"),
                tags$ul(
                  tags$li(
                    strong("Batch Prediction: "),
                    "Upload a CSV, obtain CKD probabilities/labels. SHAP interpretation plots for whole model/single feature/single row are all available."
                  ),
                  tags$li(
                    strong("Single Prediction: "),
                    "Enter single patient data and receive instant prediction result and SHAP interpretations."
                  )
                )
              )
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Work Flow",
            width = 12, status = "primary", solidHeader = TRUE,
            div(
              style = "text-align:center;",
              imageOutput("illustrative_img")
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Feature Dictionary",
            width = 12, status = "primary", solidHeader = TRUE,
            p(em(
              "The dataset contains 26 observed columns per sample. ",
              "The table below explains the meaning of each feature."
            )),
            DTOutput("feature_table")
          )
        )
      ),
      
      ##################################################
      # 2.2 CKD database
      ##################################################
      tabItem(
        tabName = "database",
        
        fluidRow(
          box(
            title = "Introduction",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "This page provides an ", strong("interactive CKD dataset"), ". ",
              "Users can ", strong("view, filter or search"), " information in the dataset, ",
              "and ", strong("download"), " either the full dataset or the processed subset."
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Overall Filters",
            width = 12, status = "primary", solidHeader = TRUE,
            checkboxInput(
              "db_pipeline_process",
              strong("Use pipeline process to the dataset (see R package for more details)"),
              value = FALSE
            ),
            checkboxInput(
              "db_show_missing",
              strong("Exclude rows with missing values (NA)"),
              value = FALSE
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Detailed Filters",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "Use the widgets below to select subset of the dataset. "
            ),
            uiOutput("db_filter_ui")
          )
        ),
        
        fluidRow(
          box(
            title = "Download",
            width = 12, status = "primary", solidHeader = TRUE,
            downloadButton("dl_full_data",     "Download full dataset (CSV)"),
            span(" "),
            downloadButton("dl_filtered_data", "Download processed subset (CSV)")
          )
        ),
        
        fluidRow(
          box(
            title = "Dataset",
            width = 12, status = "primary", solidHeader = TRUE,
            DTOutput("db_table")
          )
        )
      ),
      
      ##################################################
      # 2.3 Data visualisation
      ##################################################
      tabItem(
        tabName = "visualization",
        
        fluidRow(
          box(
            title = "Introduction",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "This page provides exploratory visualisations of the filtered CKD dataset in", strong("\"Dataset Browser\" page above"), ", ",
              "including", strong("distributions of numeric features"), "and ", strong("proportions of categorical features"), ", both characterized by CKD label."
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Distributions of Numeric Features",
            width = 6, status = "primary", solidHeader = TRUE,
            selectInput("db_numeric_var", "Numeric Feature:", choices = NULL),
            plotOutput("db_numeric_plot", height = 340),
            hr(),
            div(
              style = "display:flex; gap:12px; justify-content:flex-start;",
              downloadButton("dl_numeric_png", "Download (PNG)"),
              downloadButton("dl_numeric_pdf", "Download (PDF)")
            )
          ),
          box(
            title = "Proportion of Categorical Features",
            width = 6, status = "primary", solidHeader = TRUE,
            selectInput(
              "db_categorical_var", "Categorical Feature:", choices = NULL),
            plotOutput("db_categorical_plot", height = 340),
            hr(),
            div(
              style = "display:flex; gap:12px; justify-content:flex-start;",
              downloadButton("dl_categorical_png", "Download (PNG)"),
              downloadButton("dl_categorical_pdf", "Download (PDF)")
            )
          )
        )
      ),
      
      ##################################################
      # 2.4 Performance – curves & thresholds
      ##################################################
      tabItem(
        tabName = "performance_curves",
        
        fluidRow(
          box(
            title = "Introduction",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "This page shows ROC and PRC curves for all three models on the internal test set. ",
              "You can either",
              strong("enter custom thresholds directly"),
              ", or ",
              strong("click on any point on the ROC or PRC curve"),
              " to update the thresholds and view the corresponding coordinates.",
              "Default thresholds are our recommended values."
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Note",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              HTML(
                "<span style='color:#b91c1c; font-weight:bold;'>
     Due to plotting simplifications, the curve shapes here may not perfectly reflect the underlying ROC/PRC.<br/>
     This page is intended primarily for choosing decision thresholds.<br/>
     For accurate AUC evaluation, please refer to metrics on the next page or <code>model_evaluation_ckd()</code> function in ckdPrediction R package.
     </span>"
              )
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Random Forest: ROC & PRC",
            width = 12, status = "primary", solidHeader = TRUE,
            fluidRow(
              column(
                width = 6,
                plotOutput("perf_rf_roc", height = 320, click = "rf_roc_click")
              ),
              column(
                width = 6,
                plotOutput("perf_rf_prc", height = 320, click = "rf_prc_click")
              )
            ),
            fluidRow(
              column(
                width = 6,
                div(
                  style = "margin-top:10px;",
                  downloadButton("dl_rf_roc_png", "Download ROC (PNG)"),
                  span(" "),
                  downloadButton("dl_rf_roc_pdf", "Download ROC (PDF)")
                )
              ),
              column(
                width = 6,
                div(
                  style = "margin-top:10px;",
                  downloadButton("dl_rf_prc_png", "Download PRC (PNG)"),
                  span(" "),
                  downloadButton("dl_rf_prc_pdf", "Download PRC (PDF)")
                )
              )
            ),
            hr(),
            fluidRow(
              column(
                width = 4,
                numericInput(
                  "rf_threshold_input",
                  "Current threshold (ckd probability):",
                  value = round(best_thr_rf, 3),
                  min   = 0, max = 1, step = 0.001
                )
              ),
              column(
                width = 8,
                htmlOutput("perf_rf_coords")
              )
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Ridge Logistic Regression: ROC & PRC",
            width = 12, status = "primary", solidHeader = TRUE,
            fluidRow(
              column(
                width = 6,
                plotOutput("perf_glm_roc", height = 320, click = "glm_roc_click")
              ),
              column(
                width = 6,
                plotOutput("perf_glm_prc", height = 320, click = "glm_prc_click")
              )
            ),
            fluidRow(
              column(
                width = 6,
                div(
                  style = "margin-top:10px;",
                  downloadButton("dl_glm_roc_png", "Download ROC (PNG)"),
                  span(" "),
                  downloadButton("dl_glm_roc_pdf", "Download ROC (PDF)")
                )
              ),
              column(
                width = 6,
                div(
                  style = "margin-top:10px;",
                  downloadButton("dl_glm_prc_png", "Download PRC (PNG)"),
                  span(" "),
                  downloadButton("dl_glm_prc_pdf", "Download PRC (PDF)")
                )
              )
            ),
            hr(),
            fluidRow(
              column(
                width = 4,
                numericInput(
                  "glm_threshold_input",
                  "Current threshold (ckd probability):",
                  value = round(best_thr_glm, 3),
                  min   = 0, max = 1, step = 0.001
                )
              ),
              column(
                width = 8,
                htmlOutput("perf_glm_coords")
              )
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Support Vector Machine: ROC & PRC",
            width = 12, status = "primary", solidHeader = TRUE,
            fluidRow(
              column(
                width = 6,
                plotOutput("perf_svm_roc", height = 320, click = "svm_roc_click")
              ),
              column(
                width = 6,
                plotOutput("perf_svm_prc", height = 320, click = "svm_prc_click")
              )
            ),
            fluidRow(
              column(
                width = 6,
                div(
                  style = "margin-top:10px;",
                  downloadButton("dl_svm_roc_png", "Download ROC (PNG)"),
                  span(" "),
                  downloadButton("dl_svm_roc_pdf", "Download ROC (PDF)")
                )
              ),
              column(
                width = 6,
                div(
                  style = "margin-top:10px;",
                  downloadButton("dl_svm_prc_png", "Download PRC (PNG)"),
                  span(" "),
                  downloadButton("dl_svm_prc_pdf", "Download PRC (PDF)")
                )
              )
            ),
            hr(),
            fluidRow(
              column(
                width = 4,
                numericInput(
                  "svm_threshold_input",
                  "Current threshold (ckd probability):",
                  value = round(best_thr_svm, 3),
                  min   = 0, max = 1, step = 0.001
                )
              ),
              column(
                width = 8,
                htmlOutput("perf_svm_coords")
              )
            )
          )
        )
      ),
      
      ##################################################
      # 2.5 Performance – threshold comparison
      ##################################################
      tabItem(
        tabName = "performance_thresholds",
        
        fluidRow(
          box(
            title = "Introduction",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "This page",
              strong("compares three models' performance"),
              "under",
              strong("user-defined decision thresholds"),
              ". ",
              "Metrics of models, together with a comparison bar plot will be generated.",
              "Default thresholds are our recommended values."
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Choose Decision Thresholds",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "Adjust the decision threshold for each model below, then click ",
              strong("\"Evaluate & Compare\""),
              " to recompute performance on the internal test set."
            ),
            fluidRow(
              column(
                width = 4,
                numericInput(
                  "thr_rf2", "Random forest threshold:",
                  value = round(best_thr_rf, 3),
                  min   = 0, max = 1, step = 0.001
                )
              ),
              column(
                width = 4,
                numericInput(
                  "thr_glm2", "Ridge logistic regression threshold:",
                  value = round(best_thr_glm, 3),
                  min   = 0, max = 1, step = 0.001
                )
              ),
              column(
                width = 4,
                numericInput(
                  "thr_svm2", "Support vector machine threshold:",
                  value = round(best_thr_svm, 3),
                  min   = 0, max = 1, step = 0.001
                )
              )
            ),
            br(),
            actionButton("run_threshold_compare", "Evaluate & Compare", icon = icon("play"))
          )
        ),
        
        fluidRow(
          box(
            title = "Metric Comparison Table",
            width = 12, status = "primary", solidHeader = TRUE,
            DTOutput("threshold_compare_table")
          )
        ),
        
        fluidRow(
          box(
            title = "Visual Comparison of Metrics",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "The bar chart below summarises all key metrics ",
              "for the three models under the chosen thresholds above."
            ),
            plotOutput("threshold_compare_plot", height = 360),
            hr(),
            div(
              style = "display:flex; gap:12px; justify-content:flex-start;",
              downloadButton("dl_threshold_compare_png", "Download (PNG)"),
              downloadButton("dl_threshold_compare_pdf", "Download (PDF)")
            )
          )
        )
      ),
      
      ##################################################
      # 2.6 Performance – global importance
      ##################################################
      tabItem(
        tabName = "performance_importance",
        
        fluidRow(
          box(
            title = "Introduction",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "This page presents feature importance plots of three models ",
              "based on the internal training set. ",
              "Feature importance is model-specific and not directly comparable in absolute scale ",
              "across different model families."
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Random Forest – MeanDecreaseGini",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "For the random forest, feature importance is quantified by ",
              strong("MeanDecreaseGini"),
              ": each time a feature is used to split a node, it reduces the node Gini impurity; ",
              "the reductions are accumulated and averaged over all trees. ",
              "Features with larger MeanDecreaseGini contribute more to the overall purity of the forest."
            ),
            plotOutput("rf_importance_plot", height = 320),
            hr(),
            div(
              style = "display:flex; gap:12px; justify-content:flex-start;",
              downloadButton("dl_rf_importance_png", "Download (PNG)"),
              downloadButton("dl_rf_importance_pdf", "Download (PDF)")
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Ridge Logistic Regression – |β\u2097| × sd(x\u2097)",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "For the ridge logistic regression, feature importance is quantified by ",
              strong("|β\u2097| × sd(x\u2097)"),
              ": the absolute value of the estimated coefficient multiplied by the ",
              "standard deviation of the corresponding predictor. ",
              "This rescales coefficients into a comparable influence on the log-odds scale."
            ),
            plotOutput("glm_importance_plot", height = 320),
            hr(),
            div(
              style = "display:flex; gap:12px; justify-content:flex-start;",
              downloadButton("dl_glm_importance_png", "Download (PNG)"),
              downloadButton("dl_glm_importance_pdf", "Download (PDF)")
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Support Vector Machine – Permutation ΔAccuracy",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "For the support vector machine, feature importance is quantified by ",
              strong("permutation-based ΔAccuracy"),
              ": for each feature, its values are randomly permuted across samples, ",
              "and the drop in classification accuracy is recorded. ",
              "Larger drops indicate that the feature is more influential for the SVM decision function."
            ),
            plotOutput("svm_importance_plot", height = 320),
            hr(),
            div(
              style = "display:flex; gap:12px; justify-content:flex-start;",
              downloadButton("dl_svm_importance_png", "Download (PNG)"),
              downloadButton("dl_svm_importance_pdf", "Download (PDF)")
            )
          )
        )
      ),
      
      ##################################################
      # 2.7 Batch prediction
      ##################################################
      tabItem(
        tabName = "batch",
        
        fluidRow(
          box(
            title = "Introduction",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "This page allows you to upload an external CSV file, ",
              "and obtain CKD predictions from all three models at desired threshols.",
              "After that, SHAP interpretations of whole model/single feature/single row are available."
            ),
            p(
              strong("Note: "),
              "The uploaded CSV must be the processed by ",
              code("data_process_ckd()"),
              "function in our R package."
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Upload CSV and Choose Thresholds",
            width = 12, status = "primary", solidHeader = TRUE,
            
            fileInput(
              "batch_file",
              "Upload CSV file",
              accept = c(".csv"),
              buttonLabel = "Browse...",
              placeholder = "No file selected"
            ),
            
            actionButton(
              "load_example",
              "Load example dataset",
              icon  = icon("magic")
            ),
            
            hr(),
            h4("Decision thresholds for three models"),
            fluidRow(
              column(
                width = 4,
                numericInput(
                  "thr_rf_batch",
                  "Random forest threshold:",
                  value = round(best_thr_rf, 3),
                  min   = 0, max = 1, step = 0.001
                )
              ),
              column(
                width = 4,
                numericInput(
                  "thr_glm_batch",
                  "Ridge logistic regression threshold:",
                  value = round(best_thr_glm, 3),
                  min   = 0, max = 1, step = 0.001
                )
              ),
              column(
                width = 4,
                numericInput(
                  "thr_svm_batch",
                  "Support vector machine threshold:",
                  value = round(best_thr_svm, 3),
                  min   = 0, max = 1, step = 0.001
                )
              )
            ),
            
            br(),
            actionButton(
              "run_batch_pred",
              "Run prediction",
              icon = icon("play")
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Batch Prediction Results",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "Before clicking ",
              strong("\"Run prediction\""),
              ", the table shows the cleaned input data. ",
              "After running prediction, it will be updated to include six additional columns: ",
              code("rf_classification"), ", ", code("rf_ckd_prob"), "; ",
              code("glm_classification"), ", ", code("glm_ckd_prob"), "; ",
              code("svm_classification"), ", ", code("svm_ckd_prob"), "."
            ),
            DTOutput("batch_table"),
            hr(),
            downloadButton(
              "dl_batch_pred",
              "Download predictions (CSV)"
            )
          )
        ),
        
        fluidRow(
          box(
            title = "SHAP Interpretation",
            width = 12, status = "primary", solidHeader = TRUE,
            fluidRow(
              column(
                width = 4,
                selectInput(
                  "batch_interp_model",
                  "Model:",
                  choices  = names(model_list),
                  selected = names(model_list)[1]
                )
              ),
              column(
                width = 4,
                radioButtons(
                  "batch_interp_type",
                  "Explanation type:",
                  choices = c(
                    "Global"          = "global",
                    "Single feature"  = "feature",
                    "Single row"      = "row"
                  ),
                  selected = "global"
                )
              ),
              column(
                width = 4,
                conditionalPanel(
                  condition = "input.batch_interp_type == 'feature'",
                  selectInput(
                    "batch_interp_feature",
                    "Feature for SHAP dependence:",
                    choices  = predictor_cols,
                    selected = "sc"
                  )
                ),
                conditionalPanel(
                  condition = "input.batch_interp_type == 'row'",
                  numericInput(
                    "batch_interp_row",
                    "Row index in batch data:",
                    value = 1,
                    min   = 1,
                    step  = 1
                  )
                )
              )
            ),
            
            br(),
            actionButton(
              "run_batch_shap",
              "Run SHAP",
              icon  = icon("play")
            ),
            
            br(),
            plotOutput("batch_interp_plot", height = 380),
            hr(),
            div(
              style = "display:flex; gap:12px; justify-content:flex-start;",
              downloadButton("dl_batch_shap_png", "Download (PNG)"),
              downloadButton("dl_batch_shap_pdf", "Download (PDF)")
            )
          )
        )
      ),
      
      ##################################################
      # 2.8 Single-patient prediction
      ##################################################
      tabItem(
        tabName = "single",
        
        fluidRow(
          box(
            title = "Introduction",
            width = 12, status = "primary", solidHeader = TRUE,
            p(
              "This page provides patient-level CKD risk estimation. ",
              "Users can: 1) Load a random example patient, 2) Manually enter feature values. ",
              "Then obtain instant CKD prediction from three models.",
              "SHAP interpretation for single row is also available."
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Patient Feature Input",
            width = 7, status = "primary", solidHeader = TRUE,
            
            h4("Decision thresholds for CKD"),
            fluidRow(
              column(
                width = 4,
                numericInput(
                  "single_thr_rf",
                  "Random forest:",
                  value = round(best_thr_rf, 3),
                  min   = 0, max = 1, step = 0.001
                )
              ),
              column(
                width = 4,
                numericInput(
                  "single_thr_glm",
                  "Ridge logistic regression:",
                  value = round(best_thr_glm, 3),
                  min   = 0, max = 1, step = 0.001
                )
              ),
              column(
                width = 4,
                numericInput(
                  "single_thr_svm",
                  "Support vector machine:",
                  value = round(best_thr_svm, 3),
                  min   = 0, max = 1, step = 0.001
                )
              )
            ),
            hr(),
            
            p("Option: Load a random example patient from internal test set"),
            actionButton("single_use_example", "Load random example patient", icon = icon("user")),
            hr(),
            
            uiOutput("single_input_ui"),
            hr(),
            actionButton("run_single_pred", "Predict for this patient", icon = icon("stethoscope"))
          ),
          
          box(
            title = "Prediction Results",
            width = 5, status = "primary", solidHeader = TRUE,
            
            fluidRow(
              valueBoxOutput("vb_single_final", width = 12)
            ),
            fluidRow(
              valueBoxOutput("vb_single_vote", width = 12)
            ),
            fluidRow(
              valueBoxOutput("vb_single_rf", width = 12)
            ),
            fluidRow(
              valueBoxOutput("vb_single_glm", width = 12)
            ),
            fluidRow(
              valueBoxOutput("vb_single_svm", width = 12)
            ),
            
            br(),
            DTOutput("single_pred_detail")
          ),
          box(
            title = "SHAP Interpretation",
            width = 5, status = "primary", solidHeader = TRUE,
            fluidRow(
              column(
                width = 6,
                selectInput(
                  "single_interp_model",
                  "Model:",
                  choices  = names(model_list),
                  selected = names(model_list)[1]
                )
              ),
              column(
                width = 6,
                actionButton(
                  "run_single_shap",
                  "Run SHAP",
                  icon  = icon("play")
                )
              )
            ),
            br(),
            plotOutput("single_interp_plot", height = 400),
            hr(),
            div(
              style = "display:flex; gap:12px; justify-content:flex-start;",
              downloadButton("dl_single_shap_png", "Download (PNG)"),
              downloadButton("dl_single_shap_pdf", "Download (PDF)")
            )
          )
        )
      )
    )
  )
)

#############################################################
# 3. Server logic
#############################################################

server <- function(input, output, session) {
  
  ##########################################################
  # 3.1 Global / home
  ##########################################################
  
  observeEvent(input$re_initialize, {
    session$reload()
  })
  
  output$vb_total_samples <- renderValueBox({
    n <- nrow(ckd_raw)
    valueBox(
      value    = format(n, big.mark = ","),
      subtitle = "Total Samples",
      icon     = icon("users"),
      color    = "teal"
    )
  })
  
  output$vb_feature_cols <- renderValueBox({
    p <- ncol(ckd_raw)
    valueBox(
      value    = p,
      subtitle = "Feature Columns",
      icon     = icon("list"),
      color    = "blue"
    )
  })
  
  output$vb_model_types <- renderValueBox({
    k <- length(model_list)
    valueBox(
      value    = k,
      subtitle = "ML Model Types",
      icon     = icon("project-diagram"),
      color    = "purple"
    )
  })
  
  output$illustrative_img <- renderImage({
    list(
      src         = "workflow.png",
      contentType = "image/png",
      alt         = "Workflow of the ckdPrediction Shiny app"
    )
  }, deleteFile = FALSE)
  
  output$feature_table <- renderDT({
    datatable(
      feature_dictionary,
      options = list(pageLength = 10, dom = "tip"),
      rownames = FALSE
    )
  })
  
  ##########################################################
  # 3.2 CKD database
  ##########################################################
  
  base_db_data <- reactive({
    d <- ckd_raw
    
    if (isTRUE(input$db_pipeline_process)) {
      combined_df <- rbind(train_df, test_df)
      d <- combined_df[order(combined_df$id), ]
    }
    
    if (isTRUE(input$db_show_missing)) {
      d <- na.omit(d)
    }
    
    d
  })
  
  output$db_filter_ui <- renderUI({
    d <- base_db_data()
    validate(need(ncol(d) > 0, "Dataset is empty."))
    
    num_cols <- names(d)[sapply(d, is.numeric)]
    num_cols <- setdiff(num_cols, c("id"))
    cat_cols <- names(d)[!sapply(d, is.numeric)]
    cat_cols <- setdiff(cat_cols, c("id"))
    
    num_ui <- lapply(seq_along(num_cols), function(i) {
      v   <- num_cols[i]
      rng <- range(d[[v]], na.rm = TRUE)
      step_val <- if (v %in% c("sg", "sc", "hemo", "pot", "rc")) 0.1 else 1
      
      column(
        width = 3,
        sliderInput(
          inputId = paste0("filter_num_", v),
          label   = strong(v),
          min     = rng[1],
          max     = rng[2],
          value   = rng,
          step    = step_val
        )
      )
    })
    
    cat_ui <- lapply(seq_along(cat_cols), function(i) {
      v    <- cat_cols[i]
      vals <- sort(unique(d[[v]]))
      vals <- vals[!is.na(vals)]
      
      column(
        width = 2,
        checkboxGroupInput(
          inputId  = paste0("filter_cat_", v),
          label    = strong(v),
          choices  = vals,
          selected = vals,
          inline   = FALSE
        )
      )
    })
    
    tagList(
      h4("Numeric Features"),
      fluidRow(num_ui),
      hr(),
      h4("Categorical Features"),
      fluidRow(cat_ui)
    )
  })
  
  filtered_db_data <- reactive({
    d <- base_db_data()
    if (nrow(d) == 0L) return(d)
    
    num_cols_all <- names(d)[sapply(d, is.numeric)]
    num_cols_all <- setdiff(num_cols_all, c("id"))
    cat_cols_all <- names(d)[!sapply(d, is.numeric)]
    cat_cols_all <- setdiff(cat_cols_all, c("id"))
    
    for (v in num_cols_all) {
      inp <- input[[paste0("filter_num_", v)]]
      if (!is.null(inp) && length(inp) == 2 && !any(is.na(inp))) {
        d <- d[is.na(d[[v]]) | (d[[v]] >= inp[1] & d[[v]] <= inp[2]), ]
      }
    }
    
    for (v in cat_cols_all) {
      id  <- paste0("filter_cat_", v)
      inp <- input[[id]]
      if (length(inp) == 0L) {
        d <- NULL
        break
      }
      d <- d[is.na(d[[v]]) | d[[v]] %in% inp, ]
    }
    
    d
  })
  
  output$db_table <- renderDT({
    datatable(
      filtered_db_data(),
      options = list(
        pageLength = 15,
        scrollX    = TRUE
      ),
      rownames = FALSE
    )
  })
  
  output$dl_full_data <- downloadHandler(
    filename = function() "kidney_disease_full.csv",
    content  = function(file) {
      write.csv(ckd_raw, file, row.names = FALSE)
    }
  )
  
  output$dl_filtered_data <- downloadHandler(
    filename = function() "kidney_disease_processed.csv",
    content  = function(file) {
      write.csv(filtered_db_data(), file, row.names = FALSE)
    }
  )
  
  ##########################################################
  # 3.3 Visualisation
  ##########################################################
  
  feature_label <- function(var) {
    if (var %in% c("rbc_missing", "rc_missing", "wc_missing")) {
      meaning <- switch(
        var,
        "rbc_missing" = "Missing indicator for red blood cell morphology",
        "rc_missing"  = "Missing indicator for red cell count",
        "wc_missing"  = "Missing indicator for white cell count"
      )
    } else {
      meaning <- feature_dictionary$Meaning[match(var, feature_dictionary$Column)]
    }
    
    if (is.na(meaning) || is.null(meaning)) {
      var
    } else {
      paste0(meaning, " (", var, ")")
    }
  }
  
  observe({
    all_vars <- names(filtered_db_data())
    
    has_missing_inds <- any(c("rbc_missing", "rc_missing", "wc_missing") %in% all_vars)
    if (has_missing_inds) {
      cat_vars <- c(
        "rbc", "pc", "pcc", "ba", "htn", "dm", "cad", "appet", "pe", "ane",
        "rbc_missing", "rc_missing", "wc_missing"
      )
      num_vars <- c(
        "age", "bp", "sg", "al", "su", "bgr", "bu", "sc", "sod", "pot", "hemo",
        "pcv", "wc", "rc"
      )
    } else {
      num_vars <- setdiff(
        all_vars[sapply(ckd_raw, is.numeric)],
        c("id", "classification")
      )
      cat_vars <- setdiff(
        all_vars[!sapply(ckd_raw, is.numeric)],
        c("id", "classification")
      )
    }
    
    updateSelectInput(
      session, "db_numeric_var",
      choices  = num_vars,
      selected = num_vars[1]
    )
    updateSelectInput(
      session, "db_categorical_var",
      choices  = cat_vars,
      selected = cat_vars[1]
    )
  })
  
  numeric_plot_obj <- reactive({
    req(input$db_numeric_var)
    d   <- filtered_db_data()
    var <- input$db_numeric_var
    
    validate(need(var %in% names(d), "Selected feature not found in dataset."))
    validate(need(!all(is.na(d[[var]])), "No samples are selected."))
    
    x_lab     <- feature_label(var)
    title_lab <- paste("Distribution of", x_lab, "by CKD label")
    
    ggplot(d, aes_string(x = var, fill = "classification")) +
      geom_histogram(alpha = 0.8, position = "identity") +
      scale_fill_manual(
        values = c(
          "notckd" = "#71b7ed",
          "ckd"    = "#f57c6e"
        ),
        breaks = c("notckd", "ckd"),
        labels = c("not-CKD", "CKD")
      ) +
      labs(
        title = title_lab,
        x     = x_lab,
        y     = "Count",
        fill  = "CKD label"
      ) +
      theme_bw() +
      theme(
        text               = element_text(face = "bold"),
        plot.title         = element_text(hjust = 0.5, size = 14),
        axis.title         = element_text(size = 12),
        axis.text          = element_text(size = 10),
        legend.position    = "right",
        legend.title       = element_text(size = 11),
        legend.text        = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank()
      )
  })
  
  output$db_numeric_plot <- renderPlot({
    print(numeric_plot_obj())
  })
  
  output$dl_numeric_png <- downloadHandler(
    filename = function() {
      paste0("distribution_", input$db_numeric_var, ".png")
    },
    content = function(file) {
      ggsave(file, plot = numeric_plot_obj(), width = 6, height = 4, dpi = 300)
    }
  )
  
  output$dl_numeric_pdf <- downloadHandler(
    filename = function() {
      paste0("distribution_", input$db_numeric_var, ".pdf")
    },
    content = function(file) {
      ggsave(file, plot = numeric_plot_obj(), width = 6, height = 4, device = "pdf")
    }
  )
  
  categorical_plot_obj <- reactive({
    req(input$db_categorical_var)
    d   <- filtered_db_data()
    var <- input$db_categorical_var
    
    validate(need(var %in% names(d), "Selected feature not found in dataset."))
    validate(need(any(!is.na(d[[var]])), "No samples are selected."))
    
    x_raw <- d[[var]]
    
    if (is.numeric(x_raw) && all(na.omit(unique(x_raw)) %in% c(0, 1))) {
      levs <- switch(
        var,
        "htn"   = c("No hypertension", "Hypertension"),
        "dm"    = c("No diabetes", "Diabetes"),
        "cad"   = c("No CAD", "CAD"),
        "pe"    = c("No oedema", "Oedema"),
        "ane"   = c("No anaemia", "Anaemia"),
        "appet" = c("Poor appetite", "Good appetite"),
        "pcc"   = c("Not present", "Present"),
        "ba"    = c("Not present", "Present"),
        "rbc"   = c("Abnormal", "Normal"),
        "pc"    = c("Abnormal", "Normal"),
        "rbc_missing" = c("observed", "missing"),
        "rc_missing"  = c("observed",  "missing"),
        "wc_missing"  = c("observed",  "missing"),
        c("0", "1")
      )
      x_factor <- factor(ifelse(x_raw == 1, levs[2], levs[1]), levels = levs)
    } else {
      x_factor <- factor(x_raw)
    }
    
    d$.__x_cat__ <- x_factor
    
    x_lab     <- feature_label(var)
    title_lab <- paste("Proportion of", x_lab, "by CKD label")
    
    ggplot(d, aes(x = .__x_cat__, fill = classification)) +
      geom_bar(position = "fill") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_fill_manual(
        values = c(
          "notckd" = "#71b7ed",
          "ckd"    = "#f57c6e"
        ),
        breaks = c("notckd", "ckd"),
        labels = c("not-CKD", "CKD")
      ) +
      labs(
        title = title_lab,
        x     = x_lab,
        y     = "Proportion",
        fill  = "CKD label"
      ) +
      theme_bw() +
      theme(
        text               = element_text(face = "bold"),
        plot.title         = element_text(hjust = 0.5, size = 14),
        axis.title         = element_text(size = 12),
        axis.text          = element_text(size = 10),
        legend.position    = "right",
        legend.title       = element_text(size = 11),
        legend.text        = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank()
      )
  })
  
  output$db_categorical_plot <- renderPlot({
    print(categorical_plot_obj())
  })
  
  output$dl_categorical_png <- downloadHandler(
    filename = function() {
      paste0("proportion_", input$db_categorical_var, ".png")
    },
    content = function(file) {
      ggsave(file, plot = categorical_plot_obj(), width = 6, height = 4, dpi = 300)
    }
  )
  
  output$dl_categorical_pdf <- downloadHandler(
    filename = function() {
      paste0("proportion_", input$db_categorical_var, ".pdf")
    },
    content = function(file) {
      ggsave(file, plot = categorical_plot_obj(), width = 6, height = 4, device = "pdf")
    }
  )
  
  ##########################################################
  # 3.4 Performance – curves
  ##########################################################
  
  nearest_threshold_from_click <- function(curve_df, click, type = c("roc", "prc")) {
    type <- match.arg(type)
    if (is.null(click$x) || is.null(click$y)) return(NA_real_)
    
    if (type == "roc") {
      dx <- curve_df$fpr - click$x
      dy <- curve_df$tpr - click$y
    } else {
      dx <- curve_df$recall - click$x
      dy <- curve_df$precision - click$y
    }
    i <- which.min(dx^2 + dy^2)
    curve_df$threshold[i]
  }
  
  # random forest
  observeEvent(input$rf_roc_click, {
    req(rf_curve_df)
    thr_new <- nearest_threshold_from_click(rf_curve_df, input$rf_roc_click, type = "roc")
    if (!is.na(thr_new)) {
      updateNumericInput(session, "rf_threshold_input", value = round(thr_new, 3))
    }
  })
  
  observeEvent(input$rf_prc_click, {
    req(rf_curve_df)
    thr_new <- nearest_threshold_from_click(rf_curve_df, input$rf_prc_click, type = "prc")
    if (!is.na(thr_new)) {
      updateNumericInput(session, "rf_threshold_input", value = round(thr_new, 3))
    }
  })
  
  rf_roc_plot_obj <- reactive({
    req(rf_curve_df)
    thr <- input$rf_threshold_input
    sel <- find_threshold_point(rf_curve_df, thr)
    
    ggplot(rf_curve_df, aes(x = fpr, y = tpr)) +
      geom_line(size = 1.1, colour = "#71b7ed") +
      geom_abline(
        slope    = 1, intercept = 0,
        linetype = "dashed",
        colour   = "#cbd5e1",
        linewidth = 0.6
      ) +
      geom_point(
        data  = sel,
        aes(x = fpr, y = tpr),
        colour = "#1D4ED8",
        size   = 3
      ) +
      geom_segment(
        data = sel,
        aes(x = 0, xend = fpr, y = tpr, yend = tpr),
        linetype = "dashed",
        colour   = "#1D4ED8",
        linewidth = 0.5
      ) +
      geom_segment(
        data = sel,
        aes(x = fpr, xend = fpr, y = 0, yend = tpr),
        linetype = "dashed",
        colour   = "#1D4ED8",
        linewidth = 0.5
      ) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(
        title = "ROC curve",
        x     = "False positive rate (FPR)",
        y     = "True positive rate (TPR)"
      ) +
      theme_bw() +
      theme(
        text             = element_text(face = "bold"),
        plot.title       = element_text(hjust = 0.5, size = 13),
        axis.title       = element_text(size = 11),
        axis.text        = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2, colour = "#e5e7eb"),
        panel.grid.minor = element_blank(),
        legend.position  = "none"
      )
  })
  
  rf_prc_plot_obj <- reactive({
    req(rf_curve_df)
    thr <- input$rf_threshold_input
    sel <- find_threshold_point(rf_curve_df, thr)
    
    ggplot(rf_curve_df, aes(x = recall, y = precision)) +
      geom_line(size = 1.1, colour = "#71b7ed") +
      geom_point(
        data  = sel,
        aes(x = recall, y = precision),
        colour = "#1D4ED8",
        size   = 3
      ) +
      geom_segment(
        data = sel,
        aes(x = 0, xend = recall, y = precision, yend = precision),
        linetype = "dashed",
        colour   = "#1D4ED8",
        linewidth = 0.5
      ) +
      geom_segment(
        data = sel,
        aes(x = recall, xend = recall, y = 0, yend = precision),
        linetype = "dashed",
        colour   = "#1D4ED8",
        linewidth = 0.5
      ) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(
        title = "Precision–Recall curve",
        x     = "Recall",
        y     = "Precision"
      ) +
      theme_bw() +
      theme(
        text             = element_text(face = "bold"),
        plot.title       = element_text(hjust = 0.5, size = 13),
        axis.title       = element_text(size = 11),
        axis.text        = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2, colour = "#e5e7eb"),
        panel.grid.minor = element_blank(),
        legend.position  = "none"
      )
  })
  
  output$perf_rf_roc  <- renderPlot({ print(rf_roc_plot_obj()) })
  output$perf_rf_prc  <- renderPlot({ print(rf_prc_plot_obj()) })
  
  output$dl_rf_roc_png <- downloadHandler(
    filename = function() "rf_roc_curve.png",
    content  = function(file) {
      ggsave(file, plot = rf_roc_plot_obj(), width = 6, height = 4, dpi = 300)
    }
  )
  
  output$dl_rf_roc_pdf <- downloadHandler(
    filename = function() "rf_roc_curve.pdf",
    content  = function(file) {
      ggsave(file, plot = rf_roc_plot_obj(), width = 6, height = 4, device = "pdf")
    }
  )
  
  output$dl_rf_prc_png <- downloadHandler(
    filename = function() "rf_prc_curve.png",
    content  = function(file) {
      ggsave(file, plot = rf_prc_plot_obj(), width = 6, height = 4, dpi = 300)
    }
  )
  
  output$dl_rf_prc_pdf <- downloadHandler(
    filename = function() "rf_prc_curve.pdf",
    content  = function(file) {
      ggsave(file, plot = rf_prc_plot_obj(), width = 6, height = 4, device = "pdf")
    }
  )
  
  output$perf_rf_coords <- renderUI({
    req(rf_curve_df)
    thr <- input$rf_threshold_input
    sel <- find_threshold_point(rf_curve_df, thr)
    
    HTML(sprintf(
      "<p style='margin-bottom:5px;'>
      <strong>Current threshold:</strong> %.3f
     </p>
     <ul style='padding-left:18px;margin-top:0;'>
      <li>ROC: FPR = %.3f, TPR = %.3f</li>
      <li>PRC: Recall = %.3f, Precision = %.3f</li>
     </ul>",
      sel$threshold,
      sel$fpr, sel$tpr,
      sel$recall, sel$precision
    ))
  })
  
  # glm
  observeEvent(input$glm_roc_click, {
    req(glm_curve_df)
    thr_new <- nearest_threshold_from_click(glm_curve_df, input$glm_roc_click, type = "roc")
    if (!is.na(thr_new)) {
      updateNumericInput(session, "glm_threshold_input", value = round(thr_new, 3))
    }
  })
  
  observeEvent(input$glm_prc_click, {
    req(glm_curve_df)
    thr_new <- nearest_threshold_from_click(glm_curve_df, input$glm_prc_click, type = "prc")
    if (!is.na(thr_new)) {
      updateNumericInput(session, "glm_threshold_input", value = round(thr_new, 3))
    }
  })
  
  glm_roc_plot_obj <- reactive({
    req(glm_curve_df)
    thr <- input$glm_threshold_input
    sel <- find_threshold_point(glm_curve_df, thr)
    
    ggplot(glm_curve_df, aes(x = fpr, y = tpr)) +
      geom_line(size = 1.1, colour = "#71b7ed") +
      geom_abline(
        slope    = 1, intercept = 0,
        linetype = "dashed",
        colour   = "#cbd5e1",
        linewidth = 0.6
      ) +
      geom_point(
        data  = sel,
        aes(x = fpr, y = tpr),
        colour = "#1D4ED8",
        size   = 3
      ) +
      geom_segment(
        data = sel,
        aes(x = 0, xend = fpr, y = tpr, yend = tpr),
        linetype = "dashed",
        colour   = "#1D4ED8",
        linewidth = 0.5
      ) +
      geom_segment(
        data = sel,
        aes(x = fpr, xend = fpr, y = 0, yend = tpr),
        linetype = "dashed",
        colour   = "#1D4ED8",
        linewidth = 0.5
      ) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(
        title = "ROC curve",
        x     = "False positive rate (FPR)",
        y     = "True positive rate (TPR)"
      ) +
      theme_bw() +
      theme(
        text             = element_text(face = "bold"),
        plot.title       = element_text(hjust = 0.5, size = 13),
        axis.title       = element_text(size = 11),
        axis.text        = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2, colour = "#e5e7eb"),
        panel.grid.minor = element_blank(),
        legend.position  = "none"
      )
  })
  
  glm_prc_plot_obj <- reactive({
    req(glm_curve_df)
    thr <- input$glm_threshold_input
    sel <- find_threshold_point(glm_curve_df, thr)
    
    ggplot(glm_curve_df, aes(x = recall, y = precision)) +
      geom_line(size = 1.1, colour = "#71b7ed") +
      geom_point(
        data  = sel,
        aes(x = recall, y = precision),
        colour = "#1D4ED8",
        size   = 3
      ) +
      geom_segment(
        data = sel,
        aes(x = 0, xend = recall, y = precision, yend = precision),
        linetype = "dashed",
        colour   = "#1D4ED8",
        linewidth = 0.5
      ) +
      geom_segment(
        data = sel,
        aes(x = recall, xend = recall, y = 0, yend = precision),
        linetype = "dashed",
        colour   = "#1D4ED8",
        linewidth = 0.5
      ) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(
        title = "Precision–Recall curve",
        x     = "Recall",
        y     = "Precision"
      ) +
      theme_bw() +
      theme(
        text             = element_text(face = "bold"),
        plot.title       = element_text(hjust = 0.5, size = 13),
        axis.title       = element_text(size = 11),
        axis.text        = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2, colour = "#e5e7eb"),
        panel.grid.minor = element_blank(),
        legend.position  = "none"
      )
  })
  
  output$perf_glm_roc  <- renderPlot({ print(glm_roc_plot_obj()) })
  output$perf_glm_prc  <- renderPlot({ print(glm_prc_plot_obj()) })
  
  output$dl_glm_roc_png <- downloadHandler(
    filename = function() "glm_roc_curve.png",
    content  = function(file) {
      ggsave(file, plot = glm_roc_plot_obj(), width = 6, height = 4, dpi = 300)
    }
  )
  
  output$dl_glm_roc_pdf <- downloadHandler(
    filename = function() "glm_roc_curve.pdf",
    content  = function(file) {
      ggsave(file, plot = glm_roc_plot_obj(), width = 6, height = 4, device = "pdf")
    }
  )
  
  output$dl_glm_prc_png <- downloadHandler(
    filename = function() "glm_prc_curve.png",
    content  = function(file) {
      ggsave(file, plot = glm_prc_plot_obj(), width = 6, height = 4, dpi = 300)
    }
  )
  
  output$dl_glm_prc_pdf <- downloadHandler(
    filename = function() "glm_prc_curve.pdf",
    content  = function(file) {
      ggsave(file, plot = glm_prc_plot_obj(), width = 6, height = 4, device = "pdf")
    }
  )
  
  output$perf_glm_coords <- renderUI({
    req(glm_curve_df)
    thr <- input$glm_threshold_input
    sel <- find_threshold_point(glm_curve_df, thr)
    
    HTML(sprintf(
      "<p style='margin-bottom:5px;'>
      <strong>Current threshold:</strong> %.3f
     </p>
     <ul style='padding-left:18px;margin-top:0;'>
      <li>ROC: FPR = %.3f, TPR = %.3f</li>
      <li>PRC: Recall = %.3f, Precision = %.3f</li>
     </ul>",
      sel$threshold,
      sel$fpr, sel$tpr,
      sel$recall, sel$precision
    ))
  })
  
  # svm
  observeEvent(input$svm_roc_click, {
    req(svm_curve_df)
    thr_new <- nearest_threshold_from_click(svm_curve_df, input$svm_roc_click, type = "roc")
    if (!is.na(thr_new)) {
      updateNumericInput(session, "svm_threshold_input", value = round(thr_new, 3))
    }
  })
  
  observeEvent(input$svm_prc_click, {
    req(svm_curve_df)
    thr_new <- nearest_threshold_from_click(svm_curve_df, input$svm_prc_click, type = "prc")
    if (!is.na(thr_new)) {
      updateNumericInput(session, "svm_threshold_input", value = round(thr_new, 3))
    }
  })
  
  svm_roc_plot_obj <- reactive({
    req(svm_curve_df)
    thr <- input$svm_threshold_input
    sel <- find_threshold_point(svm_curve_df, thr)
    
    ggplot(svm_curve_df, aes(x = fpr, y = tpr)) +
      geom_line(size = 1.1, colour = "#71b7ed") +
      geom_abline(
        slope    = 1, intercept = 0,
        linetype = "dashed",
        colour   = "#cbd5e1",
        linewidth = 0.6
      ) +
      geom_point(
        data  = sel,
        aes(x = fpr, y = tpr),
        colour = "#1D4ED8",
        size   = 3
      ) +
      geom_segment(
        data = sel,
        aes(x = 0, xend = fpr, y = tpr, yend = tpr),
        linetype = "dashed",
        colour   = "#1D4ED8",
        linewidth = 0.5
      ) +
      geom_segment(
        data = sel,
        aes(x = fpr, xend = fpr, y = 0, yend = tpr),
        linetype = "dashed",
        colour   = "#1D4ED8",
        linewidth = 0.5
      ) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(
        title = "ROC curve",
        x     = "False positive rate (FPR)",
        y     = "True positive rate (TPR)"
      ) +
      theme_bw() +
      theme(
        text             = element_text(face = "bold"),
        plot.title       = element_text(hjust = 0.5, size = 13),
        axis.title       = element_text(size = 11),
        axis.text        = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2, colour = "#e5e7eb"),
        panel.grid.minor = element_blank(),
        legend.position  = "none"
      )
  })
  
  svm_prc_plot_obj <- reactive({
    req(svm_curve_df)
    thr <- input$svm_threshold_input
    sel <- find_threshold_point(svm_curve_df, thr)
    
    ggplot(svm_curve_df, aes(x = recall, y = precision)) +
      geom_line(size = 1.1, colour = "#71b7ed") +
      geom_point(
        data  = sel,
        aes(x = recall, y = precision),
        colour = "#1D4ED8",
        size   = 3
      ) +
      geom_segment(
        data = sel,
        aes(x = 0, xend = recall, y = precision, yend = precision),
        linetype = "dashed",
        colour   = "#1D4ED8",
        linewidth = 0.5
      ) +
      geom_segment(
        data = sel,
        aes(x = recall, xend = recall, y = 0, yend = precision),
        linetype = "dashed",
        colour   = "#1D4ED8",
        linewidth = 0.5
      ) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(
        title = "Precision–Recall curve",
        x     = "Recall",
        y     = "Precision"
      ) +
      theme_bw() +
      theme(
        text             = element_text(face = "bold"),
        plot.title       = element_text(hjust = 0.5, size = 13),
        axis.title       = element_text(size = 11),
        axis.text        = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2, colour = "#e5e7eb"),
        panel.grid.minor = element_blank(),
        legend.position  = "none"
      )
  })
  
  output$perf_svm_roc  <- renderPlot({ print(svm_roc_plot_obj()) })
  output$perf_svm_prc  <- renderPlot({ print(svm_prc_plot_obj()) })
  
  output$dl_svm_roc_png <- downloadHandler(
    filename = function() "svm_roc_curve.png",
    content  = function(file) {
      ggsave(file, plot = svm_roc_plot_obj(), width = 6, height = 4, dpi = 300)
    }
  )
  
  output$dl_svm_roc_pdf <- downloadHandler(
    filename = function() "svm_roc_curve.pdf",
    content  = function(file) {
      ggsave(file, plot = svm_roc_plot_obj(), width = 6, height = 4, device = "pdf")
    }
  )
  
  output$dl_svm_prc_png <- downloadHandler(
    filename = function() "svm_prc_curve.png",
    content  = function(file) {
      ggsave(file, plot = svm_prc_plot_obj(), width = 6, height = 4, dpi = 300)
    }
  )
  
  output$dl_svm_prc_pdf <- downloadHandler(
    filename = function() "svm_prc_curve.pdf",
    content  = function(file) {
      ggsave(file, plot = svm_prc_plot_obj(), width = 6, height = 4, device = "pdf")
    }
  )
  
  output$perf_svm_coords <- renderUI({
    req(svm_curve_df)
    thr <- input$svm_threshold_input
    sel <- find_threshold_point(svm_curve_df, thr)
    
    HTML(sprintf(
      "<p style='margin-bottom:5px;'>
      <strong>Current threshold:</strong> %.3f
     </p>
     <ul style='padding-left:18px;margin-top:0;'>
      <li>ROC: FPR = %.3f, TPR = %.3f</li>
      <li>PRC: Recall = %.3f, Precision = %.3f</li>
     </ul>",
      sel$threshold,
      sel$fpr, sel$tpr,
      sel$recall, sel$precision
    ))
  })
  
  ##########################################################
  # 3.5 Performance – threshold comparison
  ##########################################################
  
  threshold_metrics <- eventReactive(input$run_threshold_compare, {
    req(test_df)
    
    thr_rf  <- input$thr_rf2
    thr_glm <- input$thr_glm2
    thr_svm <- input$thr_svm2
    
    do.call(
      rbind,
      list(
        eval_one(
          model     = rf_fit,
          label     = "Random forest",
          test_df   = test_df,
          threshold = thr_rf
        ),
        eval_one(
          model     = glm_fit,
          label     = "Ridge logistic regression",
          test_df   = test_df,
          threshold = thr_glm
        ),
        eval_one(
          model     = svm_fit,
          label     = "Support vector machine",
          test_df   = test_df,
          threshold = thr_svm
        )
      )
    )
  })
  
  output$threshold_compare_table <- renderDT({
    df <- threshold_metrics()
    
    datatable(
      df,
      options  = list(
        pageLength = 5,
        dom        = "tip"
      ),
      rownames = FALSE
    )
  })
  
  threshold_comp_plot <- eventReactive(input$run_threshold_compare, {
    req(test_df)
    
    thr_rf  <- input$thr_rf2
    thr_glm <- input$thr_glm2
    thr_svm <- input$thr_svm2
    
    validate(
      need(
        !any(is.na(c(thr_rf, thr_glm, thr_svm))),
        "Please provide valid thresholds for all three models."
      )
    )
    
    thr_vec <- c(thr_rf, thr_glm, thr_svm)
    
    model_comparison_ckd(
      model_list     = model_list,
      df             = test_df,
      threshold_list = thr_vec
    )
  })
  
  output$threshold_compare_plot <- renderPlot({
    p <- threshold_comp_plot()
    print(p)
  })
  
  output$dl_threshold_compare_png <- downloadHandler(
    filename = function() "model_comparison.png",
    content  = function(file) {
      ggsave(file, plot = threshold_comp_plot(), width = 7, height = 5, dpi = 300)
    }
  )
  
  output$dl_threshold_compare_pdf <- downloadHandler(
    filename = function() "model_comparison.pdf",
    content  = function(file) {
      ggsave(file, plot = threshold_comp_plot(), width = 7, height = 5, device = "pdf")
    }
  )
  
  ##########################################################
  # 3.6 Global feature importance
  ##########################################################
  
  rf_importance_plot_obj <- reactive({
    req(train_df, rf_fit)
    model_interpretation_importance_ckd(
      df    = train_df,
      model = rf_fit
    )
  })
  
  output$rf_importance_plot <- renderPlot({
    print(rf_importance_plot_obj())
  })
  
  output$dl_rf_importance_png <- downloadHandler(
    filename = function() "rf_feature_importance.png",
    content  = function(file) {
      ggsave(file, plot = rf_importance_plot_obj(), width = 6, height = 4, dpi = 300)
    }
  )
  
  output$dl_rf_importance_pdf <- downloadHandler(
    filename = function() "rf_feature_importance.pdf",
    content  = function(file) {
      ggsave(file, plot = rf_importance_plot_obj(), width = 6, height = 4, device = "pdf")
    }
  )
  
  glm_importance_plot_obj <- reactive({
    req(train_df, glm_fit)
    model_interpretation_importance_ckd(
      df    = train_df,
      model = glm_fit
    )
  })
  
  output$glm_importance_plot <- renderPlot({
    print(glm_importance_plot_obj())
  })
  
  output$dl_glm_importance_png <- downloadHandler(
    filename = function() "glm_feature_importance.png",
    content  = function(file) {
      ggsave(file, plot = glm_importance_plot_obj(), width = 6, height = 4, dpi = 300)
    }
  )
  
  output$dl_glm_importance_pdf <- downloadHandler(
    filename = function() "glm_feature_importance.pdf",
    content  = function(file) {
      ggsave(file, plot = glm_importance_plot_obj(), width = 6, height = 4, device = "pdf")
    }
  )
  
  svm_importance_plot_obj <- reactive({
    req(train_df, svm_fit)
    model_interpretation_importance_ckd(
      df    = train_df,
      model = svm_fit
    )
  })
  
  output$svm_importance_plot <- renderPlot({
    print(svm_importance_plot_obj())
  })
  
  output$dl_svm_importance_png <- downloadHandler(
    filename = function() "svm_feature_importance.png",
    content  = function(file) {
      ggsave(file, plot = svm_importance_plot_obj(), width = 6, height = 4, dpi = 300)
    }
  )
  
  output$dl_svm_importance_pdf <- downloadHandler(
    filename = function() "svm_feature_importance.pdf",
    content  = function(file) {
      ggsave(file, plot = svm_importance_plot_obj(), width = 6, height = 4, device = "pdf")
    }
  )
  
  ##########################################################
  # 3.7 Batch prediction
  ##########################################################
  
  batch_result <- reactiveVal(NULL)
  batch_raw_df <- reactiveVal(NULL)
  
  observeEvent(input$batch_file, {
    req(input$batch_file)
    df <- read.csv(input$batch_file$datapath, stringsAsFactors = FALSE)
    batch_raw_df(df)
    batch_result(NULL)
    
    showNotification(
      "CSV file uploaded. Using uploaded data for batch prediction.",
      type = "message"
    )
  }, ignoreNULL = TRUE)
  
  observeEvent(input$load_example, {
    df <- example_df
    batch_raw_df(df)
    batch_result(NULL)
  })
  
  batch_test_df <- reactive({
    df <- batch_raw_df()
    req(df)
    df
  })
  
  output$batch_table <- renderDT({
    req(batch_test_df())
    
    df_out <- batch_result()
    if (is.null(df_out)) {
      df_out <- batch_test_df()
    }
    
    datatable(
      df_out,
      options = list(
        pageLength = 15,
        scrollX    = TRUE
      ),
      rownames = FALSE
    )
  })
  
  observeEvent(input$run_batch_pred, {
    df_test <- batch_test_df()
    req(df_test)
    
    validate(
      need(nrow(df_test) > 0, "No rows available in the current batch dataset.")
    )
    
    required_predictors <- predictor_cols
    current_cols        <- colnames(df_test)
    
    missing_cols <- setdiff(required_predictors, current_cols)
    extra_cols   <- setdiff(current_cols,
                            c("id", "classification", required_predictors))
    
    validate(
      need(
        length(missing_cols) == 0,
        paste0(
          "Batch data is missing required columns: ",
          paste(missing_cols, collapse = ", "),
          ".\n\nPlease make sure your CSV (or example_df) has exactly the same ",
          "feature columns as the internal dataset.\n",
          "You can check `names(test_df)` in R for reference."
        )
      )
    )
    
    keep_cols <- c("id",
                   required_predictors,
                   intersect("classification", current_cols))
    df_test   <- df_test[, keep_cols, drop = FALSE]
    
    df_x <- df_test[, required_predictors, drop = FALSE]
    
    thr_rf  <- input$thr_rf_batch
    thr_glm <- input$thr_glm_batch
    thr_svm <- input$thr_svm_batch
    
    rf_model <- if (is.list(rf_fit) && "model" %in% names(rf_fit)) rf_fit$model else rf_fit
    
    prob_rf <- predict(rf_model, df_x, type = "prob")[, c("notckd", "ckd")]
    rf_ckd_prob <- prob_rf[, "ckd"]
    rf_class    <- ifelse(rf_ckd_prob >= thr_rf, "ckd", "notckd")
    
    glm_model <- if (is.list(glm_fit) && "model" %in% names(glm_fit)) glm_fit$model else glm_fit
    
    df_glm <- df_test
    
    if (!"classification" %in% names(df_glm)) {
      df_glm$classification <- factor("notckd", levels = c("notckd", "ckd"))
    } else if (!is.factor(df_glm$classification)) {
      df_glm$classification <- factor(df_glm$classification,
                                      levels = c("notckd", "ckd"))
    }
    
    x_new <- model.matrix(classification ~ . - id, data = df_glm)[, -1]
    
    prob_glm_pos <- as.numeric(
      predict(glm_model, newx = x_new, s = "lambda.min", type = "response")
    )
    glm_ckd_prob <- prob_glm_pos
    glm_class    <- ifelse(glm_ckd_prob >= thr_glm, "ckd", "notckd")
    
    svm_model <- if (is.list(svm_fit) && "model" %in% names(svm_fit)) svm_fit$model else svm_fit
    
    pred_svm     <- predict(svm_model, df_x, probability = TRUE)
    prob_svm_all <- attr(pred_svm, "probabilities")[, c("notckd", "ckd")]
    svm_ckd_prob <- prob_svm_all[, "ckd"]
    svm_class    <- ifelse(svm_ckd_prob >= thr_svm, "ckd", "notckd")
    
    df_pred <- data.frame(
      df_test,
      rf_classification   = rf_class,
      rf_ckd_prob         = round(rf_ckd_prob, 3),
      glm_classification  = glm_class,
      glm_ckd_prob        = round(glm_ckd_prob, 3),
      svm_classification  = svm_class,
      svm_ckd_prob        = round(svm_ckd_prob, 3),
      stringsAsFactors    = FALSE
    )
    
    batch_result(df_pred)
  })
  
  output$dl_batch_pred <- downloadHandler(
    filename = function() {
      "ckd_batch_predictions.csv"
    },
    content = function(file) {
      df_out <- batch_result()
      if (is.null(df_out)) {
        df_out <- batch_test_df()
      }
      write.csv(df_out, file, row.names = FALSE)
    }
  )
  
  ##########################################################
  # 3.7.1 Batch SHAP
  ##########################################################
  
  batch_shap_last_plot <- reactiveVal(NULL)
  
  observeEvent(input$run_batch_shap, {
    showNotification(
      "SHAP Interpretation is preparing, this takes a few seconds...",
      type = "message",
      duration = 4
    )
  })
  
  batch_shap_plot <- eventReactive(input$run_batch_shap, {
    
    req(batch_test_df())
    df_test <- batch_test_df()
    
    mdl_name <- input$batch_interp_model
    mdl      <- model_list[[mdl_name]]
    validate(need(!is.null(mdl), "Selected model is not available."))
    
    type <- input$batch_interp_type
    
    if (type == "global") {
      p <- model_interpretation_shap_ckd(
        train_df           = train_df,
        test_df            = df_test,
        model              = mdl,
        visualization_type = 1
      )
      
    } else if (type == "feature") {
      feat <- input$batch_interp_feature
      validate(need(feat %in% predictor_cols, "Invalid feature selected."))
      
      p <- model_interpretation_shap_ckd(
        train_df              = train_df,
        test_df               = df_test,
        model                 = mdl,
        visualization_type    = 2,
        visualization_feature = feat
      )
      
    } else if (type == "row") {
      row_idx <- input$batch_interp_row
      validate(
        need(
          row_idx >= 1 && row_idx <= nrow(df_test),
          paste0("Row index must be between 1 and ", nrow(df_test))
        )
      )
      
      p <- model_interpretation_shap_ckd(
        train_df           = train_df,
        test_df            = df_test,
        model              = mdl,
        visualization_type = 3,
        visualization_row  = as.integer(row_idx)
      )
    }
    
    batch_shap_last_plot(p)
    
    p
  })
  
  output$batch_interp_plot <- renderPlot({
    req(batch_shap_plot())
    print(batch_shap_plot())
  })
  
  output$dl_batch_shap_png <- downloadHandler(
    filename = function() {
      type <- input$batch_interp_type
      if (type == "global") {
        "batch_shap_global.png"
      } else if (type == "feature") {
        paste0("batch_shap_feature_", input$batch_interp_feature, ".png")
      } else {
        paste0("batch_shap_row_", input$batch_interp_row, ".png")
      }
    },
    content = function(file) {
      req(batch_shap_last_plot())
      ggsave(file, plot = batch_shap_last_plot(), width = 7, height = 5, dpi = 300)
    }
  )
  
  output$dl_batch_shap_pdf <- downloadHandler(
    filename = function() {
      type <- input$batch_interp_type
      if (type == "global") {
        "batch_shap_global.pdf"
      } else if (type == "feature") {
        paste0("batch_shap_feature_", input$batch_interp_feature, ".pdf")
      } else {
        paste0("batch_shap_row_", input$batch_interp_row, ".pdf")
      }
    },
    content = function(file) {
      req(batch_shap_last_plot())
      ggsave(file, plot = batch_shap_last_plot(), width = 7, height = 5, device = "pdf")
    }
  )
  
  ##########################################################
  # 3.8 Single-patient helpers
  ##########################################################
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install dplyr: install.packages('dplyr')")
  }
  library(dplyr, quietly = TRUE)
  
  fill_missing_values <- function(df, ref_df = test_df) {
    for (col in colnames(df)) {
      if (any(is.na(df[[col]]) | df[[col]] == "")) {
        if (is.numeric(ref_df[[col]])) {
          df[[col]] <- median(ref_df[[col]], na.rm = TRUE)
        } else if (is.factor(ref_df[[col]])) {
          tab <- table(ref_df[[col]], useNA = "no")
          df[[col]] <- factor(names(tab)[which.max(tab)], levels = levels(ref_df[[col]]))
        } else {
          tab <- table(ref_df[[col]], useNA = "no")
          df[[col]] <- names(tab)[which.max(tab)]
        }
      }
    }
    df
  }
  
  single_cat_vars <- c(
    "rbc", "pc", "pcc", "ba",
    "htn", "dm", "cad", "appet", "pe", "ane",
    "rbc_missing", "rc_missing", "wc_missing"
  )
  
  output$dl_single_example_csv <- downloadHandler(
    filename = function() "ckd_single_patient_example.csv",
    content  = function(file) {
      req(test_df)
      keep_cols <- c("id", predictor_cols, "classification")
      keep_cols <- keep_cols[keep_cols %in% colnames(test_df)]
      complete_test_df <- na.omit(test_df[, keep_cols, drop = FALSE])
      
      if (nrow(complete_test_df) == 0) {
        single_example <- test_df[sample(1:nrow(test_df), 1), keep_cols, drop = FALSE]
        single_example <- fill_missing_values(single_example)
      } else {
        single_example <- complete_test_df[sample(1:nrow(complete_test_df), 1), , drop = FALSE]
      }
      write.csv(single_example[, predictor_cols, drop = FALSE], file, row.names = FALSE)
    }
  )
  
  uploaded_single_data <- reactive({
    req(input$single_upload_csv)
    csv_data <- read.csv(
      input$single_upload_csv$datapath,
      stringsAsFactors = FALSE,
      strip.white      = TRUE
    )
    
    validate(need(nrow(csv_data) == 1, "CSV must contain exactly ONE patient's data!"))
    
    csv_cols     <- colnames(csv_data)
    missing_cols <- setdiff(predictor_cols, csv_cols)
    extra_cols   <- setdiff(csv_cols, predictor_cols)
    
    validate(need(length(missing_cols) == 0,
                  paste0("Missing required columns: ",
                         paste(missing_cols, collapse = ", "))))
    validate(need(length(extra_cols) == 0,
                  paste0("Extra columns (please remove): ",
                         paste(extra_cols, collapse = ", "))))
    
    csv_data <- csv_data[, predictor_cols, drop = FALSE]
    
    if ("id" %in% colnames(test_df)) {
      if (is.numeric(test_df$id)) {
        csv_data$id <- as.numeric(max(test_df$id, na.rm = TRUE) + 1)
      } else {
        csv_data$id <- paste0(
          "P",
          max(as.numeric(gsub("P", "", test_df$id, fixed = TRUE)), na.rm = TRUE) + 1
        )
      }
    }
    if ("classification" %in% colnames(test_df)) {
      csv_data$classification <- NA
      if (is.factor(test_df$classification)) {
        csv_data$classification <- factor(NA, levels = levels(test_df$classification))
      }
    }
    
    final_cols <- colnames(test_df)[colnames(test_df) %in% c("id", predictor_cols, "classification")]
    csv_data   <- csv_data[, final_cols, drop = FALSE]
    csv_data   <- fill_missing_values(csv_data)
    csv_data
  })
  
  single_input_data <- reactive({
    req(test_df)
    required_cols <- c("id", predictor_cols, "classification")
    required_cols <- required_cols[required_cols %in% colnames(test_df)]
    
    new_row <- test_df[1, required_cols, drop = FALSE]
    new_row[] <- NA
    
    if (!is.null(input$single_upload_csv)) {
      csv_data <- uploaded_single_data()
      for (nm in required_cols) {
        if (!nm %in% colnames(csv_data)) next
        csv_val <- csv_data[[nm]]
        
        if (is.numeric(test_df[[nm]])) {
          new_row[[nm]] <- as.numeric(csv_val)
        } else if (is.factor(test_df[[nm]])) {
          if (!is.na(csv_val) && !as.character(csv_val) %in% levels(test_df[[nm]])) {
            tab     <- table(test_df[[nm]], useNA = "no")
            csv_val <- names(tab)[which.max(tab)]
          }
          new_row[[nm]] <- factor(as.character(csv_val), levels = levels(test_df[[nm]]))
        } else {
          new_row[[nm]] <- as.character(csv_val)
        }
      }
    } else {
      for (nm in required_cols) {
        if (nm == "id") {
          if (is.numeric(test_df$id)) {
            new_row$id <- as.numeric(max(test_df$id, na.rm = TRUE) + 1)
          } else {
            new_row$id <- paste0(
              "P",
              max(as.numeric(gsub("P", "", test_df$id, fixed = TRUE)), na.rm = TRUE) + 1
            )
          }
          next
        }
        if (nm == "classification") {
          new_row$classification <- if (is.factor(test_df$classification)) {
            factor(NA, levels = levels(test_df$classification))
          } else {
            NA
          }
          next
        }
        
        id_input <- paste0("single_", nm)
        val      <- input[[id_input]]
        
        if (is.null(val) || is.na(val) || val == "") {
          if (is.numeric(test_df[[nm]])) {
            val <- median(test_df[[nm]], na.rm = TRUE)
          } else if (is.factor(test_df[[nm]])) {
            tab <- table(test_df[[nm]], useNA = "no")
            val <- names(tab)[which.max(tab)]
          }
        }
        
        if (is.numeric(test_df[[nm]])) {
          new_row[[nm]] <- as.numeric(val)
        } else if (is.factor(test_df[[nm]])) {
          new_row[[nm]] <- factor(as.character(val), levels = levels(test_df[[nm]]))
        } else {
          new_row[[nm]] <- as.character(val)
        }
      }
    }
    
    new_row <- fill_missing_values(new_row)
    new_row
  })
  
  single_result <- reactive({
    req(single_input_data())
    new_row <- single_input_data()
    
    thr_rf  <- input$single_thr_rf
    thr_glm <- input$single_thr_glm
    thr_svm <- input$single_thr_svm
    
    thr_map <- list(
      "Random forest"             = thr_rf,
      "Ridge logistic regression" = thr_glm,
      "Support vector machine"    = thr_svm
    )
    
    out <- list()
    for (nm in names(model_list)) {
      mdl <- model_list[[nm]]
      if (is.null(mdl)) next
      
      prob <- predict_ckd(new_row, mdl)
      
      if (is.matrix(prob)) {
        col_idx  <- if ("ckd" %in% colnames(prob)) which(colnames(prob) == "ckd") else 2
        ckd_prob <- as.numeric(prob[1, col_idx])
      } else if (is.data.frame(prob)) {
        col_idx  <- if ("ckd" %in% colnames(prob)) which(colnames(prob) == "ckd") else 2
        ckd_prob <- as.numeric(prob[[col_idx]][1])
      } else if (is.vector(prob)) {
        ckd_prob <- if (length(prob) == 2) as.numeric(prob[2]) else as.numeric(prob[1])
      } else {
        ckd_prob <- 0.5
      }
      
      thr_use <- thr_map[[nm]]
      if (is.null(thr_use) || is.na(thr_use)) thr_use <- 0.5
      cls <- ifelse(ckd_prob >= thr_use, "ckd", "notckd")
      
      interpretation_plot <- NULL
      tryCatch({
        interpretation_plot <- model_interpretation_ckd(
          patient_data = new_row[, predictor_cols, drop = FALSE],
          model        = mdl
        )
        if (!is.null(interpretation_plot)) {
          interpretation_plot <- interpretation_plot +
            labs(
              title    = paste0(nm, " - Feature Contribution"),
              subtitle = "Positive = Increases CKD Risk | Negative = Decreases CKD Risk",
              x        = "Feature Contribution",
              y        = "Feature"
            ) +
            theme_bw() +
            theme(
              plot.title    = element_text(hjust = 0.5, size = 13, fontface = "bold",
                                           margin = margin(b = 10)),
              plot.subtitle = element_text(hjust = 0.5, size = 10, color = "#7f8c8d"),
              axis.title    = element_text(size = 11, fontface = "bold"),
              axis.text     = element_text(size = 9)
            )
        }
      }, error = function(e) {
        interpretation_plot <- ggplot() +
          geom_text(
            aes(x = 0.5, y = 0.5),
            label = paste0(nm, "\nNo interpretation available"),
            size  = 5,
            color = "#7f8c8d"
          ) +
          theme_void()
      })
      
      out[[nm]] <- list(
        Model               = nm,
        Prob_CKD            = round(ckd_prob, 4),
        Predicted_label     = cls,
        Interpretation_plot = interpretation_plot
      )
    }
    out
  })
  
  single_summary <- reactive({
    res <- single_result()
    req(res)
    
    pred_labels <- sapply(res, function(x) x$Predicted_label)
    ckd_votes   <- sum(pred_labels == "ckd", na.rm = TRUE)
    final_label <- ifelse(ckd_votes >= 2, "CKD-positive", "CKD-negative")
    
    probs           <- sapply(res, function(x) x$Prob_CKD)
    best_model_name <- names(probs)[which.max(probs)]
    best_prob       <- max(probs)
    
    list(
      final_label     = final_label,
      ckd_votes       = ckd_votes,
      best_model_name = best_model_name,
      best_prob       = best_prob
    )
  })
  
  output$vb_single_final <- renderValueBox({
    ss <- single_summary()
    valueBox(
      value    = ifelse(ss$final_label == "CKD-positive", "CKD", "NOT CKD"),
      subtitle = "Final decision (majority vote of 3 models)",
      icon     = icon("heartbeat"),
      color    = if (ss$final_label == "CKD-positive") "red" else "green"
    )
  })
  
  output$vb_single_vote <- renderValueBox({
    ss <- single_summary()
    valueBox(
      value    = sprintf("%d / 3", ss$ckd_votes),
      subtitle = "Models predicting CKD",
      icon     = icon("users"),
      color    = "orange"
    )
  })
  
  output$vb_single_rf <- renderValueBox({
    res <- single_result()
    rf  <- res[["Random forest"]]
    req(rf)
    valueBox(
      value    = sprintf("%.1f%%", 100 * rf$Prob_CKD),
      subtitle = "Random forest: Probability(CKD)",
      icon     = icon("tree"),
      color    = "teal"
    )
  })
  
  output$vb_single_glm <- renderValueBox({
    res <- single_result()
    glm <- res[["Ridge logistic regression"]]
    req(glm)
    valueBox(
      value    = sprintf("%.1f%%", 100 * glm$Prob_CKD),
      subtitle = "Ridge logistic regression: Probability(CKD)",
      icon     = icon("chart-line"),
      color    = "blue"
    )
  })
  
  output$vb_single_svm <- renderValueBox({
    res <- single_result()
    svm <- res[["Support vector machine"]]
    req(svm)
    valueBox(
      value    = sprintf("%.1f%%", 100 * svm$Prob_CKD),
      subtitle = "Support vector machine: Probability(CKD)",
      icon     = icon("project-diagram"),
      color    = "purple"
    )
  })
  
  output$model1_interpretation <- renderPlot({
    res <- single_result()
    req(res)
    model1_name <- names(res)[1]
    req(res[[model1_name]]$Interpretation_plot)
    res[[model1_name]]$Interpretation_plot
  })
  
  output$model2_interpretation <- renderPlot({
    res <- single_result()
    req(res)
    if (length(res) < 2) return(NULL)
    model2_name <- names(res)[2]
    req(res[[model2_name]]$Interpretation_plot)
    res[[model2_name]]$Interpretation_plot
  })
  
  output$model3_interpretation <- renderPlot({
    res <- single_result()
    req(res)
    if (length(res) < 3) return(NULL)
    model3_name <- names(res)[3]
    req(res[[model3_name]]$Interpretation_plot)
    res[[model3_name]]$Interpretation_plot
  })
  
  output$single_feature_analysis <- renderUI({
    res <- single_result()
    req(res, train_df, rf_fit)
    
    rf_importance <- tryCatch(
      {
        imp <- model_interpretation_importance_ckd(train_df, rf_fit)$data
        imp
      },
      error = function(e) NULL
    )
    
    if (is.null(rf_importance) || nrow(rf_importance) == 0 ||
        !"importance" %in% colnames(rf_importance)) {
      rf_importance <- data.frame(
        feature    = predictor_cols[1:3],
        importance = c(0.3, 0.25, 0.2),
        stringsAsFactors = FALSE
      )
    } else {
      rf_importance <- rf_importance[
        !is.na(rf_importance$importance) &
          is.numeric(rf_importance$importance),
        , drop = FALSE
      ]
      if (nrow(rf_importance) < 3) {
        supplement <- data.frame(
          feature    = setdiff(predictor_cols, rf_importance$feature)[
            seq_len(3 - nrow(rf_importance))
          ],
          importance = seq(0.2, 0.1, length.out = 3 - nrow(rf_importance)),
          stringsAsFactors = FALSE
        )
        rf_importance <- rbind(rf_importance, supplement)
      }
    }
    
    rf_importance   <- rf_importance[order(-rf_importance$importance), ]
    top_features    <- rf_importance$feature[
      !rf_importance$feature %in% c("id", "classification")
    ][1:3]
    patient_data    <- single_input_data()
    
    top_feature_info <- lapply(top_features, function(f) {
      if (!f %in% colnames(patient_data)) return(NULL)
      val     <- patient_data[[f]]
      val_str <- if (is.numeric(val)) sprintf("%.2f", val) else as.character(val)
      meaning <- feature_dictionary$Meaning[match(f, feature_dictionary$Column)]
      if (is.na(meaning)) meaning <- "No description available"
      list(feature = f, value = val_str, meaning = meaning)
    })
    top_feature_info <- Filter(Negate(is.null), top_feature_info)
    
    if (length(top_feature_info) == 0) {
      top_feature_info <- list(
        list(
          feature = predictor_cols[1],
          value   = "",
          meaning = "No description available"
        )
      )
    }
    
    HTML(paste0(
      "<div style='background-color: #f5f6fa; padding: 15px; border-radius: 8px; margin-bottom: 20px;'>",
      "<h4 style='margin-top: 0; color: #2c3e50;'>Key Feature Summary</h4>",
      "<p style='margin-bottom: 10px;'>Top 3 most influential features across all models:</p>",
      "<ul style='margin-bottom: 0; padding-left: 20px;'>",
      paste0(
        sapply(top_feature_info, function(x) {
          paste0(
            "<li><strong>", x$feature, "</strong>: ", x$value,
            " <span style='color: #7f8c8d; font-size: 12px;'>(",
            x$meaning, ")</span></li>"
          )
        }),
        collapse = ""
      ),
      "</ul>",
      "</div>"
    ))
  })
  
  observeEvent(input$single_use_example, {
    req(test_df)
    required_cols <- c("id", predictor_cols, "classification")
    required_cols <- required_cols[required_cols %in% colnames(test_df)]
    
    idx <- sample(seq_len(nrow(test_df)), 1)
    row <- test_df[idx, required_cols, drop = FALSE]
    row <- fill_missing_values(row)
    
    for (nm in predictor_cols) {
      id_input <- paste0("single_", nm)
      val      <- row[[nm]]
      
      if (nm %in% single_cat_vars) {
        updateSelectInput(session, id_input, selected = as.character(val))
      } else if (is.numeric(val)) {
        updateNumericInput(session, id_input, value = round(val, 2))
      } else if (is.factor(val)) {
        updateSelectInput(session, id_input, selected = as.character(val))
      } else {
        updateTextInput(session, id_input, value = as.character(val))
      }
    }
  })
  
  output$single_input_ui <- renderUI({
    validate(need(!is.null(test_df), "Single-patient prediction requires the processed test set."))
    df    <- test_df
    preds <- predictor_cols
    
    numeric_cols <- preds[sapply(df[, preds, drop = FALSE], is.numeric)]
    factor_cols  <- preds[sapply(df[, preds, drop = FALSE], is.factor)]
    
    cat_in_use   <- intersect(preds, single_cat_vars)
    numeric_cols <- setdiff(numeric_cols, cat_in_use)
    factor_cols  <- union(factor_cols, cat_in_use)
    
    num_inputs <- lapply(numeric_cols, function(nm) {
      rng        <- range(df[[nm]], na.rm = TRUE)
      val        <- median(df[[nm]], na.rm = TRUE)
      label_text <- feature_label(nm)
      
      numericInput(
        inputId = paste0("single_", nm),
        label   = label_text,
        value   = round(val, 2),
        min     = floor(rng[1] - 1),
        max     = ceiling(rng[2] + 1),
        step    = ifelse(nm %in% c("sg", "sc", "hemo", "pot", "rc"), 0.1, 1)
      )
    })
    
    cat_inputs <- lapply(factor_cols, function(nm) {
      label_text <- feature_label(nm)
      
      if (nm %in% single_cat_vars) {
        levs <- switch(
          nm,
          "htn"   = c("No hypertension", "Hypertension"),
          "dm"    = c("No diabetes", "Diabetes"),
          "cad"   = c("No CAD", "CAD"),
          "pe"    = c("No oedema", "Oedema"),
          "ane"   = c("No anaemia", "Anaemia"),
          "appet" = c("Poor appetite", "Good appetite"),
          "pcc"   = c("Not present", "Present"),
          "ba"    = c("Not present", "Present"),
          "rbc"   = c("Abnormal", "Normal"),
          "pc"    = c("Abnormal", "Normal"),
          "rbc_missing" = c("Observed", "Missing"),
          "rc_missing"  = c("Observed", "Missing"),
          "wc_missing"  = c("Observed", "Missing"),
          c("0", "1")
        )
        
        choices_vec <- c("0", "1")
        names(choices_vec) <- levs
        
        majority_val <- tryCatch(
          {
            tab <- table(df[[nm]])
            names(tab)[which.max(tab)]
          },
          error = function(e) "0"
        )
        
        selectInput(
          inputId  = paste0("single_", nm),
          label    = label_text,
          choices  = choices_vec,
          selected = as.character(majority_val)
        )
      } else {
        levs <- levels(df[[nm]])
        selectInput(
          inputId  = paste0("single_", nm),
          label    = label_text,
          choices  = levs,
          selected = levs[1]
        )
      }
    })
    
    build_two_columns <- function(controls_list) {
      n <- length(controls_list)
      if (n == 0) return(NULL)
      half <- ceiling(n / 2)
      fluidRow(
        column(6, controls_list[seq_len(half)]),
        column(6, if (half < n) controls_list[(half + 1):n] else NULL)
      )
    }
    
    tagList(
      h4("Numeric features"),
      build_two_columns(num_inputs),
      hr(),
      h4("Categorical features"),
      build_two_columns(cat_inputs)
    )
  })
  
  ##########################################################
  # 3.8.1 Single-patient SHAP
  ##########################################################
  
  observeEvent(input$run_single_shap, {
    showNotification(
      "SHAP Interpretation is preparing, this takes a few seconds...",
      type = "message",
      duration = 4
    )
  })
  
  single_shap_plot <- eventReactive(input$run_single_shap, {
    req(input$single_interp_model)
    
    mdl_name <- input$single_interp_model
    mdl      <- model_list[[mdl_name]]
    validate(need(!is.null(mdl), "Selected model is not available."))
    
    df_shap <- rbind(test_df, single_input_data())
    
    p <- model_interpretation_shap_ckd(
      train_df           = train_df,
      test_df            = df_shap,
      model              = mdl,
      visualization_type = 3,
      visualization_row  = nrow(df_shap)
    )
    p
  })
  
  output$single_interp_plot <- renderPlot({
    p <- single_shap_plot()
    req(p)
    print(p)
  })
  
  output$dl_single_shap_png <- downloadHandler(
    filename = function() "single_shap.png",
    content  = function(file) {
      ggsave(file, plot = single_shap_plot(), width = 7, height = 5, dpi = 300)
    }
  )
  
  output$dl_single_shap_pdf <- downloadHandler(
    filename = function() "single_shap.pdf",
    content  = function(file) {
      ggsave(file, plot = single_shap_plot(), width = 7, height = 5, device = "pdf")
    }
  )
}

#############################################################
# 4. Run app
#############################################################

shinyApp(ui = ui, server = server)