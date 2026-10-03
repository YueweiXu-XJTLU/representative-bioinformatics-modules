# ---- app.R (single-file style) ----
library(shiny)
library(shinydashboard)
library(DT)
library(shinyjs)
library(ggplot2)
suppressWarnings(library(randomForest))

# --- load fitted random forest model object ---
rf_fit <- readRDS("rf_fit.rds")

# --- load data CSV ---
m6A_data    <- read.csv("../../data/m6A_imbalanced.csv")
example_data <- read.csv("m6A_input_example.csv")

# --- define necessary functions for ML prediction ---
dna_encoding <- function(dna_strings) {
  nn    <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  stopifnot(
    all(c(
      "gc_content", "RNA_type", "RNA_region", "exon_length",
      "distance_to_junction", "evolutionary_conservation", "DNA_5mer"
    ) %in% colnames(feature_df))
  ) # Check errors if incorrect column names of input data.frame
  
  # 1. process input data
  feature_df <- cbind(feature_df, dna_encoding(feature_df$DNA_5mer))
  feature_df$RNA_type   <- factor(feature_df$RNA_type,   levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  
  # 2. predict and classify
  prob        <- predict(ml_fit, feature_df, type = "prob")
  prob_binary <- ifelse(prob[, 2] <= positive_threshold, "Negative", "Positive")
  
  # 3. add rows for output
  feature_df$predicted_m6A_prob   <- prob[, 2]  # predicted positive-class probability
  feature_df$predicted_m6A_status <- prob_binary # either "Positive" or "Negative"
  feature_df$predicted_m6A_status <- factor(feature_df$predicted_m6A_status, levels = c("Negative", "Positive"))
  
  # 4. delete expanded one-hot DNA-5mer columns
  feature_df <- subset(feature_df, select = -c(nt_pos1, nt_pos2, nt_pos3, nt_pos4, nt_pos5))
  
  return(feature_df) # return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}

prediction_single <- function(
    ml_fit, gc_content, RNA_type, RNA_region, exon_length,
    distance_to_junction, evolutionary_conservation, DNA_5mer,
    positive_threshold = 0.5
) {
  # 1. construct single-row data frame
  single_row_df <- data.frame(
    gc_content               = as.numeric(gc_content),
    RNA_type                 = RNA_type,
    RNA_region               = RNA_region,
    exon_length              = as.numeric(exon_length),
    distance_to_junction     = as.numeric(distance_to_junction),
    evolutionary_conservation = as.numeric(evolutionary_conservation),
    DNA_5mer                 = toupper(trimws(DNA_5mer))
  )
  
  # 2. evoke prediction_multiple() (factor convertion will be conducted automatically here)
  results_full <- prediction_multiple(ml_fit, single_row_df, positive_threshold)
  
  # 3. extract prediction results
  returned_vector <- c(
    results_full[1, "predicted_m6A_prob"],
    as.character(results_full[1, "predicted_m6A_status"])
  )
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  
  return(returned_vector) # return a named vector with values for predicted m6A prob and predicted m6A status
}

# --- create shiny code bellow ---
####################################################################################################
# ui
####################################################################################################
ui <- dashboardPage(
  title = "m6A Bioinformatics - P6",
  skin  = "blue",
  
  # page header
  dashboardHeader(title = "m6A Bioinformatics"),
  
  # sidebar
  dashboardSidebar(
    # a button for re-initialization
    actionButton("re_initialization", "Reinitialize Everything", icon = icon("refresh")),
    
    sidebarMenu(
      ##################################################
      # mini bioinformatics database
      ##################################################
      menuItem("Home",          tabName = "home",          icon = icon("home")),
      menuItem("Data",          tabName = "data",          icon = icon("table")),
      menuItem("Visualization", tabName = "visualization", icon = icon("chart-simple")),
      menuItem("Download",      tabName = "download",      icon = icon("download")),
      ##################################################
      # m6A prediction
      ##################################################
      menuItem(
        "ML Prediction", icon = icon("flask"), startExpanded = TRUE,
        menuSubItem("m6A Prediction - Multiple", tabName = "prediction_multiple", icon = icon("layer-group")),
        menuSubItem("m6A Prediction - Single",   tabName = "prediction_single",   icon = icon("keyboard"))
      )
    )
  ),
  
  # page body
  dashboardBody(
    useShinyjs(),
    
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
      # --- Home Page ---
      tabItem(
        tabName = "home",
        
        # static valueBox
        fluidRow(
          valueBox(value = "48,055", subtitle = "Total rows",               icon = icon("database"),  color = "light-blue", width = 3),
          valueBox(value = "23.4%",  subtitle = "Positive class proportion", icon = icon("chart-pie"), color = "green",      width = 3),
          valueBox(value = "16",     subtitle = "Feature columns",           icon = icon("list-check"), color = "teal",       width = 3),
          valueBox(value = "4",      subtitle = "Distinct RNA types",        icon = icon("dna"),        color = "yellow",     width = 3)
        ),
        
        # overview & how-to-use
        fluidRow(
          box(
            title = "Overview",
            width = 6, status = "primary", solidHeader = TRUE, height = 320,
            HTML("
    <p><b>Objective:</b> This mini bioinformatics database aggregates 48,055 candidate m<sup>6</sup>A sites with multi-modal feature annotations for downstream machine-learning and exploratory analysis.</p>
    <p><b>Intended Use:</b> The application provides a full view of the dataset, PNG/PDF-exportable graphics filtered by categorical variables, and CSV downloads for the whole table or stratified subsets. A pretrained Random Forest is included to predict \"m6A_status\" given user inputs and a user-defined decision threshold, with downloadable prediction outputs.</p>
  ")
          ),
          box(
            title = "How to Use",
            width = 6, status = "primary", solidHeader = TRUE, height = 320,
            tags$ol(
              tags$li(HTML("<b>Data:</b> Browse and search the full table; sort by any column.")),
              tags$li(HTML("<b>Visualization:</b> Render category pie charts and a filtered histogram of \"m6A_level\"; export figures as PNG/PDF.")),
              tags$li(HTML("<b>Download:</b> Export the full dataset or RNA-type specific subsets as CSV.")),
              tags$li(HTML("<b>ML Prediction:</b> Batch or single-sample inference with a pretrained Random Forest; adjust the decision threshold and export predictions as CSV."))
            )
          )
        ),
        
        # --- illustrative image ---
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
        
        # key features
        fluidRow(
          box(
            title = "Key Feature Columns",
            width = 12, status = "primary", solidHeader = TRUE,
            tags$dl(
              tags$dt("m6A_status"),               tags$dd("Ground-truth label (Positive/Negative)"),
              tags$dt("m6A_level"),                tags$dd("m6A enrichment score per site"),
              tags$dt("DNA_5mer"),                 tags$dd("Central 5-mer sequence"),
              tags$dt("gc_content"),               tags$dd("Local GC proportion"),
              tags$dt("RNA_type"),                 tags$dd("mRNA / lincRNA / lncRNA / pseudogene"),
              tags$dt("RNA_region"),               tags$dd("CDS / intron / 3'UTR / 5'UTR"),
              tags$dt("exon_length"),              tags$dd("Host exon length"),
              tags$dt("distance_to_junction"),     tags$dd("Distance to splice junction"),
              tags$dt("evolutionary_conservation"), tags$dd("Cross-species conservation score")
            )
          )
        )
      ),
      
      # --- Data Page ---
      tabItem(
        tabName = "data",
        fluidRow(
          box(
            title = "Full m6A Datasets",
            width = 12, status = "primary", solidHeader = TRUE,
            p(HTML("48,055 entities of m6A with 16 features per row."))
          )
        ),
        fluidRow(
          box(
            title = "m6A Dataset Table",
            width  = 12, status = "primary", solidHeader = TRUE,
            DTOutput("m6A_table")
          )
        )
      ),
      
      # --- Visualization Page ---
      tabItem(
        tabName = "visualization",
        
        # description
        fluidRow(
          box(
            title = "Visualize m6A Datasets",
            width = 12, status = "primary", solidHeader = TRUE,
            p("Browse the m6A dataset in graphical ways:"),
            tags$ul(
              tags$li(tags$b("Pie Chart:"),  " Select one of \"m6A_status\" , \"RNA_region\" , or \"RNA_type\", and view a pie chart of the category distribution."),
              tags$li(tags$b("Histogram:"), " Use the three selection bars to filter by any combination of \"m6A_status\" , \"RNA_region\" , and \"RNA_type\" or a \"No filter\" option. Then plots a histogram of \"m6A_level\" with a density overlay.")
            ),
            p("Figures are available for downloading in both PNG and PDF formats.")
          )
        ),
        
        # plots
        fluidRow(
          # pie chart
          box(
            title = "Pie Chart",
            width = 6, status = "primary", solidHeader = TRUE,
            selectInput(
              "category_col", "Select the column to plot:",
              choices = c("m6A_status", "RNA_region", "RNA_type"),
              selected = "m6A_status"
            ),
            plotOutput("pie", height = "300px"),
            div(
              style = "margin-top:8px;",
              downloadButton("dl_pie_png", "Download (PNG)"),
              downloadButton("dl_pie_pdf", "Download (PDF)")
            )
          ),
          
          # histogram
          box(
            title = "Histogram",
            width = 6, status = "primary", solidHeader = TRUE,
            selectInput(
              "filter_m6A_status", "m6A_status:",
              choices = c("No filter", sort(unique(m6A_data$m6A_status)))
            ),
            selectInput(
              "filter_RNA_region", "RNA_region:",
              choices = c("No filter", sort(unique(m6A_data$RNA_region)))
            ),
            selectInput(
              "filter_RNA_type", "RNA_type:",
              choices = c("No filter", sort(unique(m6A_data$RNA_type)))
            ),
            plotOutput("hist", height = "300px"),
            div(
              style = "margin-top:8px;",
              downloadButton("dl_hist_png", "Download (PNG)"),
              downloadButton("dl_hist_pdf", "Download (PDF)")
            )
          )
        )
      ),
      
      # --- Download Page ---
      tabItem(
        tabName = "download",
        fluidRow(
          box(
            title = "Download m6A Datasets",
            width = 12, status = "primary", solidHeader = TRUE,
            p("Download the following CSV tables:"),
            tags$ol(
              tags$li("Total m6A dataset (unfiltered)"),
              tags$li("m6A dataset on mRNA"),
              tags$li("m6A dataset on lincRNA"),
              tags$li("m6A dataset on lncRNA"),
              tags$li("m6A dataset on pseudogene")
            )
          )
        ),
        fluidRow(
          box(
            title = "Download Buttons",
            width = 12, status = "primary", solidHeader = TRUE,
            div(
              downloadButton("dl_all",        "Download: ALL"),
              downloadButton("dl_mrna",       "Download: mRNA"),
              downloadButton("dl_lincrna",    "Download: lincRNA"),
              downloadButton("dl_lncrna",     "Download: lncRNA"),
              downloadButton("dl_pseudogene", "Download: pseudogene")
            )
          )
        )
      ),
      
      ##################################################
      # m6A prediction
      ##################################################
      # --- prediction multiple ---
      tabItem(
        tabName = "prediction_multiple",
        
        # Overview
        fluidRow(
          box(
            title = "m6A Prediction - Multiple",
            width = 12, status = "primary", solidHeader = TRUE,
            p("Upload a CSV with the required feature columns, choose a decision threshold, then generate and download predictions. An example input table is provided.")
          )
        ),
        
        # Step 1
        fluidRow(
          box(
            title = "Step 1: Upload Input Table",
            width = 12, status = "primary", solidHeader = TRUE,
            fileInput("upload_file", "Choose CSV", accept = ".csv"),
            actionButton("upload_example", "Load Example Data")
          )
        ),
        
        # Step 2
        fluidRow(
          box(
            title = "Step 2: Decision Threshold & Run",
            width = 12, status = "primary", solidHeader = TRUE,
            sliderInput(
              "positive_threshold",
              "Positive Threshold (click “Predict” again after changing this value)",
              min = 0, max = 1, value = 0.5, step = 0.01
            ),
            actionButton("predict", "Predict")
          )
        ),
        
        # Step 3
        fluidRow(
          box(
            title = "Step 3: Download Predictions",
            width = 12, status = "primary", solidHeader = TRUE,
            downloadButton("download", "Download CSV")
          )
        ),
        
        # Step 4
        fluidRow(
          box(
            title = "Step 4: Preview (predicted table after running)",
            width = 12, status = "primary", solidHeader = TRUE,
            DTOutput("table")
          )
        )
      ),
      
      # --- prediction single ---
      tabItem(
        tabName = "prediction_single",
        
        # Overview
        fluidRow(
          box(
            title = "m6A Prediction - Single",
            width = 12, status = "primary", solidHeader = TRUE,
            p("Interactively configure feature values for a single observation and obtain immediate predictions from the ML model.")
          )
        ),
        
        # Inputs / Outputs
        fluidRow(
          box(
            title = "Parameter Input",
            width = 6, status = "primary", solidHeader = TRUE,
            helpText("Set the feature values below. The model predicts automatically."),
            helpText("Predicted m6A probability (0–1) reflects support for the positive class."),
            helpText("Predicted m6A status reports Positive/Negative."),
            selectInput("RNA_type_s",   "RNA_type",   choices = c("mRNA", "lincRNA", "lncRNA", "pseudogene"), selected = "mRNA"),
            selectInput("RNA_region_s", "RNA_region", choices = c("CDS", "intron", "3'UTR", "5'UTR"),         selected = "CDS"),
            sliderInput("gc_content_s",               "gc_content",               min = 0, max = 1,  value = 0.5, step = 0.01),
            sliderInput("exon_length_s",              "exon_length",              min = 0, max = 15, value = 10,  step = 1),
            sliderInput("distance_to_junction_s",     "distance_to_junction",     min = 0, max = 15, value = 8,   step = 1),
            sliderInput("evolutionary_conservation_s","evolutionary_conservation",min = 0, max = 1,  value = 0.5, step = 0.01),
            textInput("DNA_5mer_s", "DNA_5mer (A/T/C/G, length is 5)", value = "GGACA"),
            sliderInput("positive_threshold_single",  "Decision threshold",       min = 0, max = 1,  value = 0.5, step = 0.01)
          ),
          box(
            title = "Results",
            width = 6, status = "primary", solidHeader = TRUE,
            valueBoxOutput("single_prob_box",   width = 12),
            valueBoxOutput("single_status_box", width = 12)
          )
        )
      )
    )
  )
)

####################################################################################################
# server
####################################################################################################
server <- function(input, output, session) {
  
  ##################################################
  # mini bioinformatics database
  ##################################################
  
  # --- illustrative image ---
  output$illustrative_img <- renderImage({
    list(
      src = "img_ai.png",
      contentType = "image/png"
    )
  }, deleteFile = FALSE)
  
  # --- Data Table ---
  output$m6A_table <- renderDT({
    datatable(
      m6A_data,
      options = list(scrollX = TRUE, pageLength = 15),
      rownames = FALSE
    )
  })
  
  # --- Pie Chart ---
  pie_plot <- reactive({
    pie_col <- input$category_col
    ggplot(m6A_data) +
      aes_string(1, fill = pie_col) +
      geom_bar(stat = "count") +
      coord_polar(theta = "y") +
      theme_void() +
      theme(legend.position = "right", text = element_text(face = "bold")) +
      geom_text(
        stat = "count",
        aes(label = scales::percent(after_stat(count / sum(count)), accuracy = 0.1)),
        position = position_stack(vjust = 0.5),
        size = 4,
        fontface = "bold"
      )
  })
  
  output$pie <- renderPlot(pie_plot())
  
  # download plot (png)
  output$dl_pie_png <- downloadHandler(
    filename = function() paste0("pie_chart_", input$category_col, ".png"),
    content  = function(file) {
      ggplot2::ggsave(
        filename = file, plot = pie_plot(),
        width = 6, height = 6, units = "in", dpi = 300, bg = "white"
      )
    }
  )
  
  # download plot (pdf)
  output$dl_pie_pdf <- downloadHandler(
    filename = function() paste0("pie_chart_", input$category_col, ".pdf"),
    content  = function(file) {
      ggplot2::ggsave(
        filename = file, plot = pie_plot(),
        width = 6, height = 6, units = "in", dpi = 300, bg = "white"
      )
    }
  )
  
  # --- Histogram ---
  filtered_data <- reactive({
    d <- m6A_data
    if (input$filter_m6A_status != "No filter")  d <- d[d$m6A_status == input$filter_m6A_status, ]
    if (input$filter_RNA_region  != "No filter") d <- d[d$RNA_region  == input$filter_RNA_region, ]
    if (input$filter_RNA_type    != "No filter") d <- d[d$RNA_type    == input$filter_RNA_type, ]
    d
  })
  
  hist_plot <- reactive({
    ggplot(filtered_data(), aes(m6A_level)) +
      geom_histogram(aes(y = after_stat(density)), fill = "#528bb8", alpha = 0.5, color = "white") +
      geom_density(color = "blue") +
      theme_classic() +
      labs(x = "m6A level", y = "Density") +
      theme(text = element_text(face = "bold"))
  })
  
  output$hist <- renderPlot(hist_plot())
  
  # download histogram (png)
  output$dl_hist_png <- downloadHandler(
    filename = function() {
      paste0(
        "histogram_",
        ifelse(input$filter_RNA_type == "No filter", "all", input$filter_RNA_type),
        ".png"
      )
    },
    content = function(file) {
      ggplot2::ggsave(
        filename = file, plot = hist_plot(),
        width = 7, height = 5, units = "in", dpi = 300, bg = "white"
      )
    }
  )
  
  # download histogram (pdf)
  output$dl_hist_pdf <- downloadHandler(
    filename = function() {
      paste0(
        "histogram_",
        ifelse(input$filter_RNA_type == "No filter", "all", input$filter_RNA_type),
        ".pdf"
      )
    },
    content = function(file) {
      ggplot2::ggsave(
        filename = file, plot = hist_plot(),
        width = 7, height = 5, units = "in", dpi = 300, bg = "white"
      )
    }
  )
  
  # --- Download buttons ---
  output$dl_all <- downloadHandler(
    filename = function() "m6A_all.csv",
    content  = function(file) write.csv(m6A_data, file, row.names = FALSE)
  )
  
  output$dl_mrna <- downloadHandler(
    filename = function() "m6A_mRNA.csv",
    content  = function(file) write.csv(m6A_data[m6A_data$RNA_type == "mRNA", ], file, row.names = FALSE)
  )
  
  output$dl_lincrna <- downloadHandler(
    filename = function() "m6A_lincRNA.csv",
    content  = function(file) write.csv(m6A_data[m6A_data$RNA_type == "lincRNA", ], file, row.names = FALSE)
  )
  
  output$dl_lncrna <- downloadHandler(
    filename = function() "m6A_lncRNA.csv",
    content  = function(file) write.csv(m6A_data[m6A_data$RNA_type == "lncRNA", ], file, row.names = FALSE)
  )
  
  output$dl_pseudogene <- downloadHandler(
    filename = function() "m6A_pseudogene.csv",
    content  = function(file) write.csv(m6A_data[m6A_data$RNA_type == "pseudogene", ], file, row.names = FALSE)
  )
  
  ##################################################
  # m6A prediction
  ##################################################
  
  # --- prediction multiple ---
  input_df <- reactiveVal(NULL)
  
  observeEvent(input$upload_file, {
    input_df(read.csv(input$upload_file$datapath))
  })
  
  observeEvent(input$upload_example, {
    input_df(example_data)
  })
  
  output$table <- renderDT({
    if (!is.null(predicted_df())) {
      datatable(predicted_df(), options = list(scrollX = TRUE, pageLength = 15), rownames = FALSE)
    } else {
      req(input_df())
      datatable(input_df(), options = list(scrollX = TRUE, pageLength = 15), rownames = FALSE)
    }
  })
  
  predicted_df <- reactiveVal(NULL)
  
  observeEvent(input$predict, {
    req(input_df())
    predicted_df(prediction_multiple(rf_fit, input_df(), positive_threshold = input$positive_threshold))
  })
  
  output$download <- downloadHandler(
    filename = function() "predicted.csv",
    content  = function(file) {
      write.csv(predicted_df(), file, row.names = FALSE)
    }
  )
  
  # --- prediction single ---
  single_pred <- reactive({
    prediction_single(
      rf_fit,
      gc_content                = input$gc_content_s,
      RNA_type                  = input$RNA_type_s,
      RNA_region                = input$RNA_region_s,
      exon_length               = input$exon_length_s,
      distance_to_junction      = input$distance_to_junction_s,
      evolutionary_conservation = input$evolutionary_conservation_s,
      DNA_5mer                  = toupper(trimws(input$DNA_5mer_s)),
      positive_threshold        = input$positive_threshold_single
    )
  })
  
  output$single_prob_box <- renderValueBox({
    v <- single_pred()
    valueBox(
      value    = sprintf("%.3f", as.numeric(v["predicted_m6A_prob"])),
      subtitle = "Predicted m6A probability for positive class (0 - 1)",
      icon     = icon("percentage"),
      color    = "light-blue"
    )
  })
  
  output$single_status_box <- renderValueBox({
    v <- single_pred()
    valueBox(
      value    = as.character(v["predicted_m6A_status"]),
      subtitle = "Predicted m6A status (Positive/Negative)",
      icon     = icon("check-circle"),
      color    = ifelse(as.character(v["predicted_m6A_status"]) == "Positive", "green", "yellow")
    )
  })
  
  # re-initialization
  observeEvent(input$re_initialization, {
    # --- visualization ---
    updateSelectInput(session, "category_col", selected = "m6A_status")
    updateSelectInput(session, "filter_m6A_status", selected = "No filter")
    updateSelectInput(session, "filter_RNA_region",  selected = "No filter")
    updateSelectInput(session, "filter_RNA_type",    selected = "No filter")
    
    # --- multiple ---
    input_df(NULL)
    predicted_df(NULL)
    shinyjs::reset("upload_file")
    updateSliderInput(session, "positive_threshold", value = 0.5)
    
    # --- single ---
    updateSelectInput(session, "RNA_type_s",                selected = "mRNA")
    updateSelectInput(session, "RNA_region_s",              selected = "CDS")
    updateSliderInput(session, "gc_content_s",              value = 0.5)
    updateSliderInput(session, "exon_length_s",             value = 10)
    updateSliderInput(session, "distance_to_junction_s",    value = 8)
    updateSliderInput(session, "evolutionary_conservation_s", value = 0.5)
    updateTextInput(session,  "DNA_5mer_s",                 value = "GGACA")
    updateSliderInput(session, "positive_threshold_single", value = 0.5)
    
    showNotification("Reinitialization Succeeded!", type = "message")
  })
}

# --- run app ---
shinyApp(ui, server)