# ---- app.R (single-file style) ----
library(shiny)
library(shinydashboard)
library(DT)
library(shinyjs)
suppressWarnings( library( randomForest ) )

# --- load fitted random forest model object --- 
rf_fit <- readRDS("rf_fit.rds")

# --- load examplar data CSV --- 
example_data <- read.csv("m6A_input_example.csv")

# --- define necessary functions for ML prediction --- 
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame
  
  # Complete this function by writing code at the `___`
  # 1. process input data
  feature_df <- cbind(feature_df, dna_encoding(feature_df$DNA_5mer))
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  
  # 2. predict and classify
  prob <- predict(ml_fit, feature_df, type = "prob")
  prob_binary <- ifelse(prob[, 2] <= positive_threshold, "Negative", "Positive")
  
  # 3. add rows for output
  feature_df$predicted_m6A_prob <- prob[, 2] # predicted positive-class probability
  feature_df$predicted_m6A_status <- prob_binary # either "Positive" or "Negative"
  feature_df$predicted_m6A_status <- factor(feature_df$predicted_m6A_status, levels = c("Negative", "Positive"))
  
  # 4. delete expanded one-hot DNA-5mer columns
  feature_df <- subset(feature_df, select = -c(nt_pos1, nt_pos2, nt_pos3, nt_pos4, nt_pos5))
  
  return(feature_df) #return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}

prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  
  # Complete this function by writing code at the `___`
  # 1. construct single-row data frame
  single_row_df <- data.frame(
    gc_content = as.numeric(gc_content),
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = as.numeric(exon_length),
    distance_to_junction = as.numeric(distance_to_junction),
    evolutionary_conservation = as.numeric(evolutionary_conservation),
    DNA_5mer = toupper(trimws(DNA_5mer))
  )
  
  # 2. evoke prediction_multiple() (factor convertion will be conducted automatically here)
  results_full <- prediction_multiple(ml_fit, single_row_df, positive_threshold)
  
  # 3. extract prediction results
  returned_vector <- c(results_full[1, "predicted_m6A_prob"], as.character(results_full[1, "predicted_m6A_status"]))
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  
  return(returned_vector) #return a named vector with values for predicted m6A prob and predicted m6A status
}

# --- create shiny code bellow ---
# --- ui ---
ui <- dashboardPage(
  # page header
  dashboardHeader(title = "m6A Prediction"),
  
  # sidebar
  dashboardSidebar(
    # a button for re-initialization
    actionButton("re_initialization", "Reinitialize Everything", icon = icon("refresh")),
    
    sidebarMenu(
      menuItem("m6A Prediction Multiple", tabName = "prediction_multiple", icon = icon("inbox")),
      menuItem("m6A Prediction Single",   tabName = "prediction_single",     icon = icon("comment"))
      )
  ),
  
  # page body
  dashboardBody(
    useShinyjs(),
    tabItems(
      # --- prediction multiple ---
      tabItem(
        tabName = "prediction_multiple",
        fluidRow(
          # Section 0 - header
          box(width = 12, title = "Table Input (multiple observation)"),
          
          # Section 1 - files
          box(width = 12, title = "Step 1. Upload a CSV file containing input features / Use the example data as a demo",
              # upload data (filer up-loader)
              fileInput("upload_file", "Choose CSV", accept = ".csv"),
              
              # example data (button)
              actionButton("upload_example", "Load Example Data")),
          
          # Section 2 - predict
          box(width = 12, title = "Step 2. Select a threshold for the positive class probability & Start prediction",
              # select threshold (slider)
              sliderInput("positive_threshold", "Positive Threshold (please click \"Predict\" button again after changing the threshold)", min = 0, max = 1, value = 0.5, step = 0.01),
              
              # predict (button)
              actionButton("predict", "Predict")),
          
          # Section 3 - download
          box(width = 12, title = "Step 3. Download the predicted data into a CSV file",
              # download (button)
              downloadButton("download", "Download")),
          
          # Section 4 - tables
          box(width = 12, title = "Step 4. View the data (after prediction, the original table will be altered to predicted table)",
              #  show predicted table
              DTOutput("table"))
        )
      ),
      
      # --- prediction single ---
      tabItem(
        tabName = "prediction_single",
        fluidRow(
          # Section 0 - header
          box(width = 12, title = "Manual Input (single observation)"),
          
          # left side - input
          box(width = 6, title = "Parameters Input",
              helpText("Please set the feature values below."),
              helpText("The model will automatically and instantly conduct prediction."),
              helpText("Predicted m6A probability (on a 0 - 1 scale) indicates the model's support for calssifying the entity into positive class."),
              helpText("Predicted m6A status indicates whether this site is m6A (Positive/Negative)."),
              
              selectInput("RNA_type_s",  "RNA_type",  choices = c("mRNA","lincRNA","lncRNA","pseudogene"), selected = "mRNA"),
              selectInput("RNA_region_s","RNA_region",choices = c("CDS","intron","3'UTR","5'UTR"),          selected = "CDS"),
              sliderInput("gc_content_s",             "gc_content",               min = 0,  max = 1,  value = 0.5, step = 0.01),
              sliderInput("exon_length_s",            "exon_length",              min = 0,  max = 15, value = 10,  step = 1),
              sliderInput("distance_to_junction_s",   "distance_to_junction",     min = 0,  max = 15, value = 8,   step = 1),
              sliderInput("evolutionary_conservation_s","evolutionary_conservation", min = 0, max = 1, value = 0.5, step = 0.01),
              textInput("DNA_5mer_s", "DNA_5mer (A/T/C/G, length is 5)", value = "GGACA"),
              sliderInput("positive_threshold_single","Decision threshold", min = 0, max = 1, value = 0.5, step = 0.01)
          ),
          
          # right side - output
          box(width = 6, title = "Results Output (auto-freshed after changing input)",
              valueBoxOutput("single_prob_box", width = 12),
              valueBoxOutput("single_status_box", width = 12)
          )
        )
      )
    )
  )
)

# --- server ---
server <- function(input, output, session) {
  # --- prediction multiple ---
  # decide input data
  input_df <- reactiveVal(NULL)
  
  # upload data
  observeEvent(input$upload_file, {
    input_df(read.csv(input$upload_file$datapath))
  })
  
  # upload example data
  observeEvent(input$upload_example, {
    input_df(example_data)
  })
  
  # show table
  output$table <- renderDT({
    if (!is.null(predicted_df())) {
      datatable(predicted_df(),
                options = list(scrollX = TRUE, pageLength = 15),
                rownames = FALSE)
    } else {
      req(input_df())
      datatable(input_df(),
                options = list(scrollX = TRUE, pageLength = 15),
                rownames = FALSE)
    }
  })
  
  # predict
  predicted_df <- reactiveVal(NULL)
  observeEvent(input$predict, {
    req(input_df())
    predicted_df(
      prediction_multiple(rf_fit, input_df(), positive_threshold = input$positive_threshold)
    )
  })
  
  # download table
  output$download <- downloadHandler(
    filename = function() "predicted.csv",
    content = function(file) {
      write.csv(predicted_df(), file, row.names = FALSE)
  })
  
  # --- prediction single ---
  single_pred <- reactive({
    prediction_single(
      rf_fit,
      gc_content = input$gc_content_s,
      RNA_type = input$RNA_type_s,
      RNA_region = input$RNA_region_s,
      exon_length = input$exon_length_s,
      distance_to_junction = input$distance_to_junction_s,
      evolutionary_conservation = input$evolutionary_conservation_s,
      DNA_5mer = toupper(trimws(input$DNA_5mer_s)),
      positive_threshold = input$positive_threshold_single
    )
  })
  
  output$single_prob_box <- renderValueBox({
    v <- single_pred()
    valueBox(
      value = sprintf("%.3f", as.numeric(v["predicted_m6A_prob"])),
      subtitle = "Predicted m6A probability for positive class (0 - 1)",
      icon = icon("percentage"),
      color = "light-blue"
    )
  })
  
  output$single_status_box <- renderValueBox({
    v <- single_pred()
    valueBox(
      value = as.character(v["predicted_m6A_status"]),
      subtitle = "Predicted m6A status (Positive/Negative)",
      icon = icon("check-circle"),
      color = ifelse(as.character(v["predicted_m6A_status"]) == "Positive", "green", "yellow")
    )
  })
  
  # re-initialization
  # reinitialize everything
  observeEvent(input$re_initialization, {
    # --- multiple ---
    input_df(NULL)
    predicted_df(NULL)
    shinyjs::reset("upload_file")
    updateSliderInput(session, "positive_threshold", value = 0.5)
    
    # --- single ---
    updateSelectInput(session, "RNA_type_s",  selected = "mRNA")
    updateSelectInput(session, "RNA_region_s", selected = "CDS")
    updateSliderInput(session, "gc_content_s", value = 0.5)
    updateSliderInput(session, "exon_length_s", value = 10)
    updateSliderInput(session, "distance_to_junction_s", value = 8)
    updateSliderInput(session, "evolutionary_conservation_s", value = 0.5)
    updateTextInput(session, "DNA_5mer_s", value = "GGACA")
    updateSliderInput(session, "positive_threshold_single", value = 0.5)
    showNotification("Reinitialization Succeeded!", type = "message")
  })
}

# --- run app ---
shinyApp(ui, server)
