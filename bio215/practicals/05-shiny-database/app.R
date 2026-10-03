library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)

m6A_data <- read.csv("../../data/m6A_imbalanced.csv")

ui <- dashboardPage(
  dashboardHeader(title = "m6A Bioinformatics Database"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Data", tabName = "data", icon = icon("table")),
      menuItem("Visualization", tabName = "visualization", icon = icon("chart-simple")),
      menuItem("Download", tabName = "download", icon = icon("download"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # --- Home Page ---
      tabItem(tabName = "home",
              fluidRow(
                box(
                  title = "About This Database",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  p("This bioinformatics database contains 48,055 m6A RNA modification entries 
                    used for computational modeling and feature-based analysis."),
                  p("Each record corresponds to one potential m6A site on a transcript, 
                    annotated with genomic coordinates, sequence composition, GC content, 
                    conservation score, and other biologically relevant attributes."),
                  p("This mini bioinformatics database allows you to browse the m6A dataset used for ML modeling."),
                  p("Use the tabs on the left to view the data table, generate simple plots, and download CSV subsets."),
                  
                  tags$h4("Included Feature Columns:"),
                  tags$ul(
                    tags$li("Genomic information: chromosome, start, end, strand, width"),
                    tags$li("m6A-related values: m6A_status, m6A_level, m6A_coverage, total_coverage"),
                    tags$li("Sequence & structure: DNA_5mer, gc_content, exon_length, distance_to_junction"),
                    tags$li("Functional annotation: RNA_type, RNA_region, evolutionary_conservation")
                  ),
                  
                  p("Use the sidebar to explore this dataset:"),
                  tags$ul(
                    tags$li("View the data in the 'Data' tab"),
                    tags$li("Generate summary plots in 'Visualization'"),
                    tags$li("Download specific subsets in 'Download'")
                  )
                )
              )
      ),
      
      # --- Data Page ---
      tabItem(tabName = "data",
              fluidRow(
                box(
                  title = "m6A Dataset Table",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  DTOutput("m6A_table")
                )
              )
      ),
      
      # --- Visualization Page ---
      tabItem(tabName = "visualization",
              fluidRow(
                box(
                  title = "Pie Chart",
                  width = 6,
                  status = "primary",
                  solidHeader = TRUE,
                  selectInput(
                    "category_col",
                    "Select the column to plot:",
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
                
                box(
                  title = "Histogram",
                  width = 6,
                  status = "primary",
                  solidHeader = TRUE,
                  selectInput(
                    "filter_m6A_status", 
                    "m6A_status:",
                    choices = c("No filter", sort(unique(m6A_data$m6A_status)))
                  ),
                  selectInput(
                    "filter_RNA_region", 
                    "RNA_region:",
                    choices = c("No filter", sort(unique(m6A_data$RNA_region)))
                  ),
                  selectInput(
                    "filter_RNA_type", 
                    "RNA_type:", 
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
            ),
            div(
              style = "display:flex; gap:10px; flex-wrap:wrap;",
              downloadButton("dl_all", "Download: ALL"),
              downloadButton("dl_mrna", "Download: mRNA"),
              downloadButton("dl_lincrna", "Download: lincRNA"),
              downloadButton("dl_lncrna", "Download: lncRNA"),
              downloadButton("dl_pseudogene", "Download: pseudogene")
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # --- Data Table ---
  output$m6A_table <- renderDT({
    datatable(m6A_data,
              options = list(scrollX = TRUE, pageLength = 15),
              rownames = FALSE)
  })
  
  # --- Pie Chart ---
  pie_plot <- reactive({
    pie_col <- input$category_col
    ggplot(m6A_data) +
      aes_string(1, fill = pie_col) +
      geom_bar(stat = "count") +
      coord_polar(theta = "y") +
      theme_void() +
      theme(
        legend.position = "right",
        text = element_text(face = "bold")
      ) +
      geom_text(
        stat = "count",
        aes(label = scales::percent(after_stat(count / sum(count)), accuracy = 0.1)),
        position = position_stack(vjust = 0.5),
        size = 4,
        fontface = "bold"
      )
  })
  
  output$pie <- renderPlot(pie_plot())
  
  output$dl_pie_png <- downloadHandler(
    filename = function() paste0("pie_chart_", input$category_col, ".png"),
    content  = function(file) {
      ggplot2::ggsave(
        filename = file, plot = pie_plot(),
        width = 6, height = 6, units = "in", dpi = 300, bg = "white"
      )
    }
  )
  
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
    if (input$filter_m6A_status != "No filter") {
      d <- d[d$m6A_status == input$filter_m6A_status, ]
    }
    if (input$filter_RNA_region != "No filter") {
      d <- d[d$RNA_region == input$filter_RNA_region, ]
    }
    if (input$filter_RNA_type != "No filter") {
      d <- d[d$RNA_type == input$filter_RNA_type, ]
    }
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
  
  output$dl_hist_png <- downloadHandler(
    filename = function() {
      paste0(
        "histogram_",
        ifelse(input$filter_RNA_type == "No filter", "all", input$filter_RNA_type),
        ".png"
      )
    },
    content  = function(file) {
      ggplot2::ggsave(
        filename = file, plot = hist_plot(),
        width = 7, height = 5, units = "in", dpi = 300, bg = "white"
      )
    }
  )
  
  output$dl_hist_pdf <- downloadHandler(
    filename = function() {
      paste0(
        "histogram_",
        ifelse(input$filter_RNA_type == "No filter", "all", input$filter_RNA_type),
        ".pdf"
      )
    },
    content  = function(file) {
      ggplot2::ggsave(
        filename = file, plot = hist_plot(),
        width = 7, height = 5, units = "in", dpi = 300, bg = "white"
      )
    }
  )
  
  # --- Download Handler ---
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
}

shinyApp(ui, server)
