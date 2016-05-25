library(datasets)
# http://shiny.rstudio.com/articles/basics.html

# this shiny app will take DNA sequences as input and output their reverse compliments
# need to take text input, need to sanitize text input
# make sure to do reverse compliments line wise, not entire string-wise 

# example code used as a placeholder
server <- function(input, output) {
  
  # Return the requested dataset
  datasetInput <- reactive({
    switch(input$dataset,
           "rock" = rock,
           "pressure" = pressure,
           "cars" = cars)
  })
  
  # Generate a summary of the dataset
  output$summary <- renderPrint({
    dataset <- datasetInput()
    summary(dataset)
  })
  
  # Show the first "n" observations
  output$view <- renderTable({
    head(datasetInput(), n = input$obs)
  })
}

ui <- fluidPage(
  
  # Application title
  titlePanel("Shiny Text"),
  
  # Sidebar with controls to select a dataset and specify the
  # number of observations to view
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Choose a dataset:", 
                  choices = c("rock", "pressure", "cars")),
      
      numericInput("obs", "Number of observations to view:", 10)
    ),
    
    # Show a summary of the dataset and an HTML table with the 
    # requested number of observations
    mainPanel(
      verbatimTextOutput("summary"),
      
      tableOutput("view")
    )
  )
)

shinyApp(ui = ui, server = server)
