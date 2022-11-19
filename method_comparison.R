library(shiny)
library(shinythemes)
library(ggplot2)
library(tidyverse)

options(shiny.port = 5005)
data <- read.csv("data/data.csv")
data_long <- data %>%
  pivot_longer(cols=c("Plink_F_value",
                      "GenABEL_F_value",
                      "Plink_with_X_F_value",
                      "Freq_F_value",
                      "Freq_rsIDs_F_value"),
               names_to="Method",
               values_to="F_value")

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  xmin <- min(x)
  xmax <- max(x)
  return(c(y=m,ymin=ymin,ymax=ymax,xmin=xmin,xmax=xmax))
}

eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b~italic(x)*","~~italic(R)^2~"="~r2,
                 list(a = as.character(format(coef(m)[1], digits = 4)),
                      b = as.character(format(coef(m)[2], digits = 4)),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}

# Define UI for application that draws a histogram
ui <- fluidPage(

   # Application title
   titlePanel("Method Comparison"),

   theme = shinytheme("slate"),

   # Sidebar with a slider input for number of bins
   sidebarLayout(
      sidebarPanel(

        radioButtons("colour","Colour of histogram",choices=c("red","green","steelblue"),selected="red"),
        radioButtons("stats", "Turn on line of best fit", choices=c("Yes", "No"), selected="No")
      ),

      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(
          tabPanel("Data", dataTableOutput("data")),
          tabPanel("Comparison Plot",plotOutput("plot")),
          tabPanel("Comparison including X-Linked SNPs", plotOutput("xplot")),
          tabPanel("Summary", verbatimTextOutput("summary"))


        )

      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    x    <- data[, 2]
   output$plot <- renderPlot({

      g <- ggplot(data, aes(x = Plink_F_value, y = GenABEL_F_value)) +
        geom_point(data = data, aes(x = Plink_F_value, y = GenABEL_F_value), color = input$colour, size = 0.2, show.legend = T) +
        xlab("Plink F Value") +
        ylab("GenABEL F Value") +
        labs(title = "Plink Inbreeding Score through Refset\nvs GenABEL Inbreeding Score") +
        theme_minimal(base_size = 20) +
        theme(#plot.title = element_text(size = 20, hjust=0.5),
              legend.text = element_text(size=10),
              legend.key.size = unit(0.5, "cm"))
     if(input$stats == "Yes"){
       g <- g +
         geom_smooth(method = "lm", color = "#440154FF") +
         annotate("text",
                 x = -0.2,
                 y = 0.225,
                 label = eq(data$Plink_F_value, data$GenABEL_F_value),
                 parse = T,
                 size = 7.5)
       g
     } else{
       g
     }
   })

   output$xplot <- renderPlot({

     g <- ggplot(data, aes(x = Plink_with_X_F_value, y = GenABEL_F_value)) +
       geom_point(data = data, aes(x = Plink_with_X_F_value, y = GenABEL_F_value), color = input$colour, size = 0.2, show.legend = T) +
       xlab("Plink with X F Value") +
       ylab("GenABEL F Value") +
       labs(title = "Plink Inbreeding Score through Refset Including\nX Chromosome vs GenABEL Inbreeding Score") +
       theme_minimal(base_size = 20) +
       theme(#plot.title = element_text(size = 20, hjust = 0.5),
             legend.text = element_text(size=10),
             legend.key.size = unit(0.5, "cm"))

     if(input$stats == "Yes"){
       g <- g +
       geom_smooth(method = "lm", color = "#440154FF") +
       annotate("text",
           x = -0.2,
           y = 0.225,
           label = eq(data$Plink_with_X_F_value, data$GenABEL_F_value),
           parse = T,
           size = 7.5)
       g
     } else{
       g
     }



   })

   output$summary <- renderPrint({
     summary(x)
   })


   output$data <- renderDataTable(data)
}

# Run the application
shinyApp(ui = ui, server = server)