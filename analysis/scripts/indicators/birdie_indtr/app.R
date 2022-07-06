library(shiny)
library(dplyr)
library(ggplot2)
library(DiagrammeR)

source("createGraphDf.R")
source("createGraph.R")
source("createConditionGraph.R")

elem_list <- list(taxon = c("species", "sp_type"),
                  geo = c("site", "site_type", "global"),
                  time = "time",
                  cond = c("at", "across"),
                  indicator = c("distr", "abund", "div"))


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Indicator queries"),

    selectInput(inputId = "target", label = "What are you interested in:",
                choices = list("species", "sp_type", "site", "site_type"),
                multiple = FALSE),

    column(width = 3,
           selectInput(inputId = "cond_sp", label = "Select species condition:",
                       choices = list("at", "across", "NA"),
                       selected = "NA",
                       multiple = FALSE),
           selectInput(inputId = "cond_st", label = "Select site condition:",
                       choices = list("at", "across", "NA"),
                       selected = "NA",
                       multiple = FALSE),
           selectInput(inputId = "cond_tm", label = "Select time condition:",
                       choices = list("at", "across", "NA"),
                       selected = "NA",
                       multiple = FALSE)),

    column(width = 3,
           selectInput(inputId = "filter_sp", label = "Select species filter:",
                       choices = list("species", "sp_type", "NA"),
                       selected = "NA",
                       multiple = FALSE),
           selectInput(inputId = "filter_st", label = "Select site filter:",
                       choices = list("site", "site_type", "NA"),
                       selected = "NA",
                       multiple = FALSE),
           selectInput(inputId = "filter_tm", label = "Select time filter:",
                       choices = list("time", "NA"),
                       selected = "NA",
                       multiple = FALSE)),

    actionButton("go", "Go"),

    htmlOutput("indtrPlot")

)

# Define server logic required to draw a histogram
server <- function(input, output) {

    c1 <- eventReactive(input$go, {input$cond_sp})
    c2 <- eventReactive(input$go, {input$cond_st})
    c3 <- eventReactive(input$go, {input$cond_tm})
    f1 <- eventReactive(input$go, {input$filter_sp})
    f2 <- eventReactive(input$go, {input$filter_st})
    f3 <- eventReactive(input$go, {input$filter_tm})

    output$indtrPlot <- renderUI({

        cond1 <- c1()
        cond2 <- c2()
        cond3 <- c3()
        filter1 <- f1()
        filter2 <- f2()
        filter3 <- f3()

        cond1[cond1 == "NA"] <- NA
        cond2[cond2 == "NA"] <- NA
        cond3[cond3 == "NA"] <- NA

        filter1[filter1 == "NA"] <- NA
        filter2[filter2 == "NA"] <- NA
        filter3[filter3 == "NA"] <- NA

        # Order conditions
        if(input$target == "species"){
            cond <- c(cond2, cond3)
            filt <- c(filter2, filter3)
        } else if(input$target == "sp_type"){
            cond <- c(cond1, cond2, cond3)
            filt <- c(filter1, filter2, filter3)
        } else if(input$target == "site"){
            cond <- c(cond1, cond3)
            filt <- c(filter1, filter3)
        } else if(input$target == "site_type"){
            cond <- c(cond2, cond1, cond3)
            filt <- c(filter2, filter1, filter3)
        }


        df <- createGraphDf(elem_list, input$target)
        # gr <- createGraph(df)
        # render_graph(gr, width = 700)
        render_graph(
            createConditionGraph(df, cond, filt, output = "graph"),
            width = 700)

    })

}

# Run the application
shinyApp(ui = ui, server = server)
