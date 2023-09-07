library(shiny)
library(BIRDIE)
library(bayesplot)
library(dplyr)
library(ggplot2)

source("diagnoseJAGSrhat.R")

ui <- fluidPage(

    titlePanel("Abundance state-space model diagnostics"),

    sidebarLayout(

        sidebarPanel(
            numericInput('sp_code', 'SAFRING code', min = 1, max = 10000, value = 87, step = 1),
            numericInput('year', 'Last modelled year', min = 2021, max = 2021, value = 2021, step = 1),

            actionButton("load", "Load data"),
            br(), p(""), br(),

            numericInput('site', 'Site index (0 = all)', min = 0, max = 800, value = 0, step = 1),
            textInput('param', 'Select parameter'),

            actionButton("conv", "Rhat"), br(),
            actionButton("trace", "Trace plots"), br(),
            actionButton("ppc", "Posterior checks"),
            width = 3),

        mainPanel(
            mainPanel(
                tabsetPanel(
                    tabPanel("Rhat", plotOutput('rhat')),
                    tabPanel("Trace plots", plotOutput('tracePlot')),
                    tabPanel("PPC", plotOutput('ppcPlot'))
                ),
                width = 15)
        )
    )
)


# Define server logic required to draw a histogram
server <- function(input, output) {

    # Prepare config file -----------------------------------------------------

    set_config <- function(input){
        config <- configPipeline(
            year = 2021,
            dur = 29,
            mod_file = "cwac_ssm_two_season_mean_rev.R",
            package = "jagsUI",
            data_dir = NULL,
            out_dir = NULL, #"analysis/downloads",
            server = FALSE
        )
    }

    # Define low-level functions -------------------------------------------

    # Load fit
    get_fit <- function(input, config, raw = FALSE){
        fit <- readRDS(paste0("../../../",setSpOutFilePath("ssm_fit", config, config$years_ch, input$sp_code, ".rds")))
        if(!raw){
            # compute stats
            fit <- BIRDIE:::processJAGSoutput(fit, DIC = FALSE, params.omit = NULL)
        }
        fit
    }


    # Load data
    get_counts <- function(input, config){
        read.csv(paste0("../../../", setSpOutFilePath("abu_model_data", config, config$years_ch, input$sp_code, ".csv")))
    }


    # Create data loading functions --------------------------------------------

    load_fit <- reactive({
        message("Loading model fit...")
        config <- set_config(input)
        get_fit(input, config, raw = TRUE)

    }) %>%
        bindEvent(input$load)

    load_fit_stats <- reactive({
        message("Loading model fit...")
        config <- set_config(input)
        get_fit(input, config, raw = FALSE)

    }) %>%
        bindEvent(input$load)

    subset_fit <- reactive({
        if(input$site != 0){
            paste0(input$param, "[", input$site)
        } else {
            input$param
        }
    }) %>%
        bindEvent(input$trace)

    load_counts <- reactive({
        message("Loading counts...")
        config <- set_config(input)
        get_counts(input, config)
    }) %>%
        bindEvent(input$load)



    # Prepare model fit -------------------------------------------------------

    # observe({
    #     # Load data
    #     raw_fit <- get_fit(input, config)
    #     fit <- load_fit()
    # }) %>%
    #     bindEvent(input$load)


    # Plot -------------------------------------------------------------------

    # Rhat plots
    observe({
        output$rhat <- renderPlot({

            # Load data
            fit <- load_fit_stats()

            # Plot rhat values for input parameters
            diagnoseJAGSrhat(fit, input$param)

        }, height = 700)
    }) %>% bindEvent(input$conv)


    # Trace plots
    observe({
        output$tracePlot <- renderPlot({

            # Load data
            fit_raw <- load_fit()

            # Subset parameters if necessary
            site_param <- subset_fit()

            message("Making trace plots...")

            # Plot mixing for those specific sites/years
            mcmc_trace(fit_raw, pars = vars(contains(site_param)))

            }, height = 700)
    }) %>% bindEvent(input$trace)

    # Posterior predictive plots
    observe({
        output$ppcPlot <- renderPlot({

            if(input$site == 0){
                stop("A site index must be specified for posterior predictive checks")
            }

            # Load data
            fit_stats <- load_fit_stats()

            # # compute stats
            # fit_stats <- BIRDIE:::processJAGSoutput(fit, DIC = FALSE, params.omit = NULL)
#
#             # Load counts
             counts <- load_counts()

            # Get posterior predictive distribution
            post_sims <- postPredDistJagsSsm(fit = fit_stats,
                                             data = counts,
                                             obs_error = TRUE,
                                             500)

            # Subset sites
            post_sims <- post_sims %>%
                filter(site_id == input$site) %>%
                mutate(obs = factor(obs, levels = c("1", "0")),
                       iter = factor(iter, levels = unique(iter))) %>%
                arrange(obs)

            post_sims %>%
                ggplot() +
                geom_line(aes(x = year, y = obs_sim, group = iter, alpha = obs, col = obs)) +
                scale_alpha_manual(name = "", values = c(1, 0.1), labels = c("data", "sims")) +
                scale_colour_manual(name = "", values = c('1' = "black", '0' = "red"), labels = c("data", "sims")) +
                xlab("Year") + ylab("log counts") +
                facet_grid("season")

        }, height = 700)
    }) %>% bindEvent(input$ppc)


}


# Run the application
shinyApp(ui = ui, server = server)
