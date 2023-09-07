library(shiny)
library(BIRDIE)
library(bayesplot)
library(dplyr)
library(ggplot2)
library(ABAP)
library(spOccupancy)

ui <- fluidPage(

    titlePanel("Occupancy model diagnostics"),

    sidebarLayout(

        sidebarPanel(
            numericInput('sp_code', 'SAFRING code', min = 1, max = 10000, value = 6, step = 1),
            numericInput('year', 'Last modelled year', min = 2008, max = 2021, value = 2008, step = 1),

            actionButton("load", "Load data"),
            br(), p(""), br(),

            textInput('param', 'Select parameter'),

            actionButton("summary", "Summary"), br(),
            actionButton("trace", "Trace plots"), br(),
            actionButton("ppc", "Posterior checks"),
            actionButton("map", "Problem pentads"),
            width = 3),

        mainPanel(
            mainPanel(
                tabsetPanel(
                    tabPanel("Summary", verbatimTextOutput('summary')),
                    tabPanel("Trace plots", plotOutput('tracePlot')),
                    tabPanel("PPC", plotOutput('ppcPlot')),
                    tabPanel("Map", plotOutput('ppcMapPlot'))
                ),
                width = 15)
        )
    )
)


# Define server logic required to draw a histogram
server <- function(input, output) {

    # Prepare config file -----------------------------------------------------

    set_config <- function(input){
        config <- configPipeline(year = input$year, dur = 1,
                                 occ_mod = c("log_dist_coast", "elev", "hum.km2", "wetcon",
                                             "watrec", "watext", "log_watext", "watext:watrec",
                                             "ndvi", "prcp", "tdiff"),
                                 det_mod = c("( 1|obs_id)", "(1|site_id)", "log_hours", "prcp", "tdiff", "cwac"),
                                 fixed_vars = c("Pentad", "lon", "lat", "watocc_ever", "wetext_2018","wetcon_2018",
                                                "dist_coast", "elev"),
                                 package = "spOccupancy",
                                 out_dir = "analysis/hpc/imports",
                                 server = FALSE)
    }

    # Define low-level functions -------------------------------------------

    # Load fit
    get_fit <- function(input, config){
        readRDS(paste0("../../../",setSpOutFilePath("occu_fit", config, input$year, input$sp_code, ".rds")))
    }

    get_ppc <- function(input, config){
        readRDS(paste0("../../../",setSpOutFilePath("occu_ppc", config, input$year, input$sp_code, ".rds")))
    }


    # Create data loading functions --------------------------------------------

    load_fit <- reactive({
        message("Loading model fit...")
        config <- set_config(input)
        get_fit(input, config)

    }) %>%
        bindEvent(input$load)

    load_ppc <- reactive({
        message("Loading posterior predictive checks...")
        config <- set_config(input)
        get_ppc(input, config)

    }) %>%
        bindEvent(input$load)


    # Plot -------------------------------------------------------------------

    # Summary
    observe({
        output$summary <- renderPrint({

            # Load data
            fit <- load_fit()

            # Get summary
            summary(fit)

        })
    }) %>% bindEvent(input$summary)


    # Trace plots
    observe({
        output$tracePlot <- renderPlot({

            # Load fit
            fit <- load_fit()

            message("Making trace plots...")

            # All spOccupancy sample object have the trailing ".samples"
            param <- paste0(input$param, ".samples")

            # Plot mixing for those specific sites/years
            mcmc_trace(fit[param])

            }, height = 700)
    }) %>% bindEvent(input$trace)

    # Posterior predictive plots
    observe({
        output$ppcPlot <- renderPlot({

            # Load data
            ppc_out <- load_ppc()

            ppc.df <- data.frame(fit = ppc_out$fit.y,
                                 fit.rep = ppc_out$fit.y.rep,
                                 color = 'lightskyblue1')

            ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'

            plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21,
                 ylab = 'Fit', xlab = 'True')
            lines(ppc.df$fit, ppc.df$fit, col = 'black')

        }, height = 700)
    }) %>% bindEvent(input$ppc)

    # Problem pentad plots
    observe({
        output$ppcMapPlot <- renderPlot({

            # Load data
            fit <- load_fit()
            ppc_out <- load_ppc()

            ppc_df <- rbind(as.data.frame(t(ppc_out$fit.y.group.quants[c(1, 3, 5),])) %>%
                                mutate(obs_data = 1,
                                       pentad = rownames(fit$y)),
                            as.data.frame(t(ppc_out$fit.y.rep.group.quants[c(1, 3, 5),])) %>%
                                mutate(obs_data = 0,
                                       pentad = rownames(fit$y))) %>%
                tidyr::pivot_longer(cols = -c(obs_data, pentad), names_to = "quant", values_to = "gof")

            problem_pts <- ppc_df %>%
                group_by(pentad, quant) %>%
                summarise(diff_gof = diff(gof)) %>%
                filter(quant == "50%") %>%
                filter(abs(diff_gof) > 0.1) %>%
                pull(pentad)

            # Download South African ABAP pentads
            sites <- readRDS("pentads_sa.rds")
            sa_map <- readRDS("provinces_sa.rds")

            sites %>%
                filter(pentad %in% unique(rownames(fit$y))) %>%
                mutate(problem = if_else(pentad %in% problem_pts, 1, 0),
                       pentad_id = as.numeric(factor(pentad))) %>%
                ggplot() +
                geom_sf(data = sa_map) +
                geom_sf(aes(fill = factor(problem)), lwd = 0.01) +
                scale_fill_viridis_d(name = "Problem", direction = -1)

        }, height = 700)
    }) %>% bindEvent(input$map)


}


# Run the application
shinyApp(ui = ui, server = server)
