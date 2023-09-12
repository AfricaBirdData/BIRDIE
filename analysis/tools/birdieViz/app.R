require(shiny)
require(dplyr, warn.conflicts = FALSE)
require(sf)
require(leaflet)
require(ggplot2)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("BIRDIE shiny dashboard"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(

            checkboxGroupInput(inputId = "indtr",
                               label = "Indicator",
                               choices = c("Occupancy","Abundance"),
                               inline = FALSE,
                               selected = "Occupancy"),
            selectInput(inputId = "sp",
                        label = "Species",
                        choices = list(4, 6),
                        selected = 4),
            actionButton("load", "Load data"),
            selectInput(inputId = "time",
                        label = "Year",
                        choices = list(2008, 2009, 2010),
                        selected = 2008),
            selectInput(inputId = "site",
                        label = "Site",
                        choices = list("ALL", "26352535"),
                        selected = "26352535"),
            actionButton("do", "Go!"),

            width = 2),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("showPlot", "100vh"),

            width = 10)
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {


    # Define data-loading functions -------------------------------------------

    # Occupancy
    loadOccu <- function(input){

        url <- "https://dataviz.naturalsciences.be/application/openapi/v1/occupancy"

        # Extract data
        message("Downloading your data...")

        json_data <- sf::read_sf(url)

        gm <- json_data %>%
            dplyr::select(pentad_name, geometry)

        json_data <- json_data %>%
            sf::st_drop_geometry() %>%
            readr::type_convert(col_types = readr::cols(
                .default = readr::col_double(),
                id = readr::col_integer(),
                pentad_name = readr::col_character()))

        json_data %>%
            dplyr::left_join(gm, by = "pentad_name") %>%
            sf::st_sf()

    }


    # Abundance
    loadAbund <- function(input){

        url <- paste0("https://dataviz.naturalsciences.be/application/openapi/v1/season_count/ALL/", input$sp, "/SUMMER/2010/2020" )

        # Extract data
        message("Downloading your data...")

        json_data <- sf::read_sf(url)

        json_data %>%
            sf::st_drop_geometry() %>%
            readr::type_convert(col_types = readr::cols(
                .default = readr::col_double(),
                id = readr::col_integer(),
                site_goup_name = readr::col_character(),
                site_name = readr::col_character()))

    }


    # Create data loading function --------------------------------------------

    load_data <- reactive({

        if(input$indtr == "Occupancy"){

            json_data <- loadOccu(input)

        }

        if(input$indtr == "Abundance"){

            json_data <- loadAbund(input)

        }

    }) %>%
        bindEvent(input$load)



    # Define plotting functions -----------------------------------------------

    make_occu_plot <- reactive({

        # Load data
        gen_data <- load_data()

        message("Plotting, almost there...")


        var_name_pre <- "year"
        var_name <- paste0(var_name_pre, input$time)

        gen_data %>%
            dplyr::select(id, pentad_name, layer = dplyr::all_of(var_name)) %>%
            ggplot() +
            geom_sf(aes(fill = layer), lwd = NA) +
            scale_fill_viridis_c() +
            ggtitle(paste(input$indtr, input$sp, input$time)) +
            theme(plot.title = element_text(size = 22))

    }) %>%
        bindEvent(input$do)

    make_abund_plot <- reactive({

        # Load data
        gen_data <- load_data()

        message("Plotting, almost there...")

        if(input$site != "ALL"){

            gen_data <- gen_data %>%
                dplyr::filter(id == as.integer(input$site))

        }

        gen_data %>%
            dplyr::select(-c(site_goup_name, site_name)) %>%
            tidyr::pivot_longer(-id, names_to = "year", values_to = "vals") %>%
            dplyr::mutate(year = as.integer(gsub("year_", "", year))) %>%
            dplyr::group_by(year) %>%
            dplyr::summarise(vals = sum(vals, na.rm = TRUE)) %>%
            dplyr::ungroup() %>%
            ggplot(aes(x = year, y = vals)) +
            geom_point() +
            geom_line() +
            scale_fill_viridis_c() +
            ggtitle(paste(input$indtr, input$sp)) +
            theme(plot.title = element_text(size = 22))

    }) %>%
        bindEvent(input$do)


    # Plot -------------------------------------------------------------------

    observe({
        output$showPlot <- renderPlot({

            if(input$indtr == "Occupancy"){

                make_occu_plot()

            }

            if(input$indtr == "Abundance"){

                make_abund_plot()

            }

        }, height = 700)
    }) %>% bindEvent(input$do)


}


# Run the application
shinyApp(ui = ui, server = server)
