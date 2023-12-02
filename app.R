#clear environment and load application
quickcode::clean()

libraryAll(shiny,shinythemes,dplyr,magrittr,tidyr,ggplot2)


#BF function based on Dienes & McLatchie, 2018
##modified so that H1 is represented by normal distribution (rather than t), hence there is no 'dftheory' argument
Bf<-function(sd, obtained, dfdata, meanoftheory, sdtheory, tail = 2)
{
  area <- 0
  normarea <- 0
  theta <- meanoftheory - 5 * sdtheory
  incr <- sdtheory / 200
  for (A in -1000:1000){
    theta <- theta + incr
    dist_theta <- dnorm((theta-meanoftheory)/sdtheory)
    if(identical(tail, 1)){
      if (theta <= 0){
        dist_theta <- 0
      } else {
        dist_theta <- dist_theta * 2
      }
    }
    height <- dist_theta * dt((obtained-theta)/sd, df = dfdata)
    area <- area + height * incr
    normarea <- normarea + dist_theta*incr
  }
  LikelihoodTheory <- area/normarea
  Likelihoodnull <- dt(obtained/sd, df = dfdata)
  BayesFactor <- LikelihoodTheory / Likelihoodnull
  BayesFactor
}


# Define UI for application 
ui = fluidPage(theme = shinytheme("simplex"),
               title = "The role of Bayes factors in testing interactions",
               titlePanel("The role of Bayes factors in testing interactions"),
               tabsetPanel(type = 'tabs', id = 'tabset',
                           
                           #Between groups tab
                           tabPanel('Between groups' , value = 'between',
                                    
                                    #Data input
                                    sidebarPanel(width = 6,
                                                 h3("Data"),
                                                 hr(),
                                                 h4("Group 1 (Treatment)"),
                                                 numericInput('m1_b', 'Raw effect size', 0),
                                                 numericInput('se1_b', 'SE', 1),
                                                 numericInput('n1_b', 'Sample Size', 10),
                                                 h4("Group 2 (Control)"),
                                                 numericInput('m2_b', 'Raw effect size', 0),
                                                 numericInput('se2_b', 'SE', 1),
                                                 numericInput('n2_b', 'Sample Size', 10),
                                                 h4("Model of H1 (Group 1 > Group 2)"),
                                                 numericInput('sdH1_b', 'SD of the half-normal distribution in the same units as the raw effect sizes', 1),
                                                 
                                                 actionButton("action_between", "Compute", class = "btn-primary", width = 100),
                                                 actionButton("reset_input_between", "Reset", width = 60)
                                    ),
                                    # Reporting results as a table
                                    mainPanel(width = 6,
                                              h3("Results"),
                                              hr(),
                                              tableOutput('table_between') 
                                    )
                                    ),
                           
                           #Mixed design tab
                           tabPanel('Mixed design' , value = 'mixed',
                                    
                                    #Data input
                                    sidebarPanel(width = 6,
                                                 h3("Data"),
                                                 hr(),
                                                 h4("Group 1 (Treatment)"),
                                                 numericInput('m1_m', 'Raw effect size', 0),
                                                 numericInput('se1_m', 'SE', 1),
                                                 numericInput('n1_m', 'Sample Size', 10),
                                                 h4("Group 2 (Control)"),
                                                 numericInput('m2_m', 'Raw effect size', 0),
                                                 numericInput('se2_m', 'SE', 1),
                                                 numericInput('n2_m', 'Sample Size', 10),
                                                 h4("Model of H1 (Group 1 > Group 2)"),
                                                 numericInput('sdH1_m', 'SD of the half-normal distribution in the same units as the raw effect sizes', 1),
                                                 
                                                 actionButton("action_mixed", "Compute", class = "btn-primary", width = 100),
                                                 actionButton("reset_input_mixed", "Reset", width = 60)
                                    ),
                                    # Reporting results as a table
                                    mainPanel(width = 6,
                                              h3("Results"),
                                              hr(),
                                              tableOutput('table_mixed') 
                                    )
                                    ),
                           
                           
                           # Within subjects tab
                           tabPanel('Within subjects' , value = 'within',
                                    
                                    #Data input
                                    sidebarPanel(width = 6,
                                                 h3("Data"),
                                                 hr(),
                                                 h4("Condition 1 (Treatment)"),
                                                 numericInput('m1_w', 'Raw effect size', 0),
                                                 numericInput('se1_w', 'SE', 1),
                                                 h4("Condition 2 (Control)"),
                                                 numericInput('m2_w', 'Raw effect size', 0),
                                                 numericInput('se2_w', 'SE', 1),
                                                 h4("Difference between the conditions"),
                                                 numericInput('m_diff_w', 'Raw effect size (Condition 1 - Condition 2)', 0),
                                                 numericInput('se_diff_w', 'SE', 1),
                                                 h4("Sample size"),
                                                 numericInput('n_all_w', 'Sample Size', 10),
                                                 h4("Model of H1 (Condition 1 > Condition 2)"),
                                                 numericInput('sdH1_w', 'SD of the half-normal distribution in the same units as the raw effect sizes', 1),
                                                 
                                                 actionButton("action_within", "Compute", class = "btn-primary", width = 100),
                                                 actionButton("reset_input_within", "Reset", width = 60)
                           ),
                           
                           # Reporting results as a table
                           mainPanel(width = 6,
                                     h3("Results"),
                                     hr(),
                                     tableOutput('table_within') 
                           )
                 
               ),
               
               # Reference information at the bottom of the screen
               fluidRow(
                 column(width=11, offset=.5, h5("This SinyApp is part of a tutorial by Palfi & Dienes (2019). The pre-print of the manuscript can be accessed",
                                                a('here.', 
                                                  href='https://psyarxiv.com/qjrg4'))),
                 column(width=11, offset=.5, h5("The calculation of the Bayes factor is based on the R script of Dienes & Mclatchie (2018). For more information on how to specify the predictions of H1, check out",
                                                a('Dienes (2019).', 
                                                  href='https://psyarxiv.com/yqaj4/'))),
                 column(width=11, offset=.5, h5('Application created by Bence Palfi', a('(@bence_palfi)', href="https://twitter.com/bence_palfi",
                                                                                        target="_blank"))))
               )
               )

# Define server logic 
server <- function(input, output, session) {
    
  #Calculating the Bayes factor for between design
  observeEvent(input$action_between, {
    #calculate content of the Table
    dfdata_1 <- (input$n1_b) - 2
    dfdata_2 <- (input$n2_b) - 2
    dfdata_interaction <- (input$n1_b + input$n2_b) - 4

    t_1 <- input$m1_b/input$se1_b
    t_2 <- input$m2_b/input$se2_b
    t_interaction <- (input$m1_b - input$m2_b) / sqrt(input$se1_b^2 + input$se2_b^2)

    B_1 <- Bf(sd = input$se1_b, obtained = input$m1_b, dfdata = dfdata_1, meanoftheory = 0, sdtheory = input$sdH1_b, tail = 1)
    B_2 <- Bf(sd = input$se2_b, obtained = input$m2_b, dfdata = dfdata_2, meanoftheory = 0, sdtheory = input$sdH1_b, tail = 1)
    B_interaction <- Bf(sd = sqrt(input$se1_b^2 + input$se2_b^2), obtained = (input$m1_b - input$m2_b), dfdata = dfdata_interaction, meanoftheory = 0, sdtheory = input$sdH1_b, tail = 1)

    p_1 <- 2*pt(-abs(t_1),df=dfdata_1)
    p_2 <- 2*pt(-abs(t_2),df=dfdata_2)
    p_interaction <- 2*pt(-abs(t_interaction),df=dfdata_interaction)

    #columns of the Results Table
    df     <- c(dfdata_1, dfdata_2, dfdata_interaction)
    t      <- c(t_1, t_2, t_interaction)
    B      <- c(B_1,B_2,B_interaction)
    p      <- c(p_1,p_2,p_interaction)

    #Putting the Table together
    test  <- c('Group 1', 'Group 2', 'Interaction')
    results_between<- data.frame(test, df, t, B, p)

    #Output tables
    output$table_between <- renderTable(results_between, spacing = "l", width = 300, digits = 3)

  })
  
  #Calculating the Bayes factor for mixed design
  observeEvent(input$action_mixed, {
    #calculate content of the Table
    dfdata_1 <- (input$n1_m) - 1
    dfdata_2 <- (input$n2_m) - 1
    dfdata_interaction <- (input$n1_m + input$n2_m) - 2
    
    t_1 <- input$m1_m/input$se1_m
    t_2 <- input$m2_m/input$se2_m
    t_interaction <- (input$m1_m - input$m2_m) / sqrt(input$se1_m^2 + input$se2_m^2)
    
    B_1 <- Bf(sd = input$se1_m, obtained = input$m1_m, dfdata = dfdata_1, meanoftheory = 0, sdtheory = input$sdH1_m, tail = 1)
    B_2 <- Bf(sd = input$se2_m, obtained = input$m2_m, dfdata = dfdata_2, meanoftheory = 0, sdtheory = input$sdH1_m, tail = 1)
    B_interaction <- Bf(sd = sqrt(input$se1_m^2 + input$se2_m^2), obtained = (input$m1_m - input$m2_m), dfdata = dfdata_interaction, meanoftheory = 0, sdtheory = input$sdH1_m, tail = 1)
    
    p_1 <- 2*pt(-abs(t_1),df=dfdata_1)
    p_2 <- 2*pt(-abs(t_2),df=dfdata_2)
    p_interaction <- 2*pt(-abs(t_interaction),df=dfdata_interaction)
    
    #columns of the Results Table
    df     <- c(dfdata_1, dfdata_2, dfdata_interaction)
    t      <- c(t_1, t_2, t_interaction)
    B      <- c(B_1,B_2,B_interaction)
    p      <- c(p_1,p_2,p_interaction)
    
    #Putting the Table together
    test  <- c('Group 1', 'Group 2', 'Interaction')
    results_mixed<- data.frame(test, df, t, B, p)
    
    #Output tables
    output$table_mixed <- renderTable(results_mixed, spacing = "l", width = 300, digits = 3)
    
  })
  
  #Calculating the Bayes factor for within design
  observeEvent(input$action_within, {
    #calculate content of the Table
    dfdata <- (input$n_all_w) - 1 
    
    t_1 <- input$m1_w/input$se1_w
    t_2 <- input$m2_w/input$se2_w
    t_interaction <- (input$m_diff_w) / input$se_diff_w
    
    B_1 <- Bf(sd = input$se1_w, obtained = input$m1_w, dfdata = dfdata, meanoftheory = 0, sdtheory = input$sdH1_w, tail = 1)
    B_2 <- Bf(sd = input$se2_w, obtained = input$m2_w, dfdata = dfdata, meanoftheory = 0, sdtheory = input$sdH1_w, tail = 1)
    B_interaction <- Bf(sd = input$se_diff_w, obtained = (input$m_diff_w), dfdata = dfdata, meanoftheory = 0, sdtheory = input$sdH1_w, tail = 1)
    
    p_1 <- 2*pt(-abs(t_1),df=dfdata)
    p_2 <- 2*pt(-abs(t_2),df=dfdata)
    p_interaction <- 2*pt(-abs(t_interaction),df=dfdata)
    
    
    #columns of the Results Table
    df     <- c(dfdata, dfdata, dfdata)
    t      <- c(t_1, t_2, t_interaction)
    B      <- c(B_1,B_2,B_interaction)
    p      <- c(p_1,p_2,p_interaction)
    
    
    #Putting the Table together  
    test   <- c('Condition 1', 'Condition 2', 'Interaction')
    results_within<- data.frame(test, df, t, B, p)
    
    #Output tables  
    output$table_within  <- renderTable(results_within, spacing = "l", width = 300, digits = 3)
    
  })
    
    


  #Reset the input boxes
  observeEvent(input$reset_input_between, {
    updateNumericInput(session, "m1_b", value = 0)
    updateNumericInput(session, "se1_b", value = 1)
    updateNumericInput(session, "n1_b", value = 10)
    updateNumericInput(session, "m2_b", value = 0)
    updateNumericInput(session, "se2_b", value = 1)
    updateNumericInput(session, "n2_b", value = 10)
    updateNumericInput(session, "sdH1_b", value = 1)
    output$table_between  <- renderText({
    })
  })
  
  observeEvent(input$reset_input_mixed, {
    updateNumericInput(session, "m1_m", value = 0)
    updateNumericInput(session, "se1_m", value = 1)
    updateNumericInput(session, "n1_m", value = 10)
    updateNumericInput(session, "m2_m", value = 0)
    updateNumericInput(session, "se2_m", value = 1)
    updateNumericInput(session, "n2_m", value = 10)
    updateNumericInput(session, "sdH1_m", value = 1)
    output$table_mixed  <- renderText({
    })
  })
  
  observeEvent(input$reset_input_within, {
    updateNumericInput(session, "m1_w", value = 0)
    updateNumericInput(session, "se1_w", value = 1)
    updateNumericInput(session, "m2_w", value = 0)
    updateNumericInput(session, "se2_w", value = 1)
    updateNumericInput(session, "sdH1_w", value = 1)
    updateNumericInput(session, "m_diff_w", value = 0)
    updateNumericInput(session, "se_diff_w", value = 1)
    updateNumericInput(session, "n_all_w", value = 20)
    output$table_within  <- renderText({
    })
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)