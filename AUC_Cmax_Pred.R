library(shiny)
library(ggplot2)
library(dplyr)

# Define UI
ui <- fluidPage(
  titlePanel("AUC-Cmax Prediction Tool using Dose and the Power Model"),
  # Input boxes across the top
  fluidRow(
    column(6,
           h4("Observed AUC and Cmax Information (Up to 5 dose levels)"),
           numericInput("n_doses", "Number of Input Dose Levels (1-5):", value = 3, min = 1, max = 5, step = 1),
           uiOutput("dose_inputs")
    ),
    column(6,
           h4("New Dose Levels for Prediction (Up to 5)"),
           numericInput("n_new_doses", "Number of New Dose Levels (1-5):", value = 2, min = 1, max = 5, step = 1),
           uiOutput("new_dose_inputs")
    )
  ),
  # Prediction table and proportionality table side by side
  fluidRow(
    column(6,
           h4("Predicted AUC and Cmax for New Dose Levels"),
           tableOutput("predictions")
    ),
    column(6,
           h4("Proportionality"),
           tableOutput("b1_parameters")
    )
  ),
  # Plots side by side
  fluidRow(
    column(6,
           h4("AUC vs Dose Plot"),
           plotOutput("auc_plot")
    ),
    column(6,
           h4("Cmax vs Dose Plot"),
           plotOutput("cmax_plot")
    )
  )
)


# Define server logic
server <- function(input, output) {
  
  # Dynamically generate input fields for dose levels
  output$dose_inputs <- renderUI({
    n <- input$n_doses
    lapply(1:n, function(i) {
      fluidRow(
        column(4, numericInput(paste0("dose_", i), paste("Dose", i, "(mg):"), value = 1, min = 0)),
        column(4, numericInput(paste0("auc_", i), paste("AUC", i, "(ng·h/mL):"), value = 1, min = 0)),
        column(4, numericInput(paste0("cmax_", i), paste("Cmax", i, "(ng/mL):"), value = 1, min = 0))
      )
    })
  })
  
  # Dynamically generate input fields for new dose levels
  output$new_dose_inputs <- renderUI({
    n <- input$n_new_doses
    lapply(1:n, function(i) {
      numericInput(paste0("new_dose_", i), paste("New Dose", i, "(mg):"), value = 1, min = 0)
    })
  })
  
  # Reactive function to collect input data
  input_data <- reactive({
    n <- input$n_doses
    doses <- sapply(1:n, function(i) input[[paste0("dose_", i)]])
    aucs <- sapply(1:n, function(i) input[[paste0("auc_", i)]])
    cmaxs <- sapply(1:n, function(i) input[[paste0("cmax_", i)]])
    data.frame(Dose = doses, AUC = aucs, Cmax = cmaxs)
  })
  
  # Reactive function to collect new dose levels
  new_doses <- reactive({
    n <- input$n_new_doses
    sapply(1:n, function(i) input[[paste0("new_dose_", i)]])
  })
  
  # Calculate predictions using linear interpolation
  predictions <- reactive({
    data <- input_data()
    new_dose_levels <- new_doses()
    
    # Remove any rows with NA or invalid inputs
    data <- data[complete.cases(data) & data$Dose > 0, ]
    
    if (nrow(data) < 2) {
      return(list(
        predictions = data.frame(New_Dose = numeric(), Predicted_AUC = numeric(), Predicted_Cmax = numeric()),
        b1_parameters = data.frame(Parameter = character(), b1 = numeric(), CI_Lower = numeric(), CI_Upper = numeric(), Proportionality = character())
      ))
    }
    
    # Fit power model: log(AUC) = b1 * log(Dose) + b2
    auc_model <- lm(log(AUC) ~ log(Dose), data = data)
    b1_auc <- coef(auc_model)[2]
    b2_auc <- coef(auc_model)[1]
    auc_ci <- confint(auc_model, "log(Dose)", level = 0.90)
    
    # Fit power model: log(Cmax) = b1 * log(Dose) + b2
    cmax_model <- lm(log(Cmax) ~ log(Dose), data = data)
    b1_cmax <- coef(cmax_model)[2]
    b2_cmax <- coef(cmax_model)[1]
    cmax_ci <- confint(cmax_model, "log(Dose)", level = 0.90)
    
    # Predict for new doses
    pred_auc <- exp(b1_auc * log(new_dose_levels) + b2_auc)
    pred_cmax <- exp(b1_cmax * log(new_dose_levels) + b2_cmax)
    
    # Handle any potential NaN/Inf values (e.g., if new_dose = 0)
    pred_auc[!is.finite(pred_auc)] <- 0
    pred_cmax[!is.finite(pred_cmax)] <- 0
    
    # Determine proportionality for AUC
    auc_proportionality <- if (auc_ci[1] <= 1.0 && auc_ci[2] >= 1.0) {
      "linear"
    } else if (b1_auc > 1.0) {
      "greater than proportional"
    } else {
      "less than proportional"
    }
    
    # Determine proportionality for Cmax
    cmax_proportionality <- if (cmax_ci[1] <= 1.0 && cmax_ci[2] >= 1.0) {
      "linear"
    } else if (b1_cmax > 1.0) {
      "greater than proportional"
    } else {
      "less than proportional"
    }
    
    # Create b1 parameters table
    b1_parameters <- data.frame(
      Parameter = c("AUC", "Cmax"),
      b1 = c(b1_auc, b1_cmax),
      CI_Lower = c(auc_ci[1], cmax_ci[1]),
      CI_Upper = c(auc_ci[2], cmax_ci[2]),
      Proportionality = c(auc_proportionality, cmax_proportionality)
    )
    
    list(
      predictions = data.frame(New_Dose = new_dose_levels, Predicted_AUC = pred_auc, Predicted_Cmax = pred_cmax),
      b1_parameters = b1_parameters
    )
  })
  
  # Render prediction table
  output$predictions <- renderTable({
    predictions()$predictions
  }, digits = 2)
  
  # Render b1 parameters table
  output$b1_parameters <- renderTable({
    predictions()$b1_parameters
  }, digits = 2)
  
  # AUC vs Dose Plot
  output$auc_plot <- renderPlot({
    data <- input_data()
    pred <- predictions()$predictions
    
    # Combine input and predicted data for plotting
    plot_data <- bind_rows(
      data.frame(Dose = data$Dose, Value = data$AUC, Type = "Observation"),
      data.frame(Dose = pred$New_Dose, Value = pred$Predicted_AUC, Type = "Prediction")
    )
    
    ggplot(plot_data, aes(x = Dose, y = Value, color = Type)) +
      geom_point(size = 3) +
      geom_line(aes(group = Type)) +
      labs(x = "Dose (mg)", y = "AUC (ng·h/mL)", title = "AUC vs Dose") +
      scale_color_manual(values = c("Observation" = "blue", "Prediction" = "red")) +
      theme_minimal() +
      theme(axis.title = element_text(size = 14))
  })
  
  # Cmax vs Dose Plot
  output$cmax_plot <- renderPlot({
    data <- input_data()
    pred <- predictions()$predictions
    
    # Combine input and predicted data for plotting
    plot_data <- bind_rows(
      data.frame(Dose = data$Dose, Value = data$Cmax, Type = "Observation"),
      data.frame(Dose = pred$New_Dose, Value = pred$Predicted_Cmax, Type = "Prediction")
    )
    
    ggplot(plot_data, aes(x = Dose, y = Value, color = Type)) +
      geom_point(size = 3) +
      geom_line(aes(group = Type)) +
      labs(x = "Dose (mg)", y = "Cmax (ng/mL)", title = "Cmax vs Dose") +
      scale_color_manual(values = c("Observation" = "blue", "Prediction" = "red")) +
      theme_minimal() +
      theme(axis.title = element_text(size = 14))
  })
}

# Run the application
shinyApp(ui = ui, server = server)