library(shiny)
library(dplyr)
library(MASS)
library(lavaan)
library(glue)
library(tidyr)
library(semPlot)

# A helper function that computes mediation parameters
rsquare_med <- function(data, x, m, y) {
  rxm <- cor(data[[x]], data[[m]])
  rxy <- cor(data[[x]], data[[y]])
  rmy <- cor(data[[m]], data[[y]])

  # Regression: m ~ x (for alpha)
  model1 <- lm(as.formula(paste(m, "~", x)), data = data)
  alpha <- coef(model1)[[x]]

  # Regression: y ~ x + m (for tau_prime and beta)
  model2 <- lm(as.formula(paste(y, "~", x, "+", m)), data = data)
  tau_prime <- coef(model2)[[x]]
  beta <- coef(model2)[[m]]

  total <- tau_prime + (alpha * beta)
  mediatedeffect <- alpha * beta
  rxmsquared <- rxm^2
  partialrxy_msquared <- ((rxy - rmy * rxm) / sqrt((1 - rmy^2) * (1 - rxmsquared)))^2
  partialrmy_xsquared <- ((rmy - rxy * rxm) / sqrt((1 - rxy^2) * (1 - rxmsquared)))^2
  overallrsquared <- (((rxy^2) + (rmy^2)) - (2 * rxy * rmy * rxm)) / (1 - rxmsquared)
  rsquaredmediated <- (rmy^2) - (overallrsquared - (rxy^2))

  proportion_mediated <- if (total != 0) mediatedeffect / total else NA

  results <- list(
    alpha = alpha,
    beta = beta,
    tau_prime = tau_prime,
    total = total,
    mediatedeffect = mediatedeffect,
    rxm = rxm,
    rxmsquared = rxmsquared,
    rxy = rxy,
    rmy = rmy,
    partialrmy_xsquared = partialrmy_xsquared,
    partialrxy_msquared = partialrxy_msquared,
    overallrsquared = overallrsquared,
    rsquaredmediated = rsquaredmediated,
    proportion_mediated = proportion_mediated
  )

  return(results)
}

# ORIGINAL SIMULATION FUNCTION (with fixed variances)
run_simulation_original <- function(pop_tau_prime,
                                    pop_alpha,
                                    pop_beta,
                                    sample_size,
                                    num_reps) {
  model_string <- glue("
    y ~ {pop_tau_prime} * x + {pop_beta} * m
    m ~ {pop_alpha} * x
    x ~~ 1 * x
    y ~~ 1 * y
    m ~~ 1 * m
  ")

  d_fake <- data.frame(
    x = rnorm(sample_size),
    m = rnorm(sample_size),
    y = rnorm(sample_size)
  )
  fit <- lavaan::lavaan(model = model_string, data = d_fake)
  pop_cov <- lavaan::lavInspect(fit, "cov.all")

  pop_data <- as.data.frame(MASS::mvrnorm(n = sample_size,
                                          mu = rep(0, 3),
                                          Sigma = pop_cov,
                                          empirical = TRUE))

  pop_rs <- unlist(rsquare_med(data = pop_data, x = "x", m = "m", y = "y"))
  pop_rs["proportion_mediated"] <- if (pop_rs["total"] != 0) {
    pop_rs["mediatedeffect"] / pop_rs["total"]
  } else NA

  sim_matrix <- replicate(num_reps, {
    sim_data <- as.data.frame(MASS::mvrnorm(n = sample_size,
                                            mu = rep(0, 3),
                                            Sigma = pop_cov,
                                            empirical = FALSE))
    unlist(rsquare_med(data = sim_data, x = "x", m = "m", y = "y"))
  })
  sim_means <- rowMeans(sim_matrix)
  sim_means["proportion_mediated"] <- if (sim_means["total"] != 0) {
    sim_means["mediatedeffect"] / sim_means["total"]
  } else NA

  res <- data.frame(
    parameter = names(pop_rs),
    pop_value = pop_rs,
    sim_value = sim_means,
    bias = sim_means - pop_rs,
    stringsAsFactors = FALSE
  )

  res <- filter(res, parameter %in% c("alpha", "beta", "tau_prime",
                                      "mediatedeffect",
                                      "rxmsquared", # depend on a
                                      "partialrmy_xsquared", # depend on b
                                      "partialrxy_msquared", # depend on b and a
                                      "rsquaredmediated",
                                      "proportion_mediated",
                                      "overallrsquared"))
  return(res)
}

# STANDARDISED SIMULATION FUNCTION (with computed residual variances)
run_simulation_standardized <- function(pop_tau_prime,
                                        pop_alpha,
                                        pop_beta,
                                        sample_size,
                                        num_reps) {
  model_string <- glue("
    m ~ {pop_alpha} * x
    y ~ {pop_tau_prime} * x + {pop_beta} * m
    x ~~ 1 * x
    m ~~ {1 - pop_alpha^2} * m
    y ~~ {1 - (pop_tau_prime + pop_beta * pop_alpha)^2 - pop_beta^2 * (1 - pop_alpha^2)} * y
  ")

  d_fake <- data.frame(
    x = rnorm(sample_size),
    m = rnorm(sample_size),
    y = rnorm(sample_size)
  )
  fit <- lavaan::lavaan(model = model_string, data = d_fake)
  pop_cov <- lavaan::lavInspect(fit, "cov.all")

  pop_data <- as.data.frame(MASS::mvrnorm(n = sample_size,
                                          mu = rep(0, 3),
                                          Sigma = pop_cov,
                                          empirical = TRUE))
  pop_data <- as.data.frame(scale(pop_data))

  pop_rs <- unlist(rsquare_med(data = pop_data, x = "x", m = "m", y = "y"))
  pop_rs["proportion_mediated"] <- if (pop_rs["total"] != 0) {
    pop_rs["mediatedeffect"] / pop_rs["total"]
  } else NA

  sim_matrix <- replicate(num_reps, {
    sim_data <- as.data.frame(MASS::mvrnorm(n = sample_size,
                                            mu = rep(0, 3),
                                            Sigma = pop_cov,
                                            empirical = FALSE))
    sim_data <- as.data.frame(scale(sim_data))
    unlist(rsquare_med(data = sim_data, x = "x", m = "m", y = "y"))
  })
  sim_means <- rowMeans(sim_matrix)
  sim_means["proportion_mediated"] <- if (sim_means["total"] != 0) {
    sim_means["mediatedeffect"] / sim_means["total"]
  } else NA

  res <- data.frame(
    parameter = names(pop_rs),
    pop_value = pop_rs,
    sim_value = sim_means,
    bias = sim_means - pop_rs,
    stringsAsFactors = FALSE
  )

  res <- filter(res, parameter %in% c("alpha", "beta", "tau_prime",
                                      "mediatedeffect",
                                      "rxmsquared", # depend on a
                                      "partialrmy_xsquared", # depend on b
                                      "partialrxy_msquared", # depend on b and a
                                      "rsquaredmediated",
                                      "proportion_mediated",
                                      "overallrsquared"))
  return(res)
}

# SHINY APP UI
ui <- fluidPage(
  titlePanel("Fairchild et al., 2009 Simulation"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("sample_size", "Sample Size", min = 10, max = 1000, value = 200, step = 10),

      numericInput("pop_alpha", "Population alpha", value = 0.14, min = 0, max = 1, step = 0.01),
      fluidRow(
        column(3, actionButton("alpha_0", "0", class = "btn btn-default btn-xs")),
        column(3, actionButton("alpha_14", "0.14", class = "btn btn-default btn-xs")),
        column(3, actionButton("alpha_39", "0.39", class = "btn btn-default btn-xs")),
        column(3, actionButton("alpha_59", "0.59", class = "btn btn-default btn-xs"))
      ),

      numericInput("pop_beta", "Population beta", value = 0.5, min = 0, max = 1, step = 0.01),
      fluidRow(
        column(3, actionButton("beta_0", "0", class = "btn btn-default btn-xs")),
        column(3, actionButton("beta_14", "0.14", class = "btn btn-default btn-xs")),
        column(3, actionButton("beta_39", "0.39", class = "btn btn-default btn-xs")),
        column(3, actionButton("beta_59", "0.59", class = "btn btn-default btn-xs"))
      ),

      numericInput("pop_tau_prime", "Population tau prime", value = 0.5, min = 0, max = 1, step = 0.01),
      fluidRow(
        column(3, actionButton("tau_0", "0", class = "btn btn-default btn-xs")),
        column(3, actionButton("tau_14", "0.14", class = "btn btn-default btn-xs")),
        column(3, actionButton("tau_39", "0.39", class = "btn btn-default btn-xs")),
        column(3, actionButton("tau_59", "0.59", class = "btn btn-default btn-xs"))
      ),

      numericInput("num_reps", "Number of replications", value = 100, min = 1, max = 1000),
      actionButton("run_sim", "Run Simulations")
    ),
    mainPanel(
      fluidRow(
        column(6,
               h4("Theoretical Total Effect"),
               textOutput("theoretical_total")
        ),
      ),
      br(),
      fluidRow(
        column(6,
               h3("Original Simulation Results"),
               tableOutput("orig_results")
        ),
        column(6,
               h3("Standardised Simulation Results"),
               tableOutput("std_results")
        )
      )
    )
  )
)

# SHINY APP SERVER
server <- function(input, output, session) {

  # Update pop_alpha numeric input
  observeEvent(input$alpha_0,  { updateNumericInput(session, "pop_alpha", value = 0) })
  observeEvent(input$alpha_14, { updateNumericInput(session, "pop_alpha", value = 0.14) })
  observeEvent(input$alpha_39, { updateNumericInput(session, "pop_alpha", value = 0.39) })
  observeEvent(input$alpha_59, { updateNumericInput(session, "pop_alpha", value = 0.59) })

  # Update pop_beta numeric input
  observeEvent(input$beta_0,  { updateNumericInput(session, "pop_beta", value = 0) })
  observeEvent(input$beta_14, { updateNumericInput(session, "pop_beta", value = 0.14) })
  observeEvent(input$beta_39, { updateNumericInput(session, "pop_beta", value = 0.39) })
  observeEvent(input$beta_59, { updateNumericInput(session, "pop_beta", value = 0.59) })

  # Update pop_tau_prime numeric input
  observeEvent(input$tau_0,  { updateNumericInput(session, "pop_tau_prime", value = 0) })
  observeEvent(input$tau_14, { updateNumericInput(session, "pop_tau_prime", value = 0.14) })
  observeEvent(input$tau_39, { updateNumericInput(session, "pop_tau_prime", value = 0.39) })
  observeEvent(input$tau_59, { updateNumericInput(session, "pop_tau_prime", value = 0.59) })
  output$theoretical_total <- renderText({
    total_effect <- input$pop_tau_prime + (input$pop_alpha * input$pop_beta)
    paste(total_effect)
  })

  orig_sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running Original Simulation", value = 0, {
      incProgress(0.5)
      res <- run_simulation_original(
        pop_tau_prime = input$pop_tau_prime,
        pop_alpha = input$pop_alpha,
        pop_beta = input$pop_beta,
        sample_size = input$sample_size,
        num_reps = input$num_reps
      )
      incProgress(0.5)
      res
    })
  })

  std_sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running Standardised Simulation", value = 0, {
      incProgress(0.5)
      res <- run_simulation_standardized(
        pop_tau_prime = input$pop_tau_prime,
        pop_alpha = input$pop_alpha,
        pop_beta = input$pop_beta,
        sample_size = input$sample_size,
        num_reps = input$num_reps
      )
      incProgress(0.5)
      res
    })
  })

  output$orig_results <- renderTable({
    req(orig_sim_data())
    orig_sim_data()
  }, digits = 3)

  output$std_results <- renderTable({
    req(std_sim_data())
    std_sim_data()
  }, digits = 3)
}


shinyApp(ui = ui, server = server)




#
# pop_data <- as.data.frame(scale(pop_data))
# apply(pop_data, 2, var)
# apply(pop_data, 2, mean)
