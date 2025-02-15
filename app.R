library(shiny)
library(dplyr)
library(MASS)
library(lavaan)
library(glue)
library(tidyr)
library(mediation)
library(MBESS)
library(QuantPsyc)
library(boot)
library(foreach)
library(doParallel)

set.seed(42)

# Helper function with an added do_bootstrap flag
rsquare_med <- function(data, x, m, y, R = 100, ci_level = 0.95, do_bootstrap = TRUE) {
  # Compute correlations among the variables
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

  beta_YX_M <- QuantPsyc::lm.beta(model2)[[x]]
  beta_MX <- QuantPsyc::lm.beta(model1)[[x]]
  upsilon <- (rmy - (beta_MX * beta_YX_M))^2 - (overallrsquared - (rxy^2))

  if (do_bootstrap) {
    boot_stat <- function(d, indices) {
      d_boot <- d[indices, ]
      rxm_b <- cor(d_boot[[x]], d_boot[[m]])
      rxy_b <- cor(d_boot[[x]], d_boot[[y]])
      rmy_b <- cor(d_boot[[m]], d_boot[[y]])
      model1_b <- lm(as.formula(paste(m, "~", x)), data = d_boot)
      alpha_b <- coef(model1_b)[[x]]
      model2_b <- lm(as.formula(paste(y, "~", x, "+", m)), data = d_boot)
      tau_prime_b <- coef(model2_b)[[x]]
      beta_b <- coef(model2_b)[[m]]
      rxmsquared_b <- rxm_b^2
      overallrsquared_b <- (((rxy_b^2) + (rmy_b^2)) - (2 * rxy_b * rmy_b * rxm_b)) / (1 - rxmsquared_b)
      rsquaredmediated_b <- (rmy_b^2) - (overallrsquared_b - (rxy_b^2))
      beta_MX_b <- QuantPsyc::lm.beta(model1_b)[[x]]
      beta_YX_M_b <- QuantPsyc::lm.beta(model2_b)[[x]]
      overallrsquared_b <- summary(model2_b)$r.squared
      upsilon_b <- (rmy_b - (beta_MX_b * beta_YX_M_b))^2 - (overallrsquared_b - (rxy_b^2))
      return(c(rsquaredmediated = rsquaredmediated_b,
               upsilon = upsilon_b,
               overall_rsquared = overallrsquared_b))
    }

    boot_results <- boot(data = data, statistic = boot_stat, R = R)
    boot_ci_rsquaredmediated <- boot.ci(boot_results, conf = ci_level, type = "perc", index = 1)
    boot_ci_upsilon          <- boot.ci(boot_results, conf = ci_level, type = "perc", index = 2)
    boot_ci_overallrsquared  <- boot.ci(boot_results, conf = ci_level, type = "perc", index = 3)

    rsquaredmediated_ci_lower <- boot_ci_rsquaredmediated$percent[4]
    rsquaredmediated_ci_upper <- boot_ci_rsquaredmediated$percent[5]
    upsilon_ci_lower <- boot_ci_upsilon$percent[4]
    upsilon_ci_upper <- boot_ci_upsilon$percent[5]
    overallrsquared_ci_lower <- boot_ci_overallrsquared$percent[4]
    overallrsquared_ci_upper <- boot_ci_overallrsquared$percent[5]

    upsilon_width <- upsilon_ci_upper - upsilon_ci_lower
    rsquaredmediated_width <- rsquaredmediated_ci_upper - rsquaredmediated_ci_lower
    #ll_var <- rsquaredmediated_ci_lower / overallrsquared_ci_lower
    #uu_var <- rsquaredmediated_ci_upper / overallrsquared_ci_upper
  } else {
    rsquaredmediated_ci_lower <- NA
    rsquaredmediated_ci_upper <- NA
    upsilon_ci_lower <- NA
    upsilon_ci_upper <- NA
    overallrsquared_ci_lower <- NA
    overallrsquared_ci_upper <- NA
    upsilon_width <- NA
    rsquaredmediated_width <- NA
  }

  results <- list(
    alpha = alpha,
    beta = beta,
    tau_prime = tau_prime,
    total = total,
    mediatedeffect = mediatedeffect,
    rxmsquared = rxmsquared,
    rxm = rxm,
    rxy = rxy,
    rmy = rmy,
    partialrmy_xsquared = partialrmy_xsquared,
    partialrxy_msquared = partialrxy_msquared,
    overallrsquared = overallrsquared,
    rsquaredmediated_ci_lower = rsquaredmediated_ci_lower,
    rsquaredmediated = rsquaredmediated,
    rsquaredmediated_ci_upper = rsquaredmediated_ci_upper,
    upsilon_ci_lower = upsilon_ci_lower,
    upsilon = upsilon,
    upsilon_ci_upper = upsilon_ci_upper,
    proportion_mediated = if (total != 0) mediatedeffect / total else NA,
    upsilon_width = upsilon_width,
    rsquaredmediated_width = rsquaredmediated_width
  )

  return(results)
}

# Set up parallel backend
cl <- makeCluster(detectCores() - 1)
#registerDoParallel(cl)
registerDoSEQ()

clusterExport(cl, varlist = c("rsquare_med"), envir = environment())

##############################
# ORIGINAL SIMULATION FUNCTION
##############################
run_simulation_original <- function(pop_tau_prime,
                                    pop_alpha,
                                    pop_beta,
                                    sample_size,
                                    num_reps,
                                    do_bootstrap = TRUE) {
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
  pop_cov <- lavaan::lavInspect(fit, "cor.all")

  pop_data <- as.data.frame(MASS::mvrnorm(n = sample_size,
                                          mu = rep(0, 3),
                                          Sigma = pop_cov,
                                          empirical = TRUE))

  pop_rs <- unlist(rsquare_med(data = pop_data, x = "x", m = "m", y = "y", do_bootstrap = do_bootstrap))
  pop_rs["proportion_mediated"] <- if (pop_rs["total"] != 0) {
    pop_rs["mediatedeffect"] / pop_rs["total"]
  } else NA

  sim_matrix <- foreach(i = 1:num_reps,
                        .combine = cbind,
                        .packages = c("MASS", "QuantPsyc", "lavaan")) %dopar% {
                          sim_data <- as.data.frame(MASS::mvrnorm(n = sample_size,
                                                                  mu = rep(0, 3),
                                                                  Sigma = pop_cov,
                                                                  empirical = FALSE))
                          sim_out <- rsquare_med(data = sim_data, x = "x", m = "m", y = "y", do_bootstrap = do_bootstrap)
                          rsq_cov <- as.numeric((pop_rs["rsquaredmediated"] >= sim_out["rsquaredmediated_ci_lower"]) &&
                                                  (pop_rs["rsquaredmediated"] <= sim_out["rsquaredmediated_ci_upper"]))
                          upsilon_cov <- as.numeric((pop_rs["upsilon"] >= sim_out["upsilon_ci_lower"]) &&
                                                      (pop_rs["upsilon"] <= sim_out["upsilon_ci_upper"]))
                          unlist(c(sim_out, rsq_cov = rsq_cov, upsilon_cov = upsilon_cov))
                        }

  sim_means <- rowMeans(sim_matrix[setdiff(rownames(sim_matrix),
                                           c("rsq_cov", "upsilon_cov")), , drop = FALSE])
  sim_means["proportion_mediated"] <- if (sim_means["total"] != 0) {
    sim_means["mediatedeffect"] / sim_means["total"]
  } else NA

  rsq_cov_rate <- mean(sim_matrix["rsq_cov", ])
  upsilon_cov_rate <- mean(sim_matrix["upsilon_cov", ])

  res <- data.frame(
    parameter = names(pop_rs),
    pop_value = pop_rs,
    sim_value = sim_means,
    bias = sim_means - pop_rs,
    stringsAsFactors = FALSE
  )

  res <- dplyr::filter(res, parameter %in% c("alpha", "beta", "tau_prime",
                                             "mediatedeffect",
                                             "rxmsquared",
                                             "partialrmy_xsquared",
                                             "partialrxy_msquared",
                                             "proportion_mediated",
                                             "overallrsquared",
                                             "upsilon",
                                             "rxm",
                                             "rxy",
                                             "rmy",
                                             "rsquaredmediated_ci_lower",
                                             "rsquaredmediated",
                                             "rsquaredmediated_ci_upper",
                                             "upsilon_ci_lower",
                                             "upsilon_ci_upper",
                                             "upsilon_width",
                                             "rsquaredmediated_width"))

  cov_rows <- data.frame(
    parameter = c("rsquaredmediated_miss_perc", "upsilon_miss_perc"),
    pop_value = c(0.95, 0.95),
    sim_value = c(rsq_cov_rate, upsilon_cov_rate),
    bias = c(rsq_cov_rate - 0.95, upsilon_cov_rate - 0.95),
    stringsAsFactors = FALSE
  )

  res <- rbind(res, cov_rows)
  return(res)
}

###################################
# STANDARDISED SIMULATION FUNCTION
###################################
run_simulation_standardized <- function(pop_tau_prime,
                                        pop_alpha,
                                        pop_beta,
                                        sample_size,
                                        num_reps,
                                        do_bootstrap = TRUE) {
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

  pop_rs <- unlist(rsquare_med(data = pop_data, x = "x", m = "m", y = "y", do_bootstrap = do_bootstrap))
  pop_rs["proportion_mediated"] <- if (pop_rs["total"] != 0) {
    pop_rs["mediatedeffect"] / pop_rs["total"]
  } else NA

  sim_matrix <- foreach(i = 1:num_reps,
                        .combine = cbind,
                        .packages = c("MASS", "QuantPsyc", "lavaan")) %dopar% {
                          sim_data <- as.data.frame(MASS::mvrnorm(n = sample_size,
                                                                  mu = rep(0, 3),
                                                                  Sigma = pop_cov,
                                                                  empirical = FALSE))
                          sim_out <- rsquare_med(data = sim_data, x = "x", m = "m", y = "y", do_bootstrap = do_bootstrap)
                          rsq_cov <- as.numeric((pop_rs["rsquaredmediated"] >= sim_out["rsquaredmediated_ci_lower"]) &&
                                                  (pop_rs["rsquaredmediated"] <= sim_out["rsquaredmediated_ci_upper"]))
                          upsilon_cov <- as.numeric((pop_rs["upsilon"] >= sim_out["upsilon_ci_lower"]) &&
                                                      (pop_rs["upsilon"] <= sim_out["upsilon_ci_upper"]))
                          unlist(c(sim_out, rsq_cov = rsq_cov, upsilon_cov = upsilon_cov))
                        }

  sim_means <- rowMeans(sim_matrix[setdiff(rownames(sim_matrix),
                                           c("rsq_cov", "upsilon_cov")), , drop = FALSE])
  sim_means["proportion_mediated"] <- if (sim_means["total"] != 0) {
    sim_means["mediatedeffect"] / sim_means["total"]
  } else NA

  rsq_cov_rate <- mean(sim_matrix["rsq_cov", ])
  upsilon_cov_rate <- mean(sim_matrix["upsilon_cov", ])

  res <- data.frame(
    parameter = names(pop_rs),
    pop_value = pop_rs,
    sim_value = sim_means,
    bias = sim_means - pop_rs,
    stringsAsFactors = FALSE
  )

  res <- dplyr::filter(res, parameter %in% c("alpha", "beta", "tau_prime",
                                             "mediatedeffect",
                                             "rxmsquared",
                                             "partialrmy_xsquared",
                                             "partialrxy_msquared",
                                             "proportion_mediated",
                                             "overallrsquared",
                                             "upsilon",
                                             "rxm",
                                             "rxy",
                                             "rmy",
                                             "rsquaredmediated_ci_lower",
                                             "rsquaredmediated",
                                             "rsquaredmediated_ci_upper",
                                             "upsilon_ci_lower",
                                             "upsilon_ci_upper",
                                             "upsilon_width",
                                             "rsquaredmediated_width"))

  cov_rows <- data.frame(
    parameter = c("rsquaredmediated_miss_perc", "upsilon_miss_perc"),
    pop_value = c(0.95, 0.95),
    sim_value = c(rsq_cov_rate, upsilon_cov_rate),
    bias = c(rsq_cov_rate - 0.95, upsilon_cov_rate - 0.95),
    stringsAsFactors = FALSE
  )

  res <- rbind(res, cov_rows)
  return(res)
}

# SHINY APP UI
ui <- fluidPage(
  titlePanel("Fairchild et al., 2009 Simulation"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("sample_size", "Sample Size", min = 10, max = 5000, value = 200, step = 10),
      numericInput("pop_alpha", "Population alpha", value = 0.078, min = 0, max = 1, step = 0.01),
      fluidRow(
        column(3, actionButton("alpha_0", "0", class = "btn btn-default btn-xs")),
        column(3, actionButton("alpha_14", "0.14", class = "btn btn-default btn-xs")),
        column(3, actionButton("alpha_39", "0.39", class = "btn btn-default btn-xs")),
        column(3, actionButton("alpha_59", "0.59", class = "btn btn-default btn-xs"))
      ),
      numericInput("pop_beta", "Population beta", value = 0.609, min = 0, max = 1, step = 0.01),
      fluidRow(
        column(3, actionButton("beta_0", "0", class = "btn btn-default btn-xs")),
        column(3, actionButton("beta_14", "0.14", class = "btn btn-default btn-xs")),
        column(3, actionButton("beta_39", "0.39", class = "btn btn-default btn-xs")),
        column(3, actionButton("beta_59", "0.59", class = "btn btn-default btn-xs"))
      ),
      numericInput("pop_tau_prime", "Population tau prime", value = 0.0385, min = 0, max = 1, step = 0.01),
      fluidRow(
        column(3, actionButton("tau_0", "0", class = "btn btn-default btn-xs")),
        column(3, actionButton("tau_14", "0.14", class = "btn btn-default btn-xs")),
        column(3, actionButton("tau_39", "0.39", class = "btn btn-default btn-xs")),
        column(3, actionButton("tau_59", "0.59", class = "btn btn-default btn-xs"))
      ),
      numericInput("num_reps", "Number of replications", value = 100, min = 1, max = 1000),
      checkboxInput("do_bootstrap", "Run Bootstrap (takes much longer)", value = FALSE),
      actionButton("run_sim", "Run Simulations")
    ),
    mainPanel(
      fluidRow(
        column(6,
               h4("Theoretical Total Effect"),
               textOutput("theoretical_total")
        )
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

  observeEvent(input$alpha_0,  { updateNumericInput(session, "pop_alpha", value = 0) })
  observeEvent(input$alpha_14, { updateNumericInput(session, "pop_alpha", value = 0.14) })
  observeEvent(input$alpha_39, { updateNumericInput(session, "pop_alpha", value = 0.39) })
  observeEvent(input$alpha_59, { updateNumericInput(session, "pop_alpha", value = 0.59) })

  observeEvent(input$beta_0,  { updateNumericInput(session, "pop_beta", value = 0) })
  observeEvent(input$beta_14, { updateNumericInput(session, "pop_beta", value = 0.14) })
  observeEvent(input$beta_39, { updateNumericInput(session, "pop_beta", value = 0.39) })
  observeEvent(input$beta_59, { updateNumericInput(session, "pop_beta", value = 0.59) })

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
        num_reps = input$num_reps,
        do_bootstrap = input$do_bootstrap
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
        num_reps = input$num_reps,
        do_bootstrap = input$do_bootstrap
      )
      incProgress(0.5)
      res
    })
  })

  output$orig_results <- renderTable({
    req(orig_sim_data())
    orig_sim_data()
  }, digits = 10)

  output$std_results <- renderTable({
    req(std_sim_data())
    std_sim_data()
  }, digits = 10)
}

shinyApp(ui = ui, server = server)