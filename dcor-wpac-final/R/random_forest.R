# Random Forest Population Assignment for Sea Turtles
# Load required packages
library(randomForest)
library(caret)
library(adegenet)
library(dplyr)
library(ggplot2)
library(VIM) # for missing data visualization

# ============================================================================
# 1. DATA PREPARATION
# ============================================================================

# Convert genind object to data frame (assuming your data is in genind format)
# If you have a genind object called 'turtle_genind':
convert_genind_to_df <- function(genind_obj) {
  # Extract genotype matrix
  geno_matrix <- tab(genind_obj, freq = FALSE, NA.method = "mean")
  
  # Create data frame with populations
  df <- data.frame(
    Population = pop(genind_obj),
    geno_matrix
  )
  
  return(df)
}

# Example with simulated data structure:
# turtle_df <- convert_genind_to_df(turtle_genind)

# Or if starting from a different format:
prepare_rf_data <- function(snp_matrix, pop_info) {
  # snp_matrix: individuals x SNPs matrix
  # pop_info: vector of population assignments (known origin samples)
  
  # Handle missing data - RF can handle some NAs but better to impute
  # Option 1: Simple imputation with most frequent allele
  for(i in 1:ncol(snp_matrix)) {
    snp_matrix[is.na(snp_matrix[,i]), i] <- 
      names(sort(table(snp_matrix[,i]), decreasing = TRUE))[1]
  }
  
  # Create final data frame
  rf_data <- data.frame(
    Population = as.factor(pop_info),
    snp_matrix
  )
  
  return(rf_data)
}

# ============================================================================
# 2. DATA SPLITTING AND CROSS-VALIDATION SETUP
# ============================================================================

# Split known-origin data into training and testing sets
split_data <- function(rf_data, test_proportion = 0.3) {
  set.seed(123) # for reproducibility
  
  # Stratified sampling to maintain population proportions
  train_indices <- createDataPartition(
    rf_data$Population, 
    p = 1 - test_proportion,
    list = FALSE
  )
  
  train_data <- rf_data[train_indices, ]
  test_data <- rf_data[-train_indices, ]
  
  return(list(train = train_data, test = test_data))
}

# Example usage:
# data_split <- split_data(turtle_df[!is.na(turtle_df$Population), ])

# ============================================================================
# 3. MODEL TRAINING AND OPTIMIZATION
# ============================================================================

# Train Random Forest with parameter tuning
train_rf_model <- function(train_data, tune_parameters = TRUE) {
  
  if(tune_parameters) {
    # Grid search for optimal parameters
    tune_grid <- expand.grid(
      mtry = c(sqrt(ncol(train_data)-1), 
               ncol(train_data)/3, 
               ncol(train_data)/2)
    )
    
    # 10-fold cross-validation
    ctrl <- trainControl(
      method = "cv",
      number = 10,
      classProbs = TRUE,
      savePredictions = "final"
    )
    
    # Train with parameter tuning
    rf_model <- train(
      Population ~ .,
      data = train_data,
      method = "rf",
      tuneGrid = tune_grid,
      trControl = ctrl,
      ntree = 1000,
      importance = TRUE
    )
    
  } else {
    # Simple training without tuning
    rf_model <- randomForest(
      Population ~ .,
      data = train_data,
      ntree = 1000,
      importance = TRUE,
      do.trace = 100 # progress updates
    )
  }
  
  return(rf_model)
}

# ============================================================================
# 4. MODEL EVALUATION
# ============================================================================

# Evaluate model performance
evaluate_model <- function(rf_model, test_data) {
  # Predictions
  predictions <- predict(rf_model, test_data)
  
  # Confusion matrix
  conf_matrix <- confusionMatrix(predictions, test_data$Population)
  
  # Variable importance
  importance_scores <- varImp(rf_model)
  
  # Out-of-bag error rate
  if(class(rf_model)[1] == "randomForest") {
    oob_error <- rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]
  } else {
    oob_error <- NA
  }
  
  return(list(
    confusion_matrix = conf_matrix,
    importance = importance_scores,
    oob_error = oob_error,
    predictions = predictions
  ))
}

# ============================================================================
# 5. ASSIGN UNKNOWN SAMPLES
# ============================================================================

# Assign unknown-origin samples
assign_unknowns <- function(rf_model, unknown_data, confidence_threshold = 0.7) {
  # Get prediction probabilities
  pred_probs <- predict(rf_model, unknown_data, type = "prob")
  
  # Get most likely assignment
  assignments <- predict(rf_model, unknown_data)
  
  # Calculate confidence (max probability)
  confidence <- apply(pred_probs, 1, max)
  
  # Create results data frame
  results <- data.frame(
    Sample_ID = rownames(unknown_data),
    Assigned_Population = assignments,
    Confidence = confidence,
    High_Confidence = confidence >= confidence_threshold,
    pred_probs
  )
  
  return(results)
}

# ============================================================================
# 6. VISUALIZATION FUNCTIONS
# ============================================================================

# Plot variable importance
plot_importance <- function(rf_model, top_n = 20) {
  if(class(rf_model)[1] == "train") {
    imp_data <- varImp(rf_model)$importance
    imp_data$SNP <- rownames(imp_data)
    imp_data <- imp_data[order(imp_data$Overall, decreasing = TRUE)[1:top_n], ]
    
    ggplot(imp_data, aes(x = reorder(SNP, Overall), y = Overall)) +
      geom_col() +
      coord_flip() +
      labs(title = paste("Top", top_n, "Most Important SNPs"),
           x = "SNP", y = "Importance Score") +
      theme_minimal()
  }
}

# Plot assignment results
plot_assignments <- function(assignment_results, confidence_threshold = 0.7) {
  ggplot(assignment_results, aes(x = Assigned_Population, fill = High_Confidence)) +
    geom_bar(position = "dodge") +
    labs(title = "Population Assignments for Unknown Samples",
         x = "Assigned Population", y = "Number of Samples",
         fill = paste("Confidence ≥", confidence_threshold)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Plot confidence distribution
plot_confidence <- function(assignment_results) {
  ggplot(assignment_results, aes(x = Confidence)) +
    geom_histogram(bins = 30, alpha = 0.7) +
    geom_vline(xintercept = 0.7, color = "red", linetype = "dashed") +
    labs(title = "Distribution of Assignment Confidence",
         x = "Confidence Score", y = "Number of Samples") +
    theme_minimal()
}

# ============================================================================
# 7. EXAMPLE WORKFLOW
# ============================================================================

# Complete workflow example:
run_rf_assignment <- function(known_data, unknown_data, confidence_threshold = 0.7) {
  
  cat("1. Splitting data...\n")
  data_split <- split_data(known_data, test_proportion = 0.3)
  
  cat("2. Training Random Forest model...\n")
  rf_model <- train_rf_model(data_split$train)
  
  cat("3. Evaluating model performance...\n")
  evaluation <- evaluate_model(rf_model, data_split$test)
  
  cat("Model Accuracy:", evaluation$confusion_matrix$overall['Accuracy'], "\n")
  
  cat("4. Assigning unknown samples...\n")
  assignments <- assign_unknowns(rf_model, unknown_data, confidence_threshold)
  
  cat("5. Creating visualizations...\n")
  importance_plot <- plot_importance(rf_model)
  assignment_plot <- plot_assignments(assignments, confidence_threshold)
  confidence_plot <- plot_confidence(assignments)
  
  return(list(
    model = rf_model,
    evaluation = evaluation,
    assignments = assignments,
    plots = list(
      importance = importance_plot,
      assignments = assignment_plot,
      confidence = confidence_plot
    )
  ))
}

# ============================================================================
# 8. ADDITIONAL UTILITY FUNCTIONS
# ============================================================================

# Compare RF assignments with other methods (e.g., STRUCTURE, DAPC)
compare_methods <- function(rf_assignments, other_assignments, method_name = "Other") {
  comparison <- data.frame(
    Sample_ID = rf_assignments$Sample_ID,
    RF_Assignment = rf_assignments$Assigned_Population,
    Other_Assignment = other_assignments,
    Agreement = rf_assignments$Assigned_Population == other_assignments,
    RF_Confidence = rf_assignments$Confidence
  )
  
  agreement_rate <- mean(comparison$Agreement, na.rm = TRUE)
  cat("Agreement rate between RF and", method_name, ":", agreement_rate, "\n")
  
  return(comparison)
}

# Export results for further analysis
export_results <- function(assignments, filename = "rf_assignments.csv") {
  write.csv(assignments, filename, row.names = FALSE)
  cat("Results exported to", filename, "\n")
}

# Summary statistics
summarize_assignments <- function(assignment_results, confidence_threshold = 0.7) {
  summary_stats <- assignment_results %>%
    group_by(Assigned_Population) %>%
    summarise(
      n_samples = n(),
      mean_confidence = mean(Confidence),
      high_confidence_samples = sum(High_Confidence),
      .groups = 'drop'
    )
  
  cat("Assignment Summary:\n")
  print(summary_stats)
  
  cat("\nOverall high-confidence assignments:", 
      sum(assignment_results$High_Confidence), "/", nrow(assignment_results), "\n")
  
  return(summary_stats)
}