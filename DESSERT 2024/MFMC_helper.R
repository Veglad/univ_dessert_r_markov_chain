library(markovchain)

################################### Helper functions ###################################
get_fragment_start_idx <- function(fragment_size, fragment_number) {
  return ((fragment_number - 1) * fragment_size + 1)
}

get_fragment_end_idx <- function(fragment_size, fragment_number) {
  return (fragment_number * fragment_size)
}

##### Transform functions
# Function verification
# f1_s <- function(x) { pow_transform(D, 200, x) }
# curve(expr = f1_s, from = 0, to = 200)
linear_transform <- function(D, steps_count, step) { 
  x <- (step - 1) / (steps_count - 1)
  growth <- x * (1 - D)
  return (D + growth)
}

pow_transform <- function(D, steps_count, step) { 
  x <- (step - 1) / (steps_count - 1)
  growth <- (x ^ 0.3) * (1 - D)
  return (D + growth)
}

##### Transform functions End

# transition_matrix_1fragment - matrix for 1st fragment
# basic_state_labels - state labels
# count - fragment count
# link_data - describes how to link different fragments
# transform_last_fragment - can be used to transform last fragment (linking state can be removed)
# D - varying param per fragments - reliability of detection
# transform_D - How D changes after time (fucntion)
# build_next_fragment compose next fragment
generate_MFMC <- function(transition_matrix_1fragment, basic_state_labels, 
                          count, link_data, transform_last_fragment, build_next_fragment, D, transform_D) {
  # init labels
  all_state_labels <- do.call(paste0,expand.grid(basic_state_labels, 1:count))
  
  # create empty matrix
  fragment_size <- nrow(transition_matrix_1fragment)
  matrix_size <- fragment_size * count
  final_matrix <- matrix(
    c(rep(0, matrix_size * matrix_size)),
    byrow = TRUE, 
    nrow = matrix_size, 
    dimnames = list(all_state_labels, all_state_labels)
  )
  
  # init with fragments
  for (i in 1:count) {
    start_idx <- get_fragment_start_idx(fragment_size, i)
    end_idx <- get_fragment_end_idx(fragment_size, i)
    selected_indices <- start_idx:end_idx
    
    new_fragment <- build_next_fragment(transition_matrix_1fragment, transform_D(D, count, i))
    final_matrix[selected_indices,selected_indices] <- new_fragment   
  }
  
  # init with transition between fragments
  for (i in 1:count) {
    if (i != count) { # Last fragment has no link transitions to the next fragment
      column <- get_fragment_end_idx(fragment_size, i) + link_data$to
      row <- get_fragment_start_idx(fragment_size, i) - 1 + link_data$from
      final_matrix[row, column] <- link_data$value
    }
  }
  
  # transfrom last fragment (linking state may be removed)
  last_fragment_from_idx <- get_fragment_start_idx(fragment_size, count)
  last_fragment_to_idx <- get_fragment_end_idx(fragment_size, count)
  last_fragment_indices <- last_fragment_from_idx:last_fragment_to_idx
  last_fragment <- final_matrix[last_fragment_indices, last_fragment_indices]
  last_fragment_updated <- transform_last_fragment(last_fragment, link_data)
  final_matrix[last_fragment_indices, last_fragment_indices] <- last_fragment_updated
  
  return(list(final_matrix = final_matrix, all_state_labels = all_state_labels))
}


# Description: print report for results of readiness functions of model with and without learning state
# readiness_over_time - A(t) for steps for model with learning state
# readiness_over_time_wl - A(t) for steps for model without learning state
# markov - markov chain of the model with learning step
print_report_for_rediness_functions <- function(readiness_over_time, readiness_over_time_wl) {
  lerning_min_to_without_learning_min <- min(readiness_over_time) - min(readiness_over_time_wl)
  biggest_loose_delta <- min(readiness_over_time - readiness_over_time_wl)
  
  print("*** Readiness function report ***")
  if (biggest_loose_delta < 0) {
    has_loose_in_step <- (readiness_over_time - readiness_over_time_wl) < 0
    loose_periods <- list()
    period_start <- if (has_loose_in_step[1]) 1 else NULL
    for (i in 2:length(has_loose_in_step)) {
      if (has_loose_in_step[i - 1] != has_loose_in_step[i]) {
        if (is.null(period_start)) {
          period_start <- i
        } else {
          loose_periods <- append(loose_periods, paste("[", period_start, "-", i, "]", sep=""))
          period_start <- NULL
        }
      }
      if (i == length(has_loose_in_step) && !is.null(period_start)) {
        loose_periods <- append(loose_periods, paste("[", period_start, "-", i, "]", sep=""))
      }
    }
    print(paste("Model with learning has loose periods:", paste(loose_periods, collapse = ',')))
    print(paste("Biggest loose delta :      ", biggest_loose_delta)) 
  } else {
    print("Model with learning step always wins")
  }
  print(paste("min(A(t)) with learning:   ", min(readiness_over_time)))
  print(paste("min(A(t)) without learning:", min(readiness_over_time_wl)))
  print(paste("min(A(t)) delta:           ", lerning_min_to_without_learning_min))
}

################################### Helper functions ###################################