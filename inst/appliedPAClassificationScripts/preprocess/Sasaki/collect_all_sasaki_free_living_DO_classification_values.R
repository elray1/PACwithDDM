DO_file_path_base <- "C:/Kinesiology/JefferData/Free Living data/ACE DO data/ACE"
DO_file_path_termination <- "DO.txt"

subject_numbers <- c("01", "04", "06", "08", "11", "19", "20", "21", "22", "23", "24", "27", "32", "33", "34")

behavior_values <- c()
behavior_mod_1_values <- c()
behavior_mod_2_values <- c()

for(subj in subject_numbers) {
  DO_file_path <- paste0(DO_file_path_base, subj, DO_file_path_termination)
  DO_data <- read.csv(file = file(DO_file_path), skip = 4, stringsAsFactors = FALSE)
  
  if("" %in% DO_data$Behavior)
    browser()
  
  if("Light" %in% DO_data$Behavior.Modifier.1 || "Moderate" %in% DO_data$Behavior.Modifier.1)
    browser()
  
  if("Indoor" %in% DO_data$Behavior.Modifier.2)
    browser()
  
  behavior_values <- unique(c(behavior_values, DO_data$Behavior))
  behavior_mod_1_values <- unique(c(behavior_mod_1_values, DO_data$Behavior.Modifier.1))
  behavior_mod_2_values <- unique(c(behavior_mod_2_values, DO_data$Behavior.Modifier.2))
}

combined_data_split <- lapply(subject_numbers, function(subj) {
  subj <- subject_numbers[[1]]
  preprocess_one_free_living_file(actigraph_file_path = paste0(actigraph_file_path_base, subj, actigraph_file_path_termination),
    DO_file_path = paste0(DO_file_path_base, subj, DO_file_path_termination),
    sampling_freq = 80)
})
