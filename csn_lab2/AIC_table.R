source("displaced.R")
source("plot.R")

write_summary <- function(language, file) {
  degree_sequence = read.table(file, header = FALSE)
  x <- degree_sequence$V1
  result_geometric <- calculate_geometric(x)
  result_zeta_gamma_2 <- calculate_zeta_gamma_2(x)
  result_zeta <- calculate_zeta(x)
  result_poisson <- calculate_poisson(x)
  #result_zeta_right_truncated <- calculate_zeta_right_truncated(x)

  cat(language, "&", # language
      result_poisson$param, "&", 
      result_geometric$param, "&", 
      result_zeta$param, "&",
      #result_zeta_right_truncated$param, "&", #gamma_2
      "gamma_2 &",
      "k_max", "&", #k_max
      "\\\\ \n")
}

write_summary_AIC <- function(language, file) {
  degree_sequence = read.table(file, header = FALSE)
  x <- degree_sequence$V1
  result_geometric <- calculate_geometric(x)
  result_zeta_gamma_2 <- calculate_zeta_gamma_2(x)
  result_zeta <- calculate_zeta(x)
  result_poisson <- calculate_poisson(x)
  #result_zeta_right_truncated <- calculate_zeta_right_truncated(x)
  
  cat(language, "&", # language
      result_poisson$AIC, "&", 
      result_geometric$AIC, "&", 
      result_zeta_gamma_2$AIC, "&", 
      result_zeta$AIC, 
      "\\\\ \n")
}

source = read.table("list_in.txt", 
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)

print("Language & lambda & gamma_1 & gamma_2 & k_max")

for (x in 1:nrow(source)) {
  write_summary(source$language[x], source$file[x])
}

print("AIC")

for (x in 1:nrow(source)) {
  write_summary_AIC(source$language[x], source$file[x])
}
