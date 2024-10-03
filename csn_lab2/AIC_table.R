write_summary <- function(language,file) {
  degree_sequence = read.table(file, header = FALSE)
  dg = displaced_geometric(degree_sequence)
  plot_geometric(degree_sequence, language, dg$param)
  
  cat(language, "&", # language
      dg$AIC,
      dg$param,
      "\\\\ \n")
}

source = read.table("list_geometric.txt", 
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)
for (x in 1:nrow(source)) {
  write_summary(source$language[x], source$file[x])
}

