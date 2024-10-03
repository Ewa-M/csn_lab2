plot_poisson <- function(degree_sequence, language, param) {
  df = as.data.frame(prop.table(table(degree_sequence)))
  x = strtoi(df$V1)
  y = y=df$Freq

  x_spectrum = min(x):max(x)
  
  plot(x, y, log="xy", main=cat(language, " - Poisson"))
  lines(x_spectrum, dpois(x_spectrum, lambda=param), col="red")
}

plot_geometric <- function(degree_sequence, language, param) {
  df = as.data.frame(prop.table(table(degree_sequence)))
  x = strtoi(df$V1)
  y = y=df$Freq
  
  x_spectrum = min(x):max(x)
  
  plot(x, y, log="xy", label=cat(language, " - Geometric"))
  lines(x_spectrum, dgeom(x_spectrum, param), col="red")
}


plot_zeta <- function(degree_sequence, language, param) {
  df = as.data.frame(prop.table(table(degree_sequence)))
  x = strtoi(df$V1)
  y = y=df$Freq
  
  x_spectrum = min(x):max(x)
  
  plot(x, y, log="xy", main=cat(language, " - Zeta"))
  lines(x_spectrum, dzeta(x_spectrum, param), col="red")
}