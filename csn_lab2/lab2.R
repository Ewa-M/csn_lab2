degree_sequence = read.table("./data/English_in-degree_sequence.txt",
                             header = FALSE)

require("stats4") # for MLE
require("VGAM") # for the Riemann-zeta functiom