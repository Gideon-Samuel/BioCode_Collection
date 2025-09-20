#Using R library (Biostrings) for handling strings.
library(Biostrings)

#Asking for DNA sequence as input.
x <- readline("Give DNA sequence: ")

#DNAString() wraps into a Biostrings object.
y <- DNAString(x)

#Reverse Complementing.
z <- reverseComplement(y)

#Printing results.
cat("Reverse Complement :", as.character(z), "\n")