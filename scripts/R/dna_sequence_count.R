#Asking user an input
x <- readline("Enter a DNA sequence:")

#Conversion of input to uppercase if it's lowercase
x <- toupper(x)

#Splitting into characters
dna_char <- strsplit(x, "")[[1]]

#Counting of nucleotide bases
a <- sum(dna_char == "A")
t <- sum(dna_char == "T")
g <- sum(dna_char == "G")
c <- sum(dna_char == "C")

#Printing results
cat("A=", a, "T=", t, "G=", g, "C=", c, "\n")