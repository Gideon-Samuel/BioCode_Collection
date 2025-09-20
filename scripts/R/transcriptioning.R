# Asking user for an input
x <- readline("Enter a DNA sequence: ")

#Conversion of input to uppercase if it's lowercase
x <- toupper(x)

# Replace T with U
y <- chartr("T", "U", x)

cat("Transcribed RNA sequence:", y, "\n")
