SetTextContrastColor <- function(color){
  ifelse( mean(col2rgb(color)) > 127, "black", "white")
}