raster2df <- function(x){
  df <- data.frame(rasterToPoints(x))
  t <- names(x)
  t <- sapply(t, function(x) stringr::str_split(x, 'X')[[1]][2])
  t <- as.character(gsub('\\.', '-', t))
  colnames(df)[3:length(df)] <- t
  df <- df %>% tidyr::gather("date", 'v1', 3:length(df))
  return(df)
}