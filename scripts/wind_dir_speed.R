
# Developed by Tom Evans at Lund University: tom.evans@biol.lu.se
# You are welcome to use parts of this code, but please give credit when using it extensively.
# Code available at https://github.com/thomasevans/lbbg_gps

# Wind shear -------
# Calculate wind direction and speed from u and v wind components


wind.dir.speed <- function(uwind10, vwind10){
  # This function calculates the wind speed and direction based on the u
  # v wind vectors
  
  if(is.na(uwind10) | is.na(vwind10)) return(t(c(NA,NA))) else {
    
    #Wind speed Pythagoras theorem
    wind.speed <- sqrt((uwind10 * uwind10) + (vwind10 * vwind10))
    
    # Calculate direction in radians (0 - 90 deg)
    dir <- atan(abs(uwind10/ vwind10))
    
    #   atan(1)
    #   atan(0.5)
    #   dir <- atan(0.5)
    #   ?atan
    # Direction in degrees (0 - 90)
    dir <- dir * 180 / pi
    
    # Make into bearing from North
    if(uwind10 > 0 && vwind10 < 0){
      wind.dir <- (180 - dir)
    }else if(uwind10 < 0 && vwind10 < 0){
      wind.dir <- (dir + 180)
    }else if(uwind10 < 0 && vwind10 > 0){
      wind.dir <- (360 - dir)
    }else   wind.dir <- (dir)
    
    x <- cbind(wind.speed, wind.dir)
    return(x)
  }
}