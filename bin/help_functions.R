# line tool 
# d <- attr(a2[[1]][["svg"]][["path"]],"d")
# bez2poly <- function(x){
#   #x <- gsub("[ ].*[ ]","l", x) # removes space previous version might be nessecary for something
#   x <- gsub("(?<=\\d) +(?=\\d)",",", x, perl = T) # replaces space " " with "," to separate x and y coord
#   x <- strsplit( d, split = "[a-zA-Z]")[[1]][-1]
#   x <- gsub(" ","", x)
#   x <- strsplit(x,split = ",")[]
#   x <- tibble(x = map_chr(x,~.x[1]), y = map_chr(x,~.x[2]) )
#   x <- apply(x,1,as.numeric)
#   x <- cbind(cumsum(x[,1]), cumsum(x[,2]) )
#   x <- rbind(x,x[1,])
#   colnames(x) <- c("x", "y")
#   return(x)
# }

# use this version if problems this has been working before 
bez2poly <- function(x){
  x <- strsplit( gsub("[ ].*[ ]","l", x) ,split = "[a-zA-Z]")[[1]][-1]
  x <- as.data.frame(strsplit(x,split = ","))
  x <- apply(x,1,as.numeric)
  x <- cbind(cumsum(x[,1]), cumsum(x[,2]) )
  x <- rbind(x,x[1,])
  colnames(x) <- c("x", "y")
  return(x)
}

tmat <- function(x,n=3){
  tmp <- matrix(0,3,3)
  x <- as.numeric(strsplit(gsub("[a-z()]","",x),split = ",")[[1]])
  tmp[c(1,2),c(1,1)] <- x[1:2]
  tmp[c(1,2),c(2,2)] <- x[3:4]
  tmp[c(1,2),c(3,3)] <- x[5:6]
  tmp[3,3] <- 1
  return(tmp)
}

rec2poly <- function(x){
  x <- attributes(x)
  x <- setNames( as.numeric(unlist(x)[c("x","y","width","height")]) , c("x","y","width","height") )
  df <- data.frame(
    x = x["x"] + x["width"] * c(0,0,1,1,0),
    y = x["y"] + x["height"] * c(0,1,1,0,0))
  return(df)
}

get_shape <- function(x){
  xx <- attributes(x)
  if( "d" %in% names(xx) ){
    return( bez2poly( attr(x,"d") ) )
  } else if( "width"  %in% names(xx) ){
    return(rec2poly( x ))
  } 
}

poly_area <- function(df){
  res <- t(df)
  x <- abs( sum(  ((res[1,] * c(res[2,-1],res[2,1])) - (res[2,] * c(res[1,-1],res[1,1]) )) )/2 )
  return(x)
}
