multiform<-
function(formula1,fact){

  formulas<-c();
  for(i in 1:length(formula1)){
  
    ##############################################################################
    # on first chemical formula
    formula1[i] <- gsub("D", "[2]H", formula1[i])
    ende1 <- nchar(formula1[i])
    element1 <- number1 <- c()
    ##############################################################################
    # on formula1 (deduct)
    j <- 1
    while( j <= ende1){
      if(substr(formula1[i], j, j) == "["){
            b <- j
            while(!any(substr(formula1[i], j, j) == "]")) j <- j + 1
            k <- j
            while(!any(substr(formula1[i], j, j) == c("0","1","2","3","4","5","6","7","8","9"))) j <- j + 1
            m <- j - 1
            element1 <- c(element1,substr(formula1[i], b, m))
      }
      if(!any(substr(formula1[i], j, j) == c("0","1","2","3","4","5","6","7","8","9"))){
            k <- j
            while(!any(substr(formula1[i], j, j) == c("0","1","2","3","4","5","6","7","8","9"))) j <- j + 1
            m <- j <- j - 1
            element1 <- c(element1, substr(formula1[i], k, m))
      }
      if(any(substr(formula1[i],j,j) == c("0","1","2","3","4","5","6","7","8","9"))){
            k <- j
            while(any(substr(formula1[i], j, j) == c("0","1","2","3","4","5","6","7","8","9"))) j <- j + 1
            m <- j <- j - 1
            number1 <- c(number1, as.numeric(substr(formula1[i], k, m)))
      }
    j <- j + 1
    }
    ##############################################################################
    # multiply
    number1 <- fact * number1
    formula_fin <- ""
    for(i in 1:length(element1)) formula_fin <- paste0(formula_fin, element1[i], number1[i])
    formulas <- c(formulas, formula_fin)
  }
  return(formulas)
  ##############################################################################

}

