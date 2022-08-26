multiform <-
function(formula_in, fact){

	formula_out <- c()
	for(i in 1:length(formula_in)){
  
		##############################################################################
		# on first chemical formula
		formula_in[i] <- gsub("D(?![a-z])","[2]H", formula_in[i], perl = T)
		len <- nchar(formula_in[i])
		element1 <- number1 <- c()
		##############################################################################
		# on formula_in (deduct)
		j <- 1
		while( j <= len){
	
			if(substr(formula_in[i], j, j) == "["){
				b <- j
				while(!any(substr(formula_in[i], j, j) == "]")) j <- j + 1
				k <- j
				while(!any(substr(formula_in[i], j, j) == c("0","1","2","3","4","5","6","7","8","9"))) j <- j + 1
				m <- j - 1
				element1 <- c(element1, substr(formula_in[i], b, m))
			}
	  
			if(!any(substr(formula_in[i], j, j) == c("0","1","2","3","4","5","6","7","8","9"))){
				k <- j
				while(!any(substr(formula_in[i], j, j) == c("0","1","2","3","4","5","6","7","8","9"))) j <- j + 1
				m <- j <- j - 1
				element1 <- c(element1, substr(formula_in[i], k, m))
			}
		  
			if(any(substr(formula_in[i],j,j) == c("0","1","2","3","4","5","6","7","8","9"))){
				k <- j
				while(any(substr(formula_in[i], j, j) == c("0","1","2","3","4","5","6","7","8","9"))) j <- j + 1
				m <- j <- j - 1
				number1 <- c(number1, as.numeric(substr(formula_in[i], k, m)))
			}
	  
			j <- j + 1
		}
		##############################################################################
		# multiply
		number1 <- fact * number1
		formula_fin <- ""
		for(i in 1:length(element1)) formula_fin <- paste0(formula_fin, element1[i], number1[i])
		formula_out <- c(formula_out, formula_fin)
	
	}
	
	return(formula_out)

}

