check_chemform <-
function(
	isotopes,
	chemforms,
	get_sorted = FALSE,
	get_list = FALSE
){

    ############################################################################
    # internal function definitions ############################################
    # (A) Multiplier ###########################################################
    multif <-
    function(formula1, fact, numbers){
        formulas <- c()
        ########################################################################
        # on first chemical formula ############################################
        formula1 <- gsub("D", "[2]H", formula1)
        ende1 <- nchar(formula1)
        element1 <- number1 <- c()
        ########################################################################
        # on formula1 
        j <- 1
        while(j <= ende1){
          if(substr(formula1, j, j) == c("[")){
                b<-j
                while(any(substr(formula1, j, j) == c("]")) != TRUE){
                    j <- j + 1
                }
                k <- j
                while(any(substr(formula1, j, j) == numbers) != TRUE){
                    j <- j + 1
                }
                m <- j - 1
                element1 <- c(element1, substr(formula1, b, m))
          }
          if(any(substr(formula1,j,j)==numbers) != TRUE){
                k <- j
                while(any(substr(formula1, j, j) == numbers) != TRUE){
                  j<- j + 1
                }
                m <- j -1
                j <- j - 1
                element1 <- c(element1, substr(formula1, k, m))
          };
          if(any(substr(formula1, j, j) == numbers) == TRUE){
                k <- j
                while(any(substr(formula1, j, j) == numbers) == TRUE){
                  j <- j + 1
                }
                m <- j - 1
                j <- j - 1
                number1 <- c(number1, as.numeric(substr(formula1, k, m)))
          }
        j <- j + 1
        } # end while loop
        ########################################################################
        # multiply ! ###########################################################
        number1 <- fact * number1
        formula_fin <- ""
        for(p in 1:length(element1)) formula_fin <- paste0(formula_fin, element1[p], number1[p])
        formulas <- c(formulas, formula_fin)
      return(formulas)
      ##########################################################################
    }
    ############################################################################
    
    ############################################################################
    capitals <- c("[", LETTERS)
	numbers <- as.character(0:10)
    allpossible <- c(capitals, numbers, "(", ")", "]", letters)
    masses <- warn <- c()
    elem <- unique(as.character(isotopes[, 1]))
    isotopes2 <- matrix(nrow = length(elem), ncol = 2)
    isotopes2[, 1] <- elem                                     
    for (i in 1:length(elem)) {
        intermed <- isotopes[isotopes[, 1] == elem[i], ]
        if (is.vector(intermed) == TRUE) {
            isotopes2[, 2][isotopes2[, 1] == elem[i]] <- intermed[3]
        }
        else {
            isotopes2[, 2][isotopes2[, 1] == elem[i]] <- intermed[,
                3][as.numeric(intermed[, 4]) == max(as.numeric(intermed[,
                4]))]
        }
    }
    info <- isotopes2           
	if(get_list) listed <- vector("list", length(chemforms))
    for (i in 1:length(chemforms)) {
        
		mass <- 0
        warnit <- FALSE
        # (0) empty? ###########################################################
		if(chemforms[i] == "") warnit <- TRUE
		# (1) correct any round or missing brackets for isotopologues at formula beginning or after whitespaces within formula to square brackets
		# i.e., contains white spaces? -> only allowed to separate isotopologues, e.g. [C10 15N], [C10 (15)N] or [C10 [15]N] 		
		# (2) Then remove any white spaces
		round_brackets <- c("(15)N", "(12)C", "(13)C", "(35)Cl", "(37)Cl", "(16)O", "(18)O", "(2)H", "(33)S", "(34)S", "(35)S", "(36)S")
		missing_brackets <- c("15N", "12C", "13C", "35Cl", "37Cl", "16O", "18O", "2H", "33S", "34S", "35S", "36S")
		correct_brackets <- c("[15]N", "[12]C", "[13]C", "[35]Cl", "[37]Cl", "[16]O", "[18]O", "D", "[33]S", "[34]S", "[35]S", "[36]S")
		if(!identical(length(round_brackets), length(missing_brackets), length(correct_brackets))) stop("Debug_1 check_chemform")
		chemforms[i] <- trimws(chemforms[i])										# (a) trim whitespace at formula beginning and end
		chemforms[i] <- gsub("(^\\s+)|(\\s+$)", " ", chemforms[i], fixed = TRUE)	# (b) correct tab stops etc to whitespaces
		chemform_split <- strsplit(chemforms[i], " ", fixed = TRUE)[[1]]			# (c) split by whitespaces
		for(m in 1:length(chemform_split)){
			for(n in 1:length(round_brackets)){
				has_char <- nchar(round_brackets[n])
				if(nchar(chemform_split[m]) >= has_char) if(substr(chemform_split[m], 1, has_char) == round_brackets[n]){
					chemform_split[m] <- sub(round_brackets[n], correct_brackets[n], chemform_split[m], fixed = TRUE)
				}
			}
			for(n in 1:length(missing_brackets)){
				has_char <- nchar(missing_brackets[n])
				if(nchar(chemform_split[m]) >= has_char) if(substr(chemform_split[m], 1, has_char) == missing_brackets[n]){
					chemform_split[m] <- sub(missing_brackets[n], correct_brackets[n], chemform_split[m], fixed = TRUE) 
				}
			}
		}
		chemforms[i] <- paste(chemform_split, collapse = "")						# (d) removes empty spaces inside formula
		# ambiguous round brackets in front of isotopologues still existing?
		for(check_bracket in round_brackets) if(grepl(check_bracket, chemforms[i], fixed = TRUE)) warnit <- TRUE
        ##### split string #####################################################
        formel <- as.character(chemforms[i])
        formel <- strsplit(formel, " ")[[1]]
        m <- strsplit(formel, as.character())[[1]]
        # (3) all characters plausible? ########################################
        if(!warnit) for(j in 1:length(m)) if(any(allpossible == m[j]) == FALSE) warnit <- TRUE
        # (4) do all [(bracket)] types close & [] only contain numbers? ########
        if(!warnit){        
          if( any(m == "[") || any(m == "]") || any(m == "(") || any(m == ")") ){        
            getit1 <- getit2 <- 0
            a <- 1
            while((a) <= length(m)){
              if(m[a] == "[") getit1 <- getit1 + 1
              if(m[a] == "]") getit1 <- getit1 - 1
              if(m[a] == "(") getit2 <- getit2 + 1
              if(m[a] == ")") getit2 <- getit2 - 1               
              if(getit1 & (any(numbers == m[a]) == FALSE & m[a] != "[" & m[a] != "]") ) warnit <- TRUE
              if(getit1 < 0 | getit2 < 0) warnit <- TRUE
              a <- a + 1
            }
            if(getit1 != 0 | getit2 != 0) warnit <- TRUE
            
          }
        }
        # (5) start correct? ###################################################
        if(!warnit){
          if(!m[1] %in% c(capitals, "(", ")")) warnit <- TRUE
          if(length(m) == 1){
            m <- c(m, "1")
            formel <- paste0(formel, "1")
          }
        }
        # (6) empty brackets? ##################################################
        if(!warnit) for(k in 2:length(m)) if(m[k - 1] == "(" & m[k] == ")") warnit <- TRUE
        # (7) insert 1 where missing ###########################################
        if(!warnit){
          # for closing )-brackets #############################################  
          if(any(m == "(")){
            m2 <- c()
            for(j in 1:(length(m) - 1)){
              m2 <- if(
				(m[j] == ")") & 
                !any(m[j + 1] == numbers)
              ) c(m2, m[j], "1") else c(m2, m[j])
            }
            m2 <- c(m2, m[length(m)])
            if(m[length(m)] == ")") m2 <- c(m2, "1")
            m <- m2
          }
          # for all other cases ################################################
          m2 <- m[1]
          for(j in 2:length(m)){
            m2 <- if(
              (any(m[j] == capitals) | m[j] == ")" | m[j] == "(" ) &
              all(m[j - 1] != numbers) & 
			  m[j - 1] != "(" & 
			  m[j - 1] != "]"      
            ) c(m2,"1",m[j]) else c(m2,m[j])
          }
          if(all(m[length(m)] != numbers)) m2 <- c(m2, "1")
          m <- m2
          formel <- ""
          for(k in 1:length(m)) formel <- paste0(formel, m[k])
        }
        # (8) multiply for square brackets, with nesting #######################
        if(!warnit){
          while(any(m == "(")){
            a <- getit1 <- getit2 <- 1
            while( getit1 != 0 & getit2 != 0 & a <= length(m)){              
              if(m[a] == "("){
                getit1 <- 2  
                from <- a
              }
              if(m[a] == ")"){
                getit2 <- 2  
                to <- a
              }
              if(getit1 == 2 & getit2 == 2){
                 b <- a + 1
                 count <- ""
                 while(any(m[b] == numbers & b <= length(m))){
                    count <- paste0(count, m[b])
                    b <- b + 1
                  }
                  count <- as.numeric(count)
                  m2 <- ""
                  for(k in (from + 1):(to - 1)){
                    m2 <- paste0(m2, m[k])
                  }
                  m2 <- multif(m2, count, numbers)
                  m2 <- strsplit(m2, as.character())[[1]]
                  m3 <- c()
                  doneit <- FALSE
                  for(z in 1:length(m)){
                    if( z < from || z >= b){
                      m3 <- c(m3, m[z])
                    }else{
                      if(doneit == FALSE & (z >= from | z < b)){
                        m3 <- c(m3, m2)                     
                        doneit <- TRUE
                      }
                    }
                  }
                  m <- m3
                  getit1 <- getit2 <- 0
              }
              a <- a + 1
            }
          }
          formel<- ""
          for(k in 1:length(m)) formel <- paste0(formel, m[k])
        }
        # (9) dissassemble #####################################################
        if(!warnit){
          element1 <- number1 <- c()
          ######################################################################
          j <- 1
          while(j <= nchar(formel)){
            if(substr(formel, j, j) == c("[")){
                  b <- j
                  while(
                    any(substr(formel, j, j) == c("]")) != TRUE &
                    j <= nchar(formel)
                  ){
                      j <- j + 1
                  }
                  k <- j
                  while(any(substr(formel,j,j)==numbers) != TRUE){
                      j <- j + 1
                  }
                  z <- j - 1
                  element1 <- c(element1, substr(formel, b, z))
            }
            if(any(substr(formel, j, j) == numbers) != TRUE){
                  k <- j
                  while(
                    any(substr(formel, j, j) == numbers) != TRUE &
                    j <= nchar(formel)
                  ){
                    j <- j + 1
                  }
                  z <- j - 1
                  j <- j - 1
                  element1 <- c(element1, substr(formel, k, z))
            }
            if(any(substr(formel, j, j) == numbers) == TRUE){
                  k <- j
                  while(
                    any(substr(formel, j, j) == numbers) == TRUE &
                    j <= nchar(formel)
                  ){
                    j <- j + 1
                  }
                  z <- j - 1
                  j <- j - 1
                  number1 <- c(number1, as.numeric(substr(formel, k, z)))
            }
          j <- j + 1
          }# end while
        }
        # (10) check if all elements present in isotopes list ##################
        if(!warnit){
          for(j in 1:length(element1)){ 
            if(any(element1[j] == as.character(isotopes[, 1])) == FALSE) warnit <- TRUE
          }
          if(length(element1) != length(number1)) warnit <- TRUE
        }
        # (11) merge non-unique elements #######################################
        if(!warnit){
          element2 <- number2 <- c()
          doneit <- rep(FALSE, length(element1))
          for(j in 1:length(element1)){
            if(!doneit[j]){
              doneit[element1 == element1[j]] <- TRUE
              element2 <- c(element2, element1[element1 == element1[j]][1])
              number2 <- c(number2, as.character(sum(as.numeric(number1[element1 == element1[j]]))))
            }
          }
          element1 <- element2
		  rm(element2)
          number1 <- number2
		  rm(number2)
		  if(get_sorted){ # ensure unambiguous order of elements in the formula
				this <- order(match(element1,info))
				number1 <- number1[this]
				element1 <- element1[this]
		  }
          formel <- ""
          for(k in 1:length(element1)){
            formel <- paste0(formel, element1[k], number1[k])
            mass <- mass + (
              as.numeric(info[info[, 1] == element1[k], 2][1]) * as.numeric(number1[k])
            )  
          }
        }        
        ########################################################################
        # (12) make final entry ################################################
        if(!warnit){
			if(!get_list){
				warn <- c(warn, FALSE)
				masses <- c(masses, mass)
				chemforms[i] <- formel
			}else{
				number1 <- as.numeric(number1)
				names(number1) <- element1
				listed[[i]] <- number1
				names(listed)[i] <- chemforms[i]
			}			
        }else{
			if(!get_list){
				warn <- c(warn, TRUE)
				masses <- c(masses, -9999) 
			}else{
				listed[[i]] <- numeric()
				names(listed)[i] <- "invalid formula"
			}		 
		}
        ########################################################################
		 
	}    
    ############################################################################  
	if(!get_list){	
		checked <- data.frame(warn, chemforms, masses)
		names(checked) <- c("warning", "new_formula", "monoisotopic_mass")
		checked[, 2] <- as.character(checked[, 2])
		return(checked)
	}else return(listed)
	
}        
        
        
        
           
      
        
        
        
        
        
        
        
      