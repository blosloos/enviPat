.onAttach <- function(lib, pkg)
{
	packageStartupMessage(paste("\n \n Welcome to enviPat version 1.9 \n Check www.envipat.eawag.ch for an interactive online version\n",sep=""));

}

#.onLoad<-function(lib, pkg)
#{
  	#library.dynam(iso_pattern_Call, lib, pkg);
#}