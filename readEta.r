readEta <- function(etafile = cfg$etafile) {

etadata <- read.csv(etafile, header=TRUE)

#Build eta list structure
eta <- list()
eta$aa <- substr(etadata$Parameter[grep("Delta",etadata$Parameter)],1,1)
eta$codon <- substr(etadata$Parameter[grep("Delta",etadata$Parameter)],3,5)
eta$values <- etadata$Value[grep("Delta",etadata$Parameter)]

#The delta eta values are relative to a certain codon. These are those codons.
uniquecodon <- unique(eta$aa)

defaultcodon <- c("GCT", "TGT", "GAT", "GAG", "TTT", "GGT", "CAT", "ATT", "AAG", "TTG", "AAT", "CCT", "CAG", "CGT", "TCT", "ACT", "GTT", "TAT", "AGT")


x <- 1;


pessimalMatrix <- numeric();
optimalMatrix <- numeric();



#loop through the entire genome
while( x <= length(eta$aa) ){

stop <- x;	#stop is the first element of the current group of synonyms
pessimalEta <- 0;
optimalEta <- 0;

	#loop through each group of synonyms
	while ( eta$aa[x] == eta$aa[stop] 
			&& x <= length(eta$aa) ) {
		if(eta$values[x] > pessimalEta)
		{
			pessimalEta <- eta$values[x];
			pessimalCodon <- eta$codon[x];
		}else if(eta$values[x] < optimalEta)
		{
			optimalEta <- eta$values[x];
			optimalCodon <- eta$codon[x];
		}

		x <- x+1;
	}
	#Each group has ((x-1)-stop) synonyms

	if(pessimalEta == 0){
		pessimalCodon <- defaultcodon[grep(eta$aa[stop], uniquecodon)]
	}
	if(optimalEta == 0){
		optimalCodon <- defaultcodon[grep(eta$aa[stop], uniquecodon)]
	}

	
pessimalVec <- c(eta$aa[stop], pessimalCodon, pessimalEta);
optimalVec <- c(eta$aa[stop], optimalCodon, optimalEta);

pessimalMatrix <- c(pessimalMatrix, pessimalVec);
optimalMatrix <- c(optimalMatrix, optimalVec);

}

pessimalMatrix <- matrix(pessimalMatrix, ncol=length(pessimalVec), byrow=T);
colnames(pessimalMatrix) <- c("AminoAcid","Codon","Eta");
optimalMatrix <- matrix(optimalMatrix, ncol=length(optimalVec), byrow=T);
colnames(optimalMatrix) <- c("AminoAcid","Codon","Eta");

matrices <- list();
matrices$pessimalMatrix <- pessimalMatrix;
matrices$optimalMatrix <- optimalMatrix;

	return(matrices)
}
