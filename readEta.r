readEta <- function(etafile = cfg$etafile)
{
etadata <- read.csv(etafile, header=TRUE)

#Build eta list structure
eta <- list()
eta$aa <- substr(etadata$Parameter[grep("Delta",etadata$Parameter)],1,1)
eta$codon <- substr(etadata$Parameter[grep("Delta",etadata$Parameter)],3,5)
eta$values <- etadata$Value[grep("Delta",etadata$Parameter)]

#The delta eta values are relative to a certain codon. These are those codons.
uniquecodon <- unique(eta$aa)

defaultcodon <- c("GCT", "TGT", "GAT", "GAG", "TTT", "GGT", "CAT", "ATT", "AAG", "TTG", "AAT", "CCT", "CAT", "CGT", "TCT", "ACT", "GTT", "TAT", "AGT")


x <- 1;


minmatrix <- numeric();
maxmatrix <- numeric();



#loop through the entire genome
while( x <= length(eta$aa) ){

stop <- x;	#stop is the first element of the current group of synonyms
mineta <- 0;
maxeta <- 0;

	#loop through each group of synonyms
	while ( eta$aa[x] == eta$aa[stop] 
			&& x <= length(eta$aa) ) {
		if(eta$values[x] < mineta)
		{
			mineta <- eta$values[x];
			mincodon <- eta$codon[x];
		}else if(eta$values[x] > maxeta)
		{
			maxeta <- eta$values[x];
			maxcodon <- eta$codon[x];
		}

		x <- x+1;
	}
	#Each group has ((x-1)-stop) synonyms

	if(mineta == 0){
		mincodon <- defaultcodon[grep(eta$aa[stop], uniquecodon)]
	}
	if(maxeta == 0){
		maxcodon <- defaultcodon[grep(eta$aa[stop], uniquecodon)]
	}

	
minvec <- c(eta$aa[stop], mincodon, mineta);
maxvec <- c(eta$aa[stop], maxcodon, maxeta);

minmatrix <- c(minmatrix, minvec);
maxmatrix <- c(maxmatrix, maxvec);

}

minmatrix <- matrix(minmatrix, ncol=length(minvec), byrow=T);
colnames(minmatrix) <- c("AminoAcid","Codon","ROC");
maxmatrix <- matrix(maxmatrix, ncol=length(maxvec), byrow=T);
colnames(maxmatrix) <- c("AminoAcid","Codon","ROC");

matrices <- list();
matrices$minmatrix <- minmatrix;
matrices$maxmatrix <- maxmatrix;

	return(matrices)
}
