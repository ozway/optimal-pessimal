source("config.r")
source("readEta.r")
source("writing.r")

matrices <- readEta();

worstEta <- matrices$pessimalMatrix;    #worstEta means the least fit codons

print(worstEta);

aminoacid <- worstEta[,1];

#I'm only interested in codons that have synonyms. They're the only ones I have eta values for.
synonyms <-  list(
    A = c("GCA", "GCC", "GCG", "GCT"),
    C = c("TGC", "TGT"),
    D = c("GAC", "GAT"),
    E = c("GAA", "GAG"),
    F = c("TTC", "TTT"),
    G = c("GGA", "GGC", "GGG", "GGT"),
    H = c("CAC", "CAT"),
    I = c("ATA", "ATC", "ATT"),
    K = c("AAA", "AAG"),
    L = c("CTA", "CTC", "CTG", "CTT", "TTA", "TTG"),
#    M = c("ATG"),
    N = c("AAC", "AAT"),
    P = c("CCA", "CCC", "CCG", "CCT"),
    Q = c("CAA", "CAG"),
    R = c("AGA", "AGG", "CGA", "CGC", "CGG", "CGT"),
    S = c("TCA", "TCC", "TCG", "TCT"),  # split 2 codons to Z.
    T = c("ACA", "ACC", "ACG", "ACT"),
    V = c("GTA", "GTC", "GTG", "GTT"),
#    W = c("TGG"),
#    X = c("TAA", "TAG", "TGA"),  # stop codons.
    Y = c("TAC", "TAT"),
    Z = c("AGC", "AGT")  # obtain 2 codons from S.
)



sequence <- read.seq(cfg$genome);
pessimal <- sequence;

downgradeCount <- 0;

#for(gene in 1:3){
for(gene in 1:length(sequence)){

for(index in 1:length(sequence[[gene]])/3){
	index <- index*3 + 1;

	#find out which synonym group the codon is in
	temp <- paste(sequence[[gene]][index], sequence[[gene]][index+1], sequence[[gene]][index+2], sep="");

	for(j in 1:length(synonyms)){
		if(temp %in% synonyms[[j]]){
			break;
		}else if(j == length(synonyms)){j=0; break;}
	}
	

        #If the codon has a synonym...
	if(j != 0){

if(substr(worstEta[j,2],1,3) !=
paste(sequence[[gene]][index], sequence[[gene]][index+1], sequence[[gene]][index+2], sep="")){
                        downgradeCount <- downgradeCount + 1;
        }else{didntDowngrade <- didntDowngrade + 1}


		pessimal[[gene]][index] <- substr(worstEta[j,2],1,1);
		pessimal[[gene]][index+1] <- substr(worstEta[j,2],2,2);
		pessimal[[gene]][index+2] <- substr(worstEta[j,2],3,3);
	}

}#end this gene

}#end the genome

write.seq(pessimal, cfg$pessimalfile)

print( paste(downgradeCount, "codons downgraded",
        (100 * downgradeCount) %/% (downgradeCount+didntDowngrade), "%") );

