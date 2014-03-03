
#this command can be used to generate a histogram in the command line
#Input required is a 2 column file where the second column are numbers to bin for the histogram
#cat infile | createhist.awk binsize
#where binsize is the size of the bins on the X access
#This works well for assessing distribution of genome scaffold sizes
#requires that you pipe through a sort -k 1n to get the histogram in chronological order 
awk 'BEGIN { MIN=0; MAX=10000000000; OFS="\t"; BINSIZE="'$1'"; } 
{ A=sprintf("%d", $2/BINSIZE);
  BIN[A]++; } 
END { for (X in BIN) print X, BIN[X] }'
