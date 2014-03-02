
awk 'BEGIN { MIN=0; MAX=10000000000; OFS="\t"; BINSIZE="'$1'"; } 
{ A=sprintf("%d", $2/BINSIZE);
  BIN[A]++; } 
END { for (X in BIN) print X, BIN[X] }'
