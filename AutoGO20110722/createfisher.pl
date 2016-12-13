###This script is part of perl pipeline developed by Michelle Graham, USDA-ARS (michelle.graham@ars.usda.gov).  
###Please see GO_analysis_README for usage.  
###If you use this pipeline (or modified versions) please cite: 
###van de Mortel M, Recknor JC, Graham MA, Nettleton D, Dittman JD, Nelson RT, Godoy CV, Abdelnoor RV, Almeida AM, 
### Baum TJ, Whitham SA (2007) Distinct biphasic mRNA changes in response to Asian soybean rust infection. Mol 
###Plant Microbe Interact. 20: 887-899  


#!/usr/bin/perl                                                                                                                                                        
($file1)=@ARGV;
open (DATA, $file1) || die "Cannot open file $file1 \n";
@chip='';
@exp='';
while (<DATA>) {
    @fields= split /\t/;
    $chipcount=$fields[1];
    $expressedcount=$fields[2];
    chomp $chipcount;
    chomp $expressedcount;
    push (@chip, $chipcount);
    push (@exp, $expressedcount);
    $chipcount='';
    $expressedcount='';
}
close DATA;
$sumchip = 0;
for ( @chip) {
    $sumchip += $_;
}
$sumexp=0;
for (@exp) {
    $sumexp += $_;
}
$chipcount='';
$expressedcount='';
#print "SUM CHIP $sumchip SUM EXP $sumexp\n";
open (STATS, $file1) || die "Cannot open file $file1\n";
while (<STATS>) {
    @tmp= split /\t/;
    $GOID=$tmp[0];
    $GOTOTAL=$tmp[1];
    $EXPGO=$tmp[2];
    chomp $GOID;
    chomp $GOTOTAL;
    chomp $EXPGO;
    $NOEXPGO=$GOTOTAL - $EXPGO;
    $EXPNOGO=$sumexp-$EXPGO;
    $NOEXPNOGO=$sumchip - $GOTOTAL;
    if($GOID =~ /\D{1}/) {
	print "$GOID\t$EXPGO\t$NOEXPGO\t$EXPNOGO\t$NOEXPNOGO\n";
#2x2 contingency table is consistent to the following testing if the number of on particular GO term represented in the genes that are Differentially Expressed are chosen at random what probablility that number would be found among the total GO terms 
#DE        not DE
#DE GO     not DE GO  
#DE not GO not DE not GO
    }
}
close (STATS);
exit 0;
















