###This script is part of perl pipeline developed by Michelle Graham, USDA-ARS (michelle.graham@ars.usda.gov).  
###Please see GO_analysis_README for usage.  
###If you use this pipeline (or modified versions) please cite: 
###van de Mortel M, Recknor JC, Graham MA, Nettleton D, Dittman JD, Nelson RT, Godoy CV, Abdelnoor RV, Almeida AM, 
### Baum TJ, Whitham SA (2007) Distinct biphasic mRNA changes in response to Asian soybean rust infection. Mol 
###Plant Microbe Interact. 20: 887-899  

#!/usr/bin/perl                                                                                                                                                        
($file1,$file2)=@ARGV;
open (WHOLECHIP, $file1) || die "Cannot open file $file1 \n";
while (<WHOLECHIP>) {
    @fields= split /\t/;
    $GOID=$fields[0];
    $WHOLECOUNT=$fields[1];
    chomp $GOID;
    chomp $WHOLECOUNT;
    open (DATA, $file2) || die "Cannot open file $file2\n";
    while (<DATA>) {
        if (/$GOID/) {
            @tmp=split /\t/;
            $name=$tmp[0];
            $DATACOUNT=$tmp[1];
            print "$GOID\t$WHOLECOUNT\t$DATACOUNT\n";
        }
        else {
#    print "$GOID\tNOT EXPRESSED\n";                                                                                                                                   
        }
    }
    close (DATA);
}
close (WHOLECHIP);
exit 0;



