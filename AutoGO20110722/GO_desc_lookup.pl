###This script is part of perl pipeline developed by Michelle Graham, USDA-ARS (michelle.graham@ars.usda.gov).  
###Please see GO_analysis_README for usage.  
###If you use this pipeline (or modified versions) please cite: 
###van de Mortel M, Recknor JC, Graham MA, Nettleton D, Dittman JD, Nelson RT, Godoy CV, Abdelnoor RV, Almeida AM, 
### Baum TJ, Whitham SA (2007) Distinct biphasic mRNA changes in response to Asian soybean rust infection. Mol 
###Plant Microbe Interact. 20: 887-899  


#!/usr/bin/perl
($file1, $file2)= @ARGV;
open (DATA, $file1) || die "Cannot open file $file2 \n";
while (<DATA>) {
    @tmp= split /\t/;
    $go_desc=$tmp[4];
    $go_id=$tmp[5];
    chomp $go_id;
    chomp $go_id;
    foreach $go_id (@tmp) {
	if ($go_id ne "") {
	    $data_lookup{$go_id} = $go_desc;
	}
	
    }

    $go_desc='';
    $go_id='';
}


close DATA;

open (GOIDLIST,$file2) || die "Cannot open $infile1 \n";
while (<GOIDLIST>) {
    @fields= split /\t/;
    chomp @fields;
    $requested_go_id=$fields[0];
    chomp $requested_go_id;
    if (defined($data_lookup{$requested_go_id})){
	print "$requested_go_id\t";
	print "$data_lookup{$requested_go_id}\n";

    }
    else {
	print "$requested_go_id\tERROR\n";

    }
    
}
close GOIDLIST;
exit 0;
