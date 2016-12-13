###This script is part of perl pipeline developed by Michelle Graham, USDA-ARS (michelle.graham@ars.usda.gov).  
###Please see GO_analysis_README for usage.  
###If you use this pipeline (or modified versions) please cite: 
###van de Mortel M, Recknor JC, Graham MA, Nettleton D, Dittman JD, Nelson RT, Godoy CV, Abdelnoor RV, Almeida AM, 
### Baum TJ, Whitham SA (2007) Distinct biphasic mRNA changes in response to Asian soybean rust infection. Mol 
###Plant Microbe Interact. 20: 887-899  

#!/usr/bin/perl
($infile1,$infile2)=@ARGV;
open (DATAHASH,$infile2) || die "Cannot open $infile2 \n";
while (<DATAHASH>) {
    @tmp=split /\t/;
    $gmid=$tmp[0];
    chomp $affyid;
    $GOMF= $tmp[3];
    chomp $GOMF;
    @IDs=split ('\s',$GOMF);
#    print "TESTING....\n";
#    print @IDs;
#    print "\n";
    undef %saw; 
    @uniquedata = grep(!$saw{$_}++, @IDs);
#    print @uniquedata;
#    print "\n";
    $data= join(' ',@uniquedata);
#    print "$data\n";
    foreach $affyid (@tmp) {
	if ($affyid ne "") {
	    $data_lookup{$affyid} = $data;
	}
	
    }
#    $data='';
    $affyid='';
    $GOBP='';
    @IDs='';
    @uniquedata='';
    $data='';
}
close DATAHASH;
open (STARTERFILE,$infile1) || die "Cannot open $infile1 \n";
while (<STARTERFILE>) {
    @fields= split /\t/;
    chomp @fields;
    $requested_gm_id=$fields[0];
    chomp $requested_gm_id;
    if ($requested_gm_id=~ /(Glyma\d{1,}g\d{1,})\.\S+/){
	$new_requested_gm_id=$1;
    }
    else {
	$new_requested_gm_id=$requested_gm_id;
    }
    if ($myhash{$new_requested_gm_id}) {
	# already analyzed do nothing
    }
    else {
	if (defined($data_lookup{$new_requested_gm_id})){
	    #print "$new_requested_gm_id\t";
#	    #print $data_lookup{$new_requested_gm_id};
	    #print "\n";
	    $results=$data_lookup{$new_requested_gm_id};
	    $results=~ s/\s/\n/g;
	    print $results;
	    print "\n";
	}
	else {
	    print "ERROR\n";

	}
	
	$myhash{$new_requested_gm_id}=1;
    }
    
}


close STARTERFILE;


close STARTERFILE;
exit 0;



