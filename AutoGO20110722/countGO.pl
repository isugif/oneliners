###This script is part of perl pipeline developed by Michelle Graham, USDA-ARS (michelle.graham@ars.usda.gov).  
###Please see GO_analysis_README for usage.  
###If you use this pipeline (or modified versions) please cite: 
###van de Mortel M, Recknor JC, Graham MA, Nettleton D, Dittman JD, Nelson RT, Godoy CV, Abdelnoor RV, Almeida AM, 
### Baum TJ, Whitham SA (2007) Distinct biphasic mRNA changes in response to Asian soybean rust infection. Mol 
###Plant Microbe Interact. 20: 887-899  

#!/usr/bin/perl                                                                                                                                                        
($file1)=@ARGV;
open (FILE, $file1) || die "Cannot open file $file1 \n";
my %myhash=();
my @fields=();
my $pfam='';
while (<FILE>) {
    @fields= split /\t/;
    $GO_id=$fields[0];
    chomp $GO_id;
    $GO_id=~ s/\s//g;
    if ($GO_id=~ /GO/) {
	$count++;
        if ($myhash{$GO_id}) {
            $myhash{$GO_id}++;
        }
        else {
            $myhash{$GO_id}=1;
        }
    }
}


close (FILE);
while ( my ($key, $value) = each(%myhash) ) {
    chomp $key;
    print "$key\t$value\n";
}

close (DATA);
exit 0;
