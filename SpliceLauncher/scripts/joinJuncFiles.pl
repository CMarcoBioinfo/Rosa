#!/usr/bin/perl

use strict;
use warnings;

# my @files = @ARGV;

#Code CM
my (@files, $genes);
for (my $i = 0; $i < @ARGV; $i++) {
	if ($ARGV[$i] eq '-c') {
		while ($ARGV[$i +1] !~ /^-/ && $i + 1 < @ARGV) {
			push @files, $ARGV[++$i];		
		}
	} elsif ($ARGV[$i] eq '-g') {
		$genes = $ARGV[$i + 1];
	}
}

#Verifie si la liste de genes est fournie
my $use_genes = defined $genes;

my %genes_list;
 if ($use_genes) {
	open my $gene_file, '<', $genes or die "Impossible d'ouvrir la liste de gènes: $!"; 
	 
	while (my $line = <$gene_file>) {
		chomp $line;
		$genes_list{lc $line} = 1;
		}
		close $gene_file; 
 }

my %allJuncs;
for my $f (@files){
    open(IN,$f);
    while(<IN>){
		my ($chr,$start,$end,$strand,$gene,$count) = split(/\t/);

		if (!$use_genes || grep { $gene =~ /\Q$_\E/i } keys %genes_list) {
			$allJuncs{$start."_".$end."_".$strand."_".$gene}{data} = 
			{
	    		chr => $chr,
	    		start=>$start,
	    		end=> $end,
	    		strand=>$strand,
	    		gene=>$gene
			};
			$allJuncs{$start."_".$end."_".$strand."_".$gene}{count}{$f}=$count;
    	}
	}
    close(IN);
}

##fin code CM

print "chr\tstart\tend\tstrand\tgene";
for my $f (@files){
    print "\t".$f;
}
print "\n";


for my $key (keys(%allJuncs)){
    my $data = $allJuncs{$key}{data};
    print   $data->{chr}."\t".$data->{start}."\t".$data->{end}."\t".$data->{strand}."\t".$data->{gene};
    
    for my $f (@files){
	if(defined $allJuncs{$key}{count}{$f}){
	    my $junc_count = $allJuncs{$key}{count}{$f};
		$junc_count =~s/\n//;
	    print  "\t".$junc_count;
	}else{
	    print  "\t0";
	}
    }
    print "\n";
}
