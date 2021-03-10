#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;


##Script to parse the output and produce a binary matrix for upset plots inputs

open (FILE,$ARGV[0]);

my %hash_total;
print "x\tmRNA_fragments\ttRNA\tmiRNA\tlncRNA_fragments\tpiRNA\tmisc_RNA\tsnoRNA\tothers\n";
my $total_features=0;
while (<FILE>){
	chomp $_;
	if ($_ =~ /^chr/){
	my @splitt = split (/\t/,$_);
		#print Dumper (@splitt);
		#exit;
		
##	if ($splitt[18] =~ /\"(.+)\"/){
	if ($splitt[15] =~ /(.+)/){
		my @splitt2 = split (/\;/,$1);
		my @list=(0,0,0,0,0,0,0,0);
		foreach my $types (@splitt2){
			$total_features++;
			if (exists $hash_total{$types}){
			$hash_total{$types}++;
			}else{
			$hash_total{$types}=1;
			}
			if ($types =~ "tRNA"){
				$list[1] = 1;
			}
			elsif ($types eq "protein_coding"){
				$list[0] =1;
			}
			elsif ($types eq "miRNA"){
				$list[2] =1;
			}
			elsif ($types eq "lncRNA"){
				$list[3] =1;
			}			
			elsif ($types eq "piRNA"){
				$list[4]=1;
			}
			elsif ($types eq "misc_RNA"){
				$list[5] =1;
			}
			elsif ($types eq "snoRNA"){
				$list[6] =1;
			}
			else {
			$list[7] =1;
					
		}
	
	}
	my $flag =0;
	print "$splitt[0]\t";
	foreach my $items (@list){
		if ($flag ==0 ){
			print $items;
		}else{
			print "\t$items";
		}
		$flag++;
	}
	print "\n";
	
	}




}
}
#print Dumper (%hash_total);
#print "$total_features\n";
foreach my $name (keys %hash_total){
#print "$name ",($hash_total{$name}/$total_features)*100,"%\n";
}
