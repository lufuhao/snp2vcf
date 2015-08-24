#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 --input my.fa [Options]
Version: LUFUHAO20141003

Requirements
	Perl Modules: Getopt::Long


Descriptions:
	Extract SNPs subsets from 1001 Arabidopsis genomes
	Annotate each SNP based on an given annotation file.




Options:
	--help/-h
		Print this help/usage;
	--input|-i <FILE>
		[Msg] Input fasta file to extract from;
	--region|-r 
		[Opt] Special region like Chr1:1-10000;
		no space allowed;
	--output|-o
		[Opt] File name for output;
		Will go STDOUT if not specified;
	--verbose
		[Optional] Detailed output for trouble-shooting;
	--version/-v
		[Optional] Print current SCRIPT version;
		

Example:
	perl $0 
		

Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk
EOH
###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
our ($help, $input, $region, $output, $verbose, $ver);
GetOptions(
	"help|h!" => \$help,
	"input|i:s" => \$input,
	"region|r:s" => \$region,
	"output|o:s"=> \$output,
	"verbose!" => \$verbose,
	"version|v!" => \$ver) or die USAGE;
($help or $ver) and die USAGE;



###Default and initalization#########################################
die "Error: Multiple-sample VCF not found\n" unless (defined $input and -s $input);
unlink($output) if (defined $output);
our ($chr, $chmin, $chmax);
if (defined $region) {
	if ($region=~/(\S+):(\d{1,12})-(\d{1,12})/) {
		$chr=$1;
		$chmin=$2;
		$chmax=$3;
	}
	else {
		die "Error: Unknown parameter for --region: $region\n";
	}
	print "$region\n";
	print "Selected region: $chr: $chmin - $chmax\n";
}


our @accessions=();
our %VariantID=();
#GATK
our %gen=('-' => '-',     AA => 'A',      TT => 'T',      CC => 'C',      GG => 'G',      'AG' => 'R',    'GA' => 'R',    'CT' => 'Y',    'TC' => 'Y',    'AC' => 'M',    'CA' => 'M',    'GT' => 'K',    'TG' => 'K',    'GC' => 'S',    'CG' => 'S',    'AT' => 'W',    'TA' => 'W',    'ATC' => 'H',   'ACT' => 'H',   'TAC' => 'H',   'TCA' => 'H',   'CAT' => 'H',   'CTA' => 'H',   'GTC' => 'B',   'GCT' => 'B',   CGT => 'B',     CTG => 'B',     TCG => 'B',     TGC => 'B',     GAC => 'V',     GCA => 'V',     CGA => 'V',     CAG => 'V',     AGC => 'V',     ACG => 'V',     GAT => 'D',     GTA => 'D',     AGT => 'D',     ATG => 'D',     TGA => 'D',     TAG => 'D',     ATCG => 'N',    ATGC => 'N',    ACTG => 'N',    ACGT => 'N',    AGTC => 'N',    AGCT => 'N',    TACG => 'N',    TAGC => 'N',    TCAG => 'N',    TCGA => 'N',    TGAC => 'N',    TGCA => 'N',    CATG => 'N',    CAGT => 'N',    CTAG => 'N',    CTGA => 'N',    CGAT => 'N',    CGTA => 'N',    GATC => 'N',    GACT => 'N',    GTAC => 'N',    GTCA => 'N',    GCAT => 'N',    GCTA => 'N');#this hash was use for converting the genotype to a single letter according to th IUPAC.
our %gen2=('-' => '-', 'A' => 'A', 'T' => 'T', 'C' => 'C', 'G' => 'G', 'R' => 'A/G', 'Y' => 'C/T', 'M' => 'A/C', 'K' => 'G/T', 'S' => 'C/G', 'W' => 'A/T', 'H' => 'A/T/C', 'B' => 'G/C/T', 'D' => 'A/G/T', 'V' => 'CAG', 'N' => 'ATCG');
our %gen3=();

###Convert mpileup into vcf##########################################
open (INPUT, "$input") || die "Error: Can not open $input\n";
open (GENOTYPE, ">>$output") || die "Error: can not write to genotypes\n" if (defined $output);
my $temp01=0;
while (my $line=<INPUT>){
	next if ($line=~/^##|^\@|^\*/);
	chomp $line;
#0	1	2	3	4	5	6	77	8	99	1010
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	data	ale_stenar_56_14.997.0cf.v4_data
	if ($line=~/^#CHROM/) {
		@accessions=split(/\t/, $line);
		splice @accessions,0,9;
		next;
	}
	my @arr=();
	@arr=split(/\t/, $line);
	next if (scalar(@arr)<=10);
	%gen3=();
	%gen3=('0'=>"$arr[3]", '1' => "$arr[4]");
	next if (defined $chr and ($arr[0] ne $chr));
	next if (defined $chmin and ($arr[1]<$chmin));
	last if (defined $chmax and ($arr[1]>$chmax));
	print "$arr[0]\t$arr[1]\n" if (defined $verbose);
	my $newline='';
	$newline=$arr[0]."\t".$arr[1]."\t".$arr[3]."\t".$arr[4]."\t".&RetrieveSnpEff($arr[7]);
	for (my $i=9; $i<scalar(@arr); $i++) {
		$newline.="\t".&RetrieveGenotypes($arr[$i]);
	}
	if (defined $output) {
		print GENOTYPE $newline."\n";
	}
	else {
		print $newline."\n";
	}
}

close INPUT;
close GENOTYPE  if (defined $output);


###retrieve SNP eff from SnpEFF output
#RetrieveSnpEff(EFF=)
sub RetrieveSnpEff {
	my $RSinfo=shift @_;
	my @RSarr1=split(/;/,$RSinfo);
	my $RSeff='';
	foreach (@RSarr1) {
		$RSeff=$_;
		last if ($RSeff=~s/^EFF=//);
	}
#	print $RSeff."\n";###########################################For test
#	my @RStotal_eff=("missense_variant", "synonymous_variant", "initiator_codon_variant", "stop_gained", "stop_lost", "stop_retained_variant", "intron_variant", "non_coding_exon_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "5_prime_UTR_premature_start_codon_gain_variant", "splice_region_variant", "splice_donor_variant", "splice_acceptor_variant", "nc_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", "intergenic_region");
	my @RSarr2=split(/,/, $RSeff);
	my $RSnew_eff='';
	foreach my $RSanno (@RSarr2) {
		(my $variantname=$RSanno)=~s/\(.*\)//;
#		$VariantID{$variantname}++;
		unless ($variantname=~/stream|intergenic/) {
			$RSnew_eff.=$RSanno.';';
		}
	}
	return ($RSnew_eff ne '') ? $RSnew_eff : 0;
}
#missense_variant(MODERATE|MISSENSE|Cta/Ata|p.Leu247Ile/c.739C>A|317|Exon_1_4198870_4198999|protein_coding|CODING|AT1G12350.1-Protein|9|1)
#synonymous_variant(LOW|SILENT|aaA/aaG|p.Lys249Lys/c.747A>G|317|Exon_1_4198870_4198999|protein_coding|CODING|AT1G12350.1-Protein|9|1)
#initiator_codon_variant(LOW|MISSENSE|Atg/Ttg|p.Met1?/c.1A>T|693|Exon_1_19245503_19246010|protein_coding|CODING|AT1G51830.1-Protein|1|1)
#stop_gained(HIGH|NONSENSE|Caa/Taa|p.Gln75*/c.223C>T|77|Exon_1_19292930_19293035|protein_coding|CODING|AT1G51913.1-Protein|2|1)
#stop_lost(HIGH|MISSENSE|tAg/tGg|p.Ter125Trpext*?/c.374A>G|124|Exon_1_19425821_19425928|protein_coding|CODING|AT1G52180.1-Protein|2|1)
#stop_retained_variant(LOW|SILENT|taG/taA|p.Ter315Ter/c.945G>A|314|Exon_1_21316851_21316983|protein_coding|CODING|AT1G57560.1-Protein|3|1)
#intron_variant
#non_coding_exon_variant
#3_prime_UTR_variant
#5_prime_UTR_variant
#5_prime_UTR_premature_start_codon_gain_variant
#splice_region_variant", "splice_donor_variant
#splice_acceptor_variant
#nc_transcript_variant
#upstream_gene_variant
#downstream_gene_variant
#intergenic_region



###retrieve genotype
#RetrieveSnpEff(0/1:*:*)
sub RetrieveGenotypes {
	my $RGdigit=shift @_;
	my $RGreturn='';
	return $gen2{$gen3{0}} if ($RGdigit eq '.');
	my @RGtemp01=split(/:/, $RGdigit);
	my @RGtemp02=split(/\//, $RGtemp01[0]);
	if (scalar(@RGtemp02)==2){
		$gen2{$gen3{$RGtemp02[0]}}='?' unless (exists $gen2{$gen3{$RGtemp02[0]}});
		$gen2{$gen3{$RGtemp02[1]}}='?' unless (exists $gen2{$gen3{$RGtemp02[1]}});
		my @RGtemp03=split(/\//, $gen2{$gen3{$RGtemp02[1]}});
		foreach (@RGtemp03) {
			if ($_ eq $gen2{$gen3{$RGtemp02[0]}}) {
				$RGreturn=$gen2{$gen3{$RGtemp02[1]}};
				last;
			}
		}
		if ($RGreturn eq '') {
			$RGreturn=$gen2{$gen3{$RGtemp02[0]}}.'/'.$gen2{$gen3{$RGtemp02[1]}};
		}
	}
	else {
		$RGreturn='unknown';
	}
	return $RGreturn;
}
