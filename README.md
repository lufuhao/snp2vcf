# snp2vcf



## Version: 20141003



## Descriptions:

	Comvert mpileup SNP into vcf, annotate SNPs based on snpEff
	database, and transform VCF into genotypes file.



## Requirements:

+  VCFtools

+  Tabix, bgzip

+  Perl Molules: Getopt::Long

+  corresponding snpEff database

	Configure Program root/utils/snpEFF.config
	change data.dir = /usr/users/celldev/luf/database/DB/snpeff to yours
	You may download the database using:

	#For TAIR9

	$java -jar /path/to/snpEff.jar download athalianaTair9
	
	or
	
	$java -jar /path/to/snpEff.jar download athalianaTair10


## Running:

- 1. convert mpileup into multiple-sample vcf

	$snp2vcf -i /path/to/snp/folder -o All.samples.vcf -x athalianaTair9

- 2.Transform vcf into genotypes

	$perl /path/to/vcf2geno.pl -i All.samples.vcf --region Chr1:1-20000 -o All.genotypes.txt



## Author:
	Fu-Hao Lu
	Post-Doctoral Scientist in Micheal Bevan laboratory
	Cell and Developmental Department, John Innes Centre
	Norwich NR4 7UH, United Kingdom
	E-mail: Fu-Hao.Lu\@jic.ac.uk

## Copyright

	Copyright (c) 2014-2018 Fu-Hao Lu
