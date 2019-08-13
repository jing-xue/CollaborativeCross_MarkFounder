# CollaborativeCross_MarkFounder
Compute genotypes (any snp present in dbSNP142) for any CC strain based on founder probability tables (http://csbio.unc.edu/CCstatus/index.py?run=FounderProbs) 

"GenerateFounderStrTable.R" converts the numerical probability columns in founder probabiliy tables to a single column. The columns informs any founder origins for this position. Threshold set to 0 so that any non-zero founder can make a contribution.

"Annot_snp142.R" annonates genotypes for a particular CC strain at any snp as long as it's present in dbsnp142.


# dbSNP142 vcf file cleanup:
## dbSNP mouse build142 vcf files (http://www.sanger.ac.uk/science/data/mouse-genomes-project)
### Download
$ wget 'ftp://ftp-mouse.sanger.ac.uk/current_snps/*'

##directory: /proj/folamilb/projects/Jan16MethylSeq/dbSNP142/snp142

$ gunzip -c mgp.v5.merged.snps_all.dbSNP142.vcf.gz > mgp.v5.merged.snps_all.dbSNP142.vcf
$ gunzip -c mgp.v5.merged.indels.dbSNP142.normed.vcf.gz > mgp.v5.merged.indels.
dbSNP142.normed.vcf

### keep only CC founders
$ cat mgp.v5.merged.snps_all.dbSNP142.vcf | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5}' > dbSNP142.bed

### remove header
$ egrep -v "^#" mgp.v5.merged.snps_all.dbSNP142.vcf > no_header_snp142.vcf
$ cat no_header_snp142.vcf | awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5}' > dbSNP142_noHeader.bed

### final re-format
$ cat no_header_snp142.vcf | awk '{print "chr"$1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5"\t"$14"\t"$4"\t"$11"\
t"$35"\t"$37"\t"$25"\t"$39"\t"$44}' > dbSNP142_strain_noHeader_final.bed
