# CollaborativeCross_MarkFounder
Compute genotypes (any snp present in dbSNP142) for any CC strain based on founder probability tables (http://csbio.unc.edu/CCstatus/index.py?run=FounderProbs) 

"GenerateFounderStrTable.R" converts the numerical probability columns in founder probabiliy tables to a single column. The columns informs any founder origins for this position. Threshold set to 0 so that any non-zero founder can make a contribution.

"Annot_snp142.R" annonates genotypes for a particular CC strain any any snp as long as it's present in dbsnp142.
