The version of the assembly modified by the code described below does not exist in the public domain.  This code is provided simply as an effort to provide full documentation of the software used to produce EquCab3.  This, as it is will not function in any environment other than the one for which it was written without substantial expert modification.


EquCab3, a reference genome for the domestic horse, Equus caballus was built using DNA from an inbred Thoroughbred named Twilight.  Once the processes of contigging and scaffolding were completed, our assembly was "Content Complete".  We the aligned Illumina HiSeq and miSeq Whole Genome Shotgun sequence data to this content complete version.  As anticipated, there were several thousand errant positions in the assembly where the Illumina data suggested Twilight differed homozygously from her own reference.  The two pieces of software written to identify and repair these errors the assembly are described below

These two pieces of software were written to 
1) Identify homozygous variants in Twilight in her Illumina WGS alignment vs the content complete assembly.  
   	    Four other Thoroughbreds were aligned to the content complete assembly as well.  The UnifiedGenotyper was used to discover homozygous variants in Twilight.  These variants positions were then 
	    genotyped vs the four other Thoroughbreds.  In some of the variant positions, the reference allele was identified.  This suggested that the reference allele existed in this population, and may 
	    have been undersampled in the Twilight Illumina data.  These alleles were left alone.  Those variant positions where not second allele was identified in the other animals were selected for repair.
2) Modify the positions identified in 1) with the variant allele found in the Twilight dataset.
   	    This software moved through each position in the content complete assembly, and asked the question, was this position identified as an error.  If no, write the nucleotide, then continue.  
	    If yes, then write the variant allele identified in the Twilight mapping data and continue to the next position in the assembly.


