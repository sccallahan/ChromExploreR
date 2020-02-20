################################################################################################
# @Author: Aayush Raman
# Rai Lab, MDACC 
# Date: Dec. 28th, '15
# 
# Program is used for:
# 1. Chrom States Matrix for 2 MB/100 Kb/10 Kb region
# 2. Used for generating the matrix for the Chrom State heatmap  
# 3. Ver 2. for chromStatesMatrixforClustering
#
# Comments:
# 1. Used for the Prostate Dataset
################################################################################################

#!/tools/bin/perl

use strict;
use warnings;

## Input from user
my $LearnStateNumber = $ARGV[0];
my $STATE_NUMBER = $ARGV[1];
my $bin = $ARGV[2];
my $sampleNameFile = $ARGV[3];
my $ChromHMM_folder = $ARGV[4];

print "The arguments should be given in this order: Number_of_States_used_in_calling_ChromHMM State_Number_user_is_interested_in Bin_Size Sample_Name_File ChromHMM_Folder \n"; 

## Variables
my @chr = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
		   "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", 
		   "chr21", "chr22", "chrX");
my @cellLines;
my %chromSizeHash;
my %chromGC; ## chromosome genomic coordinates 
my %chromStateMatrix;
my $count_gc = 0;
#my $bin = 10e3; # 2e6 = 2 mb, 100e3 = 100 kb
my $start = time;


## Reading sampleName Files
open my $FILE, $ChromHMM_folder."/".$sampleNameFile or die "Cannot open file with sample names $! \n";
while(my $line = <$FILE>){
	chomp $line;
	push @cellLines, $line;
	print $line,"\n";	
}
close($FILE);

## Working directory, Output File for Chrom State Matrix and Chromosome Length File 
open my $OUT, ">", $ChromHMM_folder."/CombinedMatrix-".$STATE_NUMBER."-".$bin."bps.txt" or die "Output Folder does not exists $! \n";

## Chromosome Length File
open my $IN, "ChromSize_hg19.txt" or die "Chromosome Size File does not exists $! \n";
while(my $line = <$IN>){
	chomp $line;
	my @cols = split("\t", $line);
	
	## Defining Genomic Coordinates of region
	$chromSizeHash{$cols[0]} = $cols[1];
	for(my $i = 0; $i <= $cols[1]; $i+=$bin){
		my $start = $i;
		my $end = $i+$bin-1;
		
	# if the $end is bigger than the size of the chromosome
		$end = $cols[1] if($end >= $cols[1]);
		my $gc = $cols[0].":".$start."-".$end;
		$chromGC{$cols[0]}{$gc} = 1;
		$count_gc++;
	}
}
close($IN);
print "Number of genomic bins in hg19 for binsize $bin basepairs are = ",$count_gc,"\n";

foreach my $cell (@cellLines){

	## Opening Dense Bed files
	print "Reading ".$cell."_15_dense.bed file \n";
	my $bedFile = $ChromHMM_folder."/".$cell."_".$LearnStateNumber."_dense.bed";
	print $bedFile,"\n";
	open my $denseBED, $bedFile or die "Cannot open the file $!" ;
	my $header = <$denseBED>;
	
	## Reading the Desnse Bed file
	while(my $line = <$denseBED>){
		chomp $line;
		my @cols = split("\t", $line);
		my $chrBed = $cols[0];
		my $startBed = $cols[1];
		my $endBed = $cols[2];
		my $chromStateBed = $cols[3];
		my $chromSize = $chromSizeHash{$chrBed};
		my $lengthBed = $startBed - $endBed;

		## Defining the cols of the bed file	
		if(exists $chromGC{$chrBed} && $chromStateBed == $STATE_NUMBER){
			
			## genomic cord as possible bin
			my $genomic_StartBin = int($startBed/$bin);
			my $genomic_EndBin = int($endBed/$bin);

			## check if the $binNumber_start and $binNumber_end are in the same bin or not
			if($genomic_StartBin == $genomic_EndBin){
				my $start =  $genomic_StartBin * $bin;
				my $end = $start + $bin - 1;
				$end = $chromSize if($end > $chromSize);
				my $posBin = $chrBed.":".$start."-".$end;
				$chromStateMatrix{$chrBed}{$posBin}{$chromStateBed}{$cell} += 1;
				print $posBin,"\n";	
			}
			elsif($genomic_EndBin > $genomic_StartBin){
				for(my $i = $genomic_StartBin; $i<= $genomic_EndBin; $i++){
					my $start =  $i * $bin;
					my $end = $start + $bin - 1;
					$end = $chromSize if($end > $chromSize);
					my $posBin = $chrBed.":".$start."-".$end;
					$chromStateMatrix{$chrBed}{$posBin}{$chromStateBed}{$cell} += 1;
					print $posBin,"\n";	
				}
			}
		}
	}
	close($denseBED);
	print "Calculation of the number of states in each bin of $bin for ".$cell."_dense.bed is done \n";
	my $duration = time - $start;
	print "Execution time to complete the run for $cell: $duration s\n"; sleep(2);
}

## Forming the Chrom State Matrix
foreach my $chr (@chr){
	foreach my $posBin (sort keys %{$chromStateMatrix{$chr}}){
		print $OUT $posBin;
		foreach my $cell (@cellLines){
			print "\"",$cell."-".$STATE_NUMBER,"\",";
			print $OUT "\t",$chromStateMatrix{$chr}{$posBin}{$STATE_NUMBER}{$cell} if(exists $chromStateMatrix{$chr}{$posBin}{$STATE_NUMBER}{$cell});
			print $OUT "\t0" if(!exists $chromStateMatrix{$chr}{$posBin}{$STATE_NUMBER}{$cell});
		}
		print "\n";
		print $OUT "\n";
	}
}
close($OUT);

## Total Time taken to run the script
my $duration = time - $start;
print "Execution time to complete the entire run: $duration s\n";
