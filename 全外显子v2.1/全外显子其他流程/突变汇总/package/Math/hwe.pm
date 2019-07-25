package package::Math::hwe;
use strict;
use warnings;

###
# º∆À„HWE
###
sub run{
	my ($obs_hets,$obs_hom1,$obs_hom2) = @_;
	if($obs_hom1 < 0 || $obs_hom2 < 0 || $obs_hets <0){
		return(-1);
	}
	# rare homozygotes
	my $obs_homr;
	# common homozygotes
	my $obs_homc;
	if($obs_hom1 < $obs_hom2){
		$obs_homr = $obs_hom1;
		$obs_homc = $obs_hom2;
	}
	else{
		$obs_homr = $obs_hom2;
		$obs_homc = $obs_hom1;
	}
	# number of rare allele copies
	my $rare_copies = 2 * $obs_homr + $obs_hets;
	# total number of genotypes
	my $genotypes = $obs_homr + $obs_homc + $obs_hets;
	if($genotypes <= 0){
		return(-1);
	}
	# Initialize probability array
	my @het_probs;
	for(my $i=0; $i<=$rare_copies; $i++){
		$het_probs[$i] = 0.0;
	}
	# start at midpoint
	my $mid = int($rare_copies * (2 * $genotypes - $rare_copies) / (2 * $genotypes));
	# check to ensure that midpoint and rare alleles have same parity
	if(($rare_copies & 1) ^ ($mid & 1)){
		$mid++;
	}
	my $curr_hets = $mid;
	my $curr_homr = ($rare_copies - $mid) / 2;
	my $curr_homc = $genotypes - $curr_hets - $curr_homr;
	$het_probs[$mid] = 1.0;
	my $sum = $het_probs[$mid];
	for($curr_hets = $mid; $curr_hets > 1; $curr_hets -= 2){
		$het_probs[$curr_hets - 2] = $het_probs[$curr_hets] * $curr_hets * ($curr_hets - 1.0) / (4.0 * ($curr_homr + 1.0) * ($curr_homc + 1.0));
		$sum += $het_probs[$curr_hets - 2];
		# 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
		$curr_homr++;
		$curr_homc++;
	}
	$curr_hets = $mid;
	$curr_homr = ($rare_copies - $mid) / 2;
	$curr_homc = $genotypes - $curr_hets - $curr_homr;
	for($curr_hets = $mid; $curr_hets <= $rare_copies - 2; $curr_hets += 2){
		$het_probs[$curr_hets + 2] = $het_probs[$curr_hets] * 4.0 * $curr_homr * $curr_homc / (($curr_hets + 2.0) * ($curr_hets + 1.0));
		$sum += $het_probs[$curr_hets + 2];
		# add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
		$curr_homr--;
		$curr_homc--;
	}
	for(my $i=0; $i<=$rare_copies; $i++){
		$het_probs[$i] /= $sum;
	}
	# Initialise P-value 
	my $p_hwe = 0.0;
	# P-value calculation for p_hwe
	for(my $i = 0; $i <= $rare_copies; $i++){
		if($het_probs[$i] > $het_probs[$obs_hets]){
			next;
		}
		$p_hwe += $het_probs[$i];
	}
	if($p_hwe > 1) {
		$p_hwe = 1.0;
	}
	return($p_hwe);
}



1