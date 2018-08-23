#!/usr/bin/perl -w

### Resample masks from NeuroRx processing to subject space
#
##########################################################
##########################################################
## The lesion masks from NeuroRx are in stx space, want ot bring them back to native diffusion space
## convention: 
## X = 1st inital x=2nd initial d=patient number m00=month time point
##  - X-x_00d/m00/*_ct2f.mnc.gz
##				-t2 lesions (hypo on T1p, hyper on FLAIR)
##				-in stx space
##  - X-x_00d/m00/*_gvf.mnc.gz 
##				- Gad enhancing lesions in stx space
##  - X-x_00d/m00/RegisterToStandardSpace/*/*_t1p-to-stx152lsq6.xfm
##				- transformation from t1p (pre-Gd) to stx space
##  - X-x_00d/m00/TissueClassification/*/*ANAT-na-cerebrum-wm_ISPC-stx152lsq6.mask.mnc.gz
##				- WM in stx space
##  - X-x_00d/m00/TissueClassification/*/*segmentation_ISPC-stx152lsq6.mnc.gz
##				- classification in stx space
##  - X-x_00d/m00/ResampleToStx/*/*t1gMPR-t1map_ISPC-stx152lsq6.mnc.gz
##				- MP2Rage in stx space
## - dti-to-t1b.xfm
##				- created by g-ratio_pipeline
##				- transformation from diffusion to t1p (pre-Gd) space
###########################################
#
###########################################
# Created by Ilana Leppert March 2018
###########################################


require 5.001;
use Getopt::Tabular;
use MNI::Startup qw(nocputimes);
use MNI::Spawn;
use MNI::FileUtilities qw(check_output_dirs);
use File::Basename;
use List::MoreUtils qw(zip);
use Cwd qw(cwd);

if($0 =~ /[\/A-Za-z_0-9-]+\/([A-Za-z0-9_-]+\.pl)/) {$program = $1;}	#program name without path
$Usage = <<USAGE;

Resamples selected masks to native space


Usage: $program
-help for options

USAGE
##defaults

$pwd=cwd();
if ($pwd=~/7_([a-z])_([a-z])_(\d{3})_(m\d{2})/){$sub1=$1;$sub2=$2;$num=$3;$visit=$4;} else {die "can't find subject info for $pwd\n";}

print "\n---On loki---\n";
$dir=$sub1."-".$sub2."_".$num."/".$visit;
$tar=$sub1."-".$sub2."_".$num."_".$visit.".tar";
print "tar -cvf /tmp/$tar $dir/*_ct2f.mnc* $dir/*_gvf.mnc.gz $dir/TissueClassification/*/*ANAT-na-cerebrum-wm_ISPC-stx152lsq6.mask.mnc.gz $dir/TissueClassification/*/*segmentation_ISPC-stx152lsq6.mnc.gz $dir/RegisterToStandardSpace/*/*_t1p-to-stx152lsq6.xfm $dir/ResampleToStx/*/*t1gMPR-t1map_ISPC-stx152lsq6.mnc.gz\n";
print "scp ileppert\@loki.bic.mni.mcgill.ca:/tmp/$tar .\n";
unless (-e $tar){
	die "---> first download the data from loki\n";
}

`tar xvf $tar`;
$t2les = `\\ls $dir/*_ct2f.mnc*`; chomp($t2les);
$gd = `\\ls $dir/*_gvf.mnc.gz`; chomp($gd);
$wm = `\\ls $dir/TissueClassification/*/*ANAT-na-cerebrum-wm_ISPC-stx152lsq6.mask.mnc.gz`; chomp($wm);
$seg = `\\ls $dir/TissueClassification/*/*segmentation_ISPC-stx152lsq6.mnc.gz`; chomp($seg);
$t1p2stx = `\\ls $dir/RegisterToStandardSpace/*/*_t1p-to-stx152lsq6.xfm`;chomp($t1p2stx);
$t1map = `\\ls $dir/ResampleToStx/*/*t1gMPR-t1map_ISPC-stx152lsq6.mnc.gz`; chomp($t1map);


@t1ps=`more mnclist | grep mni_gre_3D_T1w-PreGd`;
@file=split(/\s+/,$t1ps[1]);
$t1p = $file[0];#2nd should be the normalized one used by NeuroRx for registation
$diff = `ls b0_eddy_corr_mnc.mnc`; chomp($diff);
$diff2t1p  = `\\ls dti-to-t1b.xfm`; chomp($diff2t1p);


#my @args_table = (["-t1p","string",1,\$t1p,"T1w pre-contrast."],
	#["-diff","string",1,\$diff,"b0 to which everything should be registered to"],
#		  ["-t2les","string",1,\$t2les,"T2 lesions mask in stx space"],
		   #["-gd","string",1,\$gd,"Gad lesion mask in stx space"],
#);



#Getopt::Tabular::SetHelp ($Usage, '');
#my @args;
#GetOptions(\@args_table, \@ARGV, \@args) || exit 1;
#die $Usage unless $#ARGV >=0;

print "**** $program @ARGV ****\n";

### Resample masks
#transformations
print "\n------------Compute transformations----------------\n";
#stx-to-b0
$b02stx="b0-to-stx.xfm";
print "xfmconcat $diff2t1p $t1p2stx $b02stx\n";
`xfmconcat $diff2t1p $t1p2stx $b02stx`;

print "-------Depending how they were created, have to resample masks to full stx--------\n";

if (-e $t2les){
($namebase,$dir,$ext)=fileparse($t2les,'\..*');
$t2les_stx=$namebase."-stx152lsq6.mnc";
print "mincresample $t2les -near -like $wm $t2les_stx\n";
`mincresample $t2les -near -like $wm $t2les_stx`;
$t2les = $t2les_stx;
}else{die "----Can't find t2les mask\n";}

if (-e $gd){
($namebase,$dir,$ext)=fileparse($gd,'\..*');
$gdles=$namebase."-stx152lsq6.mnc";
print "mincresample $gd -near -like $wm $gdles\n";
`mincresample $gd -near -like $wm $gdles`;
}else{die "----Can't find Gdles mask\n";}

if (-e $t2les){
print "--------Create a black-hole type mask, with a threshold on T1 (looks dark on lepto)----------\n";
$blakh = "black-holes-thres1000-stx152lsq6.mnc";
	print "minccalc -expr A[0]>1000 && A[1]  $t1map $t2les $blakh\n";
	`minccalc -expr "A[0]>1000 && A[1]" $t1map $t2les $blakh` unless -e $blakh;

}else{die "----Can't find T2les mask\n";}



print "\n----------Make sure masks are mutually exclusive----------------\n";
# Gd lesions are a subset of T2 lesions
# "black holes" are a subset of T2 lesions
$t2les_only="t2les-only.mnc"; # not black holes nor WM
$gdles_only="gdles-only.mnc";
$wm_only="wm-only.mnc";

if (-e $t2les && -e $wm){
	# not black holes nor WM 
	print "--> T2 lesions only <--\n minccalc -expr A[0]!=1 && A[1]!=1 && A[2]==1  $wm $blakh $t2les $t2les_only\n";
	`minccalc -expr "A[0]!=1 && A[1]!=1 && A[2]==1"  $wm $blakh $t2les $t2les_only`;
	# not T2 lesions, only WM
	print "--> WM only <--\n minccalc -expr A[0]!=1  && A[1]==1 $t2les $wm $wm_only\n";
   `minccalc -expr "A[0]!=1 && A[1]==1" $t2les $wm $wm_only` unless -e $wm_only;
}
if (-e $gd && -e $wm){
#Gad lesions are include in T2 lesions?
	print "--> Gad Les only <--\n minccalc -expr A[0]!=1 && A[1]>0   $wm_only $gdles $gdles_only\n";
	`minccalc -expr "A[0]!=1 && A[1]>0 "  $wm_only $gdles $gdles_only`;
}



print "\n---------------Resample masks--------------\n";
$t2les_diff="t2les-like-diff.mnc";
$gdles_diff="gdles-like-diff.mnc";
$bh_diff="blackh-like-diff.mnc";
$wm_diff="wm-like-diff.mnc";
#$wm_diff="WM_MTsatmasked.mnc";
print "mincresample -near -like $diff -transformation $b02stx -invert $t2les_only $t2les_diff\n";
`mincresample -near -like $diff -transformation $b02stx -invert $t2les_only $t2les_diff`  unless -e $t2les_diff ;
if (-e $gd){
	print "mincresample -near -like $diff -transformation $b02stx -invert $gdles_only $gdles_diff\n";
	`mincresample -near -like $diff -transformation $b02stx -invert $gdles_only $gdles_diff` unless -e $gdles_diff;
}
print "mincresample -near -like $diff -transformation $b02stx -invert $blakh $bh_diff\n";
`mincresample -near -like $diff -transformation $b02stx -invert $blakh $bh_diff`  unless -e $bh_diff;
print "mincresample -near -like $diff -transformation $b02stx -invert $wm_only $wm_diff\n";
`mincresample -near -like $diff -transformation $b02stx -invert $wm_only $wm_diff`;

## WM mask should avoid PVE because artificially increases g
#print "minccalc -expr A[0]&&A[1] brainmask_MTsat.mnc $wm_diff wm-final.mnc\n";
#`minccalc -expr "A[0]&&A[1]" brainmask_MTsat.mnc $wm_diff wm-final.mnc`;
#$wm_diff= "wm-final.mnc";
## make wm mask more stringent
# to do this, we need the classfied volume in tal space, transform back to native and use Mathieu's code for majorityvoting
print "mincresample -nearest -tfm_input_sampling -transformation $b02stx -invert -step 1 1 1 -znelements 180 $seg seg_diffspace_highres.mnc\n";
`mincresample -nearest -tfm_input_sampling -transformation $b02stx -invert -step 1 1 1 -znelements 180 $seg seg_diffspace_highres.mnc`;
print "run_matlab(maskMajorityVoting('anat-to-diff_mnc.mnc', 'seg_diffspace_highres.mnc', 9, 'final-nrx')) \n";
run_matlab("maskMajorityVoting('anat-to-diff_mnc.mnc', 'seg_diffspace_highres.mnc', 9, 'final-nrx');") unless -e "final-nrx_wm_perc.mnc";
print "minccalc -expression  A[0] > 99.1? 1 : 0 final-nrx_wm_perc.mnc wm-nrx_100pc_mask.mnc\n";
`minccalc -expression  "A[0] > 99.1? 1 : 0" final-nrx_wm_perc.mnc wm-nrx_100pc_mask.mnc`;
print "mincmath -mult wm-nrx_100pc_mask.mnc brainmask_MTsat.mnc WM-100pc_MTsatmasked.mnc -clob\n";
`mincmath -mult wm-nrx_100pc_mask.mnc brainmask_MTsat.mnc WM-100pc_MTsatmasked.mnc -clob`;
$wm_diff= "WM-100pc_MTsatmasked.mnc";

print "\n----------Some stats----------------\n";
@maps=("g_MTsat_g.mnc","g_MTsat_mvf.mnc","g_MTsat_avf.mnc");
@masks=($t2les_diff,$gdles_diff,$bh_diff,$wm_diff);
for $m (@maps){
	($namebase,$dir,$ext)=fileparse($m,'\..*');
	$lst='';
	for $s (@masks){
		$lst = $lst.' '.$m.' '.$s; #map name and mask name pair
	}
	print "box_plot.py $lst $namebase-in-masks.png\n\n";
	`box_plot.py $lst $namebase-in-masks.png`;
}

sub run_matlab {
  my ($command)=@_;
  #system("(echo \"$command\"| matlab -nosplash $logging)");
  open MATLAB,"|matlab -nosplash -nojvm -nodisplay " or die "Can't start matlab!\n";
  #print MATLAB "addpath('$mydir')\n";
  #print MATLAB "addpath(genpath('$niak_dir'))\n" if $niak_dir;
  print MATLAB $command,"\n";
  close MATLAB;
  die "Error in MATLAB!\n" if $?;
}
