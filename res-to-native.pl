#!/usr/bin/perl -w

### Resample masks from NeuroRx processing to subject space
##
## Multi-time point analysis
# naming convention: type_timepoint-resampled_to
# where:
#   type: t2les, enlarging, denovo, gdles
#   timepoint: m00 or nothing if this mask is from the current timepoint
#   resample_to : space the mask has been resampled to
#e.g.
#$t2lesm00_diff="t2lesm00-like-diff.mnc"; #these are the T2 lesions from m00 only (excluding new and enlarging), resampled to the current timepoint's diffusion space
#$denovo_diff_m00="denovo-like-diff_m00.mnc"; #new lesions at this timepoint, resampled to the m00 timepoint diffusion
#$gdles_diff_m00 = "gdles-like-diff_m00.mnc"; #gd lesions at this timepoint, resampled to the m00 timepoint diffusion
#$gdles_m00_diff = "gdles_m00-like-diff.mnc"; #gd lesions at m00, resampled to the current timepoint diffusion

##########################################################
##########################################################
## The lesion masks from NeuroRx are in stx space, want ot bring them back to native diffusion space
## convention: 
## X = 1st inital x=2nd initial d=patient number mxx=month time point
## m00
##  - X-x_00d/mxx/*_ct2f_ISPC-stx152iso.mnc.gz
##				-t2 lesions (hypo on T1p, hyper on FLAIR)
##				-in stx 1mm3 space
## > m00
##  - X-x_00d/mxx/*_newt2f_TREF-mxx_ISPC-stx152iso.mnc.gz   #new lesions
##
## mxx
##  - X-x_00d/mxx/*_gvf_ISPC-stx152*.mnc.gz 
##				- Gad enhancing lesions in 1mm3 stx space
##  - X-x_00d/mxx/RegisterToStandardSpace/*/*_t1p-to-stx152lsq6.xfm
##				- transformation from t1p (pre-Gd) to stx space
##  - X-x_00d/mxx//malfTissuePriors/*prob_cerebrum-wm_ISPC-stx152iso.mnc.gz
##				- probability of WM in 1mm iso stx space 
##  - ##X-x_00d/mxx/TissueClassification/*/*segmentation_ISPC-stx152lsq6.mnc.gz
##				- classification in stx space (don't need this)
##  - ##X-x_00d/mxx/ResampleToStx/*/*t1gMPR-t1map_ISPC-stx152iso.mnc.gz 
##				- MP2Rage T1map in stx space (don't need this)
##  - X-x_00d/mxx/RegisterOtherModalities/*/*_t1gMPR-t1map-to-stx152lsq6.xfm
##				- transformation of MP2Rage T1map in stx space
## - dti-to-t1b.xfm
##				- created by g-ratio_pipeline
##				- transformation from diffusion to t1p (pre-Gd) space
## Root directory on loki is /lab2/NOVA-LEPTO/001-MNI-7
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

$thresh = 1600; #t1 threshold to separate black holes and t2 lesions

# Pattern for each file
$t2les_pat = "/*_ct2f_ISPC-stx152iso.mnc.gz"; #in 1mmiso
$gd_pat = "/*_gvf_ISPC-stx152*.mnc.gz"; #in 1x1x3mm?? argh not all subjects have it done!!
$wm_pat = "/malfTissuePriors/*prob_cerebrum-wm_ISPC-stx152iso.mnc.gz";
#$seg_pat = 
$t1p2stx_pat = "/RegisterToStandardSpace/*/*_t1p-to-stx152lsq6.xfm";
$t1map2stx_pat = "/RegisterOtherModalities/*/*_t1gMPR-t1map-to-stx152lsq6.xfm";

# If not 1st timepoint, we have other lesion masks and no t2lesion mask
$denovo_pat = "/*_denovot2f_TREF-m00_ISPC-stx152iso.mnc.gz"; #new lesions
$enlarg_pat = "/*_enlargingt2f_TREF-m00_ISPC-stx152iso.mnc.gz"; #enlarging
$res_pat = "/*_rest2f_TREF-m00_ISPC-stx152iso.mnc.gz"; #resolving


if($0 =~ /[\/A-Za-z_0-9-]+\/([A-Za-z0-9_-]+\.pl)/) {$program = $1;}	#program name without path
$Usage = <<USAGE;

Resamples selected masks to native space


Usage: $program
-help for options

USAGE
##defaults
open(ID,">id.txt");
open(VISIT,">visit.txt");

$pwd=cwd();
if ($pwd=~/7_([a-z])_([a-z])_(\d{3})_(m\d{2})/){
  $sub1=$1;$sub2=$2;$num=$3;$visit=$4;
  $dir=$sub1."-".$sub2."_".$num."/".$visit;
  $tar=$sub1."-".$sub2."_".$num."_".$visit.".tar";
  $id=$sub1."-".$sub2."_".$num;
  print ID $id;
  print VISIT $visit;
} 
elsif ($pwd=~/7_([a-z]+)_(\d{3})_(m\d{2})/){
    $sub=$1;$num=$2;$visit=$3; #sometimes the subject names is 3 chars instead of initials separated by _???
	$dir=$sub."_".$num."/".$visit;
    $tar=$sub."_".$num."_".$visit.".tar";
    $id=$sub."_".$num;
    print ID $id;
    print VISIT $visit;
} 
else {die "can't find subject info for $pwd\n";}

print "\n---On loki---\n";

#print "tar -cvf /tmp/$tar $dir/*_ct2f.mnc* $dir/*_gvf.mnc.gz $dir/TissueClassification/*/*ANAT-na-cerebrum-wm_ISPC-stx152lsq6.mask.mnc.gz $dir/TissueClassification/*/*segmentation_ISPC-stx152lsq6.mnc.gz $dir/RegisterToStandardSpace/*/*_t1p-to-stx152lsq6.xfm $dir/ResampleToStx/*/*t1gMPR-t1map_ISPC-stx152lsq6.mnc.gz\n";
if ($visit eq "m00"){
	print "tar -cvf /tmp/$tar $dir$t2les_pat $dir$gd_pat $dir$wm_pat $dir$t1p2stx_pat $dir$t1map2stx_pat\n";
}else{
	print "tar -cvf /tmp/$tar $dir$gd_pat $dir$wm_pat $dir$t1p2stx_pat $dir$t1map2stx_pat $dir$denovo_pat $dir$enlarg_pat $dir$res_pat\n";
}

print "\nscp ileppert\@loki.bic.mni.mcgill.ca:/tmp/$tar .\n\n ";
unless (-e $tar){
	die "---> first download the data from loki\n";
}

`tar xvf $tar`;

chomp($unique=`date +%y%m%d-%H:%M`); 
print "========== Date is $unique =============\n";
$gd = `\\ls $dir$gd_pat`; chomp($gd);
$wm = `\\ls $dir$wm_pat`; chomp($wm);
#$seg = `\\ls $dir/TissueClassification/*/*segmentation_ISPC-stx152lsq6.mnc.gz`; chomp($seg);
$t1p2stx = `\\ls $dir$t1p2stx_pat`;chomp($t1p2stx);
# The T1map has been N4 corrected, which totally disrupts the values, redo the resampling
#$t1map = `\\ls $dir/ResampleToStx/*/*t1gMPR-t1map_ISPC-stx152lsq6.mnc.gz`; chomp($t1map);
$t1map2stx = `\\ls $dir$t1map2stx_pat`;chomp($t1map2stx);
if ($visit eq "m00"){
	$t2les = `\\ls $dir$t2les_pat`; chomp($t2les);
}
if ($visit ne "m00"){
	$dir_m00_pat="../*".$sub1."_".$sub2."_".$num."_m00*/*/m00/"; 
	$dir_m00 = `\\ls -d $dir_m00_pat`; chomp($dir_m00);
	$t2les_m00 = `\\ls $dir_m00$t2les_pat`; chomp($t2les_m00); 
	$denovo = `\\ls $dir$denovo_pat`; chomp($denovo);
	$enlarg = `\\ls $dir$enlarg_pat`; chomp($enlarg);
	$res = `\\ls $dir$res_pat`; chomp($res);
}


@t1ps=`more mnclist | grep mni_gre_3D_T1w-PreGd`;
@file=split(/\s+/,$t1ps[1]);
$t1p = $file[0];#2nd should be the normalized one used by NeuroRx for registation
$diff = `ls b0_eddy_corr_mnc.mnc`; chomp($diff);
$diff2t1p  = `\\ls dti-to-t1b.xfm`; chomp($diff2t1p);
@t1map_nats = `more mnclist | grep MP2RAGE_1mm_T1_Images`;
@file = split(/\s+/,$t1map_nats[0]);
$t1map_nat = $file[0];

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
print "\n________________Compute transformations________________\n";
#stx-to-b0
$b02stx="b0-to-stx.xfm";
print "xfmconcat $diff2t1p $t1p2stx $b02stx\n";
`xfmconcat $diff2t1p $t1p2stx $b02stx`;

print "________________Depending how they were created, have to resample masks to full stx________________\n";

if ($visit eq "m00"){ #visit m00 has t2 lesions
	if (-e $t2les){
	## These now seem to be already in 1mm iso stx space
	#($namebase,$dir,$ext)=fileparse($t2les,'\..*');
	#$t2les_stx=$namebase."-stx152lsq6.mnc";
	#print "mincresample $t2les -near -like $wm $t2les_stx\n";
	#`minresample $t2les -near -like $wm $t2les_stx`;
	#$t2les = $t2les_stx;
	}else{die "----Can't find t2les mask\n";}
}else{ #other visits have only denovo, enlarging, resolving lesions
	print "________________Create t2lesion mask with m00 t2 + denovo + enlarging________________\n";
	if (!-e $denovo){die "----Can't find denovo mask\n";}
	if (!-e $enlarg){die "----Can't find enlarging lesion mask\n";}
	if (!-e $res){die "----Can't find resolving lesion mask\n";}
	$t2les = "ct2f-stx152iso.mnc"; #this is ALL t2 lesions
	print "minccalc -expr A[0]+A[1]+A[2] $t2les_m00 $denovo $enlarg $t2les\n";
	`minccalc -expr "A[0]+A[1]+A[2]" $t2les_m00 $denovo $enlarg $t2les`;
}

if (-e $gd){
($namebase,$dir,$ext)=fileparse($gd,'\..*');
$gdles=$namebase."iso.mnc";
print "mincresample $gd -near -like $wm $gdles\n";
`mincresample $gd -near -like $wm $gdles`;
}else{die "----Can't find Gdles mask\n";}

if (-e $t1map2stx && -e $t1map_nat){
	print "________________Create T1map in stx space________________\n";
	$t1map = "t1map-stx152lsq6iso.mnc";
	print "mincresample $t1map_nat -like $wm -transform $t1map2stx $t1map\n";
	`mincresample $t1map_nat -like $wm -transform $t1map2stx $t1map` unless -e $t1map;
}else{die "----Can't find T1map transformation or native T1map\n";}

if (-e $t2les){
print "________________Create a black-hole type mask, with a threshold on T1 (looks dark on lepto)________________\n";
$blakh = "black-holes-thres$thresh-stx152lsq6iso.mnc";
	print "minccalc -expr A[0]>$thresh && A[1]>=1  $t1map $t2les $blakh\n";
	`minccalc -expr "A[0]>$thresh && A[1]>=1" $t1map $t2les $blakh` unless -e $blakh;

}else{die "----Can't find T2les mask\n";}



print "\n________________Make sure masks are mutually exclusive________________\n";
# Gd lesions are a subset of T2 lesions
# "black holes" are a subset of T2 lesions
# the WM  mask is a probability mask, threshold at 95%?
$thresh_wm = 0.95;
$t2les_only="t2les-only.mnc"; # not black holes nor WM
$gdles_only="gdles-only.mnc";
$wm_only="wm-only.mnc";

if (-e $t2les && -e $wm){
	# not black holes nor WM 
	print "--> T2 lesions only <--\n minccalc -expr A[0]<$thresh_wm && A[1]!=1 && A[2]>=1  $wm $blakh $t2les $t2les_only\n";
	`minccalc -expr "A[0]<$thresh_wm && A[1]!=1 && A[2]>=1"  $wm $blakh $t2les $t2les_only`;
	# not T2 lesions, only WM
	print "--> WM only <--\n minccalc -expr A[0]!=1  && A[1]>=$thresh_wm $t2les $wm $wm_only\n";
   `minccalc -expr "A[0]!=1 && A[1]>=$thresh_wm" $t2les $wm $wm_only` unless -e $wm_only;
}
if (-e $gd && -e $wm){
#Gad lesions are include in T2 lesions?
	print "--> Gad Les only <--\n minccalc -expr A[0]<$thresh_wm && A[1]>0   $wm_only $gdles $gdles_only\n";
	`minccalc -expr "A[0]<$thresh_wm && A[1]>0 "  $wm_only $gdles $gdles_only`;
}



print "\n________________Resample masks to native space________________\n";
$t2les_diff="t2les-like-diff.mnc";
$gdles_diff="gdles-like-diff.mnc";
$bh_diff="blackh-like-diff.mnc";
$wm_diff="wm-like-diff.mnc";
$denovo_diff="denovo-like-diff.mnc";
$enlarging_diff="enlarging-like-diff.mnc";
$t2lesm00_diff="t2lesm00-like-diff.mnc"; #these are the old T2 lesions only (excluding new nad enlarging)
$denovo_diff_m00="denovo-like-diff_m00.mnc"; #this is resampled to the m00 timepoint diffusion
$gdles_diff_m00 = "gdles-like-diff_m00.mnc"; #this is resampled to the m00 timepoint diffusion
$gdles_m00_diff = "gdles_m00-like-diff.mnc"; #this is resampled to the current timepoint diffusion

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

if ($visit ne "m00"){
	print "mincresample -near -like $diff -transformation $b02stx -invert $denovo $denovo_diff\n";
    `mincresample -near -like $diff -transformation $b02stx -invert $denovo $denovo_diff`  unless -e $denovo_diff;
    print "mincresample -near -like $diff -transformation $b02stx -invert $enlarg $enlarging_diff\n";
    `mincresample -near -like $diff -transformation $b02stx -invert $enlarg $enlarging_diff`  unless -e $enlarging_diff;
    print "mincresample -near -like $diff -transformation $b02stx -invert $t2les_m00 $t2lesm00_diff\n";
    `mincresample -near -like $diff -transformation $b02stx -invert $t2les_m00 $t2lesm00_diff`  unless -e $t2lesm00_diff;
    ## bring the denovo mask back to m00 space to see what the params were
    $dir_m00 = "../*".$sub1."_".$sub2."_".$num."_m00*/";
    $b02stx_m00 = `\\ls $dir_m00$b02stx`;chomp($b02stx_m00);
    $diff_m00 = `\\ls $dir_m00$diff`; chomp($diff_m00);
    $gdles_m00 = `\\ls $dir_m00$gdles_only`; chomp($gdles_m00);

     print "mincresample -near -like $diff_m00 -transformation $b02stx_m00 -invert $denovo $denovo_diff_m00\n";
    `mincresample -near -like $diff_m00 -transformation $b02stx_m00 -invert $denovo $denovo_diff_m00`  unless -e $denovo_diff_m00;
    ## bring the gd mask back to m00 space to see what the params were
    ## if there is a gd lesion at this time point
     print "mincresample -near -like $diff_m00 -transformation $b02stx_m00 -invert $gdles_only $gdles_diff_m00\n";
    `mincresample -near -like $diff_m00 -transformation $b02stx_m00 -invert $gdles_only $gdles_diff_m00`  unless -e $gdles_diff_m00;
     ## if there was a gd lesion at m00, what are the params now?
     print "mincresample -near -like $diff -transformation $b02stx -invert $gdles_m00 $gdles_m00_diff\n";
    `mincresample -near -like $diff -transformation $b02stx -invert $gdles_m00 $gdles_m00_diff`  unless -e $gdles_m00_diff;
 }

## WM mask should avoid PVE because artificially increases g
#print "minccalc -expr A[0]&&A[1] brainmask_MTsat.mnc $wm_diff wm-final.mnc\n";
#`minccalc -expr "A[0]&&A[1]" brainmask_MTsat.mnc $wm_diff wm-final.mnc`;
#$wm_diff= "wm-final.mnc";
## make wm mask more stringent
# to do this, we need the classfied volume in tal space, transform back to native and use Mathieu's code for majorityvoting
$old=0; #this is the older way, not using prob maps
if ($old==1) {
	print "mincresample -nearest -tfm_input_sampling -transformation $b02stx -invert -step 1 1 1 -znelements 180 $seg seg_diffspace_highres.mnc\n";
	`mincresample -nearest -tfm_input_sampling -transformation $b02stx -invert -step 1 1 1 -znelements 180 $seg seg_diffspace_highres.mnc`;
	print "run_matlab(maskMajorityVoting('anat-to-diff_mnc.mnc', 'seg_diffspace_highres.mnc', 9, 'final-nrx')) \n";
	run_matlab("maskMajorityVoting('anat-to-diff_mnc.mnc', 'seg_diffspace_highres.mnc', 9, 'final-nrx');") unless -e "final-nrx_wm_perc.mnc";
	print "minccalc -expression  A[0] > 99.1? 1 : 0 final-nrx_wm_perc.mnc wm-nrx_100pc_mask.mnc\n";
	`minccalc -expression  "A[0] > 99.1? 1 : 0" final-nrx_wm_perc.mnc wm-nrx_100pc_mask.mnc`;
	$wm_stringent = "wm-nrx_100pc_mask.mnc";
}else{
	$wm_stringent = $wm_diff;
}

print "mincmath -mult $wm_stringent brainmask_MTsat.mnc WM-100pc_MTsatmasked.mnc -clob\n";
`mincmath -mult $wm_stringent brainmask_MTsat.mnc WM-100pc_MTsatmasked.mnc -clob`;
$wm_diff= "WM-100pc_MTsatmasked.mnc";

print "\n________________Some stats________________\n";
# Average values in rois, based on type of lesion
@maps=("g_MTsat_g.mnc","g_MTsat_mvf.mnc","g_MTsat_avf.mnc");
@maps_raw=("MTsat_images_MTsat.mnc","AMICO/NODDI/FIT_ICVF_mnc.mnc","AMICO/NODDI/FIT_ISOVF_mnc.mnc"); #HACK I extract values for these rois to write to file using the same loop, it works because there are also 3
@maps_raw_name=("mtsat","icvf","iso");

if ($visit eq "m00"){
	@masks=($t2les_diff,$gdles_diff,$bh_diff,$wm_diff);
	@masks_name=("t2les","gd","bh","wm");
}else{
	@masks=($t2lesm00_diff,$denovo_diff,$enlarging_diff,$gdles_diff,$bh_diff,$wm_diff);
	@masks_name=("t2les","denovo","enlarging","gd","bh","wm");
}
	
for($m=0;$m<scalar(@maps);$m++) {
	($namebase,$dir,$ext)=fileparse($maps[$m],'\..*');
	$lst='';
	for ($s=0;$s<scalar(@masks);$s++){
		$lst = $lst.' '.$maps[$m].' '.$masks[$s]; #map name and mask name pair
		#extract values in raw maps  and write to file
		print "mincsample $maps_raw[$m] -mask $masks[$s] > $maps_raw_name[$m]-$masks_name[$s].txt\n";
		`mincsample $maps_raw[$m] -mask $masks[$s] > $maps_raw_name[$m]-$masks_name[$s].txt\n`;
	}
	print "--> box_plot.py $lst $namebase-in-masks.png\n\n";
	`box_plot.py $lst $namebase-in-masks.png`;

}
# Average values in rois, based on time course
if ($visit ne "m00"){
	@maps=("g_MTsat_g.mnc","g_MTsat_mvf.mnc","g_MTsat_avf.mnc");
	$dir_m00 = "../*".$sub1."_".$sub2."_".$num."_m00*/";
	@maps_m00=($dir_m00.$maps[0],$dir_m00.$maps[1],$dir_m00.$maps[2]);
	
	# denovo m06 lesions and what they were at m00
	@masks = ($denovo_diff_m00,$denovo_diff);
	for($m=0;$m<scalar(@maps);$m++) {
		($namebase,$dir,$ext)=fileparse($maps[$m],'\..*');
		$lst = ' '.$maps_m00[$m].' '.$masks[0].' '.$maps[$m].' '.$masks[1]; #map name and mask name pair
		print "--> box_plot.py $lst $namebase-denovo-in-time.png\n\n";
		`box_plot.py $lst $namebase-denovo-in-time.png`;
	}
    # gd lesions and what they were on m00 (TO DO: gdles is included in denovo so have to modify this)
	@masks = ($gdles_diff_m00,$gdles_diff);
	for($m=0;$m<scalar(@maps);$m++) {
		($namebase,$dir,$ext)=fileparse($maps[$m],'\..*');
		$lst = ' '.$maps_m00[$m].' '.$masks[0].' '.$maps[$m].' '.$masks[1]; #map name and mask name pair
		print "--> box_plot.py $lst $namebase-gd-in-time.png\n\n";
		`box_plot.py $lst $namebase-gd-in-time.png`;
	}
	# gd lesions on m00 and what they are now (TO DO: gdles is included in denovo so have to modify this)
	$gdles_m00_diff_m00 = `\\ls $dir_m00$gdles_diff`;chomp($gdles_m00_diff_m00);#in m00 space 
	@masks = ($gdles_m00_diff_m00,$gdles_m00_diff);
	for($m=0;$m<scalar(@maps);$m++) {
		($namebase,$dir,$ext)=fileparse($maps[$m],'\..*');
		$lst = ' '.$maps_m00[$m].' '.$masks[0].' '.$maps[$m].' '.$masks[1]; #map name and mask name pair
		print "--> box_plot.py $lst $namebase-m00gd-in-time.png\n\n";
		`box_plot.py $lst $namebase-m00gd-in-time.png`;
	}

}
exit(1);
### Gather data from masks and store in text files for later use (like to be read by matlab)
#####

print ID $id;
print VISIT $visit;
# extract the values in each mask
@masks=($wm_diff,,$gdles_diff);
@masks_n=('les','NAWM');
@maps=("MTsat_images_MTsat.mnc","AMICO/NODDI/FIT_ICVF_mnc.mnc","AMICO/NODDI/FIT_ISOVF_mnc.mnc");
@maps_n=('mofg','MTR','FVF');
    for ($i=0;$i<scalar(@maps);$i++){
      for ($j=0;$j<scalar(@masks);$j++){
        print "minccalc -expr A[0]*A[1] $maps[$i] $masks[$j] $tmp/$maps_n[$i]-$masks_n[$j].mnc\n";
        `minccalc -expr "A[0]*A[1]" $maps[$i] $masks[$j] $tmp/$maps_n[$i]-$masks_n[$j].mnc`;
        print "mincextract -range 0.01 10 $tmp/$maps_n[$i]-$masks_n[$j].mnc | fgrep -vx 0 >> $tmp/vals_$maps_n[$i]-$masks_n[$j]\n";
        `mincextract -range 0.01 10 $tmp/$maps_n[$i]-$masks_n[$j].mnc | fgrep -vx "0" >> $tmp/vals_$maps_n[$i]-$masks_n[$j]`;
      }
    }
#mincsample g_MTsat_g.mnc -mask t2les-like-diff.mnc

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
