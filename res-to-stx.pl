#!/usr/bin/perl -w

### Resample maps to stereotaxic space
##
##########################################################
##########################################################
## The lesion masks from NeuroRx are in stx space, want to bring the maps we have to that space
#
## X = 1st inital x=2nd initial d=patient number mxx=month time point
##
## > m00
##  - X-x_00d/mxx/*_newt2f_TREF-mxx_ISPC-stx152iso.mnc.gz   #new lesions
##
## mxx
##  - X-x_00d/mxx/*_ct2f_ISPC-stx152iso.mnc.gz
##				-t2 lesions (hypo on T1p, hyper on FLAIR)
##				-in stx 1mm3 space
##  - X-x_00d/mxx/*_gvf_ ISPC-stx152iso.mnc.gz 
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
# Modified July 2019
# - found a bug whereby the denovo were not being read properly?? double check all other cases

require 5.001;
use Getopt::Tabular;
use MNI::Startup qw(nocputimes);
use MNI::Spawn;
use MNI::FileUtilities qw(check_output_dirs);
use File::Basename;
use List::MoreUtils qw(zip);
use Cwd qw(cwd);
use List::MoreUtils qw(first_index);

$thresh = 1600; #t1 threshold to separate black holes and t2 lesions
$do_plots = 0; #flag to make plots again (requires interaction to close them) (0 means we make them :p)
# Pattern for each file
$t2les_pat = "/*_ct2f_ISPC-stx152iso.mnc.gz"; #in 1mmiso
$gd_pat = "/*_gvf_ISPC-stx152iso.mnc.gz"; #in 1x1x3mm?? argh not all subjects have it done!!
$wm_pat = "/malfTissuePriors/*prob_cerebrum-wm_ISPC-stx152iso.mnc.gz";
#$seg_pat = 
$t1p2stx_pat = "/RegisterToStandardSpace/*/*_t1p-to-stx152lsq6.xfm";
$t1map2stx_pat = "/RegisterOtherModalities/*/*_t1gMPR-t1map-to-stx152lsq6.xfm";

# If not 1st timepoint, we have other new T2 lesion mask as well, it has lesions that were not there at the previous time point
$denovo_pat = "/*_newt2f_TREF-m*_ISPC-stx152iso.mnc.gz"; #new lesions

if($0 =~ /[\/A-Za-z_0-9-]+\/([A-Za-z0-9_-]+\.pl)/) {$program = $1;}	#program name without path
$Usage = <<USAGE;

Resamples maps to Standard Space


Usage: $program
-help for options

USAGE
##defaults
open(ID,">id.txt");
open(VISIT,">visit.txt");
$threechar=0; #subject name has 3 chars
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
    $threechar=1;
} 
else {die "can't find subject info for $pwd\n";}

print "\n---On loki---\n";

#print "tar -cvf /tmp/$tar $dir/*_ct2f.mnc* $dir/*_gvf.mnc.gz $dir/TissueClassification/*/*ANAT-na-cerebrum-wm_ISPC-stx152lsq6.mask.mnc.gz $dir/TissueClassification/*/*segmentation_ISPC-stx152lsq6.mnc.gz $dir/RegisterToStandardSpace/*/*_t1p-to-stx152lsq6.xfm $dir/ResampleToStx/*/*t1gMPR-t1map_ISPC-stx152lsq6.mnc.gz\n";
if ($visit eq "m00"){
	print "tar -cvf /tmp/$tar $dir$t2les_pat $dir$gd_pat $dir$wm_pat $dir$t1p2stx_pat $dir$t1map2stx_pat\n";
}else{
	print "tar -cvf /tmp/$tar  $dir$t2les_pat $dir$gd_pat $dir$wm_pat $dir$t1p2stx_pat $dir$t1map2stx_pat $dir$denovo_pat\n";
}

print "\nscp ileppert\@loki.bic.mni.mcgill.ca:/tmp/$tar .\n\n ";
unless (-e $tar){
	die "---> first download the data from loki, there should be at least 4 files\n";
}

`tar xvf $tar`;

chomp($unique=`date +%y%m%d-%H:%M`); 
print "========== Date is $unique =============\n";
print "**** $program @ARGV ****\n";

##Get the masks
$wm = `\\ls $dir$wm_pat`; chomp($wm);
$gd = `\\ls $dir$gd_pat`; chomp($gd);
if (!-e $gd){ #make my own volume of 0s... waste of memory but simplifies computations later !NO! wait for it to be done...
	die "there is no Gd lesion mask\n";
	$gd = $dir.$visit."_gvf_ISPC-stx152iso.mnc.gz";
	print "--Make a 0 gd mask...\n";
	print "mincmath -mult -const 0 $wm $gd unless -e $gd\n";
	`mincmath -mult -const 0 $wm $gd` unless -e $gd;
}
$t1p2stx = `\\ls $dir$t1p2stx_pat`;chomp($t1p2stx);
$t1map2stx = `\\ls $dir$t1map2stx_pat`;chomp($t1map2stx);
if ($visit eq "m00"){
	$t2les = `\\ls $dir$t2les_pat`; chomp($t2les);
}else{
	$denovo = `\\ls $dir$denovo_pat`; chomp($denovo);
	$t2les = `\\ls $dir$t2les_pat`;chomp($t2les);
	if (!-e $t2les){

		die "there is no T2 lesion mask\n";
		if (0){ # I stopped doind this, should stay consistent in map making...if it's not on the server, the case is not complete
				#get t2les mask from previous timepoint and make a t2les for this timepoint (sometine it's not yet done through the pipeline so i make my own)
				@timepoints=("m00","m06","m12","m18");
				$tp=first_index { $_ eq $visit } @timepoints;
				$tp_prev=$tp-1;
				if (! $threechar){
					$dir_prev = "../*".$sub1."_".$sub2."_".$num."_".$timepoints[$tp_prev]."*/*/$timepoints[$tp_prev]";
				}else{
					$dir_prev = "../*".$sub."_".$num."_".$timepoints[$tp_prev]."*/*/$timepoints[$tp_prev]";
				}
				$t2les_prev = `\\ls $dir_prev$t2les_pat`; chomp($t2les_prev);
				if (!-e $t2les_prev){die "can't find previous t2les $t2les_prev, make sure you process previous timepoints first\n";}
				print "--Make t2les from denovo and previous t2les mask---\n";
				$t2les=$dir.$visit."_ct2f_ISPC-stx152iso.mnc.gz";
				print "minccalc -expr A[0]>0 || A[1]>0 $t2les_prev $denovo $t2les\n";
				`minccalc -expr "A[0]>0 || A[1]>0" $t2les_prev $denovo $t2les` unless -e $t2les;
		}

	}
	
}



@t1ps=`more mnclist | grep mni_gre_3D_T1w-PreGd`;
@file=split(/\s+/,$t1ps[1]);
$t1p = $file[0];#2nd should be the normalized one used by NeuroRx for registation
$diff = `ls b0_eddy_corr_mnc.mnc`; chomp($diff);
$diff2t1p  = `\\ls dti-to-t1b.xfm`; chomp($diff2t1p);
@t1map_nats = `more mnclist | grep MP2RAGE_1mm_T1_Images`;
@file = split(/\s+/,$t1map_nats[0]);
$t1map_nat = $file[0];

##Maps
$g = `\\ls g_MTsat_g.mnc`; chomp($g);
$mvf = `\\ls g_MTsat_mvf.mnc`; chomp($mvf);
$avf = `\\ls g_MTsat_avf.mnc`; chomp($avf);

#Getopt::Tabular::SetHelp ($Usage, '');
#my @args;
#GetOptions(\@args_table, \@ARGV, \@args) || exit 1;
#die $Usage unless $#ARGV >=0;


### Resample masks
#transformations
print "\n________________Compute transformations________________\n";
#stx-to-b0
$b02stx="b0-to-stx.xfm";
print "xfmconcat $diff2t1p $t1p2stx $b02stx\n";
`xfmconcat $diff2t1p $t1p2stx $b02stx`;

## Everything is registered to the diffusion, just need to apply b0-to-stx.xfm to the maps
# check that everything is here
if (!-e $t2les){die "----Can't find t2les mask\n";} #all time points now have t2les masks and >m00 have a newt2 mask
if (!-e $gd){die "----Can't find gd mask\n";} 
if ($visit ne "m00"){
	$denovo = `\\ls $dir$denovo_pat`; chomp($denovo);
	if (!-e $denovo){die "----Can't find newt2 mask\n";} 
}
if (-e $t1map2stx && -e $t1map_nat){
	print "________________Create T1map in stx space________________\n";
	$t1map = "t1map-stx152lsq6iso.mnc";
	print "mincresample $t1map_nat -like $wm -transform $t1map2stx $t1map\n";
	`mincresample $t1map_nat -like $wm -transform $t1map2stx $t1map` unless -e $t1map;
}else{die "----Can't find T1map transformation or native T1map\n";}

print "________________Create a black-hole type mask, with a threshold on T1 (looks dark on lepto)________________\n";
$blakh = "black-holes-thres$thresh-stx152lsq6iso.mnc";
	print "minccalc -expr A[0]>$thresh && A[1]>=1  $t1map $t2les $blakh\n";
	`minccalc -expr "A[0]>$thresh && A[1]>=1" $t1map $t2les $blakh` unless -e $blakh;


print "\n________________Make sure masks are mutually exclusive________________\n";
# Gd lesions are a subset of T2 lesions
# "black holes" are a subset of T2 lesions
# the WM  mask is a probability mask, threshold at 95%?
$thresh_wm = 0.95;
$t2les_only="t2les-only.mnc"; # not black holes nor WM
$gdles_only="gdles-only.mnc";
$wm_only="wm-only.mnc";
$denovo_only="denovo-only.mnc";

if (-e $t2les && -e $wm){
	if ($visit eq "m00"){
		# not black holes nor WM 
		print "--> T2 lesions only <--\n minccalc -expr A[0]<$thresh_wm && A[1]!=1 && A[2]>=1  $wm $blakh $t2les $t2les_only\n";
		`minccalc -expr "A[0]<$thresh_wm && A[1]!=1 && A[2]>=1"  $wm $blakh $t2les $t2les_only`;
	}else{ # if >m00, we also have to exclude the denovo lesions in t2les, 
		print "--> T2 lesions only <--\n minccalc -expr A[0]<$thresh_wm && A[1]<1 && A[2]!=1 && A[3]>=1  $wm $denovo $blakh $t2les $t2les_only\n";
		`minccalc -expr "A[0]<$thresh_wm && A[1]<1 && A[2]!=1 && A[3]>=1"  $wm $denovo $blakh $t2les $t2les_only`;
		# denovo should exlude gd and WM, and potentially t2les from previous timepoints?
		print "--> Denovo only <--\n minccalc -expr A[0]<$thresh_wm && A[1]<1 && A[2]>0  $wm $gd $denovo $denovo_only\n";
		`minccalc -expr "A[0]<$thresh_wm && A[1]<1 && A[2]>0" $wm $gd $denovo $denovo_only`;
	}

	# not T2 lesions, only WM
	print "--> WM only <--\n minccalc -expr A[0]!=1  && A[1]>=$thresh_wm $t2les $wm $wm_only\n";
   `minccalc -expr "A[0]!=1 && A[1]>=$thresh_wm" $t2les $wm $wm_only` unless -e $wm_only;
}
if (-e $gd && -e $wm){
#Gad lesions are include in T2 lesions?
	print "--> Gad Les only <--\n minccalc -expr A[0]<$thresh_wm && A[1]>0   $wm_only $gd $gdles_only\n";
	`minccalc -expr "A[0]<$thresh_wm && A[1]>0 "  $wm_only $gd $gdles_only`;
}


print "\n________________Register each map to stx space___________\n";
$g_stx = "g_MTsat_g-stx.mnc";
$mvf_stx ="g_MTsat_mvf-stx.mnc";
$avf_stx = "g_MTsat_avf-stx.mnc";
# for stats on reprodicibility add this redundant quantity
$vic_1_viso = "vic_1-viso.mnc";
$vic_1_viso_stx = "vic_1-viso-stx.mnc";

print "mincresample $g -like $t2les -transform $b02stx $g_stx\n";
`mincresample $g -like $t2les -transform $b02stx $g_stx` unless -e $g_stx;
print "mincresample $mvf -like $t2les -transform $b02stx $mvf_stx\n";
`mincresample $mvf -like $t2les -transform $b02stx $mvf_stx` unless -e $mvf_stx;
print "mincresample $avf -like $t2les -transform $b02stx $avf_stx\n";
`mincresample $avf -like $t2les -transform $b02stx $avf_stx` unless -e $avf_stx;

print "minccalc -expr A[0]*(1-A[1]) AMICO/NODDI/FIT_ICVF_mnc.mnc AMICO/NODDI/FIT_ISOVF_mnc.mnc $vic_1_viso\n";
`minccalc -expr "A[0]*(1-A[1])" AMICO/NODDI/FIT_ICVF_mnc.mnc AMICO/NODDI/FIT_ISOVF_mnc.mnc $vic_1_viso` unless -e $vic_1_viso;
print "mincresample $vic_1_viso -like g_MTsat_avf-stx.mnc $vic_1_viso_stx\n";
`mincresample $vic_1_viso -like g_MTsat_avf-stx.mnc $vic_1_viso_stx` unless -e $vic_1_viso_stx;



print "\n________________Some stats for this tp________________\n";
# Average values in rois, based on type of lesion
$plt="plots/";
`mkdir $plt` unless -e -d $plt;
@maps=($g_stx,$mvf_stx,$avf_stx);
#@maps_raw=("MTsat_images_MTsat.mnc","AMICO/NODDI/FIT_ICVF_mnc.mnc","AMICO/NODDI/FIT_ISOVF_mnc.mnc"); #HACK I extract values for these rois to write to file using the same loop, it works because there are also 3
#@maps_raw_name=("mtsat","icvf","iso");

if ($visit eq "m00"){
	@masks=($t2les_only,$gdles_only,$blakh,$wm_only);
	@masks_name=("t2les","gd","bh","wm");
}else{
	@masks=($t2les_only,$denovo_only,$gdles_only,$blakh,$wm_only);
	@masks_name=("t2les","denovo","gd","bh","wm");
}
	
for($m=0;$m<scalar(@maps);$m++) {
	($namebase,$dir,$ext)=fileparse($maps[$m],'\..*');
	$lst='';
	for ($s=0;$s<scalar(@masks);$s++){
		$lst = $lst.' '.$maps[$m].' '.$masks[$s]; #map name and mask name pair
		#extract values in raw maps  and write to file
		#print "mincsample $maps_raw[$m] -mask $masks[$s] > $maps_raw_name[$m]-$masks_name[$s].txt\n";
		#`mincsample $maps_raw[$m] -mask $masks[$s] > $maps_raw_name[$m]-$masks_name[$s].txt\n`;
	}
	print "--> box_plot.py $lst $plt$namebase-in-masks.png\n\n";
	`box_plot.py $lst $plt$namebase-in-masks.png` unless $do_plots;
}

if ($visit eq "m00"){exit(1);}

print "\n_______________Get the masks and maps from previous time points________\n";
@timepoints=("m00","m06","m12","m18");
$tp=first_index { $_ eq $visit } @timepoints;
## Get new lesion masks from all previous time points and create symbolic links to use later

for ($t=($tp-1);$t>=0;$t--){
	$dir_prev = "../*".$sub1."_".$sub2."_".$num."_".$timepoints[$t]."*/";
	$idx=$t-$tp;
	# make a link to the previous timepoint dir
	`ln -s $dir_prev link-to_$idx` unless -e "link-to_$idx";
	print "cp $dir_prev/type.txt . \n";
	`cp $dir_prev."/type.txt" .`;# copy the MS type
	# Get Gd lesions prev time points
	$gdles_prev = `\\ls $dir_prev$gdles_only`; chomp($gdles_prev); #all tps should have a gd mask

	if (-e $gdles_prev){`ln -s $gdles_prev gdles_from_$idx` unless -e "gdles_from_$idx";} else{die "can't find $timepoints[$t] gd les\n";} #make a consistently named link to the previous masks(easier to use with -1 -2 etc.)
	# Get  New lesions from prev time points
	if($t > 0){ # only tps >= m06 have denovo lesions
		$denovo_prev = `\\ls $dir_prev$denovo_only`; chomp($denovo_prev); 
		if (-e $denovo_prev){`ln -s $denovo_prev denovo_from_$idx` unless -e "denovo_from_$idx";} else{die "can't find $timepoints[$t] denovo les\n";} #make a consistently named link to the previous masks(easier to use with -1 -2 etc.)
	}
	# Get t2 lesions prev time points
	$t2les_prev = `\\ls $dir_prev$t2les_only`; chomp($t2les_prev); #all tps should have a t2lesion mask
	if (-e $t2les_prev){`ln -s $t2les_prev t2les_from_$idx` unless -e "t2les_from_$idx";} else{die "can't find $timepoints[$t] t2 les\n";} #make a consistently named link to the previous masks(easier to use with -1 -2 etc.)

	# Make the resolving masks
	print "--> Resolving T2 lesions\n";
	print "minccalc -expr (A[0]-A[1])>=1 t2les-only.mnc $t2les_prev resolve-t2-from_$idx.mnc\n";
	`minccalc -expr "(A[0]-A[1])>=1" t2les-only.mnc $t2les_prev resolve-t2-from_$idx.mnc`;
	print "--> Resolving Black hole lesions\n";
	$blackh_prev = `\\ls $dir_prev$blakh`; chomp($blackh_prev); #all tps should have a blackh mask
	print "minccalc -expr (A[0]-A[1])>=1 $blakh $blackh_prev resolve-blackh-from_$idx.mnc\n";
	`minccalc -expr "(A[0]-A[1])>=1" $blakh $blackh_prev resolve-blackh-from_$idx.mnc`;

	#get the maps
	$g_prev = `\\ls $dir_prev$g_stx`; chomp($g_prev); #all tps should have gratio map
    if (-e $g_prev){`ln -s $g_prev g_from_$idx` unless -e "g_from_$idx";} else{die "can't find gratio map $dir_prev$g_stx\n";} #make a consistently named link to the previous masks(easier to use with -1 -2 etc.)
	$mvf_prev = `\\ls $dir_prev$mvf_stx`; chomp($mvf_prev); #all tps should have gratio map
    if (-e $mvf_prev){`ln -s $mvf_prev mvf_from_$idx` unless -e "mvf_from_$idx";} else{die "can't find mvf map $dir_prev$mvf_stx\n";} #make a consistently named link to the previous masks(easier to use with -1 -2 etc.)
	$avf_prev = `\\ls $dir_prev$avf_stx`; chomp($avf_prev); #all tps should have avf map
    if (-e $g_prev){`ln -s $avf_prev avf_from_$idx` unless -e "avf_from_$idx";} else{die "can't find avf map $dir_prev$avf_stx\n";} #make a consistently named link to the previous masks(easier to use with -1 -2 etc.)
}


## Now let's look at some time courses
# first get the maps at each time point
for($s=-$tp;$s<=0;$s++){ #number of time points to plot
	if($s == 0){ #this time point
		push(@mg,$g_stx);
		push(@mmvf,$mvf_stx);
		push(@mavf,$avf_stx);
	}else{
		push(@mg,"g_from_$s");
		push(@mmvf,"mvf_from_$s");
		push(@mavf,"avf_from_$s");
	}		
}
# Now let's get the lesion masks we want to look at
$denovo_mask='';
for($t=-$tp;$t<=0;$t++){ #could be new at each time point
	if($t>-$tp){ #only denovos as >=m06
		if($t == 0){ #this time point
			$denovo_mask=$denovo_only;
		}else{$denovo_mask="denovo_from_$t";}
		#print "denovo_mask:$denovo_mask with these maps: @mg\n";
		print "---Denovo at time $t---\n";
		$lst = join(' '.$denovo_mask.' ',@mg);
		$lst = $lst." ".$denovo_mask;
		print "--> box_plot.py $lst $plt/g-in-denovo_$t.png\n\n";
		`box_plot.py $lst $plt/g-in-denovo_$t.png` unless $do_plots;
		
		if ($tp>=2){ #for now only plot avf and mvf when more than 2 time points
			$lst = join(' '.$gdles_mask.' ',@mmvf);
			$lst = $lst." ".$gdles_mask;
			print "--> box_plot.py $lst $plt/mvf-in-gdles_$t.png\n\n";
			`box_plot.py $lst $plt/mvf-in-gdles_$t.png` unless $do_plots;
			$lst = join(' '.$gdles_mask.' ',@mavf);
			$lst = $lst." ".$gdles_mask;
			print "--> box_plot.py $lst $plt/avf-in-gdles_$t.png\n\n";
			`box_plot.py $lst $plt/avf-in-gdles_$t.png` unless $do_plots;
		}

	}
	if($t == 0){ #this time point
		$gdles_mask=$gdles_only;
	}else{$gdles_mask="gdles_from_$t";}
	# gd les masks
	#print "gdles_mask:$gdles_mask with these maps: @mg\n";
	print "---Gd lesion at time $t---\n";
	$lst = join(' '.$gdles_mask.' ',@mg);
	$lst = $lst." ".$gdles_mask;
	print "--> box_plot.py $lst $plt/g-in-gdles_$t.png\n\n";
	`box_plot.py $lst $plt/g-in-gdles_$t.png` unless $do_plots;
	
}
# t2les at baseline, what they look like over time
$t2les_mask="t2les_from_-$tp"; #this should be m00
#print "t2les_mask:$t2les_mask with these maps: @mg\n";
print "---T2les at time -$tp---\n";
$lst = join(' '.$t2les_mask.' ',@mg);
$lst = $lst." ".$t2les_mask;
print "--> box_plot.py $lst $plt/g-in-t2les_-$tp.png\n\n";
`box_plot.py $lst $plt/g-in-t2les_-$tp.png` unless $do_plots;