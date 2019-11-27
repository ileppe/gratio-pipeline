#!/usr/bin/perl -w

### Process g-ratio data
#
##########################################################
##########################################################
##
##
###########################################
#
###########################################
# Created by Ilana Leppert Sept 2017 (copied in most part from g-ratio_pipeline)
###########################################
## CHECK INSTRUCTIONS AT BOTTOM

require 5.001;
use Getopt::Tabular;
use MNI::Startup qw(nocputimes);
use MNI::Spawn;
use MNI::FileUtilities qw(check_output_dirs);
use File::Basename;
use List::MoreUtils qw(zip);

if($0 =~ /[\/A-Za-z_0-9-]+\/([A-Za-z0-9_-]+\.pl)/) {$program = $1;}	#program name without path
$Usage = <<USAGE;

Runs the g-ratio pipeline.


Usage: $program
-help for options

USAGE
#if($ARGV[0]=~/help/i){print $Usage; exit(1);}

#defaults:
$doing_MTsat=1;
$MTW="";
#$scale_MTsat=0.089; #from comparison to F in mb110515; empirical .  by comparing to F with scaling factor 1.5!!!
$scale_MTsat=0.121;# by empirically trying to get mvf ~0.3 in WM of a control
 

#presets: DO: must create these from the inputs. these are currently hardcoded; links are for last thing processed
$acq="/data/MRS_Remyelination/ilana/repeat-lesions/acqparams.txt";

# Merged bvecs and bvals: ## bvecs has to be made for each dataset independently, because of angulation
#$bvecs="/export02/data/ilana/for-others/repeat-lesions/bvecs";# 10(b300)A-P; 10(b300)P-A;30(b700)A-P; 30(b700)P-A; 64(b2500)A-P; 64(b2500)P-A; DWIs have 1 b=0
$bvals="/data/MRS_Remyelination/ilana/repeat-lesions/bvals";


# index file that tells eddy which line/of the lines in the acqparams.txt file that are relevant for the data passed into eddy.
$index="/data/MRS_Remyelination/ilana/repeat-lesions/index.txt"; # 10(b300)A-P; 10(b300)P-A;30(b700)A-P; 30(b700)P-A; 64(b2500)A-P; 64(b2500)P-A; DWIs have 1 b=0

# NEEd this for AMICO, should be made with the the bvecs (after eddy) and bvals bvecs_bvals2camino.pl
#$scheme="NODDI.scheme"; #appended to the data dir

#totframes: DO: calculate from sum of all dim4.
$totframes=214;

# if using the NeuroRx pipeline for segmentation
$nrx=1;

my @args_table = (#["-clobber","boolean",undef,\$clobber,"clobber all output files (currently not implemented, sorry)"],#HERE not implemented
		  ["-diff","string",6,\@diff_inputs,"diffusion input gzipped nii files in order b=low AP b=mid AP b=high AP b=low PA b=mid PA b=high PA"],
    ["-bvecs","string",6,\@bvecs_input,"bvec files in order b=low AP b=mid AP b=high AP b=low PA b=mid PA b=high PA"], 
["-anat","string",2,\@anatfiles,"anatomical (T1w pre-contrast) nii file and minc file."],
    ["-MTsat","string",3,\@MTsat,"MTsat input files in order MTW PDW T1W.  **minc** files."],
    ["-B1","string",2,\@B1files,"ep_seg_se_b1 60 and 120 for AFI B1 calculation. minc files. needed for all MT"],
   #-acqparams","string",1,\$acq,"acq params for topup"],
);



Getopt::Tabular::SetHelp ($Usage, '');
my @args;
GetOptions(\@args_table, \@ARGV, \@args) || exit 1;
die $Usage unless $#ARGV >=0;

#my $subdir=$args[0]; #the input dir
chomp($unique=`date +%y%m%d-%H:%M`); 
print "--> Date is $unique\n";
print "**** $program @ARGV ****\n";

#if just processing qMT data:
#goto qMT1;

$diff_10_AP=$diff_inputs[0];
$diff_30_AP=$diff_inputs[1];
$diff_64_AP=$diff_inputs[2];
$diff_10_PA=$diff_inputs[3];
$diff_30_PA=$diff_inputs[4];
$diff_64_PA=$diff_inputs[5];

#make bvecs
print "--make bvecs--\n";
$bvecs="bvecs";
@x,@y,@z;
open(BVEC,">$bvecs");
foreach $b (@bvecs_input){
  $xvec=`awk 'NR==1' $b`;
  $xvec =~ s/\s*$//g; #remove whitespace at end (chmop isn't working?)
  push(@x,$xvec);
  $yvec=`awk 'NR==2' $b`;
  $yvec =~ s/\s*$//g;
  push(@y,$yvec);
  $zvec=`awk 'NR==3' $b`;
  $zvec =~ s/\s*$//g;
  push(@z,$zvec);
}
print BVEC "@x\n@y\n@z\n";
$anat=$anatfiles[0];
$T1w_im_mnc=$anatfiles[1];



if($MTsat[0]){$doing_MTsat=1;}
#$doing_MTsat=$MTsat[0]; #IRL why was this here?
if ($doing_MTsat)
{
    $MTW=$MTsat[0];
    $PDW=$MTsat[1];
    $T1W=$MTsat[2];
}

#HERE AFI B1

$B1_60=$B1files[0];
$B1_120=$B1files[1];


#get the number of slices and if odd, add one.
@sliceinfo=`fslinfo $diff_30_AP | grep dim3`;
chomp(@sliceinfo);
($field, $Nslices)=split /\s+/,$sliceinfo[0];

if ($Nslices % 2 == 1)
{
    $Nslices=$Nslices+1;
}


$A2P_P2A_b0="A2P_P2A_b0.nii.gz";
$diff_b0_10AP="b0-10AP.nii.gz";
$diff_b0_30AP="b0-30AP.nii.gz";
$diff_b0_64AP="b0-64AP.nii.gz";
$diff_b0_10PA="b0-10PA.nii.gz";
$diff_b0_30PA="b0-30PA.nii.gz";
$diff_b0_64PA="b0-64PA.nii.gz";

#HERE get total number of images, note 30 or 64 ? how many b=0?
#and make acqparams from them, make bvecs and bvals

# First get the b=0s for blip-up and blip-down
print "\n----Get b=0s for blip-up and blip-down and send to topup----\n";
#IRL have to get 1st frames of 30dirAP and 30dirPA
print "fslroi $diff_10_AP $diff_b0_10AP 0 1 \n";
`fslroi $diff_10_AP $diff_b0_10AP 0 1` unless -e $diff_b0_10AP;
print "fslroi $diff_30_AP $diff_b0_30AP 0 1 \n";
`fslroi $diff_30_AP $diff_b0_30AP 0 1` unless -e $diff_b0_30AP;
print "fslroi $diff_64_AP $diff_b0_64AP 0 1 \n";
`fslroi $diff_64_AP $diff_b0_64AP 0 1` unless -e $diff_b0_64PA;

print "fslroi $diff_10_PA $diff_b0_10PA 0 1 \n";
`fslroi $diff_10_PA $diff_b0_10PA 0 1` unless -e $diff_b0_10PA;
print "fslroi $diff_30_PA $diff_b0_30PA 0 1 \n";
`fslroi $diff_30_PA $diff_b0_30PA 0 1` unless -e $diff_b0_30PA;
print "fslroi $diff_64_PA $diff_b0_64PA 0 1 \n";
`fslroi $diff_64_PA $diff_b0_64PA 0 1` unless -e $diff_b0_64PA;

print "fslmerge -t $A2P_P2A_b0 $diff_b0_10AP $diff_b0_30AP $diff_b0_64AP $diff_b0_10PA $diff_b0_30PA $diff_b0_64PA unless -e $A2P_P2A_b0\n";
`fslmerge -t $A2P_P2A_b0 $diff_b0_10AP $diff_b0_30AP $diff_b0_64AP $diff_b0_10PA $diff_b0_30PA $diff_b0_64PA ` unless -e $A2P_P2A_b0;

$hifi_b0="my_hifi_b0.nii.gz";
$topup_results="topup_results";

  ## For some stupid reason topup does not like an odd # of slices!!!
$A2P_P2A_b0_even="A2P_P2A_b0_even.nii.gz";
print "\n--Make number of slices even or else topup complains--\n\n";
print "fslroi $A2P_P2A_b0 $A2P_P2A_b0_even 0 128 0 128 0 $Nslices unless -e $A2P_P2A_b0_even\n";
`fslroi $A2P_P2A_b0 $A2P_P2A_b0_even 0 128 0 128 0 $Nslices` unless -e $A2P_P2A_b0_even;

print "topup --imain=$A2P_P2A_b0_even --datain=$acq --config=b02b0.cnf --out=$topup_results --iout=$hifi_b0\n";
`topup --imain=$A2P_P2A_b0_even --datain=$acq --config=b02b0.cnf --out=$topup_results --iout=$hifi_b0` unless -e $hifi_b0;


#make a brain mask.  We use the undistorted b=0
print "\n----Make brain mask----\n";
$hifi_b0m="my_hifi_b0m";
$hifi_b0_brain="my_hifi_b0_brain";
$hifi_b0_brain_mask="my_hifi_b0_brain_mask.nii.gz";
$hifi_b0_brain_mask_unz="my_hifi_b0_brain_mask.nii";

print "fslmaths $hifi_b0 -Tmean $hifi_b0m\n";
`fslmaths $hifi_b0 -Tmean $hifi_b0m`;


print "bet $hifi_b0m $hifi_b0_brain -m -f .3\n";
`bet $hifi_b0m $hifi_b0_brain -m -f .3` unless (-e $hifi_b0_brain_mask || -e $hifi_b0_brain_mask_unz);

if (-e $hifi_b0_brain_mask)
{
    `gunzip -f $hifi_b0_brain_mask`;
}


# Merge the datasets
print "\n----Merge the datasets----\n";
$diff="diff-comb.nii.gz";

print "fslmerge -t $diff $diff_10_AP $diff_30_AP $diff_64_AP $diff_10_PA $diff_30_PA $diff_64_PA unless -e $diff\n";
`fslmerge -t $diff  $diff_10_AP $diff_30_AP $diff_64_AP $diff_10_PA $diff_30_PA $diff_64_PA` unless -e $diff;


## Make even # of slices...
$diff_even="diff-comb_even.nii.gz";
print "fslroi $diff $diff_even 0 128 0 128 0 $Nslices 0 $totframes\n";
`fslroi $diff $diff_even 0 128 0 128 0 $Nslices 0 $totframes` unless -e $diff_even;


print "\n----Run eddy----\n";


$eddy_corrected_data="diff_eddy-corr";
$eddy_corrected_data_unz="diff_eddy-corr.nii";
$eddy_corrected_data_z="diff_eddy-corr.nii.gz";

$done=0;
if (-e $eddy_corrected_data_z || -e $eddy_corrected_data_unz)
{
    $done=1;
    print "eddy_openmp --data_is_shelled --imain=$diff_even --mask=$hifi_b0_brain_mask --acqp=$acq --index=$index --bvecs=$bvecs --bvals=$bvals --topup=$topup_results --out=$eddy_corrected_data\n\n";
}

if ($done==0)
{
    print "eddy_openmp --data_is_shelled --imain=$diff_even --mask=$hifi_b0_brain_mask --acqp=$acq --index=$index --bvecs=$bvecs --bvals=$bvals --topup=$topup_results --out=$eddy_corrected_data\n\n";
    `eddy_openmp --data_is_shelled --imain=$diff_even --mask=$hifi_b0_brain_mask --acqp=$acq --index=$index --bvecs=$bvecs --bvals=$bvals --topup=$topup_results --out=$eddy_corrected_data`;
}

#HERE figure out how to rotate the bvecs when do eddy (fdt_rotate_bvecs?) and will then use that as input to NODDI

print "\n----Register anatomical to diffusion series---\n";

#$xfm="diff2str.mat";
#$xfm_in="str2diff.mat";
$b0_eddy="b0_eddy-corr.nii.gz";
$b0_eddy_unz="b0_eddy-corr.nii";
$b0_eddy_corr_mnc="b0_eddy_corr_mnc.mnc";

# get 1st frame of eddy-corrected data
print "fslroi $eddy_corrected_data $b0_eddy 0 1 unless (-e $b0_eddy || -e $b0_eddy_unz)\n";
`fslroi $eddy_corrected_data $b0_eddy 0 1 ` unless (-e $b0_eddy || -e $b0_eddy_unz);

## Aug 24th: Do the registration in mnc, it's better than flirt (and less blurring during resampling)
#make b0 in minc:
`gunzip $b0_eddy` unless -e $b0_eddy_unz;
`nii2mnc $b0_eddy_unz $b0_eddy_corr_mnc` unless -e $b0_eddy_corr_mnc;

#print "flirt -in $b0_eddy -ref $anat -omat $xfm -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6 -cost mutualinfo unless -e $xfm\n";
#`flirt -in $b0_eddy -ref $anat -omat $xfm -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6 -cost mutualinfo` unless -e $xfm;

#print "convert_xfm -omat $xfm_in -inverse $xfm\n";
#`convert_xfm -omat $xfm_in -inverse $xfm`;

#get xfm from t1 to diffusion:

#DO: do I ned to n3 the anat first for this to work?  does flirt?
print "minctracc -identity -mi $b0_eddy_corr_mnc $T1w_im_mnc -lsq6 -debug -threshold 30 5 dti-to-t1.xfm -simplex 1 \n";
`minctracc  -identity -mi $b0_eddy_corr_mnc $T1w_im_mnc -lsq6 -debug -threshold 30 5 dti-to-t1.xfm -simplex 1`;
print "minctracc  -mi $b0_eddy_corr_mnc $T1w_im_mnc -lsq6 -debug -threshold 30 5 -transformation dti-to-t1.xfm dti-to-t1b.xfm -simplex 1 -step 2 2 2\n";
`minctracc -mi $b0_eddy_corr_mnc $T1w_im_mnc -lsq6 -debug -threshold 30 5 -transformation dti-to-t1.xfm dti-to-t1b.xfm -simplex 1 -step 2 2 2`;

`xfminvert dti-to-t1b.xfm str2diff.xfm`;

`rm -f dti-to-t1.xfm`;

#$anat_to_diff="anat-to-diff.nii.gz";
#$anat_to_diff_unz="anat-to-diff.nii";

qMT2:

$anat_to_diff_mnc="anat-to-diff_mnc.mnc";

#goto qMT3;
#print "flirt -in $anat -ref $b0_eddy -applyxfm -init $xfm_in -out  $anat_to_diff\n";
#`flirt -in $anat -ref $b0_eddy -applyxfm -init $xfm_in -out  $anat_to_diff` unless -e $anat_to_diff;
print "mincresample $T1w_im_mnc -like $b0_eddy_corr_mnc -transformation str2diff.xfm $anat_to_diff_mnc\n";
`mincresample $T1w_im_mnc -like $b0_eddy_corr_mnc -transformation str2diff.xfm $anat_to_diff_mnc`;


## Run NODDI on eddy-corrected data
print "\n-----Run NODDI-----\n";


#print "gunzip $eddy_corrected_data.nii.gz\n";
`gunzip $eddy_corrected_data_z`  unless -e $eddy_corrected_data_unz;

#it is all in working dir right now.
($name,$dir,$ext)=fileparse("$eddy_corrected_data_unz", qr/\.[^.]*/);
$basename = $name;

#print "gunzip $hifi_b0_brain_mask\n";
`gunzip $hifi_b0_brain_mask` unless -e $hifi_b0_brain_mask_unz;


($maskbase,$dir,$ext)=fileparse("$hifi_b0_brain_mask", qr/\.[^.]*/);
$mask = $maskbase;
# $noddi_ficvf = $basename."_NODDI_out_ficvf.nii";
# $noddi_fiso = $basename."_NODDI_out_fiso.nii";
# $noddi_ficvf_mnc = $basename."_NODDI_out_ficvf_mnc.mnc";
# $noddi_fiso_mnc= $basename."_NODDI_out_fiso_mnc.mnc";
## Using AMICO is much faster...
$noddi_ficvf = "AMICO/NODDI/FIT_ICVF.nii";
$noddi_fiso = "AMICO/NODDI/FIT_ISOVF.nii";
$noddi_ficvf_mnc = "AMICO/NODDI/FIT_ICVF_mnc.mnc";
$noddi_fiso_mnc= "AMICO/NODDI/FIT_ISOVF_mnc.mnc";

print "--Use the corrected bvecs to run NODDI--\n";
print "bvecs_bvals2camino.pl -vec diff_eddy-corr.eddy_rotated_bvecs -val $bvals -o NODDI.scheme\n";
`bvecs_bvals2camino.pl -vec diff_eddy-corr.eddy_rotated_bvecs -val $bvals -o NODDI.scheme`;
$scheme = "NODDI.scheme";

#print "NODDI_process_JC('$basename','$mask')";
#print "NODDI_process('$basename','$mask')";

#it appears that for the parallel computing, the java virtual machin eis necessary, and there is probably a workaround to run it through this script, but for now, do it manually:
  #run_matlab12b("NODDI_process_JC('$basename','$mask')") unless -e $noddi_ficvf;
#  print "\n!!!! Have to run this within MATLAB in the subject directory...!!!\n\nUse matlab12b until further notice, and remove extra directory slashes:\n\nNODDI_process('$name','$maskbase')\n\nPress Enter when finished to run the rest of this script\n\n";
#  <STDIN>;
## USING AMICO
$tmp=`pwd`; # really annoying dir structure
@blah=split("/",$tmp);#strip the last dir
$here=pop(@blah); chomp ($here);
$diramico=join("/",@blah);

print "AMICO_process('$diramico','$here','$eddy_corrected_data_unz','$hifi_b0_brain_mask_unz','$scheme')\n";
#print "$diramico \n";
run_matlab("AMICO_process('$diramico','$here','$eddy_corrected_data_unz','$hifi_b0_brain_mask_unz','$scheme')") unless -e $noddi_ficvf;
print "...This takes some time...\n";

#afterNODDI:

#use minc from now on:
print "---make minc files from brain mask and NODDI outputs:\n\n";

`nii2mnc my_hifi_b0_brain_mask.nii my_hifi_b0_brain_mask_mnc.mnc` unless -e "my_hifi_b0_brain_mask_mnc.mnc";

`nii2mnc $noddi_ficvf $noddi_ficvf_mnc` unless -e $noddi_ficvf_mnc;
`nii2mnc $noddi_fiso $noddi_fiso_mnc` unless -e $noddi_fiso_mnc;


#the brain mask and the NODDI stuff all have the wrong starts and steps; we write this from anat-to-diff:
#convert anat in diff space to minc:

#`gunzip $anat_to_diff` unless -e $anat_to_diff_unz;
#`nii2mnc $anat_to_diff_unz $anat_to_diff_mnc` unless -e $anat_to_diff_mnc;


$xstart=`mincinfo -attvalue xspace:start $anat_to_diff_mnc`;
chomp($xstart);
$ystart=`mincinfo -attvalue yspace:start $anat_to_diff_mnc`;
chomp($ystart);
$zstart=`mincinfo -attvalue zspace:start $anat_to_diff_mnc`;
chomp($zstart);
$xstep=`mincinfo -attvalue xspace:step $anat_to_diff_mnc`;
chomp($xstep);
$ystep=`mincinfo -attvalue yspace:step $anat_to_diff_mnc`;
chomp($ystep);
$zstep=`mincinfo -attvalue zspace:step $anat_to_diff_mnc`;
chomp($zstep);


print "minc_modify_header -dinsert xspace:start=$xstart -dinsert yspace:start=$ystart -dinsert zspace:start=$zstart -dinsert xspace:step=$xstep -dinsert yspace:step=$ystep -dinsert zspace:step=$zstep <new minc files>\n\n";

  `minc_modify_header -dinsert xspace:start=$xstart -dinsert yspace:start=$ystart -dinsert zspace:start=$zstart -dinsert xspace:step=$xstep -dinsert yspace:step=$ystep -dinsert zspace:step=$zstep my_hifi_b0_brain_mask_mnc.mnc`;

  `minc_modify_header -dinsert xspace:start=$xstart -dinsert yspace:start=$ystart -dinsert zspace:start=$zstart -dinsert xspace:step=$xstep -dinsert yspace:step=$ystep -dinsert zspace:step=$zstep $noddi_ficvf_mnc`;

  `minc_modify_header -dinsert xspace:start=$xstart -dinsert yspace:start=$ystart -dinsert zspace:start=$zstart -dinsert xspace:step=$xstep -dinsert yspace:step=$ystep -dinsert zspace:step=$zstep $noddi_fiso_mnc`;



print "----make the tissue masks:---------\n\n";

$done=0;
if ($nrx ==0){ # not doing NeuroRx segmentation
    if (-d "segmentation")
  {
      $done=1;
      print "standard_pipeline.pl 1 1 $T1w_im_mnc --prefix segmentation/ --model_dir /ipl/quarantine/models/icbm152_model_09c/ --beastlib /ipl/quarantine/models/beast/\n";
  }


  if ($done==0)
  {
    print "--Doing segmentation--\n";
    print "standard_pipeline.pl 1 1 $T1w_im_mnc --prefix segmentation/ --model_dir /ipl/quarantine/models/icbm152_model_09c/ --beastlib /ipl/quarantine/models/beast/\n";
    `standard_pipeline.pl 1 1 $T1w_im_mnc --prefix segmentation/ --model_dir /ipl/quarantine/models/icbm152_model_09c/ --beastlib /ipl/quarantine/models/beast/`;
  }
}



#WM for calculations: use WM.mnc
#WM for vis:WM_blur_bin.mnc
if ($nrx ==0){ # not doing NeuroRx segmentation
  print "xfminvert segmentation/1/1/tal/tal_xfm_1_1_t1w.xfm taltot1w.xfm\n";
  `xfminvert segmentation/1/1/tal/tal_xfm_1_1_t1w.xfm taltot1w.xfm`;
  print "xfmconcat taltot1w.xfm str2diff.xfm taltodiff.xfm\n";
  `xfmconcat taltot1w.xfm str2diff.xfm taltodiff.xfm`;
  print "mincmath -eq -const 3 segmentation/1/1/tal_cls/tal_clean_1_1.mnc WM_tal.mnc\n";
  `mincmath -eq -const 3 segmentation/1/1/tal_cls/tal_clean_1_1.mnc WM_tal.mnc`;
  print "mincresample -nearest -transformation taltodiff.xfm -tfm_input_sampling -like $anat_to_diff_mnc WM_tal.mnc WM.mnc\n";
  `mincresample -nearest -transformation taltodiff.xfm -tfm_input_sampling -like $anat_to_diff_mnc WM_tal.mnc WM.mnc`;
  print "mincblur -fwhm 2 WM.mnc WM\n";
  `mincblur -fwhm 2 WM.mnc WM`;
  print "mincmath -gt -const 0.5 WM_blur.mnc WM_blur_bin.mnc\n";
  `mincmath -gt -const 0.5 WM_blur.mnc WM_blur_bin.mnc`;

  #More stringent masks using Mat's/Sivlain's script
  print "\n----New more stringent tissue masks--\n";
  print "mincresample -nearest -tfm_input_sampling -transformation taltodiff.xfm -step 1 1 3 segmentation/1/1/tal_cls/tal_clean_1_1.mnc cls_diffspace_highres.mnc\n";
  `mincresample -nearest -tfm_input_sampling -transformation taltodiff.xfm -step 1 1 3 segmentation/1/1/tal_cls/tal_clean_1_1.mnc cls_diffspace_highres.mnc`;
  #percentageTissue = maskMajorityVoting(inputMincVolume, inputMincClassified, numberOfClasses, outputBaseName)
  print "run_matlab(maskMajorityVoting('$anat_to_diff_mnc', 'cls_diffspace_highres.mnc', 3, 'final')) \n";
  run_matlab("maskMajorityVoting('$anat_to_diff_mnc', 'cls_diffspace_highres.mnc', 3, 'final');") unless -e "final_wm_perc.mnc";
  print "minccalc -expression  A[0] > 99.1? 1 : 0 final_wm_perc.mnc wm_100pc_mask.mnc\n";
  `minccalc -expression  "A[0] > 99.1? 1 : 0" final_wm_perc.mnc wm_100pc_mask.mnc`;

  #parenchyma, you won't want for g-ratio
  `mincmath -gt -const 1  segmentation/1/1/tal_cls/tal_clean_1_1.mnc parenchyma_tal.mnc`;
  `mincresample -transformation taltot1w.xfm -tfm_input_sampling -like anat-to-diff_mnc.mnc parenchyma_tal.mnc parenchyma.mnc`;
  `mincmath -gt -const 0.95 parenchyma.mnc parenchyma_bin.mnc`;
  print "minccalc -expression  (A[0]> 99.1 || A[1]> 99.1) == 1 ? 1 : 0 final_wm_perc.mnc final_gm_perc.mnc parenchyma_100pc_mask.mnc\n";
  `minccalc -expression  "(A[0]> 99.1 || A[1]> 99.1) == 1 ? 1 : 0" final_wm_perc.mnc final_gm_perc.mnc parenchyma_100pc_mask.mnc`;
}
#} #segmentation


print "--------B1map---------\n\n";

#here grab date&time or something to make tmp dir unique
$tmpdir="tmp0001";
#print "mkdir $tmpdir\n\n";
`mkdir $tmpdir` unless -d $tmpdir;


#b1, need for all MT:

$B1_60_no_ext=rm_ext($B1_60);
$B1_120_no_ext=rm_ext($B1_120);


$B1_60_diffspace=$B1_60_no_ext."_diffspace.mnc";
$B1_120_diffspace=$B1_120_no_ext."_diffspace.mnc";

#IRL these seem to have the wrong order for the registration I changed them
#`minctracc -identity -mi $anat_to_diff_mnc $B1_60 -lsq6  -simplex 1 $tmpdir/xfm.xfm -clob`;
print "minctracc -identity -mi  $B1_60 $anat_to_diff_mnc -lsq6  -simplex 1 $tmpdir/xfm.xfm -clob\n";
`minctracc -identity -mi $B1_60 $anat_to_diff_mnc -lsq6  -simplex 1 $tmpdir/xfm.xfm -clob`;
print "mincresample -transformation $tmpdir/xfm.xfm -tfm_input_sampling -like $anat_to_diff_mnc $B1_60 $B1_60_diffspace \n";
`mincresample -transformation $tmpdir/xfm.xfm -tfm_input_sampling -like $anat_to_diff_mnc $B1_60 $B1_60_diffspace`;
print "minctracc -identity -mi $B1_120 $anat_to_diff_mnc -lsq6  -simplex 1 $tmpdir/xfm.xfm -clob\n";
`minctracc -identity -mi $B1_120 $anat_to_diff_mnc -lsq6  -simplex 1 $tmpdir/xfm.xfm -clob`;
print "mincresample -transformation $tmpdir/xfm.xfm -tfm_input_sampling -like $anat_to_diff_mnc $B1_120 $B1_120_diffspace\n";
`mincresample -transformation $tmpdir/xfm.xfm -tfm_input_sampling -like $anat_to_diff_mnc $B1_120 $B1_120_diffspace`;


#make B1 map regardless of MT chosen.
`mkdir b1` unless -d "b1";

#alpha/lowerFA=arccos(mag(I2/I1))/lowerFA; mathieu filters out >1 b4 cos-1

#note that you can change this to do 'AFI' if that is what was acquired
print "\nfitDataB1_IRL('$B1_60_diffspace', '$B1_120_diffspace', 'b1/b1.mnc', 'my_hifi_b0_brain_mask_mnc.mnc','DA' )\n";
run_matlab("fitDataB1_IRL('$B1_60_diffspace', '$B1_120_diffspace', 'b1/b1.mnc', 'my_hifi_b0_brain_mask_mnc.mnc','DA' )") unless -e "b1/b1.mnc";


if ($doing_MTsat==1) {
  print "-------------MTsat processing-----------\n\n";

  print "------register MTsat base images with the diffusion series-----\n\n";



  $xfm_MTW_to_b0="MTW_to_b0.xfm";
  $xfm_PDW_to_b0="PDW_to_b0.xfm";
  $xfm_T1W_to_b0="T1W_to_b0.xfm";


  #$MTW_like_b0="MTW_like_b0.mnc";
  #$PDW_like_b0="PDW_like_b0.mnc";
  #$T1W_like_b0="T1W_like_b0.mnc";

  $MTW_reg_anat="MTW_reg-to-anat.mnc";
  $PDW_reg_anat="PDW_reg-to-anat.mnc";
  $T1W_reg_anat="T1W_reg-to-anat.mnc";

  #!!!!note that the final protocol should have sufficient coverage to not need to do this, although we can do it anyway, just should be left with a full brain
  #cut off the bad slices from slab profile: hardcoded here to take off 5 slices:



  $MTW_crop=rm_ext($MTW)."_crop.mnc";
  $PDW_crop=rm_ext($PDW)."_crop.mnc";
  $T1W_crop=rm_ext($T1W)."_crop.mnc";



  $xsize=`mincinfo -dimlength xspace $MTW`;
  chomp($xsize);
  $ysize=`mincinfo -dimlength yspace $MTW`;
  chomp($ysize);
  $zsize=`mincinfo -dimlength zspace $MTW`;
  chomp($zsize);

  $new_zsize=$zsize-10;

  print "mincreshape -start 5,0,0 -count $new_zsize,$ysize,$xsize <MT_images>\n\n";

  `mincreshape -start 5,0,0 -count $new_zsize,$ysize,$xsize $MTW $MTW_crop` unless -e "$MTW_crop";
  `mincreshape -start 5,0,0 -count $new_zsize,$ysize,$xsize $PDW $PDW_crop` unless -e "$PDW_crop";
  `mincreshape -start 5,0,0 -count $new_zsize,$ysize,$xsize $T1W $T1W_crop` unless -e "$T1W_crop";


 print "\n---  Compute MTsat  ---\n";
  #register MT images to anat (registered with diffusion) and resample to diffspace sampling:
  #HERE possibly register to the original T1w instead of the downsampled T1w
  print "minctracc -identity -mi $MTW_crop $anat_to_diff_mnc -lsq6  -simplex 1 $xfm_MTW_to_b0 unless -e $xfm_MTW_to_b0\n";
  `minctracc -identity -mi $MTW_crop $anat_to_diff_mnc -lsq6  -simplex 1 $xfm_MTW_to_b0` unless -e $xfm_MTW_to_b0;
  print "minctracc -identity -mi $PDW_crop $anat_to_diff_mnc  -lsq6  -simplex 1 $xfm_PDW_to_b0 unless -e $xfm_PDW_to_b0\n";
  `minctracc -identity -mi $PDW_crop $anat_to_diff_mnc  -lsq6  -simplex 1 $xfm_PDW_to_b0` unless -e $xfm_PDW_to_b0;
  print "minctracc -identity -mi $T1W_crop $anat_to_diff_mnc  -lsq6  -simplex 1 $xfm_T1W_to_b0 unless -e $xfm_T1W_to_b0\n";
  `minctracc -identity -mi $T1W_crop $anat_to_diff_mnc  -lsq6  -simplex 1 $xfm_T1W_to_b0` unless -e $xfm_T1W_to_b0;
  #print "mincresample -tfm_input_sampling -transformation $xfm_MTW_to_b0 -like $anat_to_diff_mnc $MTW_crop $MTW_like_b0\n";
  #`mincresample -tfm_input_sampling -transformation $xfm_MTW_to_b0 -like $anat_to_diff_mnc $MTW_crop $MTW_like_b0`;
  #print "mincresample -tfm_input_sampling -transformation $xfm_PDW_to_b0 -like $anat_to_diff_mnc $PDW_crop $PDW_like_b0\n";
  #`mincresample -tfm_input_sampling -transformation $xfm_PDW_to_b0 -like $anat_to_diff_mnc $PDW_crop $PDW_like_b0`;
  #print "mincresample -tfm_input_sampling -transformation $xfm_T1W_to_b0 -like $anat_to_diff_mnc $T1W_crop $T1W_like_b0\n";
  #`mincresample -tfm_input_sampling -transformation $xfm_T1W_to_b0 -like $anat_to_diff_mnc $T1W_crop $T1W_like_b0`;

  #Keep the MT images in high-res and only downsample the result
  print "mincresample -use_input_sampling -transformation $xfm_MTW_to_b0 $MTW_crop $MTW_reg_anat\n";
  `mincresample -use_input_sampling -transformation $xfm_MTW_to_b0 $MTW_crop $MTW_reg_anat`;
  print "mincresample -use_input_sampling -transformation $xfm_PDW_to_b0 $PDW_crop $PDW_reg_anat\n";
  `mincresample -use_input_sampling -transformation $xfm_PDW_to_b0 $PDW_crop $PDW_reg_anat`;
  print "mincresample -use_input_sampling -transformation $xfm_T1W_to_b0 $T1W_crop $T1W_reg_anat\n";
  `mincresample -use_input_sampling -transformation $xfm_T1W_to_b0 $T1W_crop $T1W_reg_anat`;
  # B1map has to sampled the same way as well
  print "mincresample b1/b1.mnc -like $T1W_reg_anat b1/b1-hig-res.mnc unless -e b1/b1-hig-res.mnc\n";
  `mincresample b1/b1.mnc -like $T1W_reg_anat b1/b1-hig-res.mnc ` unless -e "b1/b1-hig-res.mnc";

  #compute MTsat images

  $MTsat_im_base="MTsat_images-hr";
  $MTsat_im="MTsat_images_MTsat.mnc";
  $MTsat_im_masked="MTsat_images_MTsat_masked.mnc";


  #change b1 map name if nec. note currently not optional (input a mask of all one if you don't have one). b1/b1.mnc is spline smoothed b1
 
  print "\nMTsat-T1_JC.pl $MTW_reg_anat $PDW_reg_anat $T1W_reg_anat b1/b1-hig-res.mnc $MTsat_im_base\n\n";
  `MTsat-T1_JC.pl $MTW_reg_anat $PDW_reg_anat $T1W_reg_anat b1/b1-hig-res.mnc $MTsat_im_base`;

  # now resample the hig res MTsat iamges to lower res to match diffusion
  print "mincresample -tfm_input_sampling MTsat_images-hr_MTsat.mnc -like $anat_to_diff_mnc $MTsat_im\n";
  `mincresample -tfm_input_sampling MTsat_images-hr_MTsat.mnc -like $anat_to_diff_mnc $MTsat_im`; 

  #DO: possibly clamp to positive values? prob not nec with mask (below); and can make mask mask those out if nec

  #for final analysis,  make a mask with a few more slices cut off to deal with misregistration:

  $new_count=$new_zsize-4;

  print "mincreshape -start 7,0,0 -count $new_count,$ysize,$xsize $MTW mask_cut_off_edge_slices.mnc\n";
  `mincreshape -start 7,0,0 -count $new_count,$ysize,$xsize $MTW mask_cut_off_edge_slices.mnc`;


  $MTsat_good_slices="MTsat_good_slices.mnc";

  print "mincresample -tfm_input_sampling -transformation $xfm_MTW_to_b0 -like $anat_to_diff_mnc mask_cut_off_edge_slices.mnc mask_cut_off_edge_slices_like_b0.mnc\n";
  `mincresample -tfm_input_sampling -transformation $xfm_MTW_to_b0 -like $anat_to_diff_mnc mask_cut_off_edge_slices.mnc mask_cut_off_edge_slices_like_b0.mnc`;
  print "mincmath -gt -const 0.0001 mask_cut_off_edge_slices_like_b0.mnc $MTsat_good_slices\n";
  `mincmath -gt -const 0.0001 mask_cut_off_edge_slices_like_b0.mnc $MTsat_good_slices`;

  `rm -f mask_cut_off_edge_slices_like_b0.mnc mask_cut_off_edge_slices.mnc`;

  $brainmask_MTsat="brainmask_MTsat.mnc";

  print "mincmath -float -mult $MTsat_good_slices my_hifi_b0_brain_mask_mnc.mnc $brainmask_MTsat -nocheck_dimensions\n";
  `mincmath -float -mult $MTsat_good_slices my_hifi_b0_brain_mask_mnc.mnc $brainmask_MTsat -nocheck_dimensions`;

  #for better vis: the OK voxels:
  print "mincmath -float -mult $MTsat_im $brainmask_MTsat $MTsat_im_masked\n";
  `mincmath -float -mult $MTsat_im $brainmask_MTsat $MTsat_im_masked`;

  #get rid of unmasked because it's confusing:

  `mv -f $MTsat_im_masked $MTsat_im`;



  print "\n\n------Compute g ratio-----\n\n";


  #compute g-ratio:
  $outbase="g_MTsat";
  $g_name=$outbase."_g.mnc";


  print "NODDI_g $noddi_ficvf_mnc $noddi_fiso_mnc $MTsat_im $outbase $scale_MTsat $brainmask_MTsat\n\n";
  `NODDI_g $noddi_ficvf_mnc $noddi_fiso_mnc $MTsat_im $outbase $scale_MTsat $brainmask_MTsat`;
  print "minc_modify_header -sinsert :history=$program @ARGV $g_name\n";
  `minc_modify_header -sinsert :history="$program @ARGV" $g_name`;

if ($nrx==0){
    print "mincmath -mult WM.mnc brainmask_MTsat.mnc WM_MTsatmasked.mnc\n";
  `mincmath -mult WM.mnc brainmask_MTsat.mnc WM_MTsatmasked.mnc` unless -e "WM_MTsatmasked.mnc";
  print "mincmath -mult parenchyma_bin.mnc brainmask_MTsat.mnc parenchyma_MTsatmasked.mnc\n";
  `mincmath -mult parenchyma_bin.mnc brainmask_MTsat.mnc parenchyma_MTsatmasked.mnc` unless -e "parenchyma_MTsatmasked.mnc";
  ## More stringent masks
  print "mincmath -mult wm_100pc_mask.mnc brainmask_MTsat.mnc WM-100pc_MTsatmasked.mnc\n";
  `mincmath -mult wm_100pc_mask.mnc brainmask_MTsat.mnc WM-100pc_MTsatmasked.mnc`;
  print "mincmath -mult parenchyma_100pc_mask.mnc brainmask_MTsat.mnc parenchyma-100pc_MTsatmasked.mnc\n";
  `mincmath -mult parenchyma_100pc_mask.mnc brainmask_MTsat.mnc parenchyma-100pc_MTsatmasked.mnc`;
}


}#doing_MTsat




#HERE
#`rm -f $tmpdir/*`;




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

sub run_matlab12b {
  my ($command)=@_;
  #system("(echo \"$command\"| matlab -nosplash $logging)");
  open MATLAB,"|matlab12b -nosplash -nojvm -nodisplay " or die "Can't start matlab!\n";
  #print MATLAB "addpath('$mydir')\n";
  #print MATLAB "addpath(genpath('$niak_dir'))\n" if $niak_dir;
  print MATLAB $command,"\n";
  close MATLAB;
  die "Error in MATLAB!\n" if $?;
}

sub rm_ext {
    my $in=$_[0];

my @parts= split /\./,$in;

pop @parts;
my $file_no_ext = join '.', @parts;
return $file_no_ext;
}


# Instructions for processing g-ratio base images


# 1. get all of the scripts needed and read through all of this first.

# /data/mril/mril12/jcampbel/data1/source/g-ratio_pipeline/g-ratio_pipeline.pl

# other files needed:

# there are links to most necessary files in /data/mril/mril12/jcampbel/data1/source/g-ratio_pipeline. I'm sorry if I'm missing some.


# NODDI toolbox is available online (or use my version: /data/mril/mril12/jcampbel/data1/source/NODDI_toolbox2)
# qMTLab is available on gitHub (or use my version: /data/mril/mril12/jcampbel/data1/source/g-ratio_pipeline/)

# you also need:
# minc tools
# access to BIC system (for e.g. mincbeast)

# 2. configuration: setup your paths and run: source /ipl/quarantine/experimental/2013-08-14/minc-toolkit-config.csh; or have it in your .cshrc

# 3. prepare the files: this isn’t in the script yet...  You need to create both nii and mnc files:

# subject has a dicom dir $dcm_sub
# mkdir $nii_dir

# dcm2nii -o $nii_dir $dcm_sub/*

# mv $dcm_sub/*nii.gz  $dcm_sub/*bvec* $dcm_sub/*bval* $nii_dir

# also convert the images to minc:

# dcm2mnc $dcm_sub .

# if doing qMT **also** use: dcm2mnc -splitdynamic $dcm_sub .  (need this for the B0 map too)

# minc files must be unzipped


# 4. run the g-ratio_pipeline.pl script:

# The script will process the files, but assumes a standard acquisition. See the presets in g-ratio_pipeline.pl for exactly what is assumed about the dataset (diffusion encoding, B1 mapping protocol, qMT protocol).  For the diffusion acquisition, see FSL documentation on how to make the required files if you didn't do this standard acquisition.


# Here is a sample input command line:

# g-ratio_pipeline.pl -MTsat mb11052015_20150511_144100_12_mri.mnc mb11052015_20150511_144100_13_mri.mnc mb11052015_20150511_144100_14_mri.mnc -diff *b0AP* *b0PA* *30AP* *30PA* *64AP* *64PA* -anat 20150511_144100T1MPRAGEADNIs003a1001.nii mb11052015_20150511_144100_3_mri.mnc -qMT mb11052015_20150511_144100_5 -B0 mb11052015_20150511_144100_11 -B1 mb11052015_20150511_144100_6_mri.mnc mb11052015_20150511_144100_7_mri.mnc -T1 mb11052015_20150511_144100_8_mri.mnc mb11052015_20150511_144100_9_mri.mnc



# NOTES



# b) outputs currently go in working directory; should clean this up with more output directories.

# c) note constants for calculation of g-ratio are ever in flux: they are hardcoded in g-ratio_pipeline.pl

# d) there is an output file called brainmask_MTsat.mnc/brainmask_qMT.mnc that is the only slices any of the MTsat and g-ratio outputs make sense for.  if things were really misregistered, should remove even more.  The creation of this file is hardcoded for an axial acquisition as is.  MT_good_slices is unbrainmasked version.

# e) need WM / NAWM for some investigations:  A blurred (for visualization) and not blurred (for calculations) version is output.  Currently Vlad's pipeline is run: standard_pipeline.pl

# f) to do lesion segmentation (not in this script currently):

# -Ilana used Nicola's
# -original SL used Erin's
# -Sridar's group has one. give data to him.  check with Simon on whether the tissue classifier is the best;

# g) note that NODDI must be run separately when prompted.  I tried taking out the progress bar to get it to run without the jvm, but didn't succeed.  Could try again...

# h) note you can calculate the tensor separately with dtifit and view the images in fslview as a general check.  I recommend weighted linear least squares.

# example:

# /usr/share/fsl/5.0/bin/dtifit --data=/data/mril/mril13/jcampbel/g-ratio-data/normals/mb11052015/diff_eddy-corr.nii.gz --out=/data/mril/mril13/jcampbel/g-ratio-data/normals/mb11052015/dti --mask=/data/mril/mril13/jcampbel/g-ratio-data/normals/mb11052015/my_hifi_b0_brain_mask.nii.gz --bvecs=/data/mril/mril12/jcampbel/data1/source/g-ratio_pipeline/bvecs --bvals=/data/mril/mril12/jcampbel/data1/source/g-ratio_pipeline/bvals --wls --save_tensor

# i) I've been experiencing some weirdness with minc output files.  For example, some image viewers make the image values appear to be quantized into very few bins, when they really are not.  This is a common observation with register, but it seems to depend on the processing that went into creating the file.  Maybe we need to edit the image min and max.  You will probably notice this in the output files; it isn't fixed yet.
