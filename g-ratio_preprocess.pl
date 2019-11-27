#!/usr/bin/perl -w

### Preprocess g-ratio data
#
##########################################################
##########################################################
## NOVA-lepto project
##
###########################################
#
###########################################
# Created by Ilana Leppert April2018
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

Prepares the data for the g-ratio pipeline


Usage: $program
-help for options

USAGE

#my @args_table = (#["-clobber","boolean",undef,\$clobber,"clobber all output files (currently not implemented, sorry)"],#HERE not implemented);

#Getopt::Tabular::SetHelp ($Usage, '');
my @args;
#GetOptions(\@args_table, \@ARGV, \@args) || exit 1;
#die $Usage unless $#ARGV >=0;

print "---Get the list of sequences from the minc files---\n";
`lsminc *.mnc.gz > mnclist` unless -e "mnclist";

print "---Sort the dicom by sequence name and then sequence number-----\n";
$dir=`\\ls -d N*`; chomp($dir);#print "dir:$dir\n";
print "sort_dicom.pl $dir/MR* '0018,1030' $dir/\n";
`sort_dicom.pl $dir/MR* '0018,1030' $dir/`;
`sort_dicom.pl $dir/MR* '0020,0011' $dir/`;


print "\n-----Convert the diffusion to nii----------\n";
`mkdir nii`;
$done=`\\ls nii/*b2500PAN*.nii.gz`; chomp($done); 
unless (-e $done) {
    `dcm2nii -o nii $dir/cmrr_mbep2d_diff_acc6_b300/*` ;
    `dcm2nii -o nii $dir/cmrr_mbep2d_diff_acc6_b300_PA/*`;
    `dcm2nii -o nii $dir/cmrr_mbep2d_diff_acc6_b700/*`;
    `dcm2nii -o nii $dir/cmrr_mbep2d_diff_acc6_b700_PA/*`;
    `dcm2nii -o nii $dir/cmrr_mbep2d_diff_acc6_b2500/*`;
    `dcm2nii -o nii $dir/cmrr_mbep2d_diff_acc6_b2500_PA/*`;
}
print "\n---Convert the T1w-PreGd to nii-----------\n";
@preGd_all=`more mnclist | grep mni_gre_3D_T1w-PreGd`; #this should return 3 scans, the order is: unfiltered, filtered, phase (want 2nd)
@preGd = split(/\s+/,$preGd_all[1]);
@T1w = split(/\s+/,$preGd_all[0]);

if ($preGd[0]=~/_(\d+)d1_/){$preGd_s=$1 } else {die "can't find the series # for pre-Gd\n";}
#print "dcm2nii -o nii $dir/$preGd_s/*\n";
$done=`\\ls nii/201*3DT1wPreGd*.nii.gz`; chomp($done); 
`dcm2nii -o nii $dir/$preGd_s/*` unless -e $done;

print "\n--- Get the necessary inputs for the processing----\n";
@MT_all=`more mnclist | grep mni_gre_3D_MTon_helms`; #this should return 2 scans, the order is: unfiltered, filtered (want 1st)
@MT = split(/\s+/,$MT_all[0]);
@PD_all=`more mnclist | grep mni_gre_3D_MToff_helms`; #this should return 2 scans, the order is: unfiltered, filtered (want 1st)
@PD = split(/\s+/,$PD_all[0]);
@B1_60=`more mnclist | grep ep_seg_se_b1_60`; 
@b1_60 = split(/\s+/,$B1_60[0]);
@B1_120=`more mnclist | grep ep_seg_se_b1_120`; 
@b1_120 = split(/\s+/,$B1_120[0]);

print "------\n";
print "g-ratio_pipeline-repeatlesions.pl -diff nii/*{b300,b700,b2500,b300PA,b700PA,b2500PA}N*.nii.gz -bvecs nii/*{b300,b700,b2500,b300PA,b700PA,b2500PA}N*.bvec -anat nii/201*3DT1wPreGd*.nii.gz  $preGd[0] -MTsat $MT[0] $PD[0] $T1w[0] -B1 $b1_60[0] $b1_120[0]> log-reg-to-t1p\n\n";
print "------\n";