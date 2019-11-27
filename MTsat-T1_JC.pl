#!/usr/bin/perl -w

### Compute MT saturation maps
#based on Helms et al MRM 2008 and equation erratum in MRM (64) 2010
#B1 correction empirical from Weiskopf Front. Neurosci. 2013; only the MTsat map is corrected
##########################################################
# Written by Ilana Leppert
# Modified by Jennifer Campbell
####################################################
require 5.001;
use Getopt::Tabular;
use MNI::Startup qw(nocputimes);
use MNI::Spawn;
use MNI::FileUtilities qw(check_output_dirs);
use File::Basename;
use List::MoreUtils qw(zip);
use Cwd;
if($0 =~ /[\/A-Za-z_0-9-]+\/([A-Za-z0-9_-]+\.pl)/) {$program = $1;}	#program name without path
$Usage = <<USAGE;

Compute MSat maps from matched MTw, PDw and T1w images

Usage: $program <MTw> <PDw> <T1w> <B1map> <output_base>

By default, will look for FA and TR in headers


USAGE
#-help for options


#@args_table = (
#);

Getopt::Tabular::SetHelp ($Usage, '');

#GetOptions(\@args_table, \@ARGV, \@args) || exit 1;
die $Usage unless $#ARGV >=0;

if($ARGV[0]=~/help/i){print $Usage; exit(1);}

$output_base=pop(@ARGV);
$B1map=pop(@ARGV);
$t1w=pop(@ARGV);
$pdw=pop(@ARGV);
$MTw=pop(@ARGV);

$a_t1w=`mincinfo -attvalue acquisition:flip_angle $t1w`;chomp($a_t1w);
$a_pdw=`mincinfo -attvalue acquisition:flip_angle $pdw`;chomp($a_pdw);
$a_MTw=`mincinfo -attvalue acquisition:flip_angle $MTw`;chomp($a_MTw);

# Convert to radians
$a_t1w_r=$a_t1w*3.14159/180;
$a_pdw_r=$a_pdw*3.14159/180;
$a_MTw_r=$a_MTw*3.14159/180;

print "\n--FA in rads: T1W:$a_t1w_r PDW:$a_pdw_r\n\n";

$TR_t1w=`mincinfo -attvalue acquisition:repetition_time $t1w`;chomp($TR_t1w);
$TR_pdw=`mincinfo -attvalue acquisition:repetition_time $pdw`;chomp($TR_pdw);
$TR_MTw=`mincinfo -attvalue acquisition:repetition_time $MTw`;chomp($TR_MTw);

# R1app=1/2 *(St1*at1/TRt1-Spd*apd/TRpd)/(Spd/apd-St1/at1)
# Aapp=Spd*St1*(TRpd*at1/apd-TRt1*apd/at1)/(St1*TRpd*at1-Spd*TRt1*apd)
$r1app=$output_base."_R1app.mnc";
$Aapp=$output_base."_Aapp.mnc";



`minccalc -float -nocheck_dimensions -expr "clamp(0.5*(A[0]*$a_t1w_r/$TR_t1w-A[1]*$a_pdw_r/$TR_pdw)/(A[1]/$a_pdw_r-A[0]/$a_t1w_r),-10,100)" $t1w $pdw $r1app `;

#`minccalc -float -nocheck_dimensions -expr "clamp(0.5*(A[0]*$a_t1w_r/$TR_t1w-A[1]*$a_pdw_r/$TR_pdw)/(A[1]/$a_pdw_r-A[0]/$a_t1w_r),-10,10)" $t1w $pdw $r1app `;

`minccalc -float -nocheck_dimensions -expr "clamp(A[0]*A[1]*($TR_pdw*$a_t1w_r/$a_pdw_r-$TR_t1w*$a_pdw_r/$a_t1w_r)/(A[1]*$TR_pdw*$a_t1w_r-A[0]*$TR_t1w*$a_pdw_r),-10000,15000)" $pdw $t1w $Aapp `;

# Apparent MT saturation

$sapp=$output_base."_MTsat_uncorr.mnc";


`minccalc -float -nocheck_dimensions -expr "clamp(100*((A[0]*$a_MTw_r/A[1]-1)*A[2]*$TR_MTw - $a_MTw_r*$a_MTw_r/2),0,50)" $Aapp $MTw $r1app $sapp `;
#`minccalc -float -nocheck_dimensions -expr "clamp(100*((A[0]*$a_MTw_r/A[1]-1)*A[2]*$TR_MTw - $a_MTw_r*$a_MTw_r/2),0,5)" $Aapp $MTw $r1app $sapp `;


$MTsatcorr=$output_base."_MTsat.mnc";

#B1 correction:

`minccalc -float -nocheck_dimensions  -expr "clamp(A[0]*(1-0.4)/(1-0.4*A[1]),0,50)" $sapp $B1map $MTsatcorr`;

#`minccalc -float -nocheck_dimensions  -expr "clamp(A[0]*(1-0.4)/(1-0.4*A[1]),0,5)" $sapp $B1map $MTsatcorr`;
