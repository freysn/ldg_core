# provide the location of ldg_data as argument
# it can  be downloaded from here: https://www.dropbox.com/s/cp4rsq8cyjfl7qs/ldg_data.zip?dl=0
ldg_data=$1

# just runs a specified number of iterations (default low number for testing, does not yield refined grid)
niter=16

#
# STOCK1K
#

feat=$ldg_data/stock1k.config

out=/tmp/ldg_core_test/stock1k
mkdir -p $out

cmd_base="build/ldg_core -d $feat --outDir $out/ --distFuncType 1 --repAggregationType 4"

echo "cmd_base: $cmd_base"

$cmd_base  --termIterations $niter --noImgOut

$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,level:0,level:1,level:2"
$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,disparity:0.67"


#
# MCMC1K
#

feat=$ldg_data/mcmc1k.config

out=/tmp/ldg_core_test/mcmc1k
mkdir -p $out

cmd_base="build/ldg_core -d $feat --outDir $out/ --distFuncType 1 --repAggregationType 3"

echo "cmd_base: $cmd_base"

$cmd_base  --termIterations $niter --noImgOut

$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,level:0,level:1,level:2"
$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,disparity:0.67"


#
# CALTECH1K
#

feat=$ldg_data/caltech1k_feat.config
img=$ldg_data/caltech1k_img.config

out=/tmp/ldg_core_test/caltech1k
mkdir -p $out

cmd_base="build/ldg_core -d $feat --outDir $out/ --distFuncType 2 --repAggregationType 2"

echo "cmd_base: $cmd_base"

$cmd_base  --termIterations $niter --noImgOut

$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,level:0,level:1,level:2" -r $img
$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,disparity:0.67" -r $img
