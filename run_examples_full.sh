# provide the location of ldg_data as argument
# it can  be downloaded from here: https://www.dropbox.com/s/cp4rsq8cyjfl7qs/ldg_data.zip?dl=0
ldg_data=$1

# just runs a specified number of iterations (default low number for testing, does not yield refined grid)
niter=16

#
# STOCK
#

feat=$ldg_data/stock.config

out=/tmp/ldg_core_test/stock
mkdir -p $out

cmd_base="build/ldg_core -d $feat --outDir $out/ --distFuncType 1 --repAggregationType 4 --preNormData"

echo "cmd_base: $cmd_base"

$cmd_base  --termIterations $niter --noImgOut

$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,level:4,level:5,level:6"
$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,disparity:0.75"


#
# MCMC
#

feat=$ldg_data/mcmc.config

out=/tmp/ldg_core_test/mcmc
mkdir -p $out

cmd_base="build/ldg_core -d $feat --outDir $out/ --distFuncType 1 --repAggregationType 3 --nTilesAssign 131072"

echo "cmd_base: $cmd_base"

$cmd_base  --termIterations $niter --noImgOut

$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,level:3,level:4,level:5"
$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,disparity:0.75"


#
# CALTECH
#

feat=$ldg_data/caltech_feat.config
img=$ldg_data/caltech_img.config
out=/tmp/ldg_core_test/caltech
mkdir -p $out

cmd_base="build/ldg_core -d $feat --outDir $out/ --distFuncType 2 --repAggregationType 2  --nTilesAssign 10240"

echo "cmd_base: $cmd_base"

$cmd_base  --termIterations $niter --noImgOut

$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,swapXY:1,level:1,level:2" -r $img
$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,swapXY:1,disparity:0.75" -r $img
