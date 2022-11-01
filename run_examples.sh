# provide the location of ldg_data as argument
ldg_data=$1

#
# STOCK1K
#

feat=$ldg_data/stock1k.config

out=/tmp/ldg_core_test/stock1k
mkdir -p $out

cmd_base="build/ldg_core -d $feat --outDir $out/ --distFuncType 1 --repAggregationType 4"

echo "cmd_base: $cmd_base"

$cmd_base  --termIterations 8 --noImgOut

$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,level:0,level:1,level:2"


#
# MCMC1K
#

feat=$ldg_data/mcmc1k.config

out=/tmp/ldg_core_test/mcmc1k
mkdir -p $out

cmd_base="build/ldg_core -d $feat --outDir $out/ --distFuncType 1 --repAggregationType 3"

echo "cmd_base: $cmd_base"

$cmd_base  --termIterations 8 --noImgOut

$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,level:0,level:1,level:2"


#
# CALTECH1K
#

feat=$ldg_data/caltech1k_feat.config
img=$ldg_data/caltech1k_img.config

out=/tmp/ldg_core_test/caltech1k
mkdir -p $out

cmd_base="build/ldg_core -d $feat --outDir $out/ --distFuncType 2 --repAggregationType 2"

echo "cmd_base: $cmd_base"

$cmd_base  --termIterations 8 --noImgOut

$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,level:0,level:1,level:2" -r $img


#
# STOCK
#

feat=$ldg_data/stock.config

out=/tmp/ldg_core_test/stock
mkdir -p $out

cmd_base="build/ldg_core -d $feat --outDir $out/ --distFuncType 1 --repAggregationType 4"

echo "cmd_base: $cmd_base"

$cmd_base  --termIterations 8 --noImgOut

$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,level:4,level:5,level:6"


#
# MCMC
#

feat=$ldg_data/mcmc.config

out=/tmp/ldg_core_test/mcmc
mkdir -p $out

cmd_base="build/ldg_core -d $feat --outDir $out/ --distFuncType 1 --repAggregationType 3 --nTilesAssign 131072"

echo "cmd_base: $cmd_base"

$cmd_base  --termIterations 8 --noImgOut

$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,level:3,level:4,level:5"


#
# CALTECH
#

feat=$ldg_data/caltech_feat.config
img=$ldg_data/caltech_img.config
out=/tmp/ldg_core_test/caltech
mkdir -p $out

cmd_base="build/ldg_core -d $feat --outDir $out --distFuncType 2 --repAggregationType 2  --nTilesAssign 10240"

echo "cmd_base: $cmd_base"

$cmd_base  --termIterations 8 --noImgOut

$cmd_base --loadAssignments ${out}/qtLeafAssignment.raw.bz2 --termTime 0 --imgOutOpts "png:0,swapXY:1,level:0,level:1,level:2" -r $img
