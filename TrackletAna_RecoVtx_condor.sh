#! /bin/bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}
# source /opt/sphenix/core/bin/sphenix_setup.sh -n mdc2.7
source /opt/sphenix/core/bin/sphenix_setup.sh -n new

export MYINSTALL=/sphenix/u/hjheng/install
export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$MYINSTALL/include:$ROOT_INCLUDE_PATH

source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

# print the environment - needed for debugging
# printenv
./Tracklet_RecoVtx $1 $2

echo all done
