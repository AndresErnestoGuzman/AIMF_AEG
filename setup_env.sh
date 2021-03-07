unset ALMAIMF_ROOTDIR
export ALMAIMF_ROOTDIR=$PWD"/AIMF_AEG/reduction"

unset EXCLUDE_7M
export EXCLUDE_7M=True

unset DO_BSENS
export DO_BSENS=False

unset DO_BSENS_ONLY
export DO_BSENS_ONLY=False

unset FIELD_ID
#export FIELD_ID=""

printenv | \egrep "ALMAIMF_ROOTDIR|EXCLUDE_7M|DO_BSENS|FIELD_ID"
