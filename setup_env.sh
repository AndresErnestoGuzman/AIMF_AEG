unset ALMAIMF_ROOTDIR
export ALMAIMF_ROOTDIR=$PWD"/AIMF_AEG/reduction"

unset EXCLUDE_7M
export EXCLUDE_7M=True

unset DO_BSENS
export DO_BSENS=True

unset DO_BSENS_ONLY
export DO_BSENS_ONLY=True

unset FIELD_ID
export FIELD_ID="G337.40"

printenv | egrep "ALMAIMF_ROOTDIR|EXCLUDE_7M|DO_BSENS|FIELD_ID"
