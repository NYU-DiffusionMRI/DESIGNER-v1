################################################################################
####                   Execute Designer Pipeline                            ####
####                  Written by Sid on 10/22/2018                          ####
####                  Modified by ___ on XX/XX/XXXX                         ####
####    Changelog:                                                          ####
####      1. Initial creation                                               ####
################################################################################
# To run Designer, open terminal and cd to current dir, then type the command:
#./designer.sh <input dir> <output dir>

# Define input arguments
# input directory needs to contain image.nii, bvals and bvecs
# output directory can be any folder of your choice
input=$1
output=$2

# Run designer
python designer.py \
-denoise -extent 5,5,5 \
-degibbs \
-rician \
-mask \
-prealign \
-smooth 1.25 \
-rpe_header -eddy \
-fit_constraints 0,1,0 \
-median \
-DKIparams -DTIparams \
${1} ${2}
