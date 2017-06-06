#!/bin/bash

VERSION="2.0"

PROGRAM_DEPENDENCIES=( 'dwidenoise' 'unring' 'flirt' 'eddy' 'matlab')
for D in ${PROGRAM_DEPENDENCIES[@]};
  do
    if [[ ! -x $(which $D) ]];
      then
        echo "Error: cannot find $D. Please install dependency or add it to your PATH"
        exit
      fi
  done

function Usage {
    cat <<USAGE

`basename $0` performs full dwi processing with the following steps

 1. pre-check: concatenate all dwi series and make sure the all diffusion AP images and PA image have the same matrix dimensions/orientations/...
 2. generate the (distorted) brain mask
 3. extract a dwi dataset with b<=2500
 4. extract a dwi dataset with b<=1000
 5. denoise the complete dataset
 6. denoise the b<=1000 and keep sigma
 7. gibbs ringing correction on complete dataset
 8. split the gibbs corrected dataset back into original series subsets
 9.  perform n4 bias correction on each dwi subset
 10. motion correction between subsets
 11. re-concatenate the dwi subsets
 12. rician bias correction using sigma from b<=1000
 13. topup + eddy, rotate output bvecs
 14. irwlls outlier map, cwlls dki fit
 15. outlier detection and removal

Usage:

`basename $0` -d <diffusion dataset(s)>
              -o <output prefix>
              <OPTARGS>

Example:

  bash $0 -d DIFF_meso.nii DIFF_meso_research.nii -o designer_params -f SAG3DMPR.nii -k

Required arguments: 

     -d:  diffusion data	           Basename of each diffusion series input to designer. Takes multiple arguemnts.
					    the program will expect one dwi image with extention .nii or .nii.gz and corresponding 
                                            bvec/bvals with the same filename prefix and extensions of .bval and .bvec
     -o:  Outpuit directory		   All outputs will be in directory specified here

Optional arguments:

     -s <sigma>: smoothing		   Include a (csf-free) smoothing step during dwi preprocessing. <sigma> is usually 1.20 x voxel size / 2.355
     -b: bias correction                   Include a bias correction step in dwi preprocessing 0/(1)
     -e <PE img> <PE dir>: topup + eddy    Perform topup and eddy
     -t <PE img> <PE dir>                  Run topup but not eddy
     -i                                    Run eddy but not topup               
     -x <x,y,z>: extent                    Denoising extent. Default is 5,5,5 
     -j:                                   Omit denoising step (for now, this also omits rician correction step) (0)/1
     -g:                                   Omit gibbs correction step (0)/1
     -r:                                   Omit rician correction step (0)/1
     -k: keep temp files                   Keep all temproary files: useful for debugging. (0)/1
USAGE
    exit 1
}

echoParameters() {
    cat <<PARAMETERS

    Using designer with the following parameters
      diffusion image(s)     = ${DIFFUSION_IMAGES[@]}
      smooth sigma           = ${SMOOTH_SIG}
      output prefix          = ${OUTPUT_PREFIX}

PARAMETERS
}

function logCmd() {
  cmd="$*"
  echo "BEGIN >>>>>>>>>>>>>>>>>>>>"
  echo $cmd
  $cmd

  cmdExit=$?

  if [[ $cmdExit -gt 0 ]];
    then
      echo "ERROR: command exited with nonzero status $cmdExit"
      echo "Command: $cmd"
      echo
    fi

  echo "END   <<<<<<<<<<<<<<<<<<<<"
  echo
  echo

  return $cmdExit
}

################################################################################
#
# Main routine
#
################################################################################

HOSTNAME=`hostname`
DATE=`date`

CURRENT_DIR=`pwd`/
OUTPUT_DIR=${CURRENT_DIR}/tmp$RANDOM/
OUTPUT_PREFIX=${OUTPUT_DIR}/tmp
OUTPUT_SUFFIX=".nii.gz"
KEEP_TMP_IMAGES=0
BIASCORRECT=1
PREFIX=designer_params
JHU_REG=0
FREESURFER_EXEC=0
SMOOTH_SIG=0
NODENOISE=0
NOGIBBS=0
NORIC=0
DNEXTENT=5,5,5
EDDY=0
EDDY_ONLY=0
TOPUP_ONLY=0

if [[ $# -lt 3 ]] ; then
  Usage >&2
  exit 1
else
  args="$(sed -r 's/(-[A-Za-z]+ )([^-]*)( |$)/\1"\2"\3/g' <<< $@)"
  declare -a a="($args)"
  set - "${a[@]}"

  while getopts "d:o:s:b:e:k:w:f:h:j:g:r:x:t:i:" OPT
    do
      case $OPT in
          d) # diffusion images
       DIFFUSION_IMAGES=($OPTARG)
       ;;
          o) # output prefix
       PREFIX=$OPTARG
       ;;
          s) # smoothing
       SMOOTH_SIG=$OPTARG
       ;;
          h) # help
       Usage >&2
       exit 0
       ;;
          e) # eddy
       EDDY=($OPTARG)
       ;;
          k) # keep tmp images
       KEEP_TMP_IMAGES=$OPTARG
       ;;
          b) # bias correction
       BIASCORRECT=$OPTARG
       ;;
          w) # wm registration
       JHU_REG=$OPTARG
       ;;
          f) # freesurfer
       FREESURFER_EXEC=$OPTARG
       ;;
          x) # denoising extent
       DNEXTENT=$OPTARG
       ;;
 	        j) # no denoising
	     NODENOISE=$OPTARG
       ;;
 	        g) # no denoising
	     NOGIBBS=$OPTARG
       ;;
 	        r) # no denoising
	     NORIC=$OPTARG
       ;;
 	       t) # no denoising
	     TOPUP_ONLY=$OPTARG
       ;;
 	      i) # no denoising
	     EDDY_ONLY=$OPTARG
       ;;
          *) # getopts issues an error message
       echo "ERROR:  unrecognized option -$OPT $OPTARG"
       exit 1
       ;;
      esac
  done
fi

################################################################################
#
# Preliminaries:
#  1. Check existence of inputs
#  2. Figure out output directory and mkdir if necessary
#
################################################################################
N=${#DIFFUSION_IMAGES[@]}
for d in $(seq 0 $(($N - 1))); do
  if [[ ! -f ${DIFFUSION_IMAGES[$d]} ]];
    then
      echo "The specified image \"${DIFFUSION_IMAGES[$d]}\" does not exist."
      exit 1
    fi
  done

ROOT_DIR=$(dirname "/${DIFFUSION_IMAGES}")
OUTPUT_DIR=${CURRENT_DIR}${ROOT_DIR}/tmp$RANDOM/
OUTPUT_PREFIX=${OUTPUT_DIR}/tmp
OUTPUT_DIR=${OUTPUT_PREFIX%\/*}
PARAM_PREFIX=${CURRENT_DIR}${ROOT_DIR}${PREFIX}
if [[ ! -d $OUTPUT_DIR ]];
  then
    echo "The output directory \"$OUTPUT_DIR\" does not exist. Making it."
    mkdir -p $OUTPUT_DIR
  fi
if [[ ! -d $PARAM_PREFIX ]];
  then
    echo "The output directory \"$PARAM_PREFIX\" does not exist. Making it."
    mkdir -p $PARAM_PREFIX
  fi

# get series list and merge
for b in ${DIFFUSION_IMAGES[@]}; do
  IFS=. read dif ext <<< $b
  slist+=($dif)
done
BVALS=( "${slist[@]/%/.bval}" )
BVECS=( "${slist[@]/%/.bvec}" )

for n in ${DIFFUSION_IMAGES[@]}; do
  mrconvert ${CURRENT_DIR}${ROOT_DIR}$n -force -quiet -stride -1,2,3,4 ${OUTPUT_PREFIX}$n
done

fslmerge -t ${OUTPUT_PREFIX}dwi_orig${OUTPUT_SUFFIX} "${slist[@]/#/${OUTPUT_PREFIX}}"
paste -d" " ${BVALS[@]} | tr -d '\r\n' | awk '{print $0}' > ${OUTPUT_PREFIX}dwi_orig.bval
paste -d" " ${BVECS[@]} | tr -d '\r' | awk '{print $0}' > ${OUTPUT_PREFIX}dwi_orig.bvec

# generate the brain mask
dwi2mask -quiet -force -fslgrad ${OUTPUT_PREFIX}dwi_orig.bvec ${OUTPUT_PREFIX}dwi_orig.bval ${OUTPUT_PREFIX}dwi_orig${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}brainmask${OUTPUT_SUFFIX}

# if bvals go too high: extract bvals <= 3000
#echo "	...extracting images for bval <= 3000 for DKI fitting"
BVAL=( `cat ${OUTPUT_PREFIX}dwi_orig.bval | tr ' ' '\n' | sort -nu | tr '\n' ' '`)
for n in ${BVAL[@]}; do
  #if [[ "$n" -le "30000" ]]; then
    A+=$(echo $n,)
  #fi
done
A=(${A::-1})
dwiextract ${OUTPUT_PREFIX}dwi_orig${OUTPUT_SUFFIX} - -quiet -fslgrad ${OUTPUT_PREFIX}dwi_orig.bvec ${OUTPUT_PREFIX}dwi_orig.bval -shell ${A[@]} | mrconvert - ${OUTPUT_PREFIX}dwi${OUTPUT_SUFFIX} -export_grad_fsl ${OUTPUT_PREFIX}dwi.bvec ${OUTPUT_PREFIX}dwi.bval -force -quiet
unset A
# generate low b image for Rician bias correction
echo "	...extracting images for bval <= 1000 for Rician correction"
for n in ${BVAL[@]}; do
  if [[ "$n" -le "2000" ]]; then
    A+=$(echo $n,)
  fi
done
A=(${A::-1})
dwiextract ${OUTPUT_PREFIX}dwi_orig${OUTPUT_SUFFIX} - -quiet -fslgrad ${OUTPUT_PREFIX}dwi_orig.bvec ${OUTPUT_PREFIX}dwi_orig.bval -shell ${A[@]} | mrconvert - ${OUTPUT_PREFIX}dwi_lowb${OUTPUT_SUFFIX} -export_grad_fsl ${OUTPUT_PREFIX}dwi_lowb.bvec ${OUTPUT_PREFIX}dwi_lowb.bval -force -quiet
unset A
unset slist

if [ $NODENOISE -eq 1 ]; then
  cp ${OUTPUT_PREFIX}dwi${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_denoised${OUTPUT_SUFFIX}
else
  # do the denoising - generate a lowb noisemap 
  echo "	...denoising for b < 1000"
  dwidenoise -force -quiet -noise ${OUTPUT_PREFIX}noisemap_lowb${OUTPUT_SUFFIX} -extent $DNEXTENT ${OUTPUT_PREFIX}dwi_lowb${OUTPUT_SUFFIX}   ${OUTPUT_PREFIX}tmp${OUTPUT_SUFFIX}
  echo "	...denoising full dataset"
  dwidenoise -force -quiet -noise ${OUTPUT_PREFIX}noisemap${OUTPUT_SUFFIX} -extent $DNEXTENT ${OUTPUT_PREFIX}dwi${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_denoised${OUTPUT_SUFFIX}
fi

if [ $NOGIBBS -eq 1 ]; then
  cp ${OUTPUT_PREFIX}dwi_denoised${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_gibbs${OUTPUT_SUFFIX}
else
  # do the degibbsing
  unring ${OUTPUT_PREFIX}dwi_denoised${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_gibbs${OUTPUT_SUFFIX} -minW 1 -maxW 3
fi

# now that the denoising and degibbsing is done we can account for differences between diffusion series
if [ $BIASCORRECT -eq 1 ]; then 
  if [ $N -eq 1 ]; then 
    dwibiascorrect -force -fsl -fslgrad ${OUTPUT_PREFIX}dwi.bvec ${OUTPUT_PREFIX}dwi.bval -mask ${OUTPUT_PREFIX}brainmask${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_gibbs${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_gibbs${OUTPUT_SUFFIX}
  else
    echo "	...finding and correcting artifacts between series"
    ndwis=(0)
    for b in ${!BVALS[@]}; do
      bs=( `cat ${CURRENT_DIR}${ROOT_DIR}${BVALS[$b]} | tr -d '\r\n'` )
      ndwis+=(${#bs[@]})
      Sndwis=(${ndwis[@]})
      unset Sndwis[${#Sndwis[@]}-1]
      tot=0
      for n in ${Sndwis[@]}; do
	       let tot+=$n
      done	
      fslroi ${OUTPUT_PREFIX}dwi_gibbs${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi${b}_gibbs_tmp${OUTPUT_SUFFIX} $tot ${ndwis[-1]}
      echo "		...correcting B1 bias using N4 in FSL"		
      dwibiascorrect -force -fsl -fslgrad ${CURRENT_DIR}${ROOT_DIR}${BVECS[$b]} ${CURRENT_DIR}${ROOT_DIR}${BVALS[$b]} -mask ${OUTPUT_PREFIX}brainmask${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi${b}_gibbs_tmp${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi${b}_gibbsN4_tmp${OUTPUT_SUFFIX}
      fslroi ${OUTPUT_PREFIX}dwi${b}_gibbsN4_tmp${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi${b}_b0_tmp${OUTPUT_SUFFIX} 0 1
      if [ $b -ne 0 ]; then
        MI=( `flirt -in ${OUTPUT_PREFIX}dwi${b}_b0_tmp${OUTPUT_SUFFIX} -ref ${OUTPUT_PREFIX}dwi0_b0_tmp${OUTPUT_SUFFIX} -schedule $FSLDIR/etc/flirtsch/measurecost1.sch -cost mutualinfo | head -1 | cut -f1 -d' '` )
	echo "		...Mutual Information between b0 images is $MI"
	flirt -in ${OUTPUT_PREFIX}dwi${b}_b0_tmp${OUTPUT_SUFFIX} -ref ${OUTPUT_PREFIX}dwi0_b0_tmp${OUTPUT_SUFFIX} -omat ${OUTPUT_PREFIX}b02b0_${b}_tmp.mat -dof 12
	flirt -in ${OUTPUT_PREFIX}dwi${b}_gibbsN4_tmp${OUTPUT_SUFFIX} -ref ${OUTPUT_PREFIX}dwi0_b0_tmp${OUTPUT_SUFFIX} -applyxfm -init ${OUTPUT_PREFIX}b02b0_${b}_tmp.mat -out ${OUTPUT_PREFIX}dwi${b}todwi0_tmp${OUTPUT_SUFFIX} -interp spline
      fi
    done
    dnmerge=( `ls ${OUTPUT_PREFIX}*todwi0_tmp${OUTPUT_SUFFIX}` )
    fslmerge -t ${OUTPUT_PREFIX}dwi_gibbs${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi0_gibbsN4_tmp${OUTPUT_SUFFIX} ${dnmerge[@]}
  fi
  if [ $NORIC -eq 1 ]; then
    cp ${OUTPUT_PREFIX}dwi_gibbs${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_ric${OUTPUT_SUFFIX}
  else
    echo "	...rician correction = sqrt(|Gibbs_corrected^2 - noisemap_lowb^2|)"
    mrcalc ${OUTPUT_PREFIX}dwi_gibbs${OUTPUT_SUFFIX} 2 -pow ${OUTPUT_PREFIX}noisemap_lowb${OUTPUT_SUFFIX} 2 -pow -sub -abs -sqrt ${OUTPUT_PREFIX}dwi_ric${OUTPUT_SUFFIX} -force -quiet
    fslmaths ${OUTPUT_PREFIX}dwi_ric${OUTPUT_SUFFIX} -nan ${OUTPUT_PREFIX}dwi_ric${OUTPUT_SUFFIX}
  fi
  unset slist
  unset ndwis
else
  if [ $N -gt 1 ]; then
    echo "	...finding and correcting artifacts between series"
    ndwis=(0)
    for b in ${!BVALS[@]}; do
      bs=( `cat ${CURRENT_DIR}${ROOT_DIR}${BVALS[$b]}` )
      ndwis+=(${#bs[@]})
      Sndwis=(${ndwis[@]})
      unset Sndwis[${#Sndwis[@]}-1]
      tot=0
      for n in ${Sndwis[@]}; do
	let tot+=$n
      done	
      fslroi ${OUTPUT_PREFIX}dwi_gibbs${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi${b}_gibbs_tmp${OUTPUT_SUFFIX} $tot ${ndwis[-1]}
      echo "		...correcting B1 bias using N4 in FSL"		
      fslroi ${OUTPUT_PREFIX}dwi${b}_gibbs_tmp${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi${b}_b0_tmp${OUTPUT_SUFFIX} 0 1
      if [ $b -ne 0 ]; then
        MI=( `flirt -in ${OUTPUT_PREFIX}dwi${b}_b0_tmp${OUTPUT_SUFFIX} -ref ${OUTPUT_PREFIX}dwi0_b0_tmp${OUTPUT_SUFFIX} -schedule $FSLDIR/etc/flirtsch/measurecost1.sch -cost mutualinfo | head -1 | cut -f1 -d' '` )
	echo "		...Mutual Information between b0 images is $MI"
	flirt -in ${OUTPUT_PREFIX}dwi${b}_b0_tmp${OUTPUT_SUFFIX} -ref ${OUTPUT_PREFIX}dwi0_b0_tmp${OUTPUT_SUFFIX} -omat ${OUTPUT_PREFIX}b02b0_${b}_tmp.mat -dof 12
	flirt -in ${OUTPUT_PREFIX}dwi${b}_gibbs_tmp${OUTPUT_SUFFIX} -ref ${OUTPUT_PREFIX}dwi0_b0_tmp${OUTPUT_SUFFIX} -applyxfm -init ${OUTPUT_PREFIX}b02b0_${b}_tmp.mat -out ${OUTPUT_PREFIX}dwi${b}todwi0_tmp${OUTPUT_SUFFIX} -interp spline
      fi
    done
    dnmerge=( `ls ${OUTPUT_PREFIX}*todwi0_tmp${OUTPUT_SUFFIX}` )
    fslmerge -t ${OUTPUT_PREFIX}dwi_gibbs${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi0_gibbs_tmp${OUTPUT_SUFFIX} ${dnmerge[@]}
  fi
  if [ $NORIC -eq 1 ]; then
    cp ${OUTPUT_PREFIX}dwi_gibbs${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_ric${OUTPUT_SUFFIX}
  else
    echo "	...rician correction = sqrt(|Gibbs_corrected^2 - noisemap_lowb^2|)"
    mrcalc ${OUTPUT_PREFIX}dwi_gibbs${OUTPUT_SUFFIX} 2 -pow ${OUTPUT_PREFIX}noisemap_lowb${OUTPUT_SUFFIX} 2 -pow -sub -abs -sqrt ${OUTPUT_PREFIX}dwi_ric${OUTPUT_SUFFIX} -force -quiet
    fslmaths ${OUTPUT_PREFIX}dwi_ric${OUTPUT_SUFFIX} -nan ${OUTPUT_PREFIX}dwi_ric${OUTPUT_SUFFIX}
  fi
  unset slist
  unset ndwis
fi

fslroi ${OUTPUT_PREFIX}dwi_orig${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}b0_orig${OUTPUT_SUFFIX} 0 1

if [ ${#EDDY[@]} -eq 1 ]; then
  echo "	...skipping topup+eddy"	
  cp ${OUTPUT_PREFIX}dwi_ric${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_tuec${OUTPUT_SUFFIX}  
  cp ${OUTPUT_PREFIX}dwi.bvec ${OUTPUT_PREFIX}dwi_tuec.bvec
  cp ${OUTPUT_PREFIX}dwi.bval ${OUTPUT_PREFIX}dwi_tuec.bval
else
  fslroi ${CURRENT_DIR}${ROOT_DIR}${EDDY[0]} ${OUTPUT_PREFIX}rpeb0_orig${OUTPUT_SUFFIX} 0 1
  dwipreproc -force -rpe_pair ${OUTPUT_PREFIX}b0_orig${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}rpeb0_orig${OUTPUT_SUFFIX} -fslgrad ${OUTPUT_PREFIX}dwi.bvec ${OUTPUT_PREFIX}dwi.bval -export_grad_fsl ${OUTPUT_PREFIX}dwi_tuec.bvec ${OUTPUT_PREFIX}dwi_tuec.bval ${EDDY[1]} ${OUTPUT_PREFIX}dwi_ric${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_tuec${OUTPUT_SUFFIX}
  #dwipreproc -force -rpe_all ${CURRENT_DIR}${ROOT_DIR}RL_all.nii.gz -fslgrad ${OUTPUT_PREFIX}dwi.bvec ${OUTPUT_PREFIX}dwi.bval -export_grad_fsl ${OUTPUT_PREFIX}dwi_tuec.bvec ${OUTPUT_PREFIX}dwi_tuec.bval ${EDDY[1]} ${OUTPUT_PREFIX}dwi_ric${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_tuec${OUTPUT_SUFFIX}
fi

if [ ${EDDY_ONLY} -eq 1 ]; then
   dwipreproc -force -rpe_none -fslgrad ${OUTPUT_PREFIX}dwi.bvec ${OUTPUT_PREFIX}dwi.bval -export_grad_fsl ${OUTPUT_PREFIX}dwi_tuec.bvec  ${OUTPUT_PREFIX}dwi_tuec.bval AP ${OUTPUT_PREFIX}dwi_ric${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_tuec${OUTPUT_SUFFIX} 

fi

if [ ${TOPUP_ONLY} -eq 1 ]; then
    echo "TOPUP only not yet implemented"
fi

fslroi ${OUTPUT_PREFIX}dwi_tuec${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}b0_tuec${OUTPUT_SUFFIX} 0 1
dwi2mask -quiet -force -fslgrad ${OUTPUT_PREFIX}dwi_tuec.bvec ${OUTPUT_PREFIX}dwi_tuec.bval ${OUTPUT_PREFIX}dwi_tuec${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}brainmask${OUTPUT_SUFFIX}

if (( $(echo "${SMOOTH_SIG} > 0" | bc -l) )); then
  echo "        ...smoothing non-csf voxels"
  fslmaths ${OUTPUT_PREFIX}b0_tuec${OUTPUT_SUFFIX} -mas ${OUTPUT_PREFIX}brainmask${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}b0brain${OUTPUT_SUFFIX}
  fast -n 3 -t 2 -o ${OUTPUT_PREFIX}_tissue ${OUTPUT_PREFIX}b0brain${OUTPUT_SUFFIX}
  fslmaths ${OUTPUT_PREFIX}_tissue_pve_0${OUTPUT_SUFFIX} -thr 0.7 -bin ${OUTPUT_PREFIX}CSF_mask${OUTPUT_SUFFIX}
  fslmaths ${OUTPUT_PREFIX}CSF_mask${OUTPUT_SUFFIX} -mul -1 -add 1 ${OUTPUT_PREFIX}invCSF_mask${OUTPUT_SUFFIX}
  fslmaths ${OUTPUT_PREFIX}dwi_tuec${OUTPUT_SUFFIX} -s ${SMOOTH_SIG} -mas ${OUTPUT_PREFIX}invCSF_mask${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}smoothr1${OUTPUT_SUFFIX}
  fslmaths ${OUTPUT_PREFIX}invCSF_mask${OUTPUT_SUFFIX} -s ${SMOOTH_SIG} -mas ${OUTPUT_PREFIX}invCSF_mask${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}smoothr2${OUTPUT_SUFFIX}
  fslmaths ${OUTPUT_PREFIX}smoothr1${OUTPUT_SUFFIX} -mul ${OUTPUT_PREFIX}smoothr2${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}smoothr3${OUTPUT_SUFFIX}
  fslmaths ${OUTPUT_PREFIX}dwi_tuec${OUTPUT_SUFFIX} -mas ${OUTPUT_PREFIX}CSF_mask${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_csforig${OUTPUT_SUFFIX}
  fslmaths ${OUTPUT_PREFIX}dwi_csforig${OUTPUT_SUFFIX} -add ${OUTPUT_PREFIX}smoothr3${OUTPUT_SUFFIX} ${OUTPUT_PREFIX}dwi_tuec${OUTPUT_SUFFIX}
fi

gunzip -f ${OUTPUT_PREFIX}dwi_tuec${OUTPUT_SUFFIX}
gunzip -f ${OUTPUT_PREFIX}brainmask${OUTPUT_SUFFIX}

cd /data1/Hamster/Ben/scripts/DESIGNER
matlab -nodesktop -nodisplay -r "tensorfitting(); exit();"
cd $CURRENT_DIR

if [ $KEEP_TMP_IMAGES -eq 0 ]; then
  rm -r ${OUTPUT_DIR}
fi
	


