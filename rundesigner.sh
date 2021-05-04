
echo $BASH_VERSION

subjs=( M3808 M3809 )

root=/mnt/labspace/Projects/MESO_v2.0/ALLSUBJS_2.0
for i in ${subjs[@]}; do
    (
    cd $root/$i

    meso=(`ls | grep DIFF_meso.nii`)
    research=(`ls | grep DIFF_meso_research.nii`)
    pa=(`ls | grep DIFF_meso_PA.nii`)

    if [ ${#research[@]} -eq 0 ]; then
        continue
    fi
    echo $meso
    echo $research
    echo $pa
    
    /cbi05data/data1/Hamster/DESIGNER/designer/DESIGNER.py \
    -denoise \

    -nocleanup \
    -tempdir $root/$i/designer_processing_test \
    $meso,$research $root/$i/designer_out_test
    ) #&
done
wait

#     -denoise \
#     -degibbs \
#     -rician \
#     -eddy -rpe_pair $pa -pe_dir -j \
#     -smooth 1.2 \
#     -akc \
#     -mask \
#     -nocleanup \
#     -DTIparams -DKIparams -WMTIparams \
#     -tempdir $out/$i/processing \
