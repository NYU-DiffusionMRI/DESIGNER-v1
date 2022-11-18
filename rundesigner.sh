
echo $BASH_VERSION

subjs=( M9806 )

root=/mnt/labspace/Projects/MESO_v2.0/ALLSUBJS_2.0
for i in ${subjs[@]}; do
    (
    cd $root/$i

    meso=(`ls | grep DIFF_meso_1.nii`)
    research=(`ls | grep DIFF_meso_2.nii`)
    pa=(`ls | grep DIFF_meso_pa.nii`)

    if [ ${#research[@]} -eq 0 ]; then
        continue
    fi
    echo $meso
    echo $research
    echo $pa
    
    /cbi05data/data1/Hamster/DESIGNER/bin/designer \
    -rpe_pair $pa -eddy -pe_dir AP \
    -nocleanup \
    -scratch $root/$i/designer_processing_test_pe \
    $meso,$research $root/$i/designer_out_test_pe
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
