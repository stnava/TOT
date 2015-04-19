dim=3
pre=$1
seg=${pre}LapSegmentation.nii.gz
m=${pre}LapSegmentationPosteriors3.nii.gz
gm=${pre}LapSegmentationPosteriors2.nii.gz
f=${pre}sum.nii.gz
ThresholdImage $dim $seg ${pre}cort.nii.gz 2 2
ThresholdImage $dim $seg ${pre}wm.nii.gz 3 4
ImageMath 3 ${pre}wm.nii.gz GetLargestComponent ${pre}wm.nii.gz
ImageMath 3 ${pre}wm.nii.gz MD ${pre}wm.nii.gz 1
MultiplyImages 3 ${pre}wm.nii.gz ${pre}LapSegmentationPosteriors3.nii.gz ${pre}wmmod.nii.gz
ImageMath $dim $f + $gm ${pre}wmmod.nii.gz
its=[100,1.e-5,5]
smth=0vox
down=1
# if [[ ! -s  ${pre}0Warp.nii.gz  ]] ; then
antsRegistration -d $dim  \
                        -m meansquares[ $f,$m,1] \
                        -t TimeVaryingVelocityField[ 1.0, 4, 0.0,0.0, 0.5,0 ] \
                        -c $its  \
                        -s $smth  \
                        -f $down \
                       -u 1 -x ${pre}cort.nii.gz --verbose 1 \
                       -o [${pre},${pre}ex_diff.nii.gz,${pre}ex_diff_inv.nii.gz]
# fi
CreateJacobianDeterminantImage $dim ${pre}0Warp.nii.gz ${pre}jac.nii.gz 0 1
CreateJacobianDeterminantImage $dim ${pre}0InverseWarp.nii.gz ${pre}jac_inv.nii.gz 0 1
antsApplyTransforms -d $dim -i ${pre}jac_inv.nii.gz -o ${pre}thicknessj.nii.gz -t ${pre}0Warp.nii.gz -r ${pre}jac_inv.nii.gz
MultiplyImages $dim ${pre}thicknessj.nii.gz $gm ${pre}thicknessjj.nii.gz
