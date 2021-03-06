#  registration
rootdir=$PWD
templatedir=${rootdir}/template_TOT0004/
templatepriors=${templatedir}Priors/priorsj
template=${templatedir}template_JLF.nii.gz
templatebmask=${templatedir}template_JLF_brain_mask.nii.gz
##############
subjectimage=$1
outdir=$2
dti=$3
malfdir=${rootdir}/malf/
if [[ ! -s $malfdir ]] ; then echo no malfdir ; fi
atits=150
if [[ ! -s $template ]] ; then
  echo template $template does not exist. exiting.
  exit 1
fi
if [[ ! -s $subjectimage ]] ; then
  echo subject $subjectimage does not exist. exiting.
  exit 1
fi
nm=`basename $subjectimage`
nm2=` echo $nm | cut -d '.' -f 1 `
outdir=${outdir}/${nm2}
nm=${outdir}/${nm2}
#  map priors
fwdmap=" -t ${nm}1Warp.nii.gz -t ${nm}0GenericAffine.mat -r $template "
invmap=" -t [${nm}0GenericAffine.mat,1] -t ${nm}1InverseWarp.nii.gz -r $subjectimage "
if [[ ! -s $outdir ]] ; then
  echo create $outdir
  mkdir -p $outdir
fi
if [[ ! -s ${nm}1InverseWarp.nii.gz ]] ; then
  antsRegistrationSyN.sh -f $template -m $subjectimage -j 0 -o ${nm}
fi
if [[ ! -s ${nm}jacobian.nii.gz ]] ; then
  CreateJacobianDeterminantImage 3 ${nm}1Warp.nii.gz ${nm}jacobian.nii.gz 0 1
fi
for  x in 1 2 3 4 5 6  ; do # csf gm wm dgm bstem cerebellum ventricles
  antsApplyTransforms -d 3 -i ${templatepriors}${x}.nii.gz -o ${nm}_priors${x}.nii.gz $invmap
done
# initial brain mask --- doesnt have to be perfect
antsApplyTransforms -d 3 -i ${templatebmask} -o ${nm}_brainmaskt.nii.gz $invmap -n NearestNeighbor
echo got mask
####### get a better brain mask via malf ########
if [[ ${#malfdir} -gt 3 ]] && [[ ! -s ${nm}_fusionLabels.nii.gz ]] ; then
  cts=" 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 " # can cause memory problems
  echo "begin malf by building malf label list from training data $cts"
  segcmd=""
  for ct in $cts  ; do
    segcmd=" $segcmd -g ${malfdir}/${ct}_T2.nii.gz -l  ${malfdir}/${ct}_seg.nii.gz "
  done
  ImageMath 3 ${nm}_temp.nii.gz MD ${nm}_brainmaskt.nii.gz 1
  MultiplyImages 3 $subjectimage ${nm}_temp.nii.gz ${nm}_brain.nii.gz
  ImageMath 3 ${nm}_brain.nii.gz ClosestSimplifiedHeaderMatrix ${nm}_brain.nii.gz
  antsJointLabelFusion.sh -d 3 -k 0 -c 0  -j 1 -q 1 \
    -o ${nm}_fusion \
    -t ${nm}_brain.nii.gz \
    $segcmd
  echo MALF Done!
fi

if [[ -s ${nm}_fusionLabels.nii.gz ]] ; then
    ImageMath 3 ${nm}_norm.nii.gz TruncateImageIntensity $subjectimage 0.005 0.995 256
    ImageMath 3 ${nm}_norm.nii.gz Normalize ${nm}_norm.nii.gz
    CopyImageHeaderInformation $subjectimage ${nm}_fusionLabels.nii.gz ${nm}_fusionLabels.nii.gz 1 1 1
    ImageMath 3 ${nm}_temp.nii.gz MD ${nm}_brainmaskt.nii.gz 0
    ThresholdImage 3 ${nm}_fusionLabels.nii.gz ${nm}_brainmask.nii.gz 1 Inf
    ImageMath 3 ${nm}_brainmask.nii.gz ME ${nm}_brainmask.nii.gz 8
    ImageMath 3 ${nm}_brainmask.nii.gz GetLargestComponent ${nm}_brainmask.nii.gz
    ImageMath 3 ${nm}_brainmask.nii.gz MD ${nm}_brainmask.nii.gz 9
    ImageMath 3 ${nm}_brainmask.nii.gz addtozero ${nm}_brainmask.nii.gz ${nm}_temp.nii.gz
    MultiplyImages 3 ${nm}_norm.nii.gz  ${nm}_brainmask.nii.gz ${nm}_brain.nii.gz
    N3BiasFieldCorrection 3 ${nm}_brain.nii.gz  ${nm}_brain.nii.gz  4
    ThresholdImage 3 ${nm}_brain.nii.gz  ${nm}_otsu.nii.gz Otsu 3
    ThresholdImage 3 ${nm}_otsu.nii.gz ${nm}_otsu.nii.gz 2 3
    MultiplyImages 3 ${nm}_brainmask.nii.gz ${nm}_otsu.nii.gz ${nm}_brainmask.nii.gz
    MultiplyImages 3 ${nm}_norm.nii.gz  ${nm}_brainmask.nii.gz ${nm}_brain.nii.gz
    ImageMath 3 ${nm}_fusionLabels.csv LabelStats ${nm}_fusionLabels.nii.gz ${nm}_brainmask.nii.gz
else
    echo MALF failed    exit
fi
if [[ ! -s ${nm}_norm.nii.gz ]] ; then
   echo ${nm} missing norm.nii.gz
   exit
fi
# segmentation steps below - this may help cortical segmentation
if [[  -s ${nm}_norm.nii.gz ]] ; then
  antsLaplacianBoundaryCondition.R --output ${nm}_laplacian.nii.gz --mask  ${nm}_brainmask.nii.gz  --input ${nm}_norm.nii.gz
fi

# check outputs
if [[ -s ${nm}_laplacian.nii.gz ]] &&  [[ -s ${nm}_brainmask.nii.gz ]] &&  [[ -s ${nm}_norm.nii.gz ]] ; then
  echo produced ${nm}_norm.nii.gz ${nm}_brainmask.nii.gz ${nm}_laplacian.nii.gz
else
  echo should have produced ${nm}_norm.nii.gz ${nm}_brainmask.nii.gz ${nm}_laplacian.nii.gz
  exit
fi

# map jlf segmentation to probabilities
## ThresholdImage 3 ${nm}_fusionLabels.nii.gz  ${nm}_priors4.nii.gz  40 47
ThresholdImage 3 ${nm}_fusionLabels.nii.gz  ${nm}_priors6.nii.gz  17 18
## ThresholdImage 3 ${nm}_fusionLabels.nii.gz  ${nm}_priors7.nii.gz  49 50
for x in  3 4 5 6 ; do
  SmoothImage 3 ${nm}_priors${x}.nii.gz 0.5 ${nm}_priors${x}.nii.gz
  ImageMath 3 ${nm}_priors${x}.nii.gz Normalize ${nm}_priors${x}.nii.gz
done
mod2=" -a ${nm}_laplacian.nii.gz "
if [[ ${#dti} -gt 4 ]] ; then
  mod2=" $mod2 -a $dti "
fi
antsAtroposN4.sh -d 3 -u 0 -m 2 -n $atits -x ${nm}_brainmask.nii.gz -c 6 \
  -p ${nm}_priors%d.nii.gz -w 0.25 -o ${nm}Lap         \
  -r "[0.0,1x1x1]" -a ${nm}_norm.nii.gz $mod2
#  ThresholdImage 3 ${nm}_fusionLabels.nii.gz  ${nm}_ventmask.nii.gz  49 50
#  ImageMath 3 ${nm}_negventmask.nii.gz Neg ${nm}_ventmask.nii.gz
#  MultiplyImages 3 ${nm}_brainmask.nii.gz   ${nm}_negventmask.nii.gz  ${nm}_segmask.nii.gz
##  if [[ ! -s ${nm}LapSegmentation.nii.gz ]] ; then
#    antsAtroposN4.sh -d 3 -m 1 -n $atits -x ${nm}_segmask.nii.gz -c 6 -p ${nm}_priors%d.nii.gz -w 0.25 -o ${nm}Lap         -r "[0.1,1x1x1]" -a ${nm}_norm.nii.gz   -a ${nm}_laplacian.nii.gz
#  ImageMath 3  ${nm}LapSegmentation.nii.gz + ${nm}LapSegmentation.nii.gz ${nm}_ventmask.nii.gz

# use geometry to infer sulcal locations
if [[ -s ${nm}LapSegmentation.nii.gz ]] ; then
  bash ${rootdir}/thk_example.sh ${nm}
  antsApplyTransforms -d 3 -i ${nm}thicknessjj.nii.gz -o ${nm}_CortVolToTemplate.nii.gz $fwdmap
fi
if [[ -s ${nm}thicknessjj.nii.gz  ]] ; then
for  x in 1 2 3 4 5 ; do # exclude the JLF probabilities
  antsApplyTransforms -d 3 -i ${templatepriors}${x}.nii.gz -o ${nm}_priors${x}.nii.gz $invmap
done
MultiplyImages 3 ${nm}wm.nii.gz ${nm}_priors3.nii.gz  ${nm}_priors3.nii.gz
ImageMath 3 ${nm}thicknessjjsig.nii.gz  SigmoidImage ${nm}thicknessjj.nii.gz 0.01 0.6
ThresholdImage 3 ${nm}thicknessjj.nii.gz ${nm}thicknessjjcsf.nii.gz 0.05 0.5
ImageMath 3 ${nm}wmneg.nii.gz Neg ${nm}wm.nii.gz
MultiplyImages 3 ${nm}thicknessjjcsf.nii.gz ${nm}wmneg.nii.gz ${nm}thicknessjjcsf.nii.gz
SmoothImage 3 ${nm}thicknessjjcsf.nii.gz 0.5 ${nm}thicknessjjcsf.nii.gz
MultiplyImages 3 ${nm}thicknessjjsig.nii.gz  ${nm}_priors2.nii.gz ${nm}_priors2.nii.gz
ImageMath 3  ${nm}thicknessjj_neg.nii.gz  Neg  ${nm}thicknessjjsig.nii.gz
MultiplyImages 3 ${nm}thicknessjj_neg.nii.gz  ${nm}_priors1.nii.gz ${nm}_priors1.nii.gz
ImageMath 3 ${nm}_priors1.nii.gz + ${nm}thicknessjjcsf.nii.gz ${nm}_priors1.nii.gz
MultiplyImages 3 ${nm}thicknessjj_neg.nii.gz  ${nm}_priors3.nii.gz ${nm}_priors3.nii.gz
antsAtroposN4.sh -d 3 -u 0 -m 2 -n $atits -x ${nm}_brainmask.nii.gz -c 6 \
  -p ${nm}_priors%d.nii.gz -w 0.25 -o ${nm}DIFF \
  -r "[0.0,1x1x1]" -a ${nm}_norm.nii.gz
fi


########################################################################
########################################################################
echo "Finished segmentation 0 now thickness"
########################################################################
########################################################################
DIRECT=KellyKapowski
DIRECT_CONVERGENCE="[45,0.0,10]"
DIRECT_THICKNESS_PRIOR="10"
DIRECT_GRAD_STEP_SIZE="0.02"
DIRECT_SMOOTHING_SIGMA="1.5"
DIRECT_NUMBER_OF_DIFF_COMPOSITIONS="10"
dimension=3
exe_direct="$DIRECT -d $dimension -s [${nm}DIFFSegmentation.nii.gz,2,3]"
tissueprobs=" -g ${nm}DIFFSegmentationPosteriors2.nii.gz -w ${nm}DIFFSegmentationPosteriors3.nii.gz "
exe_direct="${exe_direct} $tissueprobs -o ${nm}Thickness.nii.gz"
exe_direct="${exe_direct} -c ${DIRECT_CONVERGENCE} -r ${DIRECT_GRAD_STEP_SIZE}"
exe_direct="${exe_direct} -m ${DIRECT_SMOOTHING_SIGMA} -n ${DIRECT_NUMBER_OF_DIFF_COMPOSITIONS}"
# if [[ ! -s ${nm}Thickness.nii.gz ]] ; then
  $exe_direct
# fi
antsApplyTransforms -d 3 -i ${nm}Thickness.nii.gz -o ${nm}_ThicknessToTemplate.nii.gz $fwdmap

########################################################################
########################################################################
echo " get csv files "
########################################################################
########################################################################
ImageMath 3 ${nm}_TissueSegmentationDiffeo.csv LabelStats ${nm}DIFFSegmentation.nii.gz ${nm}_brainmask.nii.gz
ImageMath 3 ${nm}_TissueSegmentation.csv LabelStats ${nm}DIFFSegmentation.nii.gz ${nm}_brainmask.nii.gz
ThresholdImage 3 ${nm}DIFFSegmentation.nii.gz ${nm}_fusionLabelsBin.nii.gz  1 Inf
ImageMath 3 ${nm}_Malf.csv LabelStats ${nm}_fusionLabels.nii.gz ${nm}_fusionLabelsBin.nii.gz
ThresholdImage 3 ${nm}DIFFSegmentation.nii.gz ${nm}_fusionLabelsGMBin.nii.gz  2 2
MultiplyImages 3 ${nm}_fusionLabelsGMBin.nii.gz ${nm}_fusionLabels.nii.gz ${nm}_fusionLabelsGM.nii.gz
ImageMath 3 ${nm}_MalfGM.csv LabelStats ${nm}_fusionLabelsGM.nii.gz ${nm}_fusionLabelsGMBin.nii.gz
ImageMath 3 ${nm}_MalfThickness.csv LabelStats ${nm}_fusionLabelsGM.nii.gz ${nm}Thickness.nii.gz
echo "Done!"
