#  registration
rootdir=$PWD
templatedir=${rootdir}/template/
template=${templatedir}local_17_template.nii.gz
templatepriors=${templatedir}PriorsX/priors
templatebmask=${templatedir}local_17_template_brain_mask.nii.gz
subjectimage=$1
outdir=$2
md=$3
t1=$4
malfdir=${rootdir}/malf/
if [[ ! -s $malfdir ]] ; then echo no malfdir ;  exit ; fi
atits=50
# if [[ ${#t1} -lt 3 ]] ; then echo no t1 ;  exit ; fi;
# if [[ ${#md} -lt 3 ]] ; then echo no md ;  exit ; fi;
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
if [[ ! -s $outdir ]] ; then 
  echo create $outdir 
  mkdir -p $outdir
fi
if [[ ! -s ${nm}1InverseWarp.nii.gz ]] ; then 
  antsRegistrationSyNQuick.sh -f $template -m $subjectimage -j 1 -o ${nm}
fi
if [[ ! -s ${nm}jacobian.nii.gz ]] ; then 
  CreateJacobianDeterminantImage 3 ${nm}1Warp.nii.gz ${nm}jacobian.nii.gz 0 1
fi
#  map priors 
fwdmap=" -t ${nm}1Warp.nii.gz -t ${nm}0GenericAffine.mat -r $template "
invmap=" -t [${nm}0GenericAffine.mat,1] -t ${nm}1InverseWarp.nii.gz -r $subjectimage "
for  x in 1 2 3 4 5 6 ; do 
  antsApplyTransforms -d 3 -i ${templatepriors}${x}.nii.gz -o ${nm}_priors${x}.nii.gz $invmap 
done
# initial brain mask --- doesnt have to be perfect
antsApplyTransforms -d 3 -i ${templatebmask} -o ${nm}_brainmaskt.nii.gz $invmap -n NearestNeighbor
####### get a better brain mask via malf ########
if [[ ${#malfdir} -gt 3 ]] && [[ ! -s ${nm}_fusionMalfLabels.nii.gz ]] ; then 
  cts=" 01 02 03 04 05 06 07 08 09 10 " # 11 12 13 14 15 16 17 18 19 20 " # can cause memory problems
  echo "begin malf by building malf label list from training data $cts"
  segcmd=""
  for ct in $cts  ; do 
    segcmd=" $segcmd -g ${malfdir}/${ct}_T2.nii.gz -l  ${malfdir}/${ct}_seg.nii.gz "
  done 
  ImageMath 3 ${nm}_temp.nii.gz MD ${nm}_brainmaskt.nii.gz 1
  MultiplyImages 3 $subjectimage ${nm}_temp.nii.gz ${nm}_brain.nii.gz 
  antsMalfLabeling.sh -k 0 -c 0 -q 1 -j $ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS \
    -d 3 \
    -o ${nm}_fusion \
    -t ${nm}_brain.nii.gz \
    $segcmd
  if [[ -s ${nm}_fusionMalfLabels.nii.gz ]] ; then 
    ImageMath 3 ${nm}_temp.nii.gz ME ${nm}_brainmaskt.nii.gz 5
    ThresholdImage 3 ${nm}_fusionMalfLabels.nii.gz ${nm}_brainmask.nii.gz 1 Inf
    ImageMath 3 ${nm}_brainmask.nii.gz GetLargestComponent ${nm}_brainmask.nii.gz
    ImageMath 3 ${nm}_brainmask.nii.gz addtozero ${nm}_brainmask.nii.gz ${nm}_temp.nii.gz
    MultiplyImages 3 $subjectimage ${nm}_brainmask.nii.gz ${nm}_brain.nii.gz 
    rm ${nm}_temp.nii.gz
    echo MALF Done!
  else 
    echo MALF failed
  fi
fi

ImageMath 3 ${nm}_norm.nii.gz TruncateImageIntensity $subjectimage 0.05 0.995 256
N3BiasFieldCorrection 3 ${nm}_norm.nii.gz ${nm}_norm.nii.gz 8
N3BiasFieldCorrection 3 ${nm}_norm.nii.gz ${nm}_norm.nii.gz 4
ImageMath 3 ${nm}_norm.nii.gz Normalize ${nm}_norm.nii.gz
if  [[ 1 == 0 ]] ; then # old brain masking 
  ThresholdImage 3 ${nm}_norm.nii.gz ${nm}_brainmask.nii.gz 0.3 1
  MultiplyImages 3 ${nm}_brainmask.nii.gz ${nm}_brainmaskt.nii.gz ${nm}_brainmask.nii.gz
  ImageMath 3 ${nm}_brainmaskt.nii.gz ME ${nm}_brainmaskt.nii.gz 10
  ImageMath 3 ${nm}_brainmask.nii.gz + ${nm}_brainmask.nii.gz ${nm}_brainmaskt.nii.gz
  ThresholdImage 3 ${nm}_brainmask.nii.gz ${nm}_brainmask.nii.gz 0.25 Inf
  ImageMath 3 ${nm}_brainmask.nii.gz FillHoles ${nm}_brainmask.nii.gz 
  ImageMath 3 ${nm}_brainmask.nii.gz ME ${nm}_brainmask.nii.gz 2 
  ImageMath 3 ${nm}_brainmask.nii.gz GetLargestComponent ${nm}_brainmask.nii.gz 2 
  ImageMath 3 ${nm}_brainmask.nii.gz MD ${nm}_brainmask.nii.gz 2 
fi

# segmentation steps below
ImageMath 3 ${nm}_norm.nii.gz Normalize ${nm}_norm.nii.gz
antsLaplacianBoundaryCondition.R --output ${nm}_laplacian.nii.gz --mask  ${nm}_brainmask.nii.gz  --input ${nm}_norm.nii.gz 
t1w=${nm}_t1_normWarped.nii.gz
echo $t1 T1
if [[ -s $t1 ]] && [[ ${#t1} -gt 3 ]] && [[ ! -s ${nm}_t1_prob1.nii.gz  ]] ; then 
  antsRegistrationSyNQuick.sh -f ${nm}_norm.nii.gz -m $t1 -j 0 -o ${nm}_t1_norm -t r
  N3BiasFieldCorrection 3 $t1w  $t1w  4
  antsLaplacianBoundaryCondition.R --output ${nm}_laplacian2.nii.gz --mask  ${nm}_brainmask.nii.gz  --input $t1w
  Atropos  -d 3 -x ${nm}_brainmask.nii.gz  -i PriorProbabilityImages[6,${nm}_priors%d.nii.gz,0.25]  -c [50,0] -o [${nm}_t1_seg.nii.gz,${nm}_t1_prob%0d.nii.gz] -m [0.1,1x1x1] -a $t1w -a ${nm}_norm.nii.gz
  MultiplyImages 3 ${nm}_t1_prob1.nii.gz 0.1 ${nm}_t1_prob1.nii.gz
fi
if [[ -s  ${nm}_t1_prob1.nii.gz ]]  ; then 
  ImageMath 3 ${nm}_norm.nii.gz + ${nm}_norm.nii.gz ${nm}_t1_prob1.nii.gz  # increase csf intensity 
fi
mdw=${nm}_md_normWarped.nii.gz
echo $md MD
if [[ -s $md ]] && [[ ${#md} -gt 3 ]]  && [[ ${#t1} -lt 3 ]] ; then 
  MultiplyImages 3 ${nm}_norm.nii.gz ${nm}_brainmask.nii.gz ${nm}_norm.nii.gz
  antsRegistrationSyNQuick.sh -f ${nm}_norm.nii.gz -m $md -j 0 -o ${nm}_md_norm -t sr
  Atropos  -d 3 -x ${nm}_brainmask.nii.gz  -i kmeans[3] -a $mdw -c [1,0] -o [${nm}_md_seg.nii.gz,${nm}_md_prob%0d.nii.gz]
  MultiplyImages 3 ${nm}_md_prob3.nii.gz 0.15 ${nm}_md_prob3.nii.gz
  ImageMath 3 ${nm}_norm.nii.gz + ${nm}_norm.nii.gz ${nm}_md_prob3.nii.gz  # increase csf intensity 
fi
if [[ -s ${nm}_laplacian.nii.gz ]] &&  [[ -s ${nm}_brainmask.nii.gz ]] &&  [[ -s ${nm}_norm.nii.gz ]] ; then
  echo produced ${nm}_norm.nii.gz ${nm}_brainmask.nii.gz ${nm}_laplacian.nii.gz
else 
  echo should have produced ${nm}_norm.nii.gz ${nm}_brainmask.nii.gz ${nm}_laplacian.nii.gz 
  exit
fi
if [[ ! -s ${nm}LapSegmentation.nii.gz ]] ; then 
#  ImageMath 3 ${nm}_norm.nii.gz PeronaMalik ${nm}_norm.nii.gz 2 0.5 
  antsAtroposN4.sh -d 3 -m 1 -n $atits -x ${nm}_brainmask.nii.gz -c 6 -p ${nm}_priors%d.nii.gz -w 0.25 -o ${nm}Lap         -r "[0.1,1x1x1]" -a ${nm}_norm.nii.gz  -a ${nm}_laplacian.nii.gz 
fi 
echo "Finished segmentation"
DIRECT=KellyKapowski
DIRECT_CONVERGENCE="[45,0.0,10]"
DIRECT_THICKNESS_PRIOR="10"
DIRECT_GRAD_STEP_SIZE="0.02"
DIRECT_SMOOTHING_SIGMA="1.5"
DIRECT_NUMBER_OF_DIFF_COMPOSITIONS="10"
dimension=3
exe_direct="$DIRECT -d $dimension -s [${nm}LapSegmentation.nii.gz,2,3]"
tissueprobs=" -g ${nm}LapSegmentationPosteriors2.nii.gz -w ${nm}LapSegmentationPosteriors3.nii.gz "
exe_direct="${exe_direct} $tissueprobs -o ${nm}Thickness.nii.gz"
exe_direct="${exe_direct} -c ${DIRECT_CONVERGENCE} -r ${DIRECT_GRAD_STEP_SIZE}"
exe_direct="${exe_direct} -m ${DIRECT_SMOOTHING_SIGMA} -n ${DIRECT_NUMBER_OF_DIFF_COMPOSITIONS}"
if [[ ! -s ${nm}Thickness.nii.gz ]] ; then 
  $exe_direct
fi
antsApplyTransforms -d 3 -i ${nm}Thickness.nii.gz -o ${nm}_ThicknessToTemplate.nii.gz $fwdmap 
echo "Finished Thickness"
# get csv files 
ImageMath 3 ${nm}_TissueSegmentation.csv LabelStats ${nm}LapSegmentation.nii.gz ${nm}_brainmask.nii.gz
ThresholdImage 3 ${nm}LapSegmentation.nii.gz ${nm}_fusionMalfLabelsBin.nii.gz  1 Inf
ImageMath 3 ${nm}_Malf.csv LabelStats ${nm}_fusionMalfLabels.nii.gz ${nm}_fusionMalfLabelsBin.nii.gz
ThresholdImage 3 ${nm}LapSegmentation.nii.gz ${nm}_fusionMalfLabelsGMBin.nii.gz  2 2 
MultiplyImages 3 ${nm}_fusionMalfLabelsGMBin.nii.gz ${nm}_fusionMalfLabels.nii.gz ${nm}_fusionMalfLabelsGM.nii.gz
ImageMath 3 ${nm}_MalfGM.csv LabelStats ${nm}_fusionMalfLabelsGM.nii.gz ${nm}_fusionMalfLabelsGMBin.nii.gz
ImageMath 3 ${nm}_MalfThickness.csv LabelStats ${nm}_fusionMalfLabelsGM.nii.gz ${nm}Thickness.nii.gz
echo "Done!"
