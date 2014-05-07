#  registration
templatedir=/Users/stnava/data/TOT/template/
template=${templatedir}local_17_template.nii.gz
templatepriors=${templatedir}PriorsX/priors
templatebmask=${templatedir}local_17_template_brain_mask.nii.gz
subjectimage=$1
outdir=$2
md=$3
t1=$4
if [[ ${#t1} -lt 3 ]] ; then echo no t1 ;  exit ; fi;
if [[ ${#md} -lt 3 ]] ; then echo no md ;  exit ; fi;
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
  antsRegistrationSyNQuick.sh -f $template -m $subjectimage -o ${nm}
fi
if [[ ! -s ${nm}jacobian.nii.gz ]] ; then 
  CreateJacobianDeterminantImage 3 ${nm}1Warp.nii.gz ${nm}jacobian.nii.gz 0 1
fi
#  map priors 
invmap=" -t [${nm}0GenericAffine.mat,1] -t ${nm}1InverseWarp.nii.gz -r $subjectimage "
# if [[ ! -s ${nm}_priors6.nii.gz ]] || [[ ! -s ${nm}_brainmask.nii.gz ]] ; then 
  for  x in 1 2 3 4 5 6 ; do 
    antsApplyTransforms -d 3 -i ${templatepriors}${x}.nii.gz -o ${nm}_priors${x}.nii.gz $invmap 
  done
  antsApplyTransforms -d 3 -i ${templatebmask} -o ${nm}_brainmaskt.nii.gz $invmap -n NearestNeighbor
  ImageMath 3 ${nm}_brainmaskt.nii.gz ME ${nm}_brainmaskt.nii.gz 2
# fi
atits=30
# segmentation 
ImageMath 3 ${nm}_norm.nii.gz TruncateImageIntensity $subjectimage 0.05 0.995 256
N3BiasFieldCorrection 3 ${nm}_norm.nii.gz ${nm}_norm.nii.gz 8
N3BiasFieldCorrection 3 ${nm}_norm.nii.gz ${nm}_norm.nii.gz 4
ImageMath 3 ${nm}_norm.nii.gz Normalize ${nm}_norm.nii.gz
ThresholdImage 3 ${nm}_norm.nii.gz ${nm}_brainmask.nii.gz 0.3 1
MultiplyImages 3 ${nm}_brainmask.nii.gz ${nm}_brainmaskt.nii.gz ${nm}_brainmask.nii.gz
ImageMath 3 ${nm}_brainmask.nii.gz + ${nm}_brainmask.nii.gz ${nm}_brainmaskt.nii.gz
ThresholdImage 3 ${nm}_brainmask.nii.gz ${nm}_brainmask.nii.gz 0.25 Inf
ImageMath 3 ${nm}_brainmask.nii.gz FillHoles ${nm}_brainmask.nii.gz 
ImageMath 3 ${nm}_brainmask.nii.gz ME ${nm}_brainmask.nii.gz 2 
ImageMath 3 ${nm}_brainmask.nii.gz GetLargestComponent ${nm}_brainmask.nii.gz 2 
ImageMath 3 ${nm}_brainmask.nii.gz MD ${nm}_brainmask.nii.gz 2 
antsLaplacianBoundaryCondition.R --output ${nm}_laplacian.nii.gz --mask  ${nm}_brainmask.nii.gz  --input ${nm}_norm.nii.gz 
t1w=${nm}_t1_normWarped.nii.gz
echo $t1 T1
if [[ -s $t1 ]] && [[ ${#t1} -gt 3 ]]  ; then 
  antsRegistrationSyNQuick.sh -f ${nm}_norm.nii.gz -m $t1 -o ${nm}_t1_norm -t r
#  Atropos  -d 3 -x ${nm}_brainmask.nii.gz  -i kmeans[3] -a $t1w -c [1,0] -o [${nm}_t1_seg.nii.gz,${nm}_t1_prob%0d.nii.gz]
#  MultiplyImages 3 ${nm}_t1_prob3.nii.gz 0.01 ${nm}_temp.nii.gz
#  MultiplyImages 3 ${nm}_t1_prob3.nii.gz 0.25 ${nm}_t1_prob3.nii.gz
#  MultiplyImages 3 ${nm}_priors3.nii.gz ${nm}_temp.nii.gz ${nm}_priors3.nii.gz  # reduce wm prob here
#  ImageMath 3 ${nm}_norm.nii.gz + ${nm}_norm.nii.gz ${nm}_t1_prob3.nii.gz  # increase csf intensity 
fi
mdw=${nm}_md_normWarped.nii.gz
echo $md MD
if [[ -s $md ]] && [[ ${#md} -gt 3 ]]  && [[ 1 == 0 ]] ; then 
  MultiplyImages 3 ${nm}_norm.nii.gz ${nm}_brainmask.nii.gz ${nm}_norm.nii.gz
  antsRegistrationSyNQuick.sh -f ${nm}_norm.nii.gz -m $md -o ${nm}_md_norm -t sr
  Atropos  -d 3 -x ${nm}_brainmask.nii.gz  -i kmeans[3] -a $mdw -c [1,0] -o [${nm}_md_seg.nii.gz,${nm}_md_prob%0d.nii.gz]
  MultiplyImages 3 ${nm}_md_prob3.nii.gz 0.01 ${nm}_temp.nii.gz
  MultiplyImages 3 ${nm}_md_prob3.nii.gz 0.25 ${nm}_md_prob3.nii.gz
  MultiplyImages 3 ${nm}_priors3.nii.gz ${nm}_temp.nii.gz ${nm}_priors3.nii.gz  # reduce wm prob here
  ImageMath 3 ${nm}_norm.nii.gz + ${nm}_norm.nii.gz ${nm}_md_prob3.nii.gz  # increase csf intensity 
fi
if [[ -s ${nm}_laplacian.nii.gz ]] &&  [[ -s ${nm}_brainmask.nii.gz ]] &&  [[ -s ${nm}_norm.nii.gz ]] ; then
  echo produced ${nm}_norm.nii.gz ${nm}_brainmask.nii.gz ${nm}_laplacian.nii.gz
else 
  echo should have produced ${nm}_norm.nii.gz ${nm}_brainmask.nii.gz ${nm}_laplacian.nii.gz 
  exit
fi
# if [[ ! -s ${nm}LapSegmentation.nii.gz ]] && [[ ! -s $mdw  ]] ; then 
# if [[ ! -s ${nm}LapSegmentation.nii.gz ]] ; then 
  antsAtroposN4.sh -d 3 -m 1 -n $atits -x ${nm}_brainmask.nii.gz -c 6 -p ${nm}_priors%d.nii.gz -w 0.25 -o ${nm}Lap -a ${nm}_norm.nii.gz -r "[0.1,1x1x1]" -a ${nm}_laplacian.nii.gz -a $t1w 
# fi 
echo "Done!"
exit 

# find regions where gm is mislabeled as csf
ThresholdImage 3 ${nm}LapSegmentation.nii.gz ${nm}LapSegmentationCSF.nii.gz 1 1
ThresholdImage 3 ${nm}LapSegmentation.nii.gz ${nm}LapSegmentationGM.nii.gz  2 2 
ThresholdImage 3 ${nm}Segmentation.nii.gz    ${nm}SegmentationGM.nii.gz  2 2 
ImageMath 3 ${nm}temp.nii.gz + ${nm}SegmentationGM.nii.gz ${nm}LapSegmentationCSF.nii.gz
ThresholdImage 3 ${nm}temp.nii.gz  ${nm}temp.nii.gz 2 2 


exit
for x in 1 2 ; do 
  let labnum=${x}+1
  ThresholdImage 3  ${nm}Segmentation.nii.gz ${nm}_cw_mask${x}.nii.gz $labnum $labnum
#  ImageMath 3 ${nm}_cw_mask${x}.nii.gz GetLargestComponent ${nm}_cw_mask${x}.nii.gz
  LabelClustersUniquely 3 ${nm}_cw_mask${x}.nii.gz ${nm}_cw_mask${x}.nii.gz 2500  # get rid of small islands of disconnected WM
  ThresholdImage 3 ${nm}_cw_mask${x}.nii.gz ${nm}_cw_mask${x}.nii.gz 1 Inf
  cp ${nm}SegmentationPosteriors${labnum}.nii.gz   ${nm}_cwpriors${x}.nii.gz
  MultiplyImages 3 ${nm}_cw_mask${x}.nii.gz ${nm}_cwpriors${x}.nii.gz ${nm}_cwpriors${x}.nii.gz
done
ImageMath 3 ${nm}_cw_mask.nii.gz + ${nm}_cw_mask1.nii.gz ${nm}_cw_mask2.nii.gz
#if [[ ! -s ${nm}_2SegmentationPosteriors1.nii.gz ]] ; then
  antsAtroposN4.sh -d 3 -m 1 -n $atits -a ${nm}_norm.nii.gz  -x ${nm}_cw_mask.nii.gz -c 2 -p ${nm}_cwpriors%d.nii.gz -w 0.25 -o ${nm}_2  -a ${nm}_laplacian.nii.gz -r "[0.1,1x1x1]"
#fi
for x in 1 2 ; do 
  let labnum=${x}+1
  cp ${nm}_2SegmentationPosteriors${x}.nii.gz ${nm}_priors${labnum}.nii.gz
  SmoothImage 3 ${nm}_priors${labnum}.nii.gz 0. ${nm}_priors${labnum}.nii.gz
  ImageMath 3 ${nm}_priors${labnum}.nii.gz Normalize ${nm}_priors${labnum}.nii.gz
done
antsAtroposN4.sh -d 3 -m 1 -n $atits -a ${nm}_norm.nii.gz -x ${nm}_brainmask.nii.gz -c 6 -p ${nm}_priors%d.nii.gz -w 0.25 -o ${nm}_3 -a ${nm}_laplacian.nii.gz -r "[0.1,1x1x1]"
# if [[ -s $t1image ]] ; then 
#  antsAtroposN4.sh -d 3 -m 1 -n 6 -a $subjectimage -a $t1image -x ${nm}_brainmask.nii.gz -c 6 -p ${nm}_priors%d.nii.gz -w 0.1 -o ${nm}_2
# fi
#  DiReCT - not done yet
#  KellyKapowski -d 3 ...
#
#
