#  registration
templatedir=/Users/stnava/data/TOT/template/
template=${templatedir}local_17_template.nii.gz
templatepriors="${templatedir}PriorsX/priors"
templatebmask=${templatedir}local_17_template_brain_mask.nii.gz
subjectimage=$1
outdir=$2
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
if [[ ! -s ${nm}_priors6.nii.gz ]] || [[ ! -s ${nm}_brainmask.nii.gz ]] ; then 
  for  x in 1 2 3 4 5 6 ; do 
    antsApplyTransforms -d 3 -i ${templatepriors}${x}.nii.gz -o ${nm}_priors${x}.nii.gz $invmap 
  done
  antsApplyTransforms -d 3 -i ${templatebmask} -o ${nm}_brainmask.nii.gz $invmap -n NearestNeighbor
fi
#  segmentation 
if [[ ! -s ${nm}Segmentation.nii.gz ]] ; then 
  antsAtroposN4.sh -d 3 -m 1 -n 6 -a $subjectimage -x ${nm}_brainmask.nii.gz -c 6 -p ${nm}_priors%d.nii.gz -w 0.25 -o ${nm}
fi 
for x in 1 2 ; do 
  let labnum=${x}+1
  ThresholdImage 3  ${nm}Segmentation.nii.gz ${nm}_cw_mask${x}.nii.gz $labnum $labnum
  LabelClustersUniquely 3 ${nm}_cw_mask${x}.nii.gz ${nm}_cw_mask${x}.nii.gz 500  # get rid of small islands of disconnected WM
  ThresholdImage 3 ${nm}_cw_mask${x}.nii.gz ${nm}_cw_mask${x}.nii.gz 1 Inf
  cp ${nm}SegmentationPosteriors${labnum}.nii.gz   ${nm}_cwpriors${x}.nii.gz
  MultiplyImages 3 ${nm}_cw_mask${x}.nii.gz ${nm}_cwpriors${x}.nii.gz ${nm}_cwpriors${x}.nii.gz
done
ImageMath 3 ${nm}_cw_mask.nii.gz + ${nm}_cw_mask1.nii.gz ${nm}_cw_mask2.nii.gz
# if [[ ! -s ${nm}_2SegmentationPosteriors1.nii.gz ]] ; then
  antsAtroposN4.sh -d 3 -m 1 -n 6 -a $subjectimage  -x ${nm}_cw_mask.nii.gz -c 2 -p ${nm}_cwpriors%d.nii.gz -w 0.0 -o ${nm}_2
# fi
for x in 1 2 ; do 
  let labnum=${x}+1
  cp ${nm}_2SegmentationPosteriors${x}.nii.gz ${nm}_priors${labnum}.nii.gz
  SmoothImage 3 ${nm}_priors${labnum}.nii.gz 1 ${nm}_priors${labnum}.nii.gz
  ImageMath 3 ${nm}_priors${labnum}.nii.gz Normalize ${nm}_priors${labnum}.nii.gz
done
antsAtroposN4.sh -d 3 -m 1 -n 6 -a $subjectimage -x ${nm}_brainmask.nii.gz -c 6 -p ${nm}_priors%d.nii.gz -w 0.25 -o ${nm}_3
# if [[ -s $t1image ]] ; then 
#  antsAtroposN4.sh -d 3 -m 1 -n 6 -a $subjectimage -a $t1image -x ${nm}_brainmask.nii.gz -c 6 -p ${nm}_priors%d.nii.gz -w 0.1 -o ${nm}_2
# fi
#  DiReCT - not done yet
#  KellyKapowski -d 3 ...
#
#