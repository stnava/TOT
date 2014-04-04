#  registration
templatedir=/Users/stnava/data/TOT/template/
template=${templatedir}local_17_template.nii.gz
templatepriors="${templatedir}PriorsUCL/priors"
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
#  map priors 
invmap=" -t [${nm}0GenericAffine.mat,1] -t ${nm}1InverseWarp.nii.gz -r $subjectimage "
if [[ ! -s ${nm}_priors6.nii.gz ]] || [[ ! -s ${nm}_brainmask.nii.gz ]] ; then 
  for  x in 1 2 3 4 5 6 ; do 
    antsApplyTransforms -d 3 -i ${templatepriors}${x}.nii.gz -o ${nm}_priors${x}.nii.gz $invmap 
  done
  antsApplyTransforms -d 3 -i ${templatebmask} -o ${nm}_brainmask.nii.gz $invmap -n NearestNeighbor
fi
#  segmentation 
antsAtroposN4.sh -d 3 -m 1 -n 6 -a $subjectimage -x ${nm}_brainmask.nii.gz -c 6 -p ${nm}_priors%d.nii.gz -w 0.25 -o ${nm}
#  DiReCT - not done yet
#  KellyKapowski -d 3 ...
#
