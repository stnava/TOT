
for x in 1 2 3 ; do
  MultiplyImages 3 priors${x}.nii.gz negpriors6.nii.gz  priors${x}.nii.gz
  MultiplyImages 3 priors${x}.nii.gz negpriors5.nii.gz  priors${x}.nii.gz
done

