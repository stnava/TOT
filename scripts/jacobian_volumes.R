##############################################################################################################
library(ANTsR)
ids<-c("7002","7003","7004","7005","7006","7009","7010","7011","7012","7013","7014","7016","7017","7019","7021","7023","7024","7026","7027","7028","7029","7030","7032","7034","7035","7036","7037","7038","7040","7041","7042","7043","7044","7046","7048","7049","7050","7053","7054","7055","7058","7059")
jfn<-c()
for ( id in ids ) {
    lj<-Sys.glob(paste("processing/",id,"*/*bian.nii.gz",sep=''))[1]
    jfn<-lappend(jfn,lj)
}
amyg<-antsImageRead("template/PriorsX/amygdala.nii.gz",3)
tcort<-antsImageRead("template/PriorsX/priors2.nii.gz",3)
tdgm<-antsImageRead("template/PriorsX/priors4.nii.gz",3)
twm<-antsImageRead("template/PriorsX/priors3.nii.gz",3)
pthr<-0.5
loginds<-( twm > pthr | tcort > pthr | tdgm > pthr )
twm[  loginds  ]<-twm[ loginds ]+tdgm[ loginds  ]+tcort[ loginds  ]
bmask<-antsImageClone( twm )
bmask[ twm > pthr ]<-1
bmask[ twm <= pthr ]<-0
labs<-sort(unique( amyg[ amyg > 0 ] ))
nsubs<-length(jfn)
vols<-data.frame(vol1=rep(0,nsubs),vol2=rep(0,nsubs),vol3=rep(0,nsubs),vol4=rep(0,nsubs), vol5=rep(0,nsubs) )
ids<-c()
for ( x in jfn )
  {
  ids<-c(ids,substr(x,12,24))
  }
wh1<-( amyg == 1 )
wh2<-( amyg == 2 )
wh3<-( amyg == 3 )
wh4<-( amyg == 4 )
wh5<-( bmask == 1 )
ct<-1
for ( x in jfn )
  {
  print(x)
  jimg<-antsImageRead(x,3)
  v1<-sum( jimg[ wh1 ] )
  v2<-sum( jimg[ wh2 ] )
  v3<-sum( jimg[ wh3 ] )
  v4<-sum( jimg[ wh4 ] )
  v5<-sum( jimg[ wh5 ] )
  vols$vol1[ct]<-v1
  vols$vol2[ct]<-v2
  vols$vol3[ct]<-v3
  vols$vol4[ct]<-v4
  vols$vol5[ct]<-v5
  ct<-ct+1
  }
vols<-cbind(ids,vols)
print( cor( data.matrix( vols ) ) )
write.csv(vols,"hipp_amyg_vols.csv",row.names=F)

########### now eanat ###########
mat<-imagesToMatrix( jfn, bmask )
lmat<-lowrankRowMatrix( mat, 10 )
eanat<-sparseDecom( lmat, bmask, sparseness=0.1, nvecs=10, its=3, cthresh=250, smooth=1, mycoption=1 )
eanatseg<-eigSeg( bmask, eanat$eigenanatomyimages )
