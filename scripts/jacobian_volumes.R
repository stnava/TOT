
library(ANTsR)
jfn<-Sys.glob("processing/*/*bian.nii.gz")
amyg<-antsImageRead("template/PriorsX/amygdala.nii.gz",3)
labs<-sort(unique( amyg[ amyg > 0 ] ))
nsubs<-length(jfn)
vols<-data.frame(vol1=rep(0,nsubs),vol2=rep(0,nsubs),vol3=rep(0,nsubs),vol4=rep(0,nsubs) )
ids<-c()
for ( x in jfn )
  {
  ids<-c(ids,substr(x,12,24))
  }
wh1<-( amyg == 1 )
wh2<-( amyg == 2 )
wh3<-( amyg == 3 )
wh4<-( amyg == 4 )
ct<-1
for ( x in jfn )
  {
  print(x)
  jimg<-antsImageRead(x,3)
  v1<-sum( jimg[ wh1 ] )
  v2<-sum( jimg[ wh2 ] )
  v3<-sum( jimg[ wh3 ] )
  v4<-sum( jimg[ wh4 ] )
  vols$vol1[ct]<-v1
  vols$vol2[ct]<-v2
  vols$vol3[ct]<-v3
  vols$vol4[ct]<-v4
  ct<-ct+1
  }
vols<-cbind(ids,vols)
write.csv(vols,"hipp_amyg_vols.csv",row.names=F)
