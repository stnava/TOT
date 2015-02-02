library(ANTsR)
library(caret)
library(randomForest)
library(RRF)
setwd("/Users/stnava/data/TOT/")
demog<-read.csv("data/TOT_raw_demog.csv")
qual<-read.csv("data/TOT_BA_Notes.csv")
date2ppdformat<-function( x )
{
temp<-unlist(strsplit(as.character(x), "-"))
return( paste( (temp) ,collapse='') )
}
ofn<-"data/TOT_demog.csv"
demog$DOB<-as.Date(as.character(demog$DOB),format="%m/%d/%Y")
ids<-unique(demog$No)
excinds<-c(grep("Date",colnames(demog)),grep("MRI",colnames(demog)))
demog4wk<-cbind(demog[1,-excinds],Age=NA,MRIDate=NA,
  AnatomicalQualityMeasure=NA,AnatomicalFileExtension=NA,
  DTQualityMeasure=NA,DTFileExtension=NA)
demog1yr<-cbind(demog[1,-excinds],Age=NA,MRIDate=NA,
  AnatomicalQualityMeasure=NA,AnatomicalFileExtension=NA,
  DTQualityMeasure=NA,DTFileExtension=NA)
computedage<-NA
if ( exists("demogreo") ) rm("demogreo")
for ( i in ids[!is.na(ids)] )
  {
  ww<-which(demog$No==i)
  demog4wk[1,1:(ncol(demog4wk)-4)]<-demog[ ww , -excinds]
  demog4wk$MRIDate=as.character(demog[ ww , ]$DateOfMRI1mo[1])
  demog4wk$MRIDate=as.Date(demog4wk$MRIDate,format="%m/%d/%Y")
  ppd1<-date2ppdformat( demog4wk$MRIDate )
  demog4wk$Age=as.numeric(demog4wk$MRIDate-demog4wk$DOB)
  fafn4wk<-Sys.glob(
    paste("data/dti_recon/*/*/*/",i,"*",ppd1,"*fa*nii.gz",sep=''))
  if ( length(fafn4wk) > 0 )
    demog4wk$DTFileExtension=unlist(strsplit(fafn4wk[1],"/"))[[6]]
  else demog4wk$DTFileExtension=NA
  if (  any( !is.na(qual[qual$ids==i,]$meanM1_M2) ) )
    {
    temp<-qual[qual$ids==i,]
    demog4wk$AnatomicalQualityMeasure=temp[which.max(temp$meanM1_M2),]$meanM1_M2
    demog4wk$AnatomicalFileExtension=temp[which.max(temp$meanM1_M2),]$Scan
    }
  else
    {
    demog4wk$AnatomicalQualityMeasure=NA
    demog4wk$AnatomicalFileExtension=NA
    }
  demog1yr[1,1:(ncol(demog4wk)-4)]<-demog[ ww , -excinds]
  demog1yr$MRIDate=as.character(demog[ ww , ]$DateOfMRI1yr[1])
  demog1yr$MRIDate=as.Date(demog1yr$MRIDate,format="%m/%d/%Y")
  demog1yr$Age=as.numeric(demog1yr$MRIDate-demog4wk$DOB)
  ppd2<-date2ppdformat( demog1yr$MRIDate )
  fafn1yr<-Sys.glob(
    paste("data/dti_recon/*/*/*/",i,"*",ppd2,"*fa*nii.gz",sep=''))
  if ( length(fafn1yr) > 0 )
    demog1yr$DTFileExtension=unlist(strsplit(fafn1yr[1],"/"))[[6]]
  else demog1yr$DTFileExtension=NA
  if ( ! exists("demogreo") )
    {
    demogreo<-rbind(demog4wk,demog1yr)
    }
  else
    {
    demogreo<-rbind(demogreo,demog4wk,demog1yr)
    }
  }
write.csv(demogreo,ofn,row.names=F)
print(
summary(lm(  AnatomicalQualityMeasure ~ MomEd + ITNEnrollment,
  data=demogreo, subset=demogreo$AnatomicalQualityMeasure>=2))
)

baddt<-c("7006_20140430_DTIfa.nii.gz","7038_20130904_DTIfa.nii.gz","7041_20140826_DTIfa.nii.gz","","","")
demogreo$DTQualityMeasure[ !is.na(demogreo$DTFileExtension) ]<-1
for ( dtfn in baddt )
  demogreo$DTQualityMeasure[ demogreo$DTFileExtension == dtfn ]<-0

# now we are set up to process dti
if ( TRUE )
  {
  library(ANTsR)
  setwd("/Users/stnava/data/TOT/tempdti")
  wfn<-na.omit(demogreo$DTFileExtension[ demogreo$DTQualityMeasure == 1])
  ref<-antsImageRead("fa_avg.nii.gz",3)
  refm<-antsImageRead("fa_ref_mask.nii.gz",3)
  for ( k in 1:length(wfn) )
    {
    fafn<-Sys.glob(paste("../data/dti_recon/*/*/*/",wfn[k],sep=''))
    onm<-paste("w",wfn[k],sep='')
    if ( ! file.exists(onm) )
      {
      fa<-antsImageRead(fafn,3)
      fam<-getMask(fa,0.02,0.8*mean(fa),2);
      ImageMath(3,fam,'MD',fam,1)
      ImageMath(3,fam,'ME',fam,1)
      ImageMath(3,fam,'FillHoles',fam)
    #  fa[fam==0]<-0
      mytx<-antsRegistration(fixed=ref , moving=fa ,
             typeofTransform = c("SyN"),
             outprefix=tempfile(),mask=refm)
      antsImageWrite( mytx$warpedmovout , onm  )
      }
    }
    if ( FALSE ) {
      ilist<-imageFileNames2ImageList( Sys.glob("w*gz") , 3 )
      rad<-rep(3,ref@dimension)
      pp<-jointIntensityFusion(ref,refm,ilist, rSearch=0, rad=rad )
    }
    sel<-demogreo$DTQualityMeasure == 1 & demogreo$Age < 100 &
      !is.na(demogreo$Age)  & !is.na(demogreo$DTQualityMeasure )
    babylist<-paste("w",na.omit( demogreo$DTFileExtension[sel] ), sep='' )
    ref[  refm == 0  ]<-0
    famask<-antsImageRead('famask_topo_skel.nii.gz',3)
    famask<-antsImageRead('famask.nii.gz',3)
    famask<-antsImageRead('famask_cerebrum.nii.gz',3)
    mat<-imagesToMatrix( babylist , famask )
    subdemog<-demogreo[sel,]
    confinds<-c(5,7,8,38)
    confoundmat<-antsrimpute(data.matrix(subdemog[,confinds]))
    rmat<-residuals(lm(mat~confoundmat))
    predmat<-antsrimpute(data.matrix(subdemog[,c(12)]))
    predmat<-antsrimpute(data.matrix(subdemog[,c(10,11,12)]))
    poorOrNot<-factor( predmat[,3]  < 1 )
    predmat[,3]<-rank(predmat[,3] )
    mat2<-antsrimpute(data.matrix(subdemog[,c(10,12,confinds)]))
    set.seed(9)
    selector<-createDataPartition( poorOrNot )$Resample1
    if ( FALSE )
    {
    matlist<-list(data.matrix( rmat[selector,]),
                  data.matrix(predmat)[selector,])
    sccan<-sparseDecom2( inmatrix=matlist,
      inmask=c(famask,NA), nvecs=2, sparseness=c( 0.05, 0.9 ),
      cthresh=c(5000,0), its=15, mycoption=1, perms=0, smooth=0 )
    avgmat<-abs(imageListToMatrix( sccan$eig1 , famask ))
    avgmat<-avgmat/rowSums(abs(avgmat))
    avgmat<- mat %*% t(avgmat)
#    cor( pp, antsrimpute(mydfc) )
    }
    eanat<-sparseDecom( inmatrix=mat, inmask=famask, nvecs=50,
      sparseness=0.005, cthresh=500, its=5, mycoption=0 )
    jeanat<-joinEigenanatomy( mat , famask, eanat$eig,
      c(1:20)/100.0 , joinMethod='multilevel' )
    useeig<-eanat$eig
    useeig<-jeanat$fusedlist
#    eanat2<-sparseDecom( inmatrix=mat, inmask=famask, nvecs=length(useeig),
#      sparseness=0, cthresh=50000, its=5, mycoption=0,
#      initializationList=useeig )
#    useeig<-eanat2$eig
    eseg<-eigSeg( famask, useeig )
    antsImageWrite( eseg, 'eseg.nii.gz' )
    avgmat<-abs(imageListToMatrix( useeig , famask ))
    avgmat<-avgmat/rowSums(abs(avgmat))
    imgmat<-(  mat %*% t(avgmat)  )
#    imgmat<-svd(mat,nv=0,nu=10)$u
    set.seed(2)
    nreps<-1000
    perfval<-rep(0,nreps)
    perfvalnoimg<-rep(0,nreps)
    perfvalperm<-rep(0,nreps) # permuted
    for ( i in 1:nreps )
      {
      selector<-createDataPartition( subdemog$EGAcalc_cr, p=4.0/5.0 )$Resample1
      mydf<-data.frame( env=factor(poorOrNot) , confoundmat )
      mdl<-randomForest( env ~ . , data=mydf[ selector,], importance = T )
      pp<-predict(mdl,newdata=mydf[-selector,]  )
      xtab <- table(pp, mydf$env[-selector]  )
      cm <- confusionMatrix(xtab)
      perfvalnoimg[i]<-cm$overall[1]
      # now use imaging
      mydf<-data.frame( env=factor(poorOrNot) , img=imgmat, confoundmat )
      mdl<-randomForest( env ~ . , data=mydf[ selector,], importance = T )
      pp<-predict(mdl,newdata=mydf[-selector,]  )
      xtab <- table(pp, mydf$env[-selector]  )
      cm <- confusionMatrix(xtab)
      perfval[i]<-cm$overall[1]
#      print(cm$overall)
      imgmatperm<-imgmat[sample(1:nrow(imgmat)), ]
      mydf<-data.frame( env=factor(poorOrNot) , img=imgmatperm, confoundmat )
      mdl<-randomForest( env ~ . , data=mydf[ selector,], importance = T )
      pp<-predict(mdl,newdata=mydf[-selector,]  )
      xtab <- table(pp, mydf$env[-selector]  )
      cm <- confusionMatrix(xtab)
      perfvalperm[i]<-cm$overall[1]
      }
    fullmdl<-randomForest( env ~ . , data=mydf, importance = T )
    print( t.test( perfval - perfvalperm  ) )
    hist(perfval-perfvalperm)
    print( t.test( perfval - perfvalnoimg  ) )
    hist(perfval-perfvalnoimg)
    print(fullmdl$importance)
    plot(density(perfval,bw=0.1))
######### ok done with prediction #########
    if ( FALSE ) {
      mylm<-(lm( avgmat ~   .  , data=subdemog[,c(12,confinds)]  ))
      mylmres<-bigLMStats(mylm)
      qv<-p.adjust( mylmres$beta.p["ITNEnrollment",] ,method='none' )
      sigct<-1
      esegq<-antsImageClone( eseg )
      for ( k in as.numeric( which(qv <= 0.05 ) ) )
        {
        esegq[ eseg != k & eseg > sigct ]<-0
        esegq[ eseg == k ]<-sigct
        sigct<-sigct+1
        }
      antsImageWrite( esegq, 'esegq.nii.gz' )
      summary( lm( avgmat[,qv<=0.05] ~   .  , data=subdemog[,c(12,confinds)]  ))
    }
  }
