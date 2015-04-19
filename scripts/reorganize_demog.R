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
#print(
#summary(lm(  AnatomicalQualityMeasure ~ MomEd + ITNEnrollment,
#  data=demogreo, subset=demogreo$AnatomicalQualityMeasure>=2))
#)

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
  ref<-antsImageRead("fa_avg.nii.gz")
  ref<-antsImageRead("fa_avg_raw.nii.gz")
  refm<-antsImageRead("fa_ref_mask.nii.gz")
  for ( k in 1:length(wfn) )
    {
    fafn<-Sys.glob(paste("../data/dti_recon/*/*/*/",wfn[k],sep=''))
    subid<-unlist( strsplit( fafn[1] , "/" ))
    subid<-paste(subid[4],subid[5],sep="_")
    srch<-paste("../segmentationsMar25_2015/processing3/*/",
      subid,"*rain.nii.gz",sep='')
    t2fn<-Sys.glob(srch)
    onm<-paste("w",wfn[k],sep='')
    onmm<-paste("m",wfn[k],sep='')
    print(fafn)
    print(t2fn)
    if ( length(t2fn) > 0 & FALSE ) # &  ! file.exists(onmm) )
      {
      t2<-antsImageRead( t2fn[1] )
      fa<-antsImageRead(fafn,3)
      fam<-getMask(fa,0.025,0.8*mean(fa),3);
      fam<-iMath(fam,'MD',1)
      fam<-iMath(fam,'ME',1)
      fam<-iMath(fam,'FillHoles')
    #  fa[fam==0]<-0
      mytx<-antsRegistration(fixed=ref , moving=fa ,
             typeofTransform = c("SyN") )
#      antsImageWrite( mytx$warpedmovout , onm  )
      fam <- antsApplyTransforms( fixed=fa, moving=refm,
       transformlist=mytx$invtransforms, interpolator="NearestNeighbor" )
      fam2 <- fa * iMath( fam, "ME", 1 )
      mytx<-antsRegistration( fixed=t2, moving=fam2,
             typeofTransform = "SyNCC") # ,mask=refm)
      antsImageWrite( mytx$warpedmovout , onmm  )
      print( onmm )
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
#    confinds<-c(5,7,38)
    confoundmat<-antsrimpute(data.matrix(subdemog[,confinds]))
    rmat<-residuals(lm(mat~confoundmat))
    predmat<-antsrimpute(data.matrix(subdemog[,c(12)]))
    predmat<-antsrimpute(data.matrix(subdemog[,c(10,11,12)]))
    subdemog$AdmHC_cr<-antsrimpute( subdemog$AdmHC_cr )
    poorOrNot<-factor( predmat[,3]  < 1 )
    tryMerfNerf <- TRUE
    if ( tryMerfNerf )
    {
      selector<-createDataPartition( poorOrNot, p=0.5 )$Resample1
      sublist<-imageFileNames2ImageList( babylist[selector] )
      sublistes<-imageFileNames2ImageList( babylist[-selector] )
      trnlist<-list()
      for ( i in 1:length(sublist) ) trnlist[[i]]<-list(  sublist[[i]] )
      predval<-as.numeric( poorOrNot )
      predval<-psych::winsor( antsrimpute(subdemog$ITNEnrollment), 0.05 )
      temp<-abs( cor(  mat[selector,] ,  predval[selector] ) )
      tmdl<-lm( mat[selector,] ~ predval[selector] + . ,
        data=subdemog[selector,confinds] )
      tmdl<-bigLMStats( tmdl )
      temp<-abs( tmdl$beta.t[1,] )
      statfamask<-makeImage( famask, temp ) %>% smoothImage(2)
      statfamask<-thresholdImage( statfamask , 0.1, Inf )
      rmat<-data.matrix( residuals( lm( mat[selector,] ~ 0 + . ,
        data=subdemog[ selector, confinds ] ) ) )
      eanat<-sparseDecom( rmat, famask, nvecs=2, its=2,
        sparseness=0.1, cthresh=50  )
      statfamask<-eigSeg(  famask, eanat$eig  ) %>% thresholdImage(1,Inf)
#      plot( famask, statfamask, window.overlay=c( mean(temp), max(temp) ) )
      statfamaskdil<-iMath( statfamask, "GD", 8 )
      cmask<-cropImage( statfamask, statfamaskdil )
      rad<-rep( 1, 3 );  mr<-c( 2, 1 )
      inds2<-1:round(length(selector)*0.75)
      crfm<-mrvnrfs( predval[selector][inds2], trnlist[inds2],
                    cmask,  rad=rad,
                    nsamples = 200,  asFactors=F,
                    ntrees = 1000 , multiResSchedule=mr )

      # now get best voxels from the training data
      crestrain<-mrvnrfs.predict( crfm$rflist, trnlist[-inds2] ,
        cmask, rad=rad, multiResSchedule=mr, asFactors=F )
      predmat<-imageListToMatrix( unlist(crestrain$probs) , cmask )
      etr<-apply( abs(predmat - predval[selector][-inds2]),
        FUN=mean, MARGIN=2 )

      # retrain on all data
      crfm<-mrvnrfs( predval[selector], trnlist,
                    cmask,  rad=rad,
                    nsamples = 200,  asFactors=F,
                    ntrees = 1000 , multiResSchedule=mr )

      # now apply to test data and make images
      teslist<-list()
      for ( i in 1:length(sublistes) )
        teslist[[i]]<-list( sublistes[[i]] )
      cres<-mrvnrfs.predict( crfm$rflist, teslist ,
        cmask, rad=rad, multiResSchedule=mr, asFactors=F )
      predmat<-imageListToMatrix( unlist(cres$probs) , cmask )
      corrmat<-antsrimpute( cor( predmat , predval[-selector] ) )
      corimg<-makeImage( cmask, corrmat )
      corimg<-decropImage( corimg, famask )
      plot( famask, corimg, window.overlay=c( 0.33, max(corimg) ) )
      besttr <- head( order(etr) , round( 0.5 * ncol(predmat) ) )
      rfpred <- apply( predmat[ , besttr ], FUN=mean, MARGIN=1 )
      print( cor.test( predval[-selector], rfpred ) )
#      plot(  predval[-selector], rfpred  )
      print( summary( lm( predval[-selector] ~  rfpred + . ,
        data = subdemog[-selector,confinds]  ) ) )
    setwd("/Users/stnava/data/TOT")
    stop("doingy-doingy!")
    }
    predmat[,3]<-rank(predmat[,3] )
    mat2<-antsrimpute(data.matrix(subdemog[,c(10,12,confinds)]))
    set.seed(9)
    selector<-createDataPartition( poorOrNot, p=0.5 )$Resample1
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
