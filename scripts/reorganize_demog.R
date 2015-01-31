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
if ( processDT )
  {
  library(ANTsR)
  setwd("/Users/stnava/data/TOT/tempdti")
  wfn<-na.omit(demogreo$DTFileExtension[ demogreo$DTQualityMeasure == 1])
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
      ref<-antsImageRead("fa_ref.nii.gz",3)
      refm<-antsImageRead("fa_ref_mask.nii.gz",3)
      mytx<-antsRegistration(fixed=ref , moving=fa ,
             typeofTransform = c("SyN"),
             outprefix=tempfile(),mask=refm)
      antsImageWrite( mytx$warpedmovout , onm  )
      }
    }
  }
