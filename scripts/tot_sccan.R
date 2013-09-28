library( robust )
library(visreg)
library(foreign)
library(psych)
library(ANTsR)
library(e1071)
ctl<-lmrob.control( "KS2011", max.it = 1000 )
labelnames<-read.csv('data/labels_surface_DKT31_fullnames.csv')
ff<-read.csv('data/TOT_Study_Preliminary_Image_Variables.csv')
demog<-data.frame( read.spss('data/TOT_Demographics_MRI_Analyses July2013.sav') )
names(demog)[1]<-"ID"
tissue<-read.csv('data/tissue_volume.csv')
ffmerge<-merge( demog, tissue )  #
inds<-c(1:nrow(ffmerge))
bvol<-apply( ffmerge[,13:15], FUN = sum, MARGIN = 1 )
ff<-ffmerge
####### some imputation ##### 
ff$BMIprepreg[ is.na( ff$BMIprepreg ) ]<-mean( ff$BMIprepreg,na.rm = T ) 
ff$BMIdel[ is.na( ff$BMIdel ) ]<-mean( ff$BMIdel,na.rm = T ) 
ff$HOMETotalScore[ is.na( ff$HOMETotalScore ) ]<-mean( ff$HOMETotalScore,na.rm = T ) 
ff$ADJ_HOMETotalScore[ is.na( ff$ADJ_HOMETotalScore ) ]<-mean( ff$ADJ_HOMETotalScore,na.rm = T ) 
######### setup study data #########
np<-5000
studyhome<-FALSE
if ( studyhome ) {
  # home variable study
  wh<-c(11)
  brain<-as.matrix( cbind( ff[ ,c(3,6,7,8,4)]  , ff[ ,13:ncol(ff) ] ) ) # with BMI
} else {
 wh<-c(3,4,6,7,8,11,12)
 brain<-as.matrix( cbind( bvol , ff[ ,13:ncol(ff) ] ) )
}
demog<-as.matrix(  ff[,wh] )  # just adj_home
colnames( demog )<-colnames( ff )[wh]
brain<-impute(brain)
######### setup analysis #########
nv<-min( 3, ncol(demog) )
sccan<-sparseDecom2( inmatrix=list( demog , brain ), inmask = c( NA , NA ) ,
  sparseness=c( 0.33, 0.1 ), nvecs=nv, its=30, smooth=0, perms=np, cthresh = c(0, 0), robust=1 )
eps<-1.e-6
for ( ind in 1:nv ) {
  print(paste("Sccanvec",ind))
  print( paste( colnames(demog)[ abs(sccan$eig1[,ind]) > eps ] ) )
  print( ( sccan$eig1[,ind] )[ abs(sccan$eig1[,ind]) > eps] )
  print("Imaging Predictors")
  myanat<-colnames(brain)[ abs(sccan$eig2[,ind]) > eps ]
  wt2<-( sccan$eig2[,ind] )[ abs(sccan$eig2[,ind]) > eps ] 
  for ( j in 1:length(myanat) )  {
    selector<-labelnames$Anatomy[  myanat[j] == labelnames$ID ] 
    if ( length( selector ) > 0 ) print(  paste( selector , wt2[j] , sum(wt2[1:j])/sum(wt2) ) )else print( paste(  myanat[j] , wt2[j]  , sum(wt2[1:j])/sum(wt2) ) )
  }
}
