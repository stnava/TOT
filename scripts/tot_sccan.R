library( robust )
library(visreg)
library(foreign)
library(psych)
library(ANTsR)
library(e1071)
ctl<-lmrob.control( "KS2011", max.it = 1000 )
ff<-read.csv('data/TOT_Study_Preliminary_Image_Variables.csv')
demog<-data.frame( read.spss('data/TOT_Demographics_MRI_Analyses July2013.sav') )
names(demog)[1]<-"ID"
tissue<-read.csv('data/tissue_volume.csv')
ffmerge<-merge( demog, tissue )  #
inds<-c(1:nrow(ffmerge))
bvol<-apply( ffmerge[,13:15], FUN = sum, MARGIN = 1 )
ff<-ffmerge
ff$BMIprepreg[ is.na( ff$BMIprepreg ) ]<-mean( ff$BMIprepreg,na.rm = T ) 
ff$BMIdel[ is.na( ff$BMIdel ) ]<-mean( ff$BMIdel,na.rm = T ) 
ff$HOMETotalScore[ is.na( ff$HOMETotalScore ) ]<-mean( ff$HOMETotalScore,na.rm = T ) 
ff$ADJ_HOMETotalScore[ is.na( ff$ADJ_HOMETotalScore ) ]<-mean( ff$ADJ_HOMETotalScore,na.rm = T ) 
#
confounds<-cbind( ff$BMIprepreg, ff$BMIdel )
confounds<-cbind( confounds, ff$MatAgeAtDel*(-1) )
confounds<-cbind( confounds, ff$EDU_MOM_HS_Status )
confounds<-cbind( confounds, ff$GestAge )
######### the study #########
np<-5000
# home vars 
demog<-as.matrix(  ff[,c(11,12)] )
brain<-as.matrix( cbind( ff[ ,c(3,6,7,8,4)]  , ff[ ,13:ncol(ff) ] ) )
# imaging vars separate
demog<-as.matrix(  ff[,c(3,4,6,7,8,11,12)] )
brain<-as.matrix( cbind( bvol , ff[ ,13:ncol(ff) ] ) )
brain<-impute(brain)
nv<-min( 3, ncol(demog) )
sccan<-sparseDecom2( inmatrix=list( demog , brain ), inmask = c( NA , NA ) ,
  sparseness=c( 0.33, 0.1 ), nvecs=nv, its=30, smooth=0, perms=np, cthresh = c(0, 0), robust=0 )
eps<-1.e-6
for ( ind in 1:nv ) {
  print(paste("Sccanvec",ind))
  print( paste( colnames(demog)[ abs(sccan$eig1[,ind]) > eps ] ) )
  print( paste( colnames(brain)[ abs(sccan$eig2[,ind]) > eps ] ) )
  print( ( sccan$eig1[,ind] )[abs(sccan$eig1[,ind]) > eps] )
  print( ( sccan$eig2[,ind] )[ abs(sccan$eig2[,ind]) > eps ] )
}
