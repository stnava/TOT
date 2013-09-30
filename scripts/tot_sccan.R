library( robust )
library( visreg )
library( foreign )
library( psych )
library( ANTsR )
library( e1071 )
eps<-1.e-6
labelnames<-read.csv('data/labels_surface_DKT31_fullnames.csv')
ff<-read.csv('data/TOT_Study_Preliminary_Image_Variables.csv')
demog<-data.frame( read.spss('data/TOT_Demographics_MRI_Analyses July2013.sav') )
names(demog)[1]<-"ID"
tissue<-read.csv('data/tissue_volume.csv')
ffmerge<-merge( demog, tissue )  #
inds<-c(1:nrow(ffmerge))
bvol<-apply( ffmerge[,13:15], FUN = sum, MARGIN = 1 )
ff<-ffmerge
ff$ADJ_HOMETotalScore<-ff$ADJ_HOMETotalScore*(-1)
####### some imputation ##### 
ff$HC[ is.na( ff$HC ) ]<-mean( ff$HC,na.rm = T ) 
ff$BMIprepreg[ is.na( ff$BMIprepreg ) ]<-mean( ff$BMIprepreg,na.rm = T ) 
ff$BMIdel[ is.na( ff$BMIdel ) ]<-mean( ff$BMIdel,na.rm = T ) 
ff$HOMETotalScore[ is.na( ff$HOMETotalScore ) ]<-mean( ff$HOMETotalScore,na.rm = T ) 
ff$ADJ_HOMETotalScore[ is.na( ff$ADJ_HOMETotalScore ) ]<-mean( ff$ADJ_HOMETotalScore,na.rm = T ) 
######### setup study data #########
if ( !exists("np") ) np<-5000
whbrain<-c(17,19:ncol(ff))
if ( !exists("studyhome") ) studyhome<-FALSE
if ( studyhome ) {
  # home variable study
  wh<-c(11)
  demog<-as.matrix(  ff[,wh] )  # just adj_home , negated 
  brain<-as.matrix( ff[ ,c(3,4,6,8,whbrain)] ) # with BMI
  nv<-1
} else {
 wh<-c( 3, 4, 6, 8 )
 demog<-as.matrix(  ff[,wh] )  # just adj_home
 rbrain<-impute( ff[ ,whbrain] )
 rbrain<-residuals( lm(   as.matrix(rbrain) ~ bvol  ) )
 brain<-as.matrix( cbind(  bvol, rbrain ) )
 nv<-3
}
colnames( demog )<-colnames( ff )[wh]
######### setup analysis #########
sccan<-sparseDecom2( inmatrix=list( demog , brain ), inmask = c( NA , NA ) ,
  sparseness=c( 0.25 ,  0.1 ), nvecs=nv, its=20, smooth=0, perms=np, cthresh = c(0, 0), robust=1 )
for ( ind in 1:nv ) {
  print(paste("Sccanvec",ind,"pvalue",sccan$ccasummary[1,ind+1],"corr",sccan$ccasummary[2,ind+1]))
  print( paste( colnames(demog)[ abs(sccan$eig1[,ind]) > eps ] ) )
  print( ( sccan$eig1[,ind] )[ abs(sccan$eig1[,ind]) > eps] / sum(  abs(sccan$eig1[,ind])  ) )
  print("Imaging Predictors")
  myanat<-colnames(brain)[ abs(sccan$eig2[,ind]) > eps ]
  wt2<-( sccan$eig2[,ind] )[ abs(sccan$eig2[,ind]) > eps ] 
  for ( j in 1:length(myanat) )  {
    selector<-labelnames$Anatomy[  myanat[j] == labelnames$ID ] 
    if ( length( selector ) > 0 ) print(  paste( selector , wt2[j]/sum(wt2)  , sum(wt2[1:j])/sum(wt2) ) )else print( paste(  myanat[j] , wt2[j]/sum(wt2)   , sum(wt2[1:j])/sum(wt2) ) )
  }
}
print( p.adjust( sccan$ccasummary[1,2:(ncol(sccan$ccasummary)-1)] , method="BH" ) )
