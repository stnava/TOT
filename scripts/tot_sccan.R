library( robust )
library( visreg )
library( foreign )
library( psych )
library( ANTsR )
library( e1071 )
eps<-1.e-6
labelnames<-read.csv('data/labels_surface_DKT31_fullnames.csv')
ff<-read.csv('data/TOT_Study_Preliminary_Image_Variables.csv')
demog<-data.frame( read.spss('data/ShortwithHOMENewbornChartReview2013-09-26.sav') )
names(demog)[1]<-"ID"
tissue<-read.csv('data/tissue_volume.csv')
ffmerge<-merge( demog, tissue )  #
inds<-c(1:nrow(ffmerge))
bvol<-apply( ffmerge[,c(19:21)], FUN = sum, MARGIN = 1 )
bvol<-c( impute(cbind(bvol,bvol))[,1]) 
ff<-ffmerge
####### some imputation ##### 
ff$BMIprepreg[ is.na( ff$BMIprepreg ) ]<-mean( ff$BMIprepreg,na.rm = T ) 
ff$BMIdel[ is.na( ff$BMIdel ) ]<-mean( ff$BMIdel,na.rm = T ) 
ff$ZTotalHOMEADJ<-ff$ZTotalHOMEADJ*(-1)  # negatively correlated with BMI 
######### setup study data #########
if ( !exists("np") ) np<-5000
whbrain<-c(23,25:ncol(ff))
rbrain<-impute( ff[ ,whbrain] )
rbrain<-residuals( lm(   as.matrix(rbrain) ~ bvol  ) )
brain<-as.matrix( cbind( bvol, rbrain ) )
if ( !exists("studyhome") ) studyhome<-FALSE
if ( studyhome ) {
  # home variable study
  wh<-c(18)
  demog<-as.matrix(  ff[,wh] )  # just adj_home , negated 
  brain<-as.matrix( cbind( ff[ ,c(5, 8, 13 )] , brain ) ) # with BMI
  nv<-1
} else {
 wh<-c( 5 , 8, 13, 18 )
 demog<-as.matrix(  ff[,wh] )  # just adj_home
 nv<-3
}
colnames( demog )<-colnames( ff )[wh]
######### setup analysis #########
sccan<-sparseDecom2( inmatrix=list( demog , brain ), inmask = c( NA , NA ) ,
  sparseness=c( 0.15 ,  0.05 ), nvecs=nv, its=15, smooth=0, perms=np, cthresh = c(0, 0), robust=1 )
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
