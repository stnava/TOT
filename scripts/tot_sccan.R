library( robust )
library( visreg )
library( foreign )
library( psych )
library( ANTsR )
library( e1071 )
library( pheatmap ) 
eps<-1.e-6
labelnames<-read.csv('data/labels_surface_DKT31_fullnames.csv')
ff<-read.csv('data/TOT_Study_Preliminary_Image_Variables.csv')
demogin<-data.frame( read.spss('data/ShortwithHOMENewbornChartReview2013-09-26.sav') )
demog<-data.frame( read.spss('data/ShortwithHOMENewbornChartReview2013-09-26.sav') )
names(demog)[1]<-"ID"
tissue<-read.csv('data/tissue_volume.csv')
ffmerge<-merge( demog, tissue )  #
inds<-c(1:nrow(ffmerge))
bvol<-apply( ffmerge[,c(19:21)], FUN = sum, MARGIN = 1 )
bvol<-c( impute(cbind(bvol,bvol))[,1])
initialN<-length(bvol)
ffmerge<-cbind( ffmerge, bvol ) 
mysel<-rep( c(1,2) , (length( bvol )/2-1) )
mysel<-rnorm( length( bvol ) ) < 1.e7
ff<-subset( ffmerge , mysel )
####### some imputation ##### 
ff$BMIprepreg[ is.na( ff$BMIprepreg ) ]<-mean( ff$BMIprepreg,na.rm = T ) 
ff$BMIdel[ is.na( ff$BMIdel ) ]<-mean( ff$BMIdel,na.rm = T ) 
######### setup study data #########
if ( !exists("np") ) np<-1000
if ( !exists("sp") ) sp<-0.1
if ( !exists("sp1") ) sp1<-( -0.2 )
whbrain<-c(25:(ncol(ff)-1))
rbrain<-impute( ff[ ,whbrain] )
rbrain<-residuals( lm(   as.matrix(rbrain) ~ ff$bvol ) )
brain<-as.matrix(  data.frame( bvol = ff$bvol, rbrain  ) )
wh<-c( 5 , 8, 12 , 18 )
wh<-c( 5 , 12 , 18, 10 )  # 10 = EGAcalc,  11 = AdmHC
wh<-c( 5 , 8, 10, 12, 18 )  # 10 = EGAcalc,  11 = AdmHC
demog<-as.matrix(  ff[,wh] )  # just adj_home
nv<-as.numeric( ncol(demog) )
colnames( demog )<-colnames( ff )[wh]
######### setup analysis #########
myrob<-1
searchsp<-(c(1:21)-11)/40
# searchsp<-(c(1:21))/40
spmat<-matrix( rep(NA, length(searchsp)*length(searchsp) ) , ncol=length(searchsp) )
spval1<-matrix( rep(NA, length(searchsp)*length(searchsp) ) , ncol=length(searchsp) )
spval2<-matrix( rep(NA, length(searchsp)*length(searchsp) ) , ncol=length(searchsp) )
ct1<-ct2<-1
for ( sp1 in searchsp ) {
  ct2<-1
  for ( sp in searchsp )  {
    sccan<-sparseDecom2( inmatrix=list( demog , brain ), inmask = c( NA , NA ) , mycoption = 0,
                        sparseness=c( sp1, sp ), nvecs=nv, its=25, smooth=0, perms=1, cthresh = c(0, 0), robust=myrob )
    spmat[ct1,ct2]<-mean(as.numeric(sccan$ccasummary[2,2:(ncol(sccan$ccasummary)-1)]),na.rm=F)
    spval1[ct1,ct2]<-sp1
    spval2[ct1,ct2]<-sp
    ct2<-ct2+1
  }
  ct1<-ct1+1
}
rownames( spmat )<-searchsp
colnames( spmat )<-searchsp
pheatmap( spmat , cluster_rows = F, cluster_cols = F)
# inspect the heat map and decide on sparsenss 
sp<-0.1 # spval1[which.max( spmat )]
sp1<-0.2 # spval2[which.max( spmat )]
sccan<-sparseDecom2( inmatrix=list( demog , brain ), inmask = c( NA , NA ) , mycoption = 0,
                    sparseness=c( sp1, sp ), nvecs=nv, its=25, smooth=0, perms=np, cthresh = c(0, 0), robust=myrob )

for ( ind in 1:nv ) {
  print(paste("Sccanvec",ind,"pvalue",sccan$ccasummary[1,ind+1],"corr",sccan$ccasummary[2,ind+1]))
  print( paste( colnames(demog)[ abs(sccan$eig1[,ind]) > eps ] ) )
  print( ( sccan$eig1[,ind] )[ abs(sccan$eig1[,ind]) > eps] / sum(  abs(sccan$eig1[,ind])  ) )
  print("Imaging Predictors")
  myanat<-colnames(brain)[ abs(sccan$eig2[,ind]) > eps ]
  wt2<-( sccan$eig2[,ind] )[ abs(sccan$eig2[,ind]) > eps ] 
  for ( j in 1:length(myanat) )  {
    selector<-labelnames$Anatomy[  myanat[j] == labelnames$ID ] 
    if ( length( selector ) > 0 ) print(  paste( selector , wt2[j]/sum(abs(wt2))  , sum(abs(wt2[1:j]))/sum(abs(wt2)) , which( colnames(brain) == myanat[j] ) ) )else print( paste(  myanat[j] , wt2[j]/sum(abs(wt2))  , sum(abs(wt2[1:j]))/sum(abs(wt2)) ) )
  }
}
print("Adjusted p-values") 
print( p.adjust( sccan$ccasummary[1,2:(ncol(sccan$ccasummary)-1)] , method="BH" ) )
print(sum(as.numeric(mysel))/initialN)
print( colnames(demog) ) 
write.csv(brain,'temp2.csv',row.names=F,quote=F)
write.csv(demog,'temp1.csv',row.names=F,quote=F)
print( mean(as.numeric(sccan$ccasummary[2,2:(ncol(sccan$ccasummary)-1)]),na.rm=F) )

