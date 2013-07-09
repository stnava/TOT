

library(visreg)
ff<-read.csv('data/TOT_Study_Preliminary_Image_Variables.csv')
mdl<-lm( ADJ_HOMETotalScore  ~ Cortex +  deepGrayMatter +  whiteMatter + csf +
GestAge + BMIdel +   MatAgeAtDel , data=ff ) 
print( summary( mdl ) )
pp<-predict( mdl ) 
plot( pp ~ ff$ADJ_HOMETotalScore[as.numeric(names(pp))]  )

x<-ff
y<-ff$ADJ_HOMETotalScore
nfold<-8
mycode<-c(1:nfold)
grouping <- sample( rep( mycode , 1000 ) )[1:nrow(ff)]
ferr<-rep(NA,nfold)
predicted<-rep(NA,nrow(ff))
for ( ii in 1:nfold )  {
  x1<-subset( x, subset = ( grouping != ii ) )
  x2<-subset( x, subset = ( grouping == ii ) )
  y1<-y[ grouping != ii ]
  y2<-y[ grouping == ii ]
  mdl<-lm( ADJ_HOMETotalScore  ~ Cortex +  deepGrayMatter +  whiteMatter + csf +
GestAge + BMIdel +   MatAgeAtDel , data=x1 ) 
  fitter<-predict( mdl , newdata=x2 )
  predicted[ grouping == ii ]<-fitter
  ferr[ii]<-mean( abs(  fitter - y2  ) )
}
print( cor.test( predicted, y )  )
print( mean( ferr , na.rm = T ) )
crossValHome<-predicted
realHome<-y 
mdl<-lm( realHome ~ crossValHome )
tit<-paste("TOT",nfold,"fold cross-validation predicted Home vs real Home")
fn<-paste(tit,'.pdf',sep='')
pdf(fn)
visreg(  mdl , main=tit )
dev.off()
print( summary( mdl ) )

