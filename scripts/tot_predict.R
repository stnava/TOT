

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
predicted<-rep(0,nrow(ff))
predictedct<-rep(0,nrow(ff))
for ( run in 1:100 ) {
  grouping <- sample( rep( mycode , 1000 ) )[1:nrow(ff)]
for ( ii in 1:nfold )  {
  x1<-subset( x, subset = ( grouping != ii ) )
  x2<-subset( x, subset = ( grouping == ii ) )
  y1<-y[ grouping != ii ]
  y2<-y[ grouping == ii ]
  mdl<-lm( ADJ_HOMETotalScore  ~ Cortex +  deepGrayMatter +  whiteMatter + csf +
GestAge + BMIdel +   MatAgeAtDel , data=x1 ) 
  fitter<-predict( mdl , newdata=x2 )
  predicted[ grouping == ii ]<-predicted[ grouping == ii ]+fitter
  fitter<-rep( 1 , sum( grouping == ii ) )
  predictedct[ grouping == ii ]<-predictedct[ grouping == ii ]+fitter
}
}
crossValHome<-predicted/predictedct
realHome<-y 
mdl<-lm( realHome ~ crossValHome )
tit<-paste("TOT",nfold,"fold cross-validation predicted Home vs real Home")
tit2<-paste("figs/TOT",nfold,"fold_cross-validation_predicted_Home_vs_real_Home",sep="_")
fn<-paste(tit2,'.pdf',sep='')
pdf(fn)
visreg(  mdl , main=tit )
dev.off()
print( summary( mdl ) )

