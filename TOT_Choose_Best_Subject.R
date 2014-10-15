
totreview<-read.csv("TOT_BA_Notes.csv")
totreview$Keep<-totreview$meanM1_M2>=2
ids<-unique( totreview$ids )
keepnotes<-rep("",nrow(totreview))
nties<-0
for ( id in ids )
{
keeper<-totreview$meanM1_M2[ totreview$ids == id ]
if ( length(keeper) == 1 )
keepnotes[ totreview$ids == id ]<-
  as.character(totreview$Keep[ totreview$ids == id ])
if ( length(keeper) > 1 )
  {
  bestReviewed<-which(  totreview$ids == id &
    totreview$meanM1_M2== max(keeper) )
  worstReviewed<-which(  totreview$ids == id &
  totreview$meanM1_M2!= max(keeper) )
  if ( length(worstReviewed) > 0 )
   keepnotes[ worstReviewed ]<-"FALSE"
  if ( length(bestReviewed) == 1  )
    {
    keepnotes[ bestReviewed ]<-"TRUE"
    }
  else
    {
    msg<-"Tie: Need to decide which to keep"
    keepnotes[ bestReviewed ]<-msg
    nties<-nties+1
    }
  }
}
print(paste("Will Keep at min:",sum(keepnotes=="TRUE")+nties))
print(totreview$id[keepnotes=="TRUE"])
print(totreview$Scan[keepnotes=="TRUE"])
l1<-length(totreview$id[keepnotes=="TRUE"])
l2<-length(unique(totreview$id[keepnotes=="TRUE"] ))
if ( l1 != l2 ) print("there is a bug")
print("Subjects with ties")
print(totreview$Scan[grep("Tie",keepnotes)])
