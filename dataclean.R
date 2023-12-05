source('NHANESDiagnostics/parameters.R')
nhanes=read.csv('Data_raw/nh3bpdat290716.csv')
nhanes$yrsfu=nhanes$yrsfuHome

whichsys=match(c('systolicA','systolicB','systolicC'),names(nhanes))
whichdias=match(c('diastolicA','diastolicB','diastolicC'),names(nhanes))  

whichsyshome=match(c('systolicAhome','systolicBhome','systolicChome'),names(nhanes))
whichdiashome=match(c('diastolicAhome','diastolicBhome','diastolicChome'),names(nhanes))
whichBP=c(whichsys,whichdias)
whichBPhome=c(whichsyshome,whichdiashome)
allBP=c(whichBP,whichBPhome)
allsys=c(whichsys,whichsyshome)
alldias=c(whichdias,whichdiashome)


### Means and variances

sys=nhanes[,whichsys]
dias=nhanes[,whichdias]

sysH=nhanes[,whichsyshome]
diasH=nhanes[,whichdiashome]

nhanes$meandias=apply(dias,1,mean)
nhanes$meansys=apply(sys,1,mean)
nhanes$sysDel=apply(sysH,1,mean)-apply(sys,1,mean)
nhanes$diasDel=apply(diasH,1,mean)-apply(dias,1,mean)
nhanes$vardias=apply(dias,1,var)
nhanes$varsys=apply(sys,1,var)
nhanes$sddias=sqrt(nhanes$vardias)
nhanes$sdsys=sqrt(nhanes$varsys)

nhanes$precdias=1/(nhanes$vardias+1/3)
nhanes$precsys=1/(nhanes$varsys+1/3)

require('survival')

race <- with(nhanes, ifelse(white==1, 2, ifelse(black==1, 1,ifelse(mexican==1,3, 0))))
nhanes$race <- as.factor(race)
#sexrace=factor(c('other','female, black','female, white','female, Mex','male, black','male, white','male, Mex'),
#   levels=c('other','female, black','female, white','female, Mex','male, black','male, white','male, Mex'),ordered=TRUE)
type=factor(1+race+3*(1-nhanes$female)*(race>0)+6*(race==0)) # other race all type 7
levels(type)=c('female, black','female, white','female, Mex','male, black','male, white','male, Mex','other')
nhanes$type=type

nhanesna=!(is.na(apply(nhanes[,allBP],1,prod)))

nhanessysrange=(apply(nhanes[,allsys],1,min)<60)|(apply(nhanes[,allsys],1,max)>250)
nhanesdiasrange=(apply(nhanes[,alldias],1,min)<40)|(apply(nhanes[,alldias],1,max)>140)
nhanesgood=(nhanesna&(nhanes$yrsfu>0)&!nhanessysrange&!nhanesdiasrange&nhanes$other==0)

######  Causes of death

nhanes$eventhrt <- with(nhanes, mrtOtherCVD+mrtHrt)
nhanes$eventother <- with(nhanes, dead-eventhrt)

frs=read.csv('Data_raw/FRS.csv')

#Choose versions of FRS
if(whichfrs==1){nhanes$FRS=frs$ATP.FRS[match(nhanes$SEQN,frs$SEQN)]} else{nhanes$FRS=frs$X1998.FRS[match(nhanes$SEQN,frs$SEQN)]}

nhanesA=nhanes[nhanesgood,]
nhanesA$type=factor(nhanesA$type,exclude=7)

#### Restrict based on FRS? 0= use everything, 1= only with FRS, -1 = only w/o FRS

restrict=function(DF, restcodename, restcode){
  if (!is.numeric(restcodename)){
    restcodename=match(restcodename,names(DF))
  }
  if (restcode==0){return(DF)}
  if (restcode==1){return(DF[!is.na(DF[,restcodename]),])}
  if (restcode==-1){return(DF[is.na(DF[,restcodename]),])}
}

################### FRS############################
nhanesC=subset(restrict(nhanesA,'FRS',frsrestrict),age<agemax&age>agemin)

sex=2-nhanesC$female
race=as.numeric(nhanesC$race)-1

age=nhanesC$age
finage=nhanesC$age+nhanesC$yrsfu
if (whichevents==2){dead=nhanesC$eventhrt}else {dead=nhanesC$dead}
n=length(age)

save(nhanesA,file='Data_cleaned/nhanesA.RData')
