require('tidyverse')
require('magrittr')

source('../NHANES_diagnostics2/parameters.R')
nhanes_old=read.csv('Data_raw/nh3bpdat290716.csv')
nhanes_update=read.csv('NHANES_ML_update.csv',colClasses = 'integer') %>%
  mutate(yrsfuExam = round(permth_exm/12,2), yrsfuHome = round(permth_int/12, 2)) %>%
  filter(eligstat == 1)

nhanes_old %<>% mutate(yrsfuExam = round(yrsfuExam, 2), yrsfuHome = round(yrsfuHome,2), permth_exm = round(yrsfuExam *12), permth_int = round(yrsfuHome*12))

# Compare old and new mortality
old_dead <- sort(nhanes_old$SEQN[nhanes_old$dead==1])
new_dead <- sort(nhanes_update$seqn[nhanes_update$mortstat==1])
conflict_dead <- old_dead[ !(old_dead %in% new_dead)] # Tested: No resurrections
missing_update <- nhanes_old$SEQN[!(nhanes_old$SEQN %in% nhanes_update$seqn)] #Tested: No one from old data set missing from followup
missing_old <- !(nhanes_update$seqn %in% nhanes_old$SEQN) # 6 subjects in the updated data set with mortality data (all alive) who were not in the old data set.

nhanes_update %<>% filter(!missing_old) # remove them
identical(nhanes_old$SEQN,nhanes_update$seqn) # Sequence of ids now identical
fu_diff <- nhanes_update$permth_int - nhanes_update$permth_exm
fu_diff_update <- nhanes_update$permth_int - nhanes_old$permth_int

identical(nhanes_old$UCOD_LEADING[!is.na(nhanes_old$UCOD_LEADING)], nhanes_update$ucod_leading[!is.na(nhanes_old$UCOD_LEADING)])
#TRUE, so cause-of-death codes for all individuals who were dead in the first data set are identical

conflict_old <- subset(nhanes_old, nhanes_old$SEQN %in% conflict_dead)
conflict_new <- subset(nhanes_update, nhanes_update$seqn %in% conflict_dead)

nhanes <- nhanes_old
nhanes$UCOD_LEADING <- nhanes_update$ucod_leading
nhanes$yrsfuExam <- nhanes_update$yrsfuExam
nhanes$yrsfuHome <- nhanes_update$yrsfuHome
nhanes$dead <- nhanes_update$mortstat
nhanes$mrtHrt <- as.integer(nhanes$UCOD_LEADING==1)  # 1 codes for Heart
nhanes$mrtNeo <- as.integer(nhanes$UCOD_LEADING==2)  # 2 codes for Neoplasm
nhanes$mrtInj <- as.integer(nhanes$UCOD_LEADING==4)  # 4 codes for Injury 
nhanes$mrtOtherCVD <- as.integer(nhanes$UCOD_LEADING==5)  # 5 codes for other CVD

# Change NAs to 0
nhanes$mrtHrt[is.na(nhanes$mrtHrt)] <- 0
nhanes$mrtNeo[is.na(nhanes$mrtNeo)] <- 0
nhanes$mrtInj[is.na(nhanes$mrtInj)] <- 0
nhanes$mrtOtherCVD[is.na(nhanes$mrtOtherCVD)] <- 0

nhanes$eventhrt <- with(nhanes, mrtOtherCVD+mrtHrt)
nhanes$eventother <- with(nhanes, dead-eventhrt)


nhanes$yrsfu=pmin(nhanes_update$yrsfuHome,nhanes$yrsfuExam, na.rm = TRUE)
# Note: The most recent examination (home or clinic) is a left truncation time,
#   so follow-up is the minimum. For some the clinic is missing.
#   We call the follow-up the home follow-up time, but it doesn't matter,
#   since those individuals will be excluded from analysis.

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
