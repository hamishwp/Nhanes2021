nho = order(nhanesA$yrsfu+nhanesA$age)
nhsh = nhanesA[,c('age','yrsfu','mrtHrt','mrtOtherCVD','sysDel','meansys','white','female','sdsys')]
nhsh$death=nhsh$mrtOtherCVD+nhsh$mrtHrt
nhsh$end = nhsh$age+nhsh$yrsfu

nhom = order(nhsh$end,1-nhsh$death)
nhsh=nhsh[nhom,]

wnh= filter(nhsh, white==1 & female ==0)
ncon = 0
ndis = 0
L = dim(wnh)[1]

for (i in seq(L-1)){
  a1= nhsh$end[i]
  a2=nhsh$end[seq(i+1,L)]
  m1= nhsh$death[i]
  m2=nhsh$death[seq(i+1,L)]
  b1= abs(nhsh$sdsys[i])
  b2 = abs(nhsh$sdsys[seq(i+1,L)])
  s1= nhsh$age[i]
  s2 = nhsh$age[seq(i+1,L)]
  if (m1==1){
    ncon = ncon + sum((b1>b2)*(a1>=s2))+ sum((b1==b2)*(a1>=s2))/2
    ndis = ndis + sum((b1<b2)*(a1>=s2))+ sum((b1==b2)*(a1>=s2))/2
  }
}

scorename='meansys'

CindCI = function(scorename, sex,race){
  ch = vector('integer',L)
  dh= vector('integer',L)
  wnh= filter(nhsh, white==race & female ==sex)
  L = dim(wnh)[1]

for (i in seq(L)){
  a1= wnh$end[i]
  a2=wnh$end[-i]
  m1= wnh$death[i]
  m2=wnh$death[-i]
  b1= abs(wnh[i,scorename])
  b2 = abs(wnh[-i,scorename])
  s1= wnh$age[i]
  s2 = wnh$age[-i]
    ch[i] = sum((b1>b2)*(a1>=s2)*(a2>a1)*m1)+ sum((b1==b2)*(a1>=s2)*(a2>a1)*m1)/2 +
      sum((b1<b2)*(a1>=s2)*(a2<a1)*m2)+ sum((b1==b2)*(a1>=s2)*(a2<=a1)*m2)/2 +
      sum((b1>b2)*(a2==a1)*m2*(1-m1))+ sum((b1<b2)*(a2==a1)*m1*(1-m2))
    dh[i] = sum((b1<b2)*(a1>=s2)*(a2>a1)*m1)+ sum((b1==b2)*(a1>=s2)*(a2>a1)*m1)/2 +
      sum((b1>b2)*(a1>=s2)*(a2<=a1)*m2)+ sum((b1==b2)*(a1>=s2)*(a2<=a1)*m2)/2 +
      sum((b1<b2)*(a2==a1)*m2*(1-m1))+ sum((b1>b2)*(a2==a1)*m1*(1-m2))
}

pc=sum(ch)/L/(L-1)
pd=sum(dh)/L/(L-1)

pcc= sum(ch*(ch-1)/L/(L-1)/(L-2))
pdd= sum(dh*(dh-1)/L/(L-1)/(L-2))
pcd= sum(ch*dh/L/(L-1)/(L-2))

Cind= pc/(pc+pd)
Cvar = 4*(pd^2*pcc - 2*pc*pd*pcd + pc^2*pdd)/(pc+pd)^4/L
Csd = sqrt(Cvar)

print(paste0(scorename,': C-Index = ',Cind, ', Conf. Int. = (',Cind-2*Csd, ' , ', Cind+2*Csd, ')'))
}

CindCI('sysDel',1,1)
CindCI('sdsys',1,1)
CindCI('meansys',1,1)
CindCI('sysDel',0,1)
CindCI('sdsys',0,1)
CindCI('meansys',0,1)
CindCI('sysDel',1,0)
CindCI('sdsys',1,0)
CindCI('meansys',1,0)
CindCI('sysDel',1,1)
CindCI('sdsys',1,1)
CindCI('meansys',1,1)




require(survival)
s=Surv(nhanesA$age,nhanesA$age+nhanesA$yrsfu,m)

cp=coxph(s~b+strata(nhanesA$female)+strata(nhanesA$race))

wnh= filter(nhanesA, white==1 & female ==1)

s=Surv(wnh$age,wnh$age+wnh$yrsfu,(wnh$mrtHrt + wnh$mrtOtherCVD ))

cp=coxph(s~wnh$meansys)
summary(cp)

allrel=NULL
adth= wnh$yrsfu +wnh$age
m = (wnh$mrtHrt + wnh$mrtOtherCVD == 1)
b = (wnh$meansys)
for (i in 1:(length(adth)-1)){
  a1= adth[i]
  a2=adth[seq(i+1,length(adth))]
  m1= m[i]
  m2=m[seq(i+1,length(adth))]
  b1= b[i]
  b2 = b[seq(i+1,length(adth))]
  keepb = (a2>=a1 & m2==1) | (a2<= a1 & m1==1)
  allrel = c(allrel, exp(0.009231*abs(b1-b2[keepb])))
}