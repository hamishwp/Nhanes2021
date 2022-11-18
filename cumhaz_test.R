B=c(16.2,3.12,26.6,42.9,5.08,52.7)*1E-5
theta=c(.072,.0906,.0625,.0641,.09,.058)

# race black=1, white=2, mex=3
cumhaz = function(b,th,T){
  b*(exp(th*T)-1)/th
}

r= as.integer(nhanesA$race)-1
dem <- (1-nhanesA$female)*3 +r

CHstart = cumhaz(B[dem],theta[dem],nhanesA$age)
CHstop = cumhaz(B[dem],theta[dem],nhanesA$yrsfu + nhanesA$age)
CHdif = CHstop-CHstart

y1= nhanesA$yrsfu
y2= 

ord = sample(1:length(CHdif),length(CHdif)) # random rearrangement
success1 = CHdif > CHdif[ord] & 