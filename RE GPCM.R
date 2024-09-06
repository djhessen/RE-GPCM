library(fastDummies)
library(calculus)
library(msm) # delta method

dat<-read.table('Dekovic_Aggr.dat',header=FALSE)

v<-5
k<-dim(dat)[2]
if (min(dat)==0){m<-max(dat)}
if (v<1){v<-1} else
if (v%%2==1){q<-(v-1)/2} else
if (v%%2==0){q<-(v-2)/2}

fr<-rep(1,dim(dat)[1])
dat1<-data.frame(dat,fr)
names(dat1)<-NULL
names(dat1)<-c(paste('i',1:k,sep=''),'fr')

D<-aggregate(fr~.,data=dat1,sum)

D2<-dummy_cols(D,select_columns=paste('i',1:k,sep=''),remove_first_dummy=TRUE,remove_selected_columns=FALSE)

X<-as.matrix(D2[,(k+2):(k+m*k+1)]) 
Y<-as.matrix(D2[,1:k])
n_is<-as.vector(colSums(X*D2$fr))
n<-sum(D2$fr)

# RE-GPC model

ay<-function(x){
A<-matrix(,nrow=dim(Y)[1],ncol=(2*q+1))
for (c in 0:(2*q)){A[,c+1]<-(Y%*%x)^(c+2)/(c^2+3*c+2)}
return(A)}

ay1<-function(x){
A1<-matrix(,nrow=dim(Y)[1],ncol=(2*q+1))
for (c in 0:(2*q)){A1[,c+1]<-(Y%*%x)^(c+1)/(c+1)}
return(A1)}

e10<-3/5 #sqrt(1/2)
e20<-4/5

hc<-function(e1,e2){
if (v%%2==0){
H<-c(e10,e1)%*%t(c(e10,e1))+c(e20,e2)%*%t(c(e20,e2))
h<-rep(0,(2*q+1))
for (c in 2:(2*q+2)){for (i in 1:(q+1)){for (j in 1:(q+1)){if (i+j==c){h[c-1]<-h[c-1]+H[i,j]}}}}} else
if (v%%2==1){
if (v==3){H<-c(e10,e1)%*%t(c(e10,e1))+c(e20,0)%*%t(c(e20,0))} else
H<-c(e10,e1)%*%t(c(e10,e1))+c(e20,e2,0)%*%t(c(e20,e2,0))
h<-rep(0,(2*q+1))
for (c in 2:(2*q+2)){for (i in 1:(q+1)){for (j in 1:(q+1)){if (i+j==c){h[c-1]<-h[c-1]+H[i,j]}}}}}
return(h)}

K<-function(x){
if (v%%2==0){
ay(x[(k*m+1):(k*m+k)])%*%hc(x[(k*m+k+1):(k*m+k+q)],x[(k*m+k+q+1):(k*m+k+2*q)])} else
if (v%%2==1){
if (v==3){ay(x[(k*m+1):(k*m+k)])%*%hc(x[(k*m+k+1)],0)} else
ay(x[(k*m+1):(k*m+k)])%*%hc(x[(k*m+k+1):(k*m+k+q)],x[(k*m+k+q+1):(k*m+k+2*q-1)])}
} # cumulant generating function

if (v%%2==0){
stv<-c(log(n_is/(n-n_is)),rep(1,k),rep(0,q),rep(0,q))} else
if (v%%2==1&v==3){
stv<-c(log(n_is/(n-n_is)),rep(1,k),rep(0,q))} else
if (v%%2==1&v>1){
stv<-c(log(n_is/(n-n_is)),rep(1,k),rep(0,q),rep(0,q-1))} else
if (v%%2==1&v==1){
stv<-log(n_is/(n-n_is))}

dK<-function(x){
if (v%%2==0){
ay1(x[(k*m+1):(k*m+k)])%*%hc(x[(k*m+k+1):(k*m+k+q)],x[(k*m+k+q+1):(k*m+k+2*q)])} else
if (v%%2==1){
ay1(x[(k*m+1):(k*m+k)])%*%hc(x[(k*m+k+1):(k*m+k+q)],x[(k*m+k+q+1):(k*m+k+2*q-1)])}}

dKa<-function(x){
if (v==2){c(ay1(x[(k*m+1):(k*m+k)]))*Y} else
c(dK(x))*Y}

dhc1<-function(e1){
E1<-matrix(rep(0,((q+1)*(2*q+1))),nrow=q+1,ncol=2*q+1)
for (j in 0:q){for (c in j:(j+q)){if (j==c){E1[j+1,c+1]<-2*e10} else if (j<c){E1[j+1,c+1]<-2*e1[c-j]}}}
return(E1)}

dhc2<-function(e2){
if (v%%2==0){
E2<-matrix(rep(0,((q+1)*(2*q+1))),nrow=q+1,ncol=2*q+1)
for (j in 0:q){for (c in j:(j+q)){if (j==c){E2[j+1,c+1]<-2*e20} else if (j<c){E2[j+1,c+1]<-2*e2[c-j]}}}} else
if (v%%2==1){
E2<-matrix(rep(0,(q*(2*q+1))),nrow=q,ncol=2*q+1)
for (j in 0:(q-1)){for (c in j:(j+q-1)){if (j==c){E2[j+1,c+1]<-2*e20} else if (j<c){E2[j+1,c+1]<-2*e2[c-j]}}}}
return(E2)}

dKe1<-function(x){dhc1(x[(k*m+k+1):(k*m+k+q)])[2:(q+1),]%*%t(ay(x[(k*m+1):(k*m+k)]))}

dKe2<-function(x){
if (v%%2==0){
dhc2(x[(k*m+k+q+1):(k*m+k+2*q)])[2:(q+1),]%*%t(ay(x[(k*m+1):(k*m+k)]))} else
if (v%%2==1){
dhc2(x[(k*m+k+q+1):(k*m+k+2*q-1)])[2:q,]%*%t(ay(x[(k*m+1):(k*m+k)]))}}

by<-function(x){
if (v==1){exp(X%*%x[1:(k*m)])} else
if (v==2){exp(X%*%x[1:(k*m)]+ay(x[(k*m+1):(k*m+k)]))} else
if (v>2){exp(X%*%x[1:(k*m)]+K(x))}
}

l<-function(x){
if (v==1){n*log(sum(by(x)))-t(x[1:(k*m)])%*%n_is} else
if (v==2){n*log(sum(by(x)))-t(x[1:(k*m)])%*%n_is-t(D2$fr)%*%ay(x[(k*m+1):(k*m+k)])} else
if (v>2){n*log(sum(by(x)))-t(x[1:(k*m)])%*%n_is-t(D2$fr)%*%K(x)}
}

db<-function(x) {-n_is+n*t(X)%*%by(x)/sum(by(x))} # derivatives beta
da<-function(x){-t(dKa(x))%*%D2$fr+(n/sum(by(x)))*t(dKa(x))%*%by(x)} # derivatives alpha
de1<-function(x){-dKe1(x)%*%D2$fr+(n/sum(by(x)))*dKe1(x)%*%by(x)} # derivatives eta1
de2<-function(x){-dKe2(x)%*%D2$fr+(n/sum(by(x)))*dKe2(x)%*%by(x)} # derivatives eta2

der<-function(x){
if (v==1){db(x)} else
if (v==2){c(db(x),da(x))} else
if (v>2){c(db(x),da(x),de1(x),de2(x))}}

if (v%%2==0){
if (v==2){stv<-c(log(n_is/(n-n_is)),exp(runif(k,0,1)))} else
if (v>2){stv<-c(log(n_is/(n-n_is)),exp(runif(k,-1,0)),rep(0,2*q))}} else
if (v%%2==1){
if (v==1){stv<-log(n_is/(n-n_is))} else
if (v==3){stv<-c(log(n_is/(n-n_is)),rep(1,k),rep(0,q))} else
if (v>3){stv<-c(log(n_is/(n-n_is)),exp(runif(k,-1,0)),rep(0,2*q-1))}}
Afit<-optim(stv,l,der,
      method = "BFGS",
      control=list(trace=0,reltol=1e-20,abstol=1e-20,maxit=1000),
      lower=-Inf,upper=Inf,
      hessian = FALSE)

par<-Afit$par
round(der(par),3)

HES<-hessian(l,par,drop=TRUE,accuracy=4) # package calculus
cbind(round(par,3),round(sqrt(diag(solve(HES))),3))
cm<-solve(HES)[1:25,1:25]

# DELTA METHOD
se<-c(deltamethod(list(~-x1/x5,~(x1-x2)/x5,~(x2-x3)/x5,~(x3-x4)/x5),par[c(1:4,21)],cm[c(1:4,21),c(1:4,21)]),
      deltamethod(list(~-x1/x5,~(x1-x2)/x5,~(x2-x3)/x5,~(x3-x4)/x5),par[c(5:8,22)],cm[c(5:8,22),c(5:8,22)]),
      deltamethod(list(~-x1/x5,~(x1-x2)/x5,~(x2-x3)/x5,~(x3-x4)/x5),par[c(9:12,23)],cm[c(9:12,23),c(9:12,23)]),
      deltamethod(list(~-x1/x5,~(x1-x2)/x5,~(x2-x3)/x5,~(x3-x4)/x5),par[c(13:16,24)],cm[c(13:16,24),c(13:16,24)]),
      deltamethod(list(~-x1/x5,~(x1-x2)/x5,~(x2-x3)/x5,~(x3-x4)/x5),par[c(17:20,25)],cm[c(17:20,25),c(17:20,25)]))
delta<-c(-par[1]/par[21],(par[1]-par[2])/par[21],(par[2]-par[3])/par[21],(par[3]-par[4])/par[21],
         -par[5]/par[22],(par[5]-par[6])/par[22],(par[6]-par[7])/par[22],(par[7]-par[8])/par[22],
         -par[9]/par[23],(par[9]-par[10])/par[23],(par[10]-par[11])/par[23],(par[11]-par[12])/par[23],
         -par[13]/par[24],(par[13]-par[14])/par[24],(par[14]-par[15])/par[24],(par[15]-par[16])/par[24],
         -par[17]/par[25],(par[17]-par[18])/par[25],(par[18]-par[19])/par[25],(par[19]-par[20])/par[25])
cbind(round(delta,3),round(se,3))

omega<-c(-par[1]/par[21],-par[2]/(2*par[21]),-par[3]/(3*par[21]),-par[4]/(4*par[21]),
         -par[5]/par[22],-par[6]/(2*par[22]),-par[7]/(3*par[22]),-par[8]/(4*par[22]),
         -par[9]/par[23],-par[10]/(2*par[23]),-par[11]/(3*par[23]),-par[12]/(4*par[23]),
         -par[13]/par[24],-par[14]/(2*par[24]),-par[15]/(3*par[24]),-par[16]/(4*par[24]),
         -par[17]/par[25],-par[18]/(2*par[25]),-par[19]/(3*par[25]),-par[20]/(4*par[25]))
se3<-c(deltamethod(list(~-x1/x5,~-x2/(2*x5),~-x3/(3*x5),~-x4/(4*x5)),par[c(1:4,21)],cm[c(1:4,21),c(1:4,21)]),
       deltamethod(list(~-x1/x5,~-x2/(2*x5),~-x3/(3*x5),~-x4/(4*x5)),par[c(5:8,22)],cm[c(5:8,22),c(5:8,22)]),
       deltamethod(list(~-x1/x5,~-x2/(2*x5),~-x3/(3*x5),~-x4/(4*x5)),par[c(9:12,23)],cm[c(9:12,23),c(9:12,23)]),
       deltamethod(list(~-x1/x5,~-x2/(2*x5),~-x3/(3*x5),~-x4/(4*x5)),par[c(13:16,24)],cm[c(13:16,24),c(13:16,24)]),
       deltamethod(list(~-x1/x5,~-x2/(2*x5),~-x3/(3*x5),~-x4/(4*x5)),par[c(17:20,25)],cm[c(17:20,25),c(17:20,25)]))
cbind(round(omega,3),round(se3,3))

# cumulants

cum<-function(x){
f<-c()
for (c in 0:(2*q)){f[c+1]<-factorial(c)}
if (v%%2==0){
cumulants<-f*hc(x[(k*m+k+1):(k*m+k+q)],x[(k*m+k+q+1):(k*m+k+2*q)])} else
if (v%%2==1){
cumulants<-f*hc(x[(k*m+k+1):(k*m+k+q)],x[(k*m+k+q+1):(k*m+k+2*q-1)])}
return(cumulants)}

round(cum(par),3)

cm2<-solve(HES)[26:28,26:28]

# standard errors for cumulants in the case v=5
se2<-c(deltamethod(list(~6*x1/5+8*x3/5,~2*x1^2+2*x3^2+12*x2/5,~12*x1*x2,~24*x2^2),par[26:28],cm2))
round(se2,3)

ay2<-function(x){
A2<-matrix(,nrow=dim(Y)[1],ncol=(2*q+1))
for (c in 0:(2*q)){A2[,c+1]<-(Y%*%x)^c}
return(A2)}

dK2<-function(x){
if (v%%2==0){
ay2(x[(k*m+1):(k*m+k)])%*%hc(x[(k*m+k+1):(k*m+k+q)],x[(k*m+k+q+1):(k*m+k+2*q)])} else
if (v%%2==1){
ay2(x[(k*m+1):(k*m+k)])%*%hc(x[(k*m+k+1):(k*m+k+q)],x[(k*m+k+q+1):(k*m+k+2*q-1)])}}

t<-Y%*%par[(k*m+1):(k*m+k)]
skew<-cum(par)[2]+cum(par)[3]*t+cum(par)[4]*(1/2)*t^2+cum(par)[5]*(1/6)*t^3
kurt<-cum(par)[3]+cum(par)[4]*t+cum(par)[5]*(1/2)*t^2

# EAP
EAP<-dK(par)
VAR<-dK2(par)
cbind(round(t,3),round(EAP,3),round(VAR,3),round(skew,3),round(kurt,3))

m1<-EAP
m2<-VAR+EAP^2
m3<-skew+3*VAR*EAP+EAP^3
m4<-kurt+4*skew*EAP+3*VAR^2+6*VAR*EAP^2+EAP^4

m1B<-t(m1)%*%by(par)/sum(by(par)) # first non-central population moment
m2B<-t(m2)%*%by(par)/sum(by(par)) # second non-central population moment
m3B<-t(m3)%*%by(par)/sum(by(par)) # third non-central population moment
m4B<-t(m4)%*%by(par)/sum(by(par)) # fourth non-cnetral population moment
round(c(m1B,m2B,m3B,m4B),3)

round(m1B,3) # mean population
round(m2B-m1B^2,3) # variance population
round(m3B-3*m2B*m1B+2*m1B^3,3) # third population cumulant
round(m4B-4*m3B*m1B-3*m2B^3+12*m2B*m1B^2-6*m1B^4,3) # fourth population cumulant

round((m3B-3*m2B*m1B+2*m1B^3)/(m2B-m1B^2)^{3/2},3) # coefficient of skewness (population)
round(3+(m4B-4*m3B*m1B-3*m2B^3+12*m2B*m1B^2-6*m1B^4)/(m2B-m1B^2)^2,3) # coefficient of kurtosis (population)

par(mfrow=c(4,1),mai=c(4/5,4/5,1/3,1/2))
plot(cbind(t,dK(par))[order(t),],type='b',las=1,xlab='weighted sum of item scores',ylab='conditional latent mean')
plot(cbind(t,dK2(par))[order(t),],type='b',las=1,xlab='weighted sum of item scores',ylab='conditional latent variance')
plot(cbind(t,skew)[order(t),],type='b',las=1,xlab='weighted sum of item scores',ylab='conditional latent skewness')
plot(cbind(t,kurt)[order(t),],type='b',las=1,xlab='weighted sum of item scores',ylab='conditional latent kurtosis')

