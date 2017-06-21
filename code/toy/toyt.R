ngrid<-50
x<-seq(0,1,len=ngrid)
y<-seq(0,1,len=ngrid)

t<-1:10

for(k in 1:10){
 
D<-NULL
Z<-NULL
for(i in 1:ngrid)
 for(j in 1:ngrid){
    
   
    u<-x[i]
    v<-y[j]
    w<-t[k]
    
    z<-2*sin(2*pi*u)*sin(2*pi*u)*sin(w)+6*cos(2*pi*v)*cos(2*pi*v)*cos(w)
    
    D<-rbind(D, c(u,v,w,z))
  }
  write(t(D), ncol=4, file=paste("toyt",k,".dat",sep=""))
 }
