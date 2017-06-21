x<-seq(0,1,len=50)
y<-seq(0,1,len=50)

Err<-NULL

pdf("function.pdf")

########### t=1 ###########################
x<-seq(0,1,len=50)
y<-seq(0,1,len=50)
z<-matrix(scan("t=1/aa.fit"), ncol=3, byrow=T)

f<-matrix(z[,2], ncol=50)
fhat<-matrix(z[,3], ncol=50)

e<-mean((f-fhat)*(f-fhat))
Err<-rbind(Err,c(1,e))

par(mfrow=c(2,2))

persp(x,y,f, main="true function with t=1")
persp(x,y,fhat, main="fitted function with t=1")

contour(x,y,f, main="true function with t=1")
contour(x,y,fhat, main="fitted function with t=1")


########### t=2 ###########################
x<-seq(0,1,len=50)
y<-seq(0,1,len=50)
z<-matrix(scan("t=2/aa.fit"), ncol=3, byrow=T)

f<-matrix(z[,2], ncol=50)
fhat<-matrix(z[,3], ncol=50)

e<-mean((f-fhat)*(f-fhat))
Err<-rbind(Err,c(2,e))

par(mfrow=c(2,2))
persp(x,y,f, main="true function with t=2")
persp(x,y,fhat, main="fitted function with t=2")

contour(x,y,f, main="true function with t=2")
contour(x,y,fhat, main="fitted function with t=2")

########### t=3 ###########################
x<-seq(0,1,len=50)
y<-seq(0,1,len=50)
z<-matrix(scan("t=3/aa.fit"), ncol=3, byrow=T)

f<-matrix(z[,2], ncol=50)
fhat<-matrix(z[,3], ncol=50)

e<-mean((f-fhat)*(f-fhat))
Err<-rbind(Err,c(3,e))

par(mfrow=c(2,2))
persp(x,y,f, main="true function with t=3")
persp(x,y,fhat, main="fitted function with t=3")

contour(x,y,f, main="true function with t=3")
contour(x,y,fhat, main="fitted function with t=3")

########### t=4 ###########################
x<-seq(0,1,len=50)
y<-seq(0,1,len=50)
z<-matrix(scan("t=4/aa.fit"), ncol=3, byrow=T)

f<-matrix(z[,2], ncol=50)
fhat<-matrix(z[,3], ncol=50)

e<-mean((f-fhat)*(f-fhat))
Err<-rbind(Err,c(4,e))

par(mfrow=c(2,2))
persp(x,y,f, main="true function with t=4")
persp(x,y,fhat, main="fitted function with t=4")

contour(x,y,f, main="true function with t=4")
contour(x,y,fhat, main="fitted function with t=4")


########### t=5 #########################
x<-seq(0,1,len=50)
y<-seq(0,1,len=50)
z<-matrix(scan("t=5/bb.fit"), ncol=3, byrow=T)

f<-matrix(z[,2], ncol=50)
fhat<-matrix(z[,3], ncol=50)

e<-mean((f-fhat)*(f-fhat))
Err<-rbind(Err,c(5,e))

par(mfrow=c(2,2))
persp(x,y,f, main="true function with t=5")
persp(x,y,fhat, main="fitted function with t=5")

contour(x,y,f, main="true function with t=5")
contour(x,y,fhat, main="fitted function with t=5")


########### t=6 #########################
x<-seq(0,1,len=50)
y<-seq(0,1,len=50)
z<-matrix(scan("t=6/bb.fit"), ncol=3, byrow=T)

f<-matrix(z[,2], ncol=50)
fhat<-matrix(z[,3], ncol=50)

e<-mean((f-fhat)*(f-fhat))
Err<-rbind(Err,c(6,e))

par(mfrow=c(2,2))
persp(x,y,f, main="true function with t=6")
persp(x,y,fhat, main="fitted function with t=6")

contour(x,y,f, main="true function with t=6")
contour(x,y,fhat, main="fitted function with t=6")


########### t=7 #########################
x<-seq(0,1,len=50)
y<-seq(0,1,len=50)
z<-matrix(scan("t=7/bb.fit"), ncol=3, byrow=T)

f<-matrix(z[,2], ncol=50)
fhat<-matrix(z[,3], ncol=50)

e<-mean((f-fhat)*(f-fhat))
Err<-rbind(Err,c(7,e))

par(mfrow=c(2,2))
persp(x,y,f, main="true function with t=7")
persp(x,y,fhat, main="fitted function with t=7")

contour(x,y,f, main="true function with t=7")
contour(x,y,fhat, main="fitted function with t=7")


########### t=8 #########################
x<-seq(0,1,len=50)
y<-seq(0,1,len=50)
z<-matrix(scan("t=8/bb.fit"), ncol=3, byrow=T)

f<-matrix(z[,2], ncol=50)
fhat<-matrix(z[,3], ncol=50)

e<-mean((f-fhat)*(f-fhat))
Err<-rbind(Err,c(8,e))

par(mfrow=c(2,2))
persp(x,y,f, main="true function with t=8")
persp(x,y,fhat, main="fitted function with t=8")

contour(x,y,f, main="true function with t=8")
contour(x,y,fhat, main="fitted function with t=8")


########### t=9 #########################
x<-seq(0,1,len=50)
y<-seq(0,1,len=50)
z<-matrix(scan("t=9/bb.fit"), ncol=3, byrow=T)

f<-matrix(z[,2], ncol=50)
fhat<-matrix(z[,3], ncol=50)

e<-mean((f-fhat)*(f-fhat))
Err<-rbind(Err,c(9,e))

par(mfrow=c(2,2))
persp(x,y,f, main="true function with t=9")
persp(x,y,fhat, main="fitted function with t=9")

contour(x,y,f, main="true function with t=9")
contour(x,y,fhat, main="fitted function with t=9")


########### t=10 #########################
x<-seq(0,1,len=50)
y<-seq(0,1,len=50)
z<-matrix(scan("t=10/bb.fit"), ncol=3, byrow=T)

f<-matrix(z[,2], ncol=50)
fhat<-matrix(z[,3], ncol=50)

e<-mean((f-fhat)*(f-fhat))
Err<-rbind(Err,c(10,e))

par(mfrow=c(2,2))
persp(x,y,f, main="true function with t=10")
persp(x,y,fhat, main="fitted function with t=10")

contour(x,y,f, main="true function with t=10")
contour(x,y,fhat, main="fitted function with t=10")

#######################################
dev.off()

