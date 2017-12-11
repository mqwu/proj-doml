########################################################################################
 ##### call the package h2o 
 library(h2o, lib="~/lib")
 h2o.init(nthreads = -1)

 ### Data extraction 
 system("./extdat0 Dim.dat XTrestart.dat FDTrestart.dat")

 x<-scan("Dim.par")
 X_dim<-x[1]
 Y_dim<-x[2]
 miniter<-x[3]
 maxiter<-x[4]
 NES<-x[5]
 
 for(generation in miniter:maxiter){

    recentgeneration<-2
    max_sample_size<-1000
    DataFormat<-c(3, X_dim, Y_dim)
    #############################################################
    ######  data reading ########################################
    DataPast<-NULL
    if(generation-recentgeneration-1>miniter){ 
       
      for(k in miniter:(generation-recentgeneration-1)){ 
          fname<-paste("BimG",k,".dat",sep="")
          D<-read.table(fname)
          DataPast<-rbind(DataPast, D)
        }
     }

    DataRecent<-NULL 
    if(generation-recentgeneration<miniter) m2<-miniter

    if(generation>=m2 && generation>=miniter){
       for(k in m2:generation){
           fname<-paste("BimG",k,".dat",sep="")
           D<-read.table(fname)
           DataRecent<-rbind(DataRecent, D)
          }
       }

     #################################################################
     #### E2 contains recent samples to be used in prediction 
     #### prepare to draw samples from E1=DataPast   
     #################################################################

     E1<-DataPast
     E2<-DataRecent
  
     #### prepare data information file 
     if(length(E1)>0 && length(E2)>0){ 

         ss<-c(nrow(E1), nrow(E2), max_sample_size)
         d<-rbind(ss,DataFormat)
         write(t(d), ncol=3, file="temp") 
  
         write(t(E1), ncol=ncol(E1), file="old.dat")
         write(t(E2), ncol=ncol(E2), file="new.dat") 
    
         system("./maxmin3 temp old.dat new.dat combine.dat") 
     
         Dtraining<-read.table("combine.dat") 
    
         system("rm temp")
         system("rm old.dat")
         system("rm new.dat")
         system("rm combine.dat")
       }
 
     if(length(E1)==0) Dtraining<-DataRecent 
     ################################################################
 
     D.hex<-as.h2o(Dtraining)

     #### Deep neural network training for each response variable #####
     pdf(paste("UncA",generation,".pdf",sep=""))
     R<-NULL
     DERI<-NULL

     ### for(k in (sum(DataFormat[1:2])+1):sum(DataFormat[1:3])){
     for(k in (sum(DataFormat[1:2])+1):sum(DataFormat)){
         print(paste("response=V",k,sep=""))

         X_train<-Dtraining[,(DataFormat[1]+1):sum(DataFormat[1:2])]
         Y_train<-Dtraining[,k]
         YX_train<-cbind(Y_train, X_train)
         data_num<-nrow(YX_train)
         # write(t(YX_train), ncol=ncol(YX_train), file=paste("YX",k,".txt",sep=""))
    
         hidden_nodes<-c(6*X_dim, 3*X_dim, 2*X_dim)
         hidden_Act<-"Rectifier"  #### either Rectifier or Tanh

         ####### DNN training ###############
         r1<-0.75
         count<-0
         while(r1<0.95 && count<5){
       
             G0.dl<-h2o.deeplearning(x=(DataFormat[1]+1):sum(DataFormat[1:2]), y=k, 
                activation=hidden_Act, hidden=hidden_nodes, epochs=500,
                export_weights_and_biases=TRUE, training_frame=D.hex)

             zfit<-h2o.predict(G0.dl,D.hex)
             zfit<-as.matrix(zfit)
             zfit<-zfit[,1]
              
             r1<-cor(Dtraining[,k], zfit)
             count<-count+1
            }
         R<-c(R,r1)

         par(mfrow=c(1,1))
         plot(Dtraining[,k],zfit)
         #####################################

        nodes<-c(hidden_nodes,1)
        for(L in 1:length(nodes)){

            w<-h2o.weights(G0.dl, matrix_id=L)
            b<-h2o.biases(G0.dl, vector_id=L)
            w<-as.matrix(w)
            b<-as.matrix(b)
            wb<-cbind(b,w)
            write(t(wb), ncol=ncol(wb), file=paste("wb",L,".txt",sep=""))
          }

        if(hidden_Act=="LOGISTIC") act<-0
        if(hidden_Act=="Tanh") act<-1
        if(hidden_Act=="Rectifier") act<-2
        Layer<-length(nodes)+1  

        if(generation<10) fname<-paste("XBASECY000",generation,".dat", sep="")
        if(generation>=10 && generation<100) fname<-paste("XBASECY00",generation,".dat", sep="")  
        if(generation>=100 && generation<1000) fname<-paste("XBASECY0",generation,".dat", sep="")
        if(generation>=1000 && generation<10000) fname<-paste("XBASECY",generation,".dat", sep="") 
        xobs<-read.table(fname) 
        NV<-nrow(xobs)/NES 
        ND<-Y_dim

        nodes<-c(NV,ND,NES,generation,data_num,k,act,Layer,X_dim,nodes)
        write(t(nodes), ncol=1, file="nodes.inf")

        ### output the derivative information saved in the file DNNk.der
        system("./backp3 nodes.inf")

        Z<-read.table(paste("XY",k,".der",sep=""))
        DERI<-rbind(DERI,Z)
      }
    dev.off() 

    R<-cbind((sum(DataFormat[1:2])+1):sum(DataFormat[1:3]), R)
    write(t(R), ncol=ncol(R), file=paste("corrA",generation,".dat",sep=""))
 
    DERI<-DERI[order(DERI[,1]),1:ncol(DERI)] 

    if(generation<10) fname<-paste("JACOBIANCY000",generation,".dat",sep="")
    if(generation>=10 && generation<100) fname<-paste("JACOBIANCY00",generation,".dat",sep="")
    if(generation>=100 && generation<1000) fname<-paste("JACOBIANCY0",generation,".dat",sep="")
    if(generation>=1000) fname<-paste("JACOBIANCY",generation,".dat",sep="")
           
    write(t(DERI[,2:ncol(DERI)]), ncol=ncol(DERI)-1, file=fname)
    ### remove median step files ########################
 } ### end for generation variable 

    system("rm nodes.inf")
    system("rm wb*.txt")
    system("rm *.der")
    system("rm BimG*.dat")
    system("rm Dim.par")
#################################################################################################

