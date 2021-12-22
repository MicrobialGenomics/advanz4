get_metabolomic_table<-function(metabolome, lvl ){
  metabolome<-as.data.frame(t(metabolome))
  
  firstMetabolite<-min(which(!nchar(metabolome[1,])==0))
  metadataHeader<-min(which(!nchar(metabolome[,1])==0))
  
  # Get the upper half of the table, which contains the different names for the compound
  
  myMetaboliteNames<-metabolome[1:metadataHeader,c(firstMetabolite:ncol(metabolome))] %>% remove_rownames() %>% 
    column_to_rownames(.,colnames(metabolome)[firstMetabolite])
  
  # Get the table of values and sample data, column order should remain the same as the metabolic names data frame
  columnsToKeep<- c(1:firstMetabolite-1, (firstMetabolite+1):ncol(metabolome))
  mySampleDataValues<-metabolome[metadataHeader:nrow(metabolome),columnsToKeep] 
  # mySampleDataValues
  
  colnames(mySampleDataValues)<-mySampleDataValues[1,]
  
  colnames(mySampleDataValues)[firstMetabolite:ncol(mySampleDataValues)]<-as.character(myMetaboliteNames[lvl,])
  
  colnames(mySampleDataValues)<-sapply(colnames(mySampleDataValues), FUN = as.character, simplify = T)
  
  mySampleDataValues<-mySampleDataValues[2:nrow(mySampleDataValues),]
  
  mySampleDataValues[,firstMetabolite:ncol(mySampleDataValues)]<-apply(mySampleDataValues[,firstMetabolite:ncol(mySampleDataValues)], 2 ,function(x){as.numeric(gsub(",",".",x))})
  return(list(mySampleDataValues, colnames(mySampleDataValues)[firstMetabolite], myMetaboliteNames))
}





## From https://www.r-bloggers.com/identify-describe-plot-and-remove-the-outliers-from-the-dataset/



myPlotPCAmixOmics<-function(mixOmicsPCAResult,variable,myPalette="Darjeeling1",paletteSeed=1234,title="Lorem Ipsum...",
                            file=paste (Sys.time(), ".pdf", sep=""),fillColor=T,lineColor=F ){
  ##TO-DO add p-values in boxplots.
  require(grid)
  require(gridExtra)
  require(ggplot2)
  require(ggthemes)
  require(ggpubr)
  require(RColorBrewer)
  require(wesanderson)
  require(mixOmics)
  class(variable)
  if(! is.data.frame(variable)){
    #If a vector is passed for a response variable we are going to assume it follows the same order as the data.frame
    variable<-data.frame(variable)
    rownames(variable)<-mixOmicsPCAResult$names$sample
  }
  variableVector<-factor(variable$variable,levels=levels(variable$variable),ordered=T)
  if(myPalette %in% names(wes_palettes)){
   availableColors<-length(wes_palette(name=myPalette))
   if(nlevels(variable) > availableColors){
     print ("Not enough colors in palette for variable ")
   }
   myColors<-wes_palette(n=max(4,nlevels(variableVector)),name=myPalette)
  }else{
    myColors<-brewer.pal(n=max(4,nlevels(variableVector)),name=myPalette)
  }
  set.seed(paletteSeed)
  myColors<-sample(myColors,nlevels(variableVector))

  # if(is.numeric(variableVector)){
  #   require(leaflet)
  #   myColors<-sample(myColors,2)
  #   mypal<-colorRampPalette(c(mycolors),length(variableVector))
  # }
 
  set.seed(paletteSeed)
  
  print("Plotting PCA 1&2 on device")
  p1<-plotIndiv(mixOmicsPCAResult, comp = c(1, 2), ind.names = F, 
                group = variableVector, 
                legend = F, ellipse=T,title=paste(title,"- PC 1&2"),col.per.group = myColors)
  # p1$graph<-p1$graph+theme_few()
  p1$graph<-p1$graph+theme(legend.position="none")
  print("Plotting PCA 1&3 on device")
  p2<-plotIndiv(mixOmicsPCAResult, comp = c(1, 3), ind.names = F, 
                group = variableVector, 
                legend = F, ellipse=T,title=paste(title, "- PC 1&3"),col.per.group = myColors)
  # p2$graph<-p2$graph+theme_few()
  p2$graph<-p2$graph+theme(legend.position="none")
  print("Plotting PCA 2&3 on device")
  p3<-plotIndiv(mixOmicsPCAResult, comp = c(2, 3), ind.names = F, 
                group = variableVector, 
                legend = F, ellipse=T,title=paste(title, "- PC 2&3"),col.per.group = myColors)
  # p3$graph<-p3$graph+theme_few()
  p3$graph<-p3$graph+theme(legend.position="none")
  print("merging PCA plots")
  mixOmicsPCAResult.b<-merge(mixOmicsPCAResult$variates,variable,by="row.names")
  mixOmicsPCAResult.b<-mixOmicsPCAResult.b[,-1]
  colnames(mixOmicsPCAResult.b)<-c("PC1","PC2","PC3","variable")
  print("Plotting PCA comp 1,2,3 boxplots")
  
  if(fillColor==T & lineColor==F){
    p4<-ggplot(mixOmicsPCAResult.b,aes(variable,PC1,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)+ theme(axis.text.x = element_text(angle = 45, size=8))
    p5<-ggplot(mixOmicsPCAResult.b,aes(variable,PC2,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)+ theme(axis.text.x = element_text(angle = 45, size=8))
    p6<-ggplot(mixOmicsPCAResult.b,aes(variable,PC3,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)+ theme(axis.text.x = element_text(angle = 45, size=8))
  }else if(fillColor==F & lineColor==T){
    p4<-ggplot(mixOmicsPCAResult.b,aes(variable,PC1,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)+theme(axis.text.x = element_text(angle = 45, size=8))
    p5<-ggplot(mixOmicsPCAResult.b,aes(variable,PC2,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)+theme(axis.text.x = element_text(angle = 45, size=8))
    p6<-ggplot(mixOmicsPCAResult.b,aes(variable,PC3,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)+theme(axis.text.x = element_text(angle = 45, size=8))
  }
  print("Arranging plots")
  pT<-ggarrange(p1$graph, p2$graph, p3$graph, ggarrange(p4, p5, p6 )+ rremove("x.text"), 
                labels = c("A", "B", "C"),
                ncol = 2, nrow = 2)
  print(paste("Saving image to ",file,sep=""))
  ggsave(file=file,height = 210, width = 297,dpi=1200, units="mm", pT,)
  
  
  return(pT)
}

myPlotiPCAmixOmics<-function(mixOmicsiPCAResult,variable,myPalette="Darjeeling1",paletteSeed=1234,title="Lorem Ipsum...",file=paste (Sys.time(), ".pdf", sep="")){
  ##TO-DO add p-values in boxplots.
  require(grid)
  require(gridExtra)
  require(ggplot2)
  require(ggthemes)
  require(ggpubr)
  require(RColorBrewer)
  require(wesanderson)
  if(myPalette %in%  c("GrandBudapest","Moonrise1","Moonrise2","Royal1","GrandBudapest2","Chevalier")){
    myColors<-wes_palette(n=4, name=myPalette)
  } else if (myPalette %in% c("Royal2","GrandBudapest2","Cavalcanti","Moonrise3","Zissou","Darjeeling","RushMore","FantasticFox")){
    myColors<-wes_palette(n=5, name=myPalette)
  } else{
    myColors<-brewer.pal(8,myPalette)
  }
  set.seed(paletteSeed)
  myColors<-sample(myColors,nlevels(variable))
  p1<-plotIndiv(mixOmicsiPCAResult, comp = c(1, 2), ind.names = F, 
                group = variable, 
                legend = F, ellipse=T,title=paste(title,"- PC 1&2"),col=myColors)
  p1$graph<-p1$graph+theme_few()
  p1$graph<-p1$graph+theme(legend.position="none")
  
  p2<-plotIndiv(mixOmicsiPCAResult, comp = c(1, 3), ind.names = F, 
                group = variable, 
                legend = F, ellipse=T,title=paste(title, "/PC 2&3"),col=myColors)
  p2$graph<-p2$graph+theme_few()
  p2$graph<-p2$graph+theme(legend.position="none")
  
  p3<-plotIndiv(mixOmicsiPCAResult, comp = c(2, 3), ind.names = F, 
                group = variable, 
                legend = F, ellipse=T,title=paste(title, "/PC 1&3"),col=myColors)
  p3$graph<-p3$graph+theme_few()
  p3$graph<-p3$graph+theme(legend.position="none")
  
  mixOmicsiPCAResult.b<-merge(mixOmicsiPCAResult$x,Y,by="row.names")
  p4<-ggplot(mixOmicsiPCAResult.b,aes(variable,IPC1,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_fill_manual(values=myColors)+scale_color_manual(values=myColors)
  p5<-ggplot(mixOmicsiPCAResult.b,aes(variable,IPC2,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_fill_manual(values=myColors)+scale_color_manual(values=myColors)
  p6<-ggplot(mixOmicsiPCAResult.b,aes(variable,IPC3,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_fill_manual(values=myColors)+scale_color_manual(values=myColors)
  
  pT<-ggarrange(p1$graph, p2$graph, p3$graph, ggarrange(p4, p5, p6 )+ rremove("x.text"), 
                labels = c("A", "B", "C"),
                ncol = 2, nrow = 2)
  ggsave(file=file,height = 210, width = 297,dpi=1200, units="mm", pT,)
  return(pT)
}

myPlotsiPCAmixOmics<-function(mixOmicssiPCAResult,variable,myPalette="Darjeeling1",paletteSeed=1234,title="Lorem Ipsum...",file=paste (Sys.time(), ".pdf", sep="")){
  ##TO-DO add p-values in boxplots.
  require(grid)
  require(gridExtra)
  require(ggplot2)
  require(ggthemes)
  require(ggpubr)
  require(RColorBrewer)
  require(wesanderson)
  colnames(mixOmicssiPCAResult$x)<-c("sIPC1","sIPC2","sIPC3")
  if(myPalette %in%  c("GrandBudapest","Moonrise1","Moonrise2","Royal1","GrandBudapest2","Chevalier")){
    myColors<-wes_palette(n=4, name=myPalette)
  } else if (myPalette %in% c("Royal2","GrandBudapest2","Cavalcanti","Moonrise3","Zissou","Darjeeling","RushMore","FantasticFox")){
    myColors<-wes_palette(n=5, name=myPalette)
  } else{
    myColors<-brewer.pal(8,myPalette)
  }
  set.seed(paletteSeed)
  myColors<-sample(myColors,nlevels(variable))
  p1<-plotIndiv(mixOmicssiPCAResult, comp = c(1, 2), ind.names = F, 
                group = variable, 
                legend = F, ellipse=T,title=paste(title,"- PC 1&2"),col=myColors)
  p1$graph<-p1$graph+theme_few()
  p1$graph<-p1$graph+theme(legend.position="none")
  
  p2<-plotIndiv(mixOmicssiPCAResult, comp = c(1, 3), ind.names = F, 
                group = variable, 
                legend = F, ellipse=T,title=paste(title, "/PC 2&3"),col=myColors)
  p2$graph<-p2$graph+theme_few()
  p2$graph<-p2$graph+theme(legend.position="none")
  
  p3<-plotIndiv(mixOmicssiPCAResult, comp = c(2, 3), ind.names = F, 
                group = variable, 
                legend = F, ellipse=T,title=paste(title, "/PC 1&3"),col=myColors)
  p3$graph<-p3$graph+theme_few()
  p3$graph<-p3$graph+theme(legend.position="none")
  
  mixOmicssiPCAResult.b<-merge(mixOmicssiPCAResult$x,Y,by="row.names")
  p4<-ggplot(mixOmicssiPCAResult.b,aes(variable,sIPC1,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_fill_manual(values=myColors)+scale_color_manual(values=myColors)
  p5<-ggplot(mixOmicssiPCAResult.b,aes(variable,sIPC2,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_fill_manual(values=myColors)+scale_color_manual(values=myColors)
  p6<-ggplot(mixOmicssiPCAResult.b,aes(variable,sIPC3,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_fill_manual(values=myColors)+scale_color_manual(values=myColors)
  
  pT<-ggarrange(p1$graph, p2$graph, p3$graph, ggarrange(p4, p5, p6 )+ rremove("x.text"), 
                labels = c("A", "B", "C"),
                ncol = 2, nrow = 2)
  ggsave(file=file,height = 210, width = 297,dpi=1200, units="mm", pT,)
  return(pT)
}

myBoxCoxTransform<-function(matrix, standardize=F){
  require(AID)
  auxDataMatrix<-matrix
  print("Running univariate col-wise BoxCox Power Transformation")
  for (column in 1:ncol(auxDataMatrix)){
    myBC<-boxcoxnc(matrix[,column],lambda2=0.00000001,verbose=F,plot=F)
    auxDataMatrix[,column]<-myBC$tf.data
  }
  print(paste("Standardizing Box-Cox Transformed Data"))
  if(standardize==T){
  auxDataMatrixStd<-vegan::decostand(auxDataMatrix,method="standardize")
  return(auxDataMatrixStd)
  }
  return(auxDataMatrix)
}
### Returns mixOmicsPLSDA object
myRunPLSDAmixOmics<-function(DataMatrix,variable,variableName="Variable",myPalette="Darjeeling1",dirForStats=NULL,
                             paletteSeed=1234,title="Lorem Ipsum...",file=paste (Sys.time(), ".pdf", sep=""),
                             fillColor=T,lineColor=F,performance=T,transformMethod=NULL,plot=F,varPlot=T,threshold=0.9,scale=T){
  require(mixOmics) 
  require(phyloseq)
  require(AID)
  require(gdata)
  require(stringr)
  require(vegan)
  myData<-DataMatrix
  if(transformMethod=="BoxCox"){
    auxDataMatrix<-DataMatrix
    print("Running univariate col-wise BoxCox Power Transformation")
    for (column in 1:ncol(DataMatrix)){
      myBC<-boxcoxnc(DataMatrix[,column],lambda2=0.00000001,verbose=F,plot=F)
      auxDataMatrix[,column]<-myBC$tf.data
    }
    myData<-auxDataM
  } else if(transformMethod=="StandardizedBoxCox"){
    auxDataMatrix<-DataMatrix
    print("Running univariate col-wise BoxCox Power Transformation")
    for (column in 1:ncol(DataMatrix)){
      myBC<-boxcoxnc(DataMatrix[,column],lambda2=0.00000001,verbose=F,plot=F)
      auxDataMatrix[,column]<-myBC$tf.data
    }
    print(paste("Standardizing Box-Cox Transformed Data"))
    auxDataMatrixStd<-vegan::decostand(auxDataMatrix,method="standardize")
    myData<-auxDataMatrixStd
  }else if(transformMethod=="Standardize"){
    print("Running Data Standarization")
    myData<-vegan::decostand(DataMatrix,method="standardize")
  }else if(transformMethod=="None" | transformMethod=="none"){
    transformMethod<-"None"
  }else if(! is.null(transformMethod ) ){
    print(paste("Runnning ",transformMethod," on data",sep=""))
    myData<-vegan::decostand(DataMatrix,method=transformMethod)
  }else{
    transformMethod<-"None"
  }
  if(! is.data.frame(variable)){
    #If a vector is passed for a response variable we are going to assume it follows the same order as the data.frame
    variable<-data.frame(variable)
    rownames(variable)<-rownames(DataMatrix)
  }
  variableVector<-factor(variable$variable,levels=levels(variable$variable),ordered=T)
  
  ##### Start PLS-DA Analysis 
  myResult <- plsda( myData, variableVector, ncomp = 15)
  # plotContrib(result)
  # plotVar(myResult)
  set.seed(paletteSeed) # for reproducibility here, only when the `cpus' argument is not used
  print(paste("Running PLS-DA Tuning using MFold Validation, 10 folds, 10 repeats"))
  perf.plsda <- perf(myResult, validation = "Mfold", folds = 10, 
                     progressBar = T, auc = TRUE, nrepeat = 10,scale=scale) 
  # myPlot<-plot(perf.plsda)
  print(paste("Saving Performance Plot in PLSDAPerformance_",transformMethod,"_",variableName,".pdf",sep=""))
  ### Keerp Performance
  # ggsave(file=paste(dirForStats,"/PLSDAPerformance_",transformMethod,"_",variableName,".pdf",sep=""),myPlot)
  pdf(paste(dirForStats,"/PLSDAPerformance_",transformMethod,"_",variableName,".pdf",sep=""))
   plot(perf.plsda)
   dev.off()
  myncomp<-max(3,perf.plsda$choice.ncomp["BER","mahalanobis.dist"])
  if(is.null(dirForStats)){
    dirForStats<-getwd()
  }
  write.csv(file=paste(dirForStats,"/PLSDAPerformance_",transformMethod,"_",variableName,".csv",sep=""),data.frame(perf.plsda$choice.ncomp))
  write.csv(file=paste(dirForStats,"/PLSDAErrorRate_",transformMethod,"_",variableName,".csv",sep=""),data.frame(perf.plsda$error.rate))
  write.csv(file=paste(dirForStats,"/PLSDAErrorRateAll_",transformMethod,"_",variableName,".csv",sep=""),data.frame(perf.plsda$error.rate.all))
  print(paste("Keeping ",myncomp," components"))
  myResult <- plsda( myData, variableVector, ncomp = myncomp,scale=scale)
  if(isTRUE(plot)){
  myPlotPLSDAmixOmicsWithPerformance(myResult,variableVector,title=title,file=file,
                                    myPalette=myPalette,paletteSeed=paletteSeed,fillColor=T,lineColor=F)
  }
  
  if(isTRUE(varPlot)){
    myLoadings<-data.frame(myResult$loadings$X)
    colnames(myLoadings)<-c("x","y","z")
    # myLoadings<-myLoadings[! grepl("^-X",rownames(myLoadings)),]
    myLoadings<-t(myLoadings)
    # colnames(myLoadings)<-sub("-X.*","",colnames(myLoadings))
    myLoadings<-data.frame(t(myLoadings))
    myPlot<-myPlotVarFunction(myLoadings[,c("x","y")],myTitle=paste("PLSDA vs ",variableName," ",transformMethod," varPlot",sep=""),size=4,threshold=NULL)
    # myPlot
    ggsave(paste(dirForStats,"/PLSDA_",variableName,"_All-Data_",transformMethod,"_Loadings_Comp1_2.pdf",sep=""),myPlot)
  }
  print(paste("Keeping Loadings in Csv file"))
  write.csv(file=paste(dirForStats,"/PLSDALoadings_",transformMethod,"_",variableName,".csv",sep=""),data.frame(myResult$loadings$X))

  return(myResult)
}

myRunsPLSDAmixOmics<-function(DataMatrix,variable,variableName="Variable",myPalette="Darjeeling1",
                             paletteSeed=1234,title="Lorem Ipsum...",file=paste (Sys.time(), ".pdf", sep=""),
                             fillColor=T,lineColor=F,performance=T,transformMethod=NULL,plot=F,varPlot=T,threshold=0.9,
                             list.keepX=NULL,ncomp=NULL,ncompTest=5,plsMethod="classic",scale=T){
  require(mixOmics) 
  require(phyloseq)
  require(AID)
  require(gdata)
  require(stringr)
  require(vegan)
  myData<-DataMatrix
  if(transformMethod=="BoxCox"){
    auxDataMatrix<-DataMatrix
    print("Running univariate col-wise BoxCox Power Transformation")
    for (column in 1:ncol(X)){
      myBC<-boxcoxnc(DataMatrix[,column],lambda2=0.00000001,verbose=F,plot=F)
      auxDataMatrix[,column]<-myBC$tf.data
    }
    myData<-auxDataM
  } else if(transformMethod=="StandardizedBoxCox"){
    auxDataMatrix<-DataMatrix
    print("Running univariate col-wise BoxCox Power Transformation")
    for (column in 1:ncol(X)){
      myBC<-boxcoxnc(DataMatrix[,column],lambda2=0.00000001,verbose=F,plot=F)
      auxDataMatrix[,column]<-myBC$tf.data
    }
    print(paste("Standardizing Box-Cox Transformed Data"))
    auxDataMatrixStd<-vegan::decostand(auxDataMatrix,method="standardize")
    myData<-auxDataMatrixStd
  }else if(transformMethod=="Standardize"){
    myData<-vegan::decostand(DataMatrix,method="standardize")
  }else if(transformMethod=="None"| transformMethod=="none"){
    transformMethod<-"None"
  }else if(! is.null(transformMethod ) ){
    myData<-vegan::decostand(DataMatrix,method=transformMethod)
  }else{
    transformMethod<-"None"
  }
  if(! is.data.frame(variable)){
    #If a vector is passed for a response variable we are going to assume it follows the same order as the data.frame
    variable<-data.frame(variable)
    rownames(variable)<-rownames(DataMatrix)
  }
  variableVector<-factor(variable$variable,levels=levels(variable$variable),ordered=T)
  if(is.null(list.keepX)){
    # list.keepX <- c(1:5, seq(10, 18, 2), seq(20,50,5))
    list.keepX<-c(1:20)
  } 
  if(ncompTest <= 2){
    print(paste("Ncomp to test is 2 or less, using a minimum of 3 components for tuning sPLS-DA model"))
    ncompTest<-3
  }
  print(paste("Running sPLS-DA Tuning using MFold Validation, 5 folds, 10 repeats, on ",ncompTest,"components",sep=""))
  tune.splsda.srbct <- tune.splsda(myData, variableVector, ncomp = ncompTest, validation = 'Mfold', folds = 10, 
                                   progressBar = T, dist = 'mahalanobis', auc=T, measure="BER",
                                   test.keepX = list.keepX, nrepeat = 10,scale = scale) #nrepeat 50-100
  print(paste("Saving Performance Plot in sPLSDATuning_",transformMethod,"_",variableName,".pdf",sep=""))
  # pdf(paste(variableName,"/sPLSDATuning_",transformMethod,"_",variableName,".pdf",sep=""))
 ggsave(paste(variableName,"/sPLSDATuning_",transformMethod,"_",variableName,".pdf",sep=""), plot(tune.splsda.srbct))
  # dev.off()
  if(is.null(ncomp)){
    ncomp<-max(3,tune.splsda.srbct$choice.ncomp$ncomp)
  }
  print(paste("Using",ncomp,"and",tune.splsda.srbct$choice.keepX,"variables"))
 
  myResult <- splsda( myData, variableVector, ncomp = ncomp ,mode=plsMethod,
                      keepX=tune.splsda.srbct$choice.keepX,scale=scale)
  perf.splsda<-perf(myResult, validation = "Mfold", folds = 10, 
       progressBar = F, auc = TRUE, nrepeat = 10,scale=scale) 
  write.csv(file=paste(variableName,"/sPLSDAPerformance_",transformMethod,"_",variableName,".csv",sep=""),data.frame(perf.splsda$choice.ncomp))
  write.csv(file=paste(variableName,"/sPLSDAErrorRate_",transformMethod,"_",variableName,".csv",sep=""),data.frame(perf.splsda$error.rate))
  write.csv(file=paste(variableName,"/sPLSDAErrorRateAll_",transformMethod,"_",variableName,".csv",sep=""),data.frame(perf.splsda$error.rate.all))
   # # Do we want to plot the resulting ordination and performance plots?
   # if(isTRUE(plot)){
   #  myPlotPLSDAmixOmicsWithPerformance(myResult,variable=variableVector,title=title,file=file,
   #                                     myPalette=myPalette,paletteSeed=paletteSeed,fillColor=T,lineColor=F)
   # }
  ### Do we need to plot loadings?
  if(isTRUE(varPlot)){
    myLoadings<-data.frame(myResult$loadings$X)
    colnames(myLoadings)<-c("x","y","z")
    # myLoadings<-myLoadings[! grepl("^-X",rownames(myLoadings)),]
    myLoadings<-t(myLoadings)
    # colnames(myLoadings)<-sub("-X.*","",colnames(myLoadings))
    myLoadings<-data.frame(t(myLoadings))
    myPlot<-myPlotVarFunction(myLoadings[!rowSums(myLoadings[,c("x","y")])==0,c("x","y")],myTitle=paste("sPLSDA vs ",variableName," ",transformMethod," varPlot",sep=""),size=4,threshold=0)
    # myPlot
    ggsave(paste(variableName,"/sPLSDA_",variableName,"_All-Data_",transformMethod,"_Loadings_Comp1_2.pdf",sep=""),myPlot)
  }
  print(paste("Keeping Loadings in Csv file"))
  write.csv(file=paste(variableName,"/sPLSDALoadings_",transformMethod,"_",variableName,".csv",sep=""),data.frame(myResult$loadings$X))
  return(myResult)
}


myPlotPLSDAmixOmics<-function(mixOmicsPLSDAResult,variable,myPalette="Darjeeling1",paletteSeed=1234,title="Lorem Ipsum...",file=paste (Sys.time(), ".pdf", sep=""),fillColor=T,lineColor=F){
  ##TO-DO add p-values in boxplots.
  require(grid)
  require(gridExtra)
  require(ggplot2)
  require(ggthemes)
  require(ggpubr)
  require(RColorBrewer)
  require(wesanderson)
  require(mixOmics)
  if(! is.data.frame(variable)){
    #If a vector is passed for a response variable we are going to assume it follows the same order as the data.frame
    variable<-data.frame(variable)
    rownames(variable)<-mixOmicsPLSDAResult$names$sample
  }
  variableVector<-factor(variable$variable,levels=levels(variable$variable),ordered=T)
  print(levels(variableVector))
  if(myPalette %in% names(wes_palettes)){
    availableColors<-length(wes_palette(name=myPalette))
    if(nlevels(variable) > availableColors){
      print ("Not enough colors in palette for variable ")
    }
    myColors<-wes_palette(n=max(4,nlevels(variableVector)),name=myPalette)
  }else{
    myColors<-brewer.pal(n=max(4,nlevels(variableVector)),name=myPalette)
  }
  set.seed(paletteSeed)
  myColors<-sample(myColors,nlevels(variable))
  p1<-plotIndiv(mixOmicsPLSDAResult, comp = c(1, 2), ind.names = F, 
                group = variableVector, 
                legend = F, ellipse=T,title=paste(title,"- PC 1&2"),col=myColors,plot=F)
  # p1$graph<-p1$graph+theme_few()
  p1$graph<-p1$graph+theme(legend.position="none")
  if ( mixOmicsPLSDAResult$ncomp>2){
    p2<-plotIndiv(mixOmicsPLSDAResult, comp = c(1, 3), ind.names = F, 
                  group = variableVector, 
                  legend = F, ellipse=T,title=paste(title, "- PC 1&3"),col=myColors,plot=F)
    # p2$graph<-p2$graph+theme_few()
    p2$graph<-p2$graph+theme(legend.position="none")
    
    p3<-plotIndiv(mixOmicsPLSDAResult, comp = c(2, 3), ind.names = F, 
                  group = variableVector, 
                  legend = F, ellipse=T,title=paste(title, "- PC 2&3"),col=myColors,plot=F)
    # p3$graph<-p3$graph+theme_few()
    p3$graph<-p3$graph+theme(legend.position="none")
  }else{

    df <- data.frame()
    p2<-ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
    p3<-ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
    p2$graph<-p2$graph+theme(legend.position="none")
    p3$graph<-p3$graph+theme(legend.position="none")
  }
  mixOmicsPLSDAResult.b<-mixOmicsPLSDAResult$variates$X
  if ( mixOmicsPLSDAResult$ncomp==2){
    colnames(mixOmicsPLSDAResult.b)<-c("Comp1","Comp2")
  }else {
    colnames(mixOmicsPLSDAResult.b)<-c("Comp1","Comp2","Comp3")
  }
  mixOmicsPLSDAResult.b<-as.data.frame(mixOmicsPLSDAResult.b)
  mixOmicsPLSDAResult.b$variable<-mixOmicsPLSDAResult$Y
  print(myColors)
  if(fillColor==T & lineColor==F){
    p4<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp1,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)+ theme(axis.text.x = element_text(angle = 45, size=8))
    p5<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp2,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)+ theme(axis.text.x = element_text(angle = 45, size=8))
    if( mixOmicsPLSDAResult$ncomp==2){
      p6<-ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
      p6$graph<-p6$graph+theme(legend.position="none")    
    }else{
      p6<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp3,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)+ theme(axis.text.x = element_text(angle = 45, size=8))
    }
  }else if(fillColor==F & lineColor==T){
    p4<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp1,col=variable))+geom_boxplot()+theme(legend.position="none")+scale_color_manual(values=myColors)+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)+ theme(axis.text.x = element_text(angle = 45, size=8))
    p5<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp2,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)+ theme(axis.text.x = element_text(angle = 45, size=8))
    if( mixOmicsPLSDAResult$ncomp==2){
      p6<-ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
      p6$graph<-p6$graph+theme(legend.position="none")   
      p6$graph<-p6$graph+theme(legend.position="none")
    }else{
      p6<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp3,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)+ theme(axis.text.x = element_text(angle = 45, size=8))
    }
  }
  # p4<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp1,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_fill_manual(values=myColors)+scale_color_manual(values=myColors)
  # p5<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp2,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_fill_manual(values=myColors)+scale_color_manual(values=myColors)
  # # if ( nlevels(variableVector)>=2){
  #   p6<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp3,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_fill_manual(values=myColors)+scale_color_manual(values=myColors)
  # }
  if ( nlevels(variableVector)>=2){
    pT<-ggarrange(p1$graph, p3$graph, p2$graph, ggarrange(p4, p5, p6 ), 
                  labels = c("A", "B", "C"),
                  ncol = 2, nrow = 2)
    ggsave(file=file,height = 210, width = 297,dpi=1200, units="mm", pT)
  }else{
    pT<-ggarrange(p1$graph,  ggarrange(p4, p5 ),
                  labels = c("A", "B", "C"),
                  ncol = 2, nrow = 2)
    ggsave(file=file,height = 210, width = 297,dpi=1200, units="mm", pT)
  }
  
  return(pT)
}

myPlotPLSDAmixOmicsWithPerformance<-function(mixOmicsPLSDAResult,variable,myPalette="Darjeeling1",paletteSeed=1234,title="Lorem Ipsum...",file=paste (Sys.time(), ".pdf", sep=""),fillColor=T,lineColor=F){
  ##TO-DO add p-values in boxplots.
  require(grid)
  require(gridExtra)
  require(ggplot2)
  require(ggthemes)
  require(ggpubr)
  require(RColorBrewer)
  require(wesanderson)
  require(mixOmics)
  if(! is.data.frame(variable)){
    #If a vector is passed for a response variable we are going to assume it follows the same order as the data.frame
    variable<-data.frame(variable)
    rownames(variable)<-mixOmicsPLSDAResult$names$sample
  }
  variableVector<-factor(variable$variable,levels=levels(variable$variable),ordered=T)
  
  if(myPalette %in% names(wes_palettes)){
    availableColors<-length(wes_palette(name=myPalette))
    if(nlevels(variableVector) > availableColors){
      print ("Not enough colors in palette for variable ")
    }
    myColors<-wes_palette(n=max(4,nlevels(variableVector)),name=myPalette)
  }else{
    # myColors<-brewer.pal(nlevels(variableVector),myPalette)[1:nlevels(variableVector)]
    myColors<-brewer.pal(n=max(4,nlevels(variableVector)),name=myPalette)
  }
  set.seed(paletteSeed)
  myColors<-sample(myColors,nlevels(variableVector))
  print(myColors)
  p1<-plotIndiv(mixOmicsPLSDAResult, comp = c(1, 2), ind.names = F, 
                group = variableVector, 
                legend = F, ellipse=T,title=paste(title,"- PC 1&2"),col=myColors)
  # p1$graph<-p1$graph+theme_few()
  p1$graph<-p1$graph+theme(legend.position="none")
  if ( nlevels(variableVector)>=2){
    p2<-plotIndiv(mixOmicsPLSDAResult, comp = c(1, 3), ind.names = F, 
                  group = variableVector, 
                  legend = F, ellipse=T,title=paste(title, "- PC 1&3"),col=myColors)
    # p2$graph<-p2$graph+theme_few()
    p2$graph<-p2$graph+theme(legend.position="none")
    
    p3<-plotIndiv(mixOmicsPLSDAResult, comp = c(2, 3), ind.names = F, 
                  group = variableVector, 
                  legend = F, ellipse=T,title=paste(title, "- PC 2&3"),col=myColors)
    # p3$graph<-p3$graph+theme_few()
    p3$graph<-p3$graph+theme(legend.position="none")
  }
  mixOmicsPLSDAResult.b<-mixOmicsPLSDAResult$variates$X
  if ( nlevels(variableVector)>=2){
    colnames(mixOmicsPLSDAResult.b)[1:3]<-c("Comp1","Comp2","Comp3")
  }else{
    colnames(mixOmicsPLSDAResult.b)[1:3]<-c("Comp1","Comp2","Comp3")
  }
  mixOmicsPLSDAResult.b<-as.data.frame(mixOmicsPLSDAResult.b)
  mixOmicsPLSDAResult.b$variable<-mixOmicsPLSDAResult$Y
  print(myColors)
  if(fillColor==T & lineColor==F){
    p4<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp1,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)
    p5<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp2,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)
    p6<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp3,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)
  }else if(fillColor==F & lineColor==T){
    p4<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp1,col=variable))+geom_boxplot()+theme(legend.position="none")+scale_color_manual(values=myColors)+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)
    p5<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp2,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)
    p6<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp3,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_color_manual(values=myColors)+scale_color_manual(values=myColors)+stat_compare_means(label.x=1.3)
  }
  # p4<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp1,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_fill_manual(values=myColors)+scale_color_manual(values=myColors)
  # p5<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp2,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_fill_manual(values=myColors)+scale_color_manual(values=myColors)
  # # if ( nlevels(variableVector)>=2){
  #   p6<-ggplot(mixOmicsPLSDAResult.b,aes(variable,Comp3,col=variable))+geom_boxplot()+theme_few()+theme(legend.position="none")+scale_fill_manual(values=myColors)+scale_color_manual(values=myColors)
  # }
  if ( nlevels(variableVector)>=2){
    pT<-ggarrange(p1$graph, p2$graph, p3$graph, ggarrange(p4, p5, p6 ), 
                  labels = c("A", "B", "C"),
                  ncol = 2, nrow = 2)
    ggsave(file=file,height = 210, width = 297,dpi=1200, units="mm", pT)
  }else{
    pT<-ggarrange(p1$graph,  ggarrange(p4, p5 ),
                  labels = c("A", "B", "C"),
                  ncol = 2, nrow = 2)
    ggsave(file=file,height = 210, width = 297,dpi=1200, units="mm", pT)
  }
  print(colnames(loadings(mixOmicsPLSDAResult)$X))
  Comp1Loadings<-as.data.frame(loadings(mixOmicsPLSDAResult)$X[,"comp1"][as.vector(loadings(mixOmicsPLSDAResult)$X[,"comp1"]!=0)])
  Comp2Loadings<-as.data.frame(loadings(mixOmicsPLSDAResult)$X[,"comp2"][as.vector(loadings(mixOmicsPLSDAResult)$X[,"comp2"]!=0)])
  Comp3Loadings<-as.data.frame(loadings(mixOmicsPLSDAResult)$X[,"comp3"][as.vector(loadings(mixOmicsPLSDAResult)$X[,"comp3"]!=0)])
  perfFile<-paste(gsub("\\..*","",file),"_Perf.pdf",sep="")
  
  p10<-(auroc(mixOmicsPLSDAResult,roc.comp=1))
  # ggsave(file="prova_p10.pdf",height=210,width=297,dpi=1200,units="mm")
  p11<-(auroc(mixOmicsPLSDAResult,roc.comp=2))
  # ggsave(file="prova_p11.pdf",height=210,width=297,dpi=1200,units="mm")
  p12<-(auroc(mixOmicsPLSDAResult,roc.comp=3))
  # ggsave(file="prova_p12.pdf",height=210,width=297,dpi=1200,units="mm")
  
  # print(p10)
  pT2<-ggarrange(p10$graph.Comp1,p11$graph.Comp2,p12$graph.Comp3,ncol=2,nrow=2, labels=c("Comp1","Comp2","Comp3"))
  ggsave(file=perfFile,height=210,width=297,dpi=1200,units="mm",pT2)
  return(pT)
}

myPlotVarFunction<-function(PC,myTitle=NULL,file=NULL,threshold=NULL,percentile=0.9,size=3){
  require(ggrepel)
  require(ggforce)
  maxCoord<-max(abs(c(PC$x,PC$y)))
  PC$label<-NULL
  if(is.null(threshold)){
    myPercentile<-quantile(sqrt((PC$x)*(PC$x)+(PC$y*PC$y)),p=percentile)
  }else{
    myPercentile<-threshold
  }
  myPlottableLabels<-sqrt((PC$x)*(PC$x)+(PC$y*PC$y))>myPercentile
  PC[myPlottableLabels,"label"]<-rownames(PC[myPlottableLabels,])
  # maxCoord<-max(abs(max(PC$x),min(PC$x),abs(max(PC$y)),abs(min(PC$y)))=
  myPlot<-ggplot(PC,aes(x=x,y=y))+ geom_point()+
    # xlim(-max(abs(PC$x)),max(abs(PC$x)))+
    # ylim(-max(abs(PC$y)),max(abs(PC$y)))
    xlim(-maxCoord,maxCoord)+ylim(-maxCoord,maxCoord)
  myPlot<-myPlot+geom_hline(aes(0), size=0.4,yintercept=0,col="darkred")
  myPlot<-myPlot+geom_vline(aes(0), size=0.4, xintercept=0,col="darkred")
  # myPlot<-myPlot+geom_text(alpha=0.8,size=9,aes(label=PC$names))
  myPlot<-myPlot+geom_segment(data=PC[myPlottableLabels,],aes(x=0,y=0,xend=x,yend=y),arrow=arrow(length=unit(0.2,"cm")),alpha=0.5,color="blue")
  myPlot<-myPlot+geom_text_repel(aes(label = PC$label),size=4,box.padding=0)
  myPlot<-myPlot+ ggplot2::annotate("path", x=0+myPercentile*cos(seq(0,2*pi,length.out=100)), y=0+myPercentile*sin(seq(0,2*pi,length.out=100)))
  myPlot<-myPlot+ggtitle(myTitle)
  return(myPlot)
}

myCollapseCorrelatingMetabolites<-function(df,threshold=0.99){
  require(Hmisc)
  myCor<-cor(df,method="spearman")
  require(gplots)
  # heatmap.2(myCor,trace="none")
  if(max(as.vector(myCor))< threshold){
    return(as.matrix(df))
  }
  myIndex<-which(myCor >=threshold, arr.ind = TRUE)
  myIndex<-myIndex[! (myIndex[,1]==myIndex[,2]),]
  if(nrow(myIndex) == 0){
    return(df)
  }
  myIndex<-cbind(myIndex,0)
  for(i in 1:nrow(myIndex)){
    myIndex[i,3]<-myCor[myIndex[i,"row"],myIndex[i,"col"]]
  }
  colnames(myIndex)<-c("row","col","value")
  myIndex<-myIndex[order(myIndex[,"value"],decreasing=T),]
  myRow<-myIndex[1,"row"]
  myCol<-myIndex[1,"col"]
  myModel<-lm(df[,myRow]~df[,myCol])
  
  plot(df[,myRow],df[,myCol])
  newVector<-df[,myRow]+df[,myCol]
  summary(myModel) 
  myNames<-colnames(df)
  df<-cbind(df,newVector)
  print(paste(myNames[myRow],myNames[myCol]))
  colnames(df)<-c(myNames,paste(myNames[myRow],myNames[myCol],sep="_"))
  df<-df[,-max(myRow,myCol)]
  df<-df[,-min(myRow,myCol)]
  newDF<-myCollapseCorrelatingMetabolites(df,threshold)
  return(newDF)
  
}

myCollapseCorrelatingMetabolitesWithAnnot<-function(df,threshold=0.99,annot){

  require(Hmisc)
  myCor<-cor(df,method="spearman")
  require(gplots)
  # heatmap.2(myCor,trace="none")
  if(max(as.vector(myCor))< threshold){
    return(as.matrix(df))
  }
  myIndex<-which(myCor >=threshold, arr.ind = TRUE)
  myIndex<-myIndex[! (myIndex[,1]==myIndex[,2]),]
  if(nrow(myIndex) == 0){
    return(df)
  }
  myIndex<-cbind(myIndex,0)
  for(i in 1:nrow(myIndex)){
    myIndex[i,3]<-myCor[myIndex[i,"row"],myIndex[i,"col"]]
  }
  colnames(myIndex)<-c("row","col","value")
  myIndex<-myIndex[order(myIndex[,"value"],decreasing=T),]
  myRow<-myIndex[1,"row"]
  myCol<-myIndex[1,"col"]
  myNames<-colnames(df)
  rowName<-annot[annot$AlignID==sub("X","",myNames[myRow]),"myCompoundName"]
  colName<-annot[annot$AlignID==sub("X","",myNames[myCol]),"myCompoundName"]
  
  myModel<-lm(df[,myRow]~df[,myCol])
  plot(df[,myRow],df[,myCol],main=paste(rowName,"vs.",colName),xlab=paste(rowName),ylab=paste(colName))
  abline(coef(myModel)[1:2])
  cf <- round(coef(myModel), 2) 
  eq <- paste0("mpg = ", cf[1],
               ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), rowName,
               ifelse(sign(cf[3])==1, " + ", " - "), abs(cf[3]), colName)
  mtext(eq, 3, line=-2)
  
  newVector<-df[,myRow]+df[,myCol]
  summary(myModel) 
 
  df<-cbind(df,newVector)
    print(paste(rowName,colName))
  colnames(df)<-c(myNames,paste(myNames[myRow],myNames[myCol],sep="_"))
  df<-df[,-max(myRow,myCol)]
  df<-df[,-min(myRow,myCol)]
  newDF<-myCollapseCorrelatingMetabolitesWithAnnot(df,threshold,annot)
  return(newDF)
  
}

#### Will extract metabolites from a mixOmics plsda result according to weight in PLS components and some criteria
#### Returns a data.frame with the metabolites that are within the threshold.
myExtractHighestLoadingsFromPLSDAMixOmics<-function(myResult,Module=T,percentile=0.10,components=3){
  myPmetabolites<-vector()
  if(Module == F ){
    for(i in 1:components){
     higher<-quantile(myResult$loadings$X[,1],(1-percentile))
     lower<-quantile(myResult$loadings$X[,1],percentile)
     higherMetabolites<-rownames(myResult$loadings$X)[myResult$loadings$X[,i]>=higher]
     lowerMetabolites<-rownames(myResult$loadings$X[,i])[myResult$loadings$X[,i]<=lower]
     myPmetabolites<-c(myPmetabolites,higherMetabolites,lowerMetabolites)
    }
  }else if(Module == T){
    ## Or extract metabolites with higher modules in the Component space
    moduleVector<-vector()
    for (i in 1:nrow(myResult$loadings$X)){
     print(i)
       moduleVector<-c(moduleVector,sqrt((myResult$loadings$X[i,1]^2)+(myResult$loadings$X[i,2]^2)+(myResult$loadings$X[i,3]^2)))
     }
     higherModule<-quantile(moduleVector,1-(percentile/2))
     lowerModule<-quantile(moduleVector,(percentile/2))
     myPmetabolites<-vector()
     higherMetabolites<-rownames(myResult$loadings$X)[moduleVector>=higherModule]
     lowerMetabolites<-rownames(myResult$loadings$X)[moduleVector<=lowerModule]
     myPmetabolites<-c(higherMetabolites,lowerMetabolites)
     myPmetabolites<-unique(myPmetabolites)
  }
  #Returns a matrix containing the metabolite values as used for PLSDA after the same transform 
  return((myResult$X[,myPmetabolites]))
}

#### Takes a numeric matrix and a categorical vector and plots a heatmap. No data transformation by now.
myPlotHeatmap2WithLabel<-function(df,variable, variableName, paletteName, paletteSeed, title,file=NULL, legend=T,hclustMethod="ward.D2",distMethod="euclidean"){
  n<-nlevels(variable)
  variableVector<-variable
  require(wesanderson)
  if(paletteName %in% names(wes_palettes)){
    availableColors<-length(wes_palette(name=paletteName))
    if(nlevels(variableVector) > availableColors){
      print ("Not enough colors in palette for variable ")
    }
    myColors<-wes_palette(n=nlevels(variableVector),name=paletteName)
  }else{
    myColors<-brewer.pal(nlevels(variableVector),paletteName)[1:nlevels(variableVector)]
  }
  # selcol <- colorRampPalette(brewer.pal(n,paletteName))
  selcol<-colorRampPalette(myColors)
  clustcol.height<-selcol(n)
  require(gplots)
  if(is.null(file)){
   pdf(paste(variableName,"/heatmap_by_",variableName,"_PLSDAmetabolites.pdf",sep=""))
  }else{
   pdf(file)
  }
  heatmap.2((as.matrix(df)),trace="none",col=colorpanel(500,"blue","white","red"),
            hclustfun=function(x) hclust(x, method=hclustMethod),
            distfun = function(x) dist(x,method=distMethod),
            key=T,RowSideColors=c(clustcol.height[variable]),margins=c(15,12))
  legend("topright",      
         legend = unique(variable),
         col = unique(clustcol.height[variable]), 
         lty= 1,             
         lwd = 5,           
         cex=.7
  )
  dev.off()
}

myMixOmicsCIMwithLabel<-function(myResult,variableVector=NULL, variableName, paletteName, paletteSeed, 
                                 title,file=NULL, legend=T,hclustMethod="ward.D2",distMethod="euclidean",scale=T){
   
   variableVector<-factor(myResult$Y,levels=(levels(myResult$Y)),ordered=T)
   n<-nlevels(variableVector)
   require(wesanderson)
   require(mixOmics)
   if(paletteName %in% names(wes_palettes)){
     availableColors<-length(wes_palette(name=paletteName))
     if(nlevels(variableVector) > availableColors){
       print ("Not enough colors in palette for variable ")
     }
     myColors<-wes_palette(n=max(4,nlevels(variableVector)),name=paletteName)
   }else{
     myColors<-brewer.pal(n=max(4,nlevels(variableVector)),name=paletteName)
   }
   # selcol <- colorRampPalette(brewer.pal(n,paletteName))
   set.seed(paletteSeed)
   myColors<-sample(myColors,nlevels(variableVector))
   selcol<-colorRampPalette(myColors)
   clustcol.height<-selcol(n)
   require(gplots)
   if(is.null(file)){
     pdf(paste(variableName,"/cim_by_",variableName,"_PLSDAmetabolites.pdf",sep=""))
   }else{
     pdf(file)
   }
   require(viridis)
   cim(myResult,color = colorRampPalette(c("red", "white", "blue"))(200),
             clust.method=c(hclustMethod,hclustMethod),
              dist.method = c(distMethod,distMethod),
             row.sideColors=c(clustcol.height[variableVector]),margins=c(15,12))
   legend("topright",      
          legend = unique(variableVector),
          col = unique(clustcol.height[variableVector]), 
          lty= 1,             
          lwd = 5,           
          cex=.7
   )
   dev.off()
 }

#### Computes correlation between two matrices, both r2, p-values, corrected p-values, and optionally
#### Keeps a csv, a network file and/or plots a heatmap or igraph. Can return igraph object.
myTwoMatrixCorrelation<-function(){
  
}

### Plot biplot from a prcomp object.
PCbiplot <- function(PC, x="PC1", y="PC2") {
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=3, aes(label=obsnames))
  plot <- plot + geom_hline(aes(0), size=.2,yintercept=0) + geom_vline(aes(0), size=.2,xintercept=0)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  plot
}

myGgplotRegression <- function (fit) {
  ### From https://susanejohnston.wordpress.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/
  require(ggplot2)
  
  myPlot<-ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
  return(myPlot)
}

myLmp <- function (modelobject) {
  ### From https://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

myFindUnivariateCorrelation<-function(mixOmicsPCAResult,variables,threshold=0.05){
  require(grid)
  require(gridExtra)
  require(ggplot2)
  require(ggthemes)
  require(ggpubr)
  require(RColorBrewer)
  require(wesanderson)
  components<-mixOmicsPCAResult$x
  variables <- variables[,sapply(variables,is.numeric)]
  components<-components[rownames(components) %in% rownames(variables),]
  variables<-variables[rownames(variables) %in% rownames(components),]
  components<-components[order(row.names(components)),]
  variables<-variables[order(row.names(variables)),]
  myVariables<-colnames(variables)
  myComponents<-colnames(components)
  myPlotList<-list()
  for(variable in myVariables){
    for (component in myComponents){
      myDF<-cbind(components[,component],variables[,variable])
      colnames(myDF)<-c(component,variable)
      myDF<-as.data.frame(myDF)
      myModel<-lm(myDF)
      if(myLmp(myModel) <= threshold){
        myPlotList[[variable]][[component]]<-myGgplotRegression(myModel)
      }
    } 
  }
  return(myPlotList)
}

### Will delete non-connected nodes from an igraph class object
delete.isolates <- function(graph, mode = 'all') {
  isolates <- which(degree(graph, mode = mode) == 0) 
  print(isolates)
  delete.vertices(graph, isolates)
}
####
# Will clean a dataframe:
#  - Replace Empty values for NAs
#  - Remove outliers from continuous variables
#  - Refactor factor variables to filter out unused levels.
myCleanDataframe<-function(my.df){
  
}
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
# Returns a matrix containing p.values
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

#### Remove outliers from a continuous variable
#### And replace it by NA
outlierKD<-function(dt, var,remove="yes",coef=1.5) {
  var_name <- eval(substitute(var),eval(dt))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name,coef=coef)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "n")
  cat("Propotion (%) of outliers:", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "n")
  cat("Mean of the outliers:", round(mo, 2), "n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "n")
  cat("Mean if we remove outliers:", round(m2, 2), "n")
  #response <- readline(prompt="Do you want to remove outliers and to replace with NA? [yes/no]: ")
  if(remove == "y" | remove == "yes"){
    dt[as.character(substitute(var))] <- invisible(var_name)
    assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    cat("Outliers successfully removed", "n")
    return(invisible(dt))
  } else{
    cat("Nothing changed", "n")
    return(invisible(var_name))
  }
}
# #### Meant for igraph objects
# delete.isolates <- function(graph, mode = 'all') {
#   isolates <- which(degree(graph) == 0) 
#   delete.vertices(graph, isolates)
# }

### Rusty function to classify Age
classifyAge<-function(age){
  if(!is.numeric(age)){
    age<-as.numeric(age)
  }
  if(age<=5){
    return("0-5")
  }else if(age<=10){
    return("5-10")
  }else if(age<=15){
    return("10-15")
  }else if(age <= 20){
    return("15-20")
  }else if(age <= 25){
    return("20-25")
  }else if(age <= 30){
    return("25-30")
  }else if(age <=35){
    return("30-35")
  }else if(age <= 40){
    return("35-40")
  }else if(age <= 45){
    return("40-45")
  }else if(age <= 50){
    return("45-50")
  }else if(age <= 55){
    return("50-55")
  }else if(age<=60){
    return("55-60")
  }else if(age<=65){
    return("60-65")
  }else if(age<=70){
    return("65-70")
  }else{
    return (">70")
  }
}

#' @title Exporting phyloseq Data in CSV Files
#' @description Writes the otu, taxonomy and metadata in csv files.
#' @param x \code{\link{phyloseq-class}} object
#' @param type 'OTU' or 'TAXONOMY' or 'METADATA'
#' @param path Path to the directory/folder where the data will be written.
#' Uses the working directory by default.
#' @return  Output file path (a string)
#' @seealso read_phyloseq
#' @export
#' @examples \dontrun{
#' data(dietswap)
#' pseq <- dietswap
#' # By default writes all info at once (ie OTU/TAXONOMY/METADATA)
#' write_phyloseq(pseq) 
#' write_phyloseq(pseq, 'OTU')
#' write_phyloseq(pseq, 'TAXONOMY')
#' write_phyloseq(pseq, 'METADATA')
#' }
#' @keywords utilities
write_phyloseq <- function(x, type="all", path=getwd()) {
  
  # TODO make read_phyloseq as well
  if (type == "OTU" || type == "all") {
    f <- paste(path, "otu_table.csv", sep="/")
    message("Writing OTU in the file ", f)
    # y <- as.data.frame(x@otu_table);
    if (f %in% dir(path)) {
      warning(paste("The file with the same name", f,
                    "exists in the given path and is overwritten."))
    }
    # Let us use abundances function here as it is guaranteed to be taxa x
    # samples always
    y <- abundances(x)
    write.csv(y, file=f, fileEncoding="UTF-16LE")
  } else if (type == "TAXONOMY" || type == "all") {
    # Renamed from TAXA to TAXONOMY as the latter is used elsewhere
    f <- paste(path, "taxonomy_table.csv", sep="/")
    message("Writing TAXONOMY in the file ", f)
    if (f %in% dir(path)) {
      warning(paste("The file with the same name", f,
                    "exists in the given path and is overwritten."))
    }
    y <- as.data.frame(tax_table(x))
    write.csv(y, file=f, fileEncoding="UTF-16LE")
  } else if (type == "METADATA" || type == "all") {
    f <- paste(path, "metadata_table.csv", sep="/")
    message("Writing METADATA in the file ", f)
    if (f %in% dir(path)) {
      warning(paste("The file with the same name", f,
                    "exists in the given path and is overwritten."))
    }
    y <- meta(x)
    write.csv(y, file=f, fileEncoding="UTF-16LE")
  }
  
  return(path)
  
}

#' @title Abundance Matrix from Phyloseq
#' @description Retrieves the taxon abundance table from
#' phyloseq-class object and ensures it is systematically returned as
#' taxa x samples matrix.
#' @inheritParams transform
#' @return Abundance matrix (OTU x samples).
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @export
#' @aliases ab, otu
#' @examples
#' data(dietswap)
#' a <- abundances(dietswap)
#' # b <- abundances(dietswap, transform='compositional')
#' @keywords utilities
abundances <- function(x, transform="identity") {
  
  # Pick the OTU data
  if (class(x) == "phyloseq") {
    
    # Pick OTU matrix
    otu <- as(otu_table(x), "matrix") # get_taxa(x)
    
    # Ensure that taxa are on the rows
    if (!taxa_are_rows(x) && ntaxa(x) > 1 && nsamples(x) > 1) {
      otu <- t(otu)
    }
    
    if (ntaxa(x) == 1) {
      otu <- matrix(otu, nrow=1)
      rownames(otu) <- taxa(x)
      colnames(otu) <- sample_names(x)
    }
    
    if (nsamples(x) == 1) {
      otu <- matrix(otu, ncol=1)
      rownames(otu) <- taxa(x)
      colnames(otu) <- sample_names(x)
    }
    
  } else if (is.vector(x)) {
    
    otu <- as.matrix(x, ncol=1)
    
  } else {
    
    # If x is not vector or phyloseq object then let us assume it is a
    # taxa x samples
    # count matrix
    otu <- x
    
  }
  
  # Apply the indicated transformation
  if (!transform == "identity") {
    otu <- transform(otu, transform)
  }
  otu
  
}

myPhyloseqToCsv<-function(psObject,filename){
  require(phyloseq)
  myTaxonDF<-data.frame(t(otu_table(psObject)))
  myTaxaDF<-data.frame((tax_table(psObject)))
  myTaxonTaxaDF<-cbind(myTaxonDF,myTaxaDF)
  write.csv(myTaxonTaxaDF,file=filename)
}
myPhyloseqToBiom<-function(psObject,filename){
  require(biomformat)
  require(phyloseq)
  data<-data.frame(t(otu_table(psObject)))
  smd<-data.frame(sample_data(psObject))
  omd<-data.frame(tax_table(psObject))
  myBiom<-make_biom(data,smd,omd)
  write_biom(myBiom,biom_file=filename)
}
myLefseFunction<-function(psObject,myRank){
  TaxonomicalLevels<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  innerPS<-tax_glom(psObject,taxrank=myRank)
  otus<-data.frame(otu_table(innerPS))
  taxons<-data.frame(tax_table(innerPS)[,1:match(myRank,colnames(tax_table(innerPS)))],check.names = F,check.rows=F)
  taxons[,(match(myRank,TaxonomicalLevels)+1):length(TaxonomicalLevels)]<-NA
  colnames(taxons)<-TaxonomicalLevels
  for (j in 1:dim(taxons)[1]){
    taxid<-taxons[j,1]
    for (i in 2:dim(taxons)[2]){
      if(!(is.na(taxons[j,i]))){
        taxid<-paste(taxid,taxons[j,i],sep="|")
      }
    }
    taxons[j,"Taxon"]<-taxid
  }
  taxons$OtuID<-rownames(taxons)
  taxons<-(taxons[,c("OtuID","Taxon")])
  otus$OtuID<-rownames(otus)
  taxons<-merge(taxons,otus,by="OtuID")
  taxons$OtuID<-NULL
  return(taxons)
}
myPhyloseqToLefse<-function(psObject,lowestTaxonomicalLevel="Genus"){
  TaxonomicalLevels<-c("Species","Genus","Family","Order","Class","Phylum","Kingdom")
  xx<- transform_sample_counts(psObject, function(x) ((( x)/sum(x))))
  xx.genus<-myLefseFunction(xx,myRank="Genus")
  xx.family<-myLefseFunction(xx,myRank="Family")
  xx.order<-myLefseFunction(xx,myRank="Order")
  xx.class<-myLefseFunction(xx,myRank="Class")
  xx.phylum<-myLefseFunction(xx,myRank="Phylum")
  xx.kingdom<-myLefseFunction(xx,myRank="Kingdom")
  xx<-rbind(xx.genus,xx.family,xx.order,xx.class,xx.phylum,xx.kingdom)
  rownames(xx)<-xx$Taxon
  xx$Taxon<-NULL
  xx<-t(xx)
  rownames(xx)<-gsub("X","",rownames(xx))
  rownames(xx)<-gsub("\\.","-",rownames(xx))
  xx<-merge(data.frame(sample_data(psObject)),xx,by="row.names")
  xx$Row.names<-NULL
  xx<-t(xx)
  colnames(xx)<-xx[1,]
  xx<-data.frame(xx[-1,])
  xx$SampleID<-rownames(xx)
  xx<-xx[,c(ncol(xx),1:(ncol(xx)-1))]
  return(xx)
}


### 
# Returns colors RGB codes using wes anderson palettes
# uses n( number of colors) to either sum up colors from different palettes
# or just only ones. Palettes to be used are in paletteNames
# returns a vector with RGB codes.
myCreateWesAndersonPalette<-function(n=1,paletteNames="Darjeeling1"){
  require(wesanderson)
  if(! (paletteNames %in% names(wes_palettes) )){
    return(get_palette(paletteNames,n))      
  }
  
  if(n<=8){possiblePalettes=names(wes_palettes[lapply(wes_palettes,length)>=n])} ### Wes Anderson palettes maximum length is limited to 7
  if(length(unlist(wes_palettes[paletteNames])) < n){ 
    print("Not enought color in palettes")
  }
  return(as.vector(unlist(wes_palettes[paletteNames]))[1:n])
  
}

myVeganOTUFromPhyloseq <- function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

myBoxPlotSingleVariable<-function(myDF,Xlabel=NULL,Ylabel=NULL,title=NULL,filename=NULL,filedir=NULL,scale=F,fontSize=14,center=F,
                                  stats=TRUE,method=NULL,paletteNames=c("Darjeeling1","Royal1","Moonrise1"),WhatToReturn="stats",savePlot=T){
  #Assumes Categorical variable is in first column and continuous variable is in second column
  if(is.null(Xlabel)){
    Xlabel=colnames(myDF)[1]
  }
  if( ! is.factor(myDF[,1])){
    myDF[,1]<-factor(myDF[,1],levels=sort(unique(myDF[,1])),ordered=T)
  }
  if(is.null(Ylabel)){
    Ylabel=colnames(myDF)[2]
  }
  if(is.null(title)){
    title=paste(Ylabel," vs ",Xlabel,sep="")
  }
  myDF[,1]<-factor(myDF[,1],levels=levels(myDF[,1]),ordered=T)
  if(isTRUE(stats) & is.null(method)){
    if(nlevels(myDF[,1]) > 2){
      method<-"kruskal.test"
    } else if(nlevels(myDF[,1] == 2)){
      method<-"wilcox.test"
    }
  }
  require(ggplot2)
  require(ggpubr)
  require(icesTAF)
  colnames(myDF)<-c("Group","value")
  # print(colnames(myDF)) 
  n<-nlevels(myDF[,"Group"])
  # print(class(n))
  myPlot<-ggplot(myDF,aes(x=Group,y=value,fill=Group))+geom_violin()+geom_boxplot(fill="white",width=0.3)+
    scale_fill_manual(values=myCreateWesAndersonPalette(as.numeric(n),paletteNames))+
    ggtitle(label = title)+xlab(Xlabel)+theme(text = element_text(size=fontSize))+
    labs(fill=Xlabel)+ylab(Ylabel)+stat_compare_means(label.x=1.3,size=fontSize-10)
  myTest<-compare_means(value~Group,data=myDF,method=method)
  if( ! is.null(filedir)){
    ##Directory for output is defined
    if (is.null(filename)){
      ### CreateFile Name
      filename=paste(Xlabel,"_vs_",Ylabel,"_",format(Sys.Date(),"%d%m%y"),"_",format(Sys.time(),"%H%M%OS"))      
    }
    print(filename)
    if(myTest$p<0.05){
      outputFile=paste(filedir,"SignificantPairwise/",myTest$p.signif,"_",filename,sep="")
      mkdir(paste(filedir,"SignificantPairwise/",sep=""))
    }else{
      outputFile=paste(filedir,"NonSignificantPairwise/",myTest$p.signif,"_",filename,sep="")
      mkdir(paste(filedir,"NonSignificantPairwise/",sep=""))
    }
    if(! grepl(".pdf",outputFile)){
      outputFile=paste(outputFile,".pdf",sep="")
    }
  }
  if ( savePlot==F){outputFile<-NULL}
  if ( is.null(outputFile)){
    print("Output File not defined, not printing output")
  }else{
    ggsave(myPlot,filename = outputFile)
  }
  
  if(WhatToReturn=="stats"){
    return(myTest)
  }else if(WhatToReturn=="plots"){
    return(myPlot)
  }else if(WhatToReturn=="both"){
    return(list(myTest,myPlot))
  }
}
####
#### Compare Two Groups using DESeq2
myCompareTwoGroups<-function(mydata=NULL,variable=NULL,category1=NULL,category2=NULL,fileForPlot=NULL,fileForTable=NULL,baseDir="",minCounts=500,maxAlpha=0.01, design=NULL, fontSizeY=16, minBaseMean=10){
  fileForPlot=paste(baseDir,"NegBin_DiffTest_",as.character(deparse(substitute(mydata))),"_",variable,"_",category1,"vs",category2,".pdf",sep="")
  print(paste("Output file for Plots: ",fileForPlot))
  fileForTable=paste(baseDir,"NegBin_DiffTest_",as.character(deparse(substitute(mydata))),"_",eval(variable),"_",eval(category1),"vs",eval(category2),".txt",sep="")
  print(paste("Output file for Table: ",fileForTable))
  stringForTitle=paste(variable,"/",category1,"vs",category2)
  colnames(tax_table(mydata))<- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  myTaxTable<-tax_table(mydata)
  kostic <- mydata
  LastRank<-rank_names(mydata)[length(rank_names(mydata))][[1]]
  myLastRankName<-sort(unique(myTaxTable[,"Phylum"]))
  print(myLastRankName)
  require(phyloseq)
  require(DESeq2)
  kostic <- prune_samples(sample_sums(kostic) > minCounts, kostic)
  diagdds = phyloseq_to_deseq2(kostic, as.formula(design))
  #print(diagdds)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(diagdds), 1, gm_mean)
  print("geoMeans")
  diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
  #colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("untreated","treated"))
  #colData(diagdds)$condition<-factor(colData(diagdds)$condition,levels=c(category1,category2))
  print("performing DEseq")
  diagdds = DESeq(diagdds, fitType="parametric",test="Wald",parallel = F)
  print(paste("DESeq done, filtering results for the following thresholds: MinCounts_", minCounts, "minBaseMean: ",minBaseMean, "maxAlpha",maxAlpha))
  res=results(diagdds,contrast = c(variable,category1,category2))
  #res=results(diagdds)
  print(res)
  res = res[order(res$padj, na.last=NA), ]
  sigtab = res[(res$padj < 1), ]
  #   sigtab = res[(res$padj < maxAlpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kostic)[rownames(sigtab), ], "matrix"))
  sigtab$OtuID<-rownames(sigtab)
  #   head(sigtab)
  write.table(sigtab,file=fileForTable,sep="\t")
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=levels(sort(sigtab$Phylum)))
  #   sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=levels(myLastRankName))
  print(dim(sigtab))
  #Cleanup for Positive enrichment in csigtabarcinoma
  posigtab=sigtab
  #posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
  #posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
  posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", rank_names(mydata))]
  library("ggplot2")
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set3", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  sigtabgen=sigtab
  sigtabgen = subset(sigtab, !is.na(Genus))
  sigtabgen = subset(sigtabgen, padj <= maxAlpha)
  sigtabgen = subset(sigtabgen, baseMean >= minBaseMean)
  sigtabgen = subset(sigtabgen, sigtabgen$Genus != "unclassified")
  # Phylum order
  #   x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
  #   x = sort(x, TRUE)
  #   sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
  #   sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=sort(as.character(sigtabgen$Phylum)))
  
  # Genus order
  if(as.character(LastRank)== "Genus"){
    sigtabgen$LastRank<-sigtabgen[,"Genus"]
  }else{
    sigtabgen$LastRank<-paste(sigtabgen[,"Genus"]," ",sigtabgen[,as.character(LastRank)])
  }
  x = tapply(sigtabgen$log2FoldChange, sigtabgen$LastRank, function(x) min(x))
  x = sort(x, TRUE)
  
  sigtabgen$LastRank = factor(as.character(sigtabgen$LastRank), levels=names(x))
  sigtabgen$log2Counts<-sigtabgen$baseMean
  sigtabgen$alpha<- 1 - sigtabgen$padj
  #pdf(fileForPlot)
  p<-ggplot(sigtabgen,aes(x=LastRank,y=log2FoldChange))
  p<-p+geom_point(aes(colour=Phylum,size=baseMean,alpha=alpha))
  p<-p+ guides(colour = guide_legend(override.aes = list(size=10)))
  p<-p+scale_size_continuous(range=c(10,20))
  #+geom_point(aes(size=sigtabgen$log2Counts))+scale_size_area()
  p<-p+theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5,size=10))
  #p<-p+theme(legend.key.size=unit(1,"cm"))
  p<-p+ ggtitle(paste(stringForTitle," Data:",as.character(deparse(substitute(mydata))))) + 
    theme(plot.title = element_text(lineheight=.7, face="bold"))+coord_flip()+
    theme(axis.text.y = element_text( size=fontSizeY)) + geom_vline(xintercept=0,colour="darkred", linetype = "longdash")
  print(p)
  #ggplot(sigtabgen, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6)+scale_size(range=c(1,5))+
  # theme(axis.text.x = element_text(angle = -90, hjust = 0,size=3, vjust=0.5), legend.key.size=unit(0.5,"cm"),legend.text=element_text(size=3))
  ggsave(filename = fileForPlot,dpi=600,width=11, height=8.5)
  return(sigtab)
  #dev.off()
}
### Returns Boxplot with anova or glm test, density and jitter
#### Numeric Boxplots
myBoxplotNumericByGroup<-function(mydata,category,variable,nbvariable,test,Rank=NULL){
  # if(is.null(Rank)){
  #    title<-paste(as.character(variable), " by ",category)
  #    fileForOutput<-paste(as.character(variable),"by",category,sep="_")
  #  }else{
  Phylum<-unique(as.vector(mydata[,Rank]))
  title<-paste(as.character(variable), " of\n ",Phylum,as.character(Rank)," by ",category)
  fileForOutput<-paste(as.character(variable),"of",Phylum,as.character(Rank),"by",category,sep="_")
  #  }
  numberOfLevels<-length(unique(mydata[,category]))
  colorsPlotRight<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
  require(gridExtra)
  if(numberOfLevels==2){
    test<-t.test(mydata[,variable]~mydata[,category])
    testString<-paste("Students t-test. p-value",statisticValue(test$p.value))
  }
  if(numberOfLevels>=3){
    require(MASS)
    glm.nb.model<-glm.nb(mydata[,nbvariable]~mydata[,category],method="glm.fit")
    glm.nb.aov<-aov(glm.nb.model)
    test<-aov(mydata[,variable]~mydata[,category])
    tukey.test<-TukeyHSD(test)
    print(tukey.test$mydata)
    mymatrix<-tukey.test$mydata
    mymatrix<-mymatrix[,c("diff","p adj")]
    mymatrix<-as.matrix(mymatrix)
    #mymatrix<-round(mymatrix,3)
    for(i in 1:nrow(mymatrix)){
      print(mymatrix[i,"p adj"])
      mymatrix[i,"p adj"]<-statisticValue(as.numeric(mymatrix[i,"p adj"]))
    }
    text2.df<-as.table(mymatrix)
    testString<-paste("ANOVA PR(>F)", statisticValue(summary(test)[[1]][["Pr(>F)"]][[1]]))
    testString<-paste(testString,"\n","NegBin ANOVA PR(>F)",statisticValue(summary(glm.nb.aov)[[1]][["Pr(>F)"]][[1]]))
  }
  #  if(test=="lm"){
  #    model<-lm(mydata[,variable]~mydata[,category])
  #  }
  mydata$xaxis<-mydata[,category]
  mydata$yvalue<-mydata[,variable]
  p<-ggplot(mydata,aes(x=xaxis,y=yvalue,fill=as.factor(xaxis)))+
    geom_boxplot()+geom_jitter(color="DarkRed")+
    #ggtitle(title)+
    xlab(category)+ylab(variable)+
    ylim(min(mydata$yvalue-1),1)+
    scale_fill_manual(values=colorsPlotRight[1:numberOfLevels])+
    theme(legend.position=c(1,1),legend.justification=c(1,1))+
    #annotate("text",x=numberOfLevels/2.5,y=0.5,label=testString,size=3)+
    annotate("text",x=1,y=0,label=testString,size=3)+
    annotate("text",y=0.4,x=1.5,label=title,size=4)
  plotRight<-ggplot(mydata,aes(yvalue,fill=xaxis))+geom_density(alpha=.5)+
    coord_flip()+scale_fill_manual(values=colorsPlotRight[1:numberOfLevels])+
    theme(legend.position="none")+
    xlim(min(mydata$yvalue-1),max(mydata$yvalue+1))
  if(numberOfLevels>=3){
    p<-p+annotation_custom(tableGrob(text2.df,gp=gpar(fontsize=8)), ymin=min(mydata$yvalue)-1, ymax=min(mydata$yvalue), xmax=numberOfLevels/1.2, xmin=numberOfLevels/2)
  }  #p2<-tableGrob(text2.df)
  #grid.arrange(p2,p,main="prova",ncol=2)
  fileForPlot <- paste(fileForOutput,".pdf")
  pdf(fileForPlot,paper="a4r")
  grid.arrange(p,plotRight,nrow=1,ncol=2,widths=c(4,1),heights=c(1,4))
  #p2<-ggplot(p2)
  #ggsave(filename = fileForPlot,dpi=600,width=11, height=8.5)
  dev.off()
  grid.arrange(p,plotRight,nrow=1,ncol=2,widths=c(4,1),heights=c(1,4))
  #print(p2)
  #return(p2)
}

##### TO - DO ##########
# myBoxplotManyVariableUnivariate
# Expects a data.frame class object, Where exploratory variables are the first to last-1 columns, and all of them are numeric
# Expect Response variable in the last column to be categoric
myBoxplotManyVariablesUnivariate<-function(myDF,VariableName=NULL,facet=F,wrap=F,stats=T,resultsDir=NULL, paletteNames="Darjeeling1"){
  
}



############### Linear Mixed Models Functions #################
############### Created by Javier Rivera-Pinto, 2018 ##########
###############################################################
#-----------------------------------------------------------------------------#
#                   LINEAR MIXED MODEL: MAKE AND PLOT IT
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# NAME: t.piecewise
# FUNCTION: creates auxiliary variables T1, . . . Ti, as many as desired
#           breakpoints for LMM
#-----------------------------------------------------------------------------#
myLMM_t.piecewise<-function(Data,Variable,t){
  
  # Data:     data where Variable is contained
  # Variable: TimePoint variable
  # t:        vector with break-points for the TimePoint variable
  
  # Variable of interest
  V<-Data[,Variable]
  # Length of breakpoints
  n<-length(t)
  # Are t values presented in Variable?
  stopifnot(sum(t%in%Data[,Variable])==n)
  # Build a matrix with the new variables
  M<-matrix(0,ncol=n,nrow=length(V))
  # Names of variables
  Names<-vector()
  # Define new variables (T1, . . . , Tn)
  for (i in 1:n){
    M[,i]<-(V-t[i])*(V>=t[i])
    Names<-c(Names,paste("T",i,sep=""))
  }
  # Colnames for T1, . . . , Tn
  colnames(M)<-Names
  # Previous data concatenated with new variables
  DataM<-cbind(Data,M)
  # Return information
  return(DataM)
}

myLMM_t.piecewise2<-function(Data,Variable,t){
  
  # Data:     data where Variable is contained
  # Variable: TimePoint variable
  # t:        vector with only one breakpoint
  
  # Variable of interest
  V<-Data[,Variable]
  # Length of breakpoints
  n<-length(t)
  # Are t values presented in Variable?
  stopifnot(sum(t%in%Data[,Variable])==n)
  # Build a matrix with the new variables
  M<-matrix(0,ncol=n,nrow=length(V))
  # Names of variables
  Names<-vector()
  # Define new variables (T1, . . . , Tn)
  for (i in 1:n){
    M[,i]<-(V-t[i])*(V>=t[i])
    Names<-c(Names,paste("T",i,sep=""))
  }
  # Colnames for T1, . . . , Tn
  colnames(M)<-Names
  # Aditional variable
  newt <- ifelse(V<t,V,t)
  # Previous data concatenated with new variables
  DataM<-cbind(Data,M,newt)
  # Return information
  return(DataM)
}
#-----------------------------------------------------------------------------#
# NAME: two.piece.GLMM
# FUNCTION: plots the LMM calculated for a certain variable and diferent
#           levels for a factor (for example Gender = (M,F))
#-----------------------------------------------------------------------------#

myLMM_two.piece.GLMM<-function(Datos,Y,Time,t=NULL,Factor=NULL,ID){
  
  # Datos    : data.frame containing the variables
  # Y        : response variable
  # Time     : variable type numeric (referring to Time)
  # t        : vector with break-points for the TimePoint variable
  # Factor   : variable type factor, a plot for each level
  # ID       : variable for patient identification (PatientID)
  
  if (is.null(Factor)){Datos[,"Factor"]<-as.factor(1);Factor<-"Factor"}
  
  # Proof that Factor is a factor class object and t is a numeric one
  stopifnot(class(Datos[,Factor])=="factor")
  
  # Number of breakpoints
  n<-length(t)     
  # Sort t values
  if(n) t<-sort(t)
  # Load library
  library(nlme)
  
  # Colour palette
  Col<-c("deepskyblue4","darkolivegreen3","darkorchid2","gold3")
  
  # ncol Datos
  colDatos<-ncol(Datos)
  # New data with T1, . . . , Tn variables (breakpoints)
  if(n==0){D<-Datos
  }else {D<-myLMM_t.piecewise(Data=Datos,Variable=Time,t=t)
  D2<-myLMM_t.piecewise2(Data=Datos,Variable=Time,t=t)}
  # Delete problems with NA
  D<-D[!is.na(D[,Y]),]
  if (n>0){D2<-D2[!is.na(D2[,Y]),]}
  # Factor levels
  lev<-levels(D[,Factor])
  # Number of Factor levels
  f<-length(lev)
  # Minimum value of t
  m<-min(D[,Time]); M<-max(D[,Time])
  T.Plot<-c(m,t,M)
  nt<-length(T.Plot)
  # Vector of increments (for plots)
  stopifnot(m!=t[1]); stopifnot(M!=t[n])
  
  # Build formula ()
  # Form: formula for fixed effects
  if (n==0){ Form<-formula(paste(Y," ~ ",Time))
  }else{
    Form<-formula(paste(Y," ~ ",Time,"+",
                        paste(colnames(D)[(colDatos + 1):(colDatos + n)],
                              collapse="+")))
  }
  # Form2: formula for random effects
  Form2<-formula(paste("~ 1|",ID))
  # Form3: formula for the second slope (new parametrization)
  Form3<-formula(paste(Y,"~ newt + T1",sep=" "))
  
  # List with the LMM results
  A<-list()
  B<-list()
  
  # For each factor level, . . .
  for (i in 1:f){
    # Model
    model<-lme(Form, random = Form2, D[D[,Factor]==lev[i],], method = "REML")
    S <- summary(model)
    if(n!=0){
      model2<-lme(Form3, random = Form2, D2[D2[,Factor]==lev[i],], method = "REML")
      S2 <- summary(model2)}
    
    # Table with coefficients and p-values
    A[[i]]<-S["tTable"][[1]]
    if (n!=0){
      B[[i]]<-S2["tTable"][[1]]
    }
  }
  
  # Elements to build a data.frame for the plot
  # Global elements (independent from n)
  # Factor
  Gen<-rep(lev,each=n+1)
  # x (initial points)
  x.ini<-rep(T.Plot[-nt],f)
  # xend (ending x positions)
  x.end<-rep(T.Plot[-1],f)
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  if (n==0){
    y.ini<-vector();y.end<-vector()
    for(j in 1:f){y.ini<-c(y.ini,A[[j]][1,1] + m*A[[j]][2,1])
    y.end<-c(y.end,A[[j]][1,1] + M*A[[j]][2,1])}
  }else if (n==1){
    y.ini<-vector();y.end<-vector()
    for(j in 1:f){
      y.ini<-c(y.ini,A[[j]][1,1] + m*A[[j]][2,1],
               A[[j]][1,1] + t*A[[j]][2,1])
      y.end<-c(y.end,A[[j]][1,1] + t*A[[j]][2,1],
               A[[j]][1,1] + M*A[[j]][2,1] + (M-t)*A[[j]][3,1])
    }
  } else{ # If n>=2
    
    # y (initial points) Deleting t_min
    Mat<-rbind(rep(1,length(T.Plot)-1),matrix(rep(T.Plot[-1],n+1),
                                              nrow=n+1,byrow=T))
    # Substracting t_i values to each row
    Mat[-c(1:2),]<-apply(Mat[-c(1:2),],2,function(x) x-t)
    # Negative values to 0
    Mat[Mat<0]<-0
    # Multiply Value colum of A[[i]] and Mat
    Val<-vector()
    for (j in 1:f){Val<-c(Val,A[[j]][,1]%*%Mat)}
    Val.ini<-vector()
    for (j in 1:f){Val.ini<-c(Val.ini,A[[j]][1,1] + m*A[[j]][2,1])}
    # New vector
    y.ini<-vector()
    for (j in 1:f){y.ini<-c(y.ini,Val.ini[j],
                            Val[((nt-1)*(j-1) + 1): ((nt-1)*(j-1) + n)])
    }
    
    # yend (ending x positions)
    y.end<-Val
  }
  
  # Geom.Data: a data.frame used for plot by Factor
  geom.data <- data.frame(x = x.ini,xend = x.end, y=y.ini,yend = y.end,
                          Factor=Gen)
  # To have different lines for Factor
  colnames(geom.data)[5]<-Factor
  
  # Save plot
  P<-ggplot(data=D,aes_string(x=Time,y=Y, color=Factor, group=ID)) +
    geom_line(size=1, alpha=0.4)+theme_bw()+
    facet_wrap(formula(paste("~", Factor)), scales="free_x", ncol=f) +
    scale_colour_manual(values=brewer.pal(f,"Set1"))+
    coord_cartesian(ylim=c(min(D[,Y])*0.95,max(D[,Y])*1.05))+
    scale_x_continuous(breaks=unique(D[,Time])) +
    theme(legend.position="none",strip.background=element_blank())+
    geom_point(alpha=0.4)+
    geom_vline(xintercept=6, linetype="dotted", color="black", size=2, inherit.aes=F)+
    geom_segment(data=geom.data,aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,size=1, color="black") +
    theme(axis.text=element_text(size=12), axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    theme(panel.border=element_rect(linetype="solid", color="black", size=0.5))+
    ggtitle(Y) 

    
  
  # List with results
  # Names for A object (Factor levels value)
  names(A)<-lev
  # If there is no timepoint, A = B
  if(n==0){B<-A}
  # Retunred object
  Ret<-list(A,P,B)
  # Return it    
  return(Ret)
}
################################################################################
# Example To Run previous LMM Functions /// Not to be executed
################################################################################
# MyData <- read.csv("C:/Users/user/Downloads/javi.csv")
# row.names(MyData) <- MyData$X  
# MyData$X <- NULL
# LL <- two.piece.GLMM(Datos=MyData, Y = "Observed",
#                      Time = "TimePoint",
#                      t = 6,
#                      Factor = NULL,
#                      ID = "Patient")
# 
# LL <- two.piece.GLMM(Datos=MyData, Y = "Observed",
#                      Time = "TimePoint",
#                      t = NULL,
#                      Factor = NULL,
#                      ID = "Patient")



################################################################################
#' Convert phyloseq data to MetagenomeSeq MRexperiment object
#' https://rdrr.io/bioc/phyloseq/src/R/extend_metagenomeSeq.R
#' No testing is performed by this function. The phyloseq data is converted
#' to the relevant \code{\link[metagenomeSeq]{MRexperiment-class}} object, which can then be
#' tested in the zero-inflated mixture model framework
#' (e.g. \code{\link[metagenomeSeq]{fitZig}}) 
#' in the metagenomeSeq package.
#' See the
#' \href{http://joey711.github.io/phyloseq-extensions}{phyloseq-extensions}
#' tutorials for more details.
#'
#' @param physeq (Required). \code{\link{phyloseq-class}}.
#' @param ... (Optional). Additional named arguments passed 
#'  to \code{\link[metagenomeSeq]{newMRexperiment}}.
#'  Most users will not need to pass any additional arguments here.
#'  
#' @return A \code{\link[metagenomeSeq]{MRexperiment-class}} object.
#' 
#' @seealso
#'
#'  \code{\link[metagenomeSeq]{fitTimeSeries}}
#'  \code{\link[metagenomeSeq]{fitLogNormal}}
#'  \code{\link[metagenomeSeq]{fitZig}}
#'  \code{\link[metagenomeSeq]{MRtable}}
#'  \code{\link[metagenomeSeq]{MRfulltable}}
#'
#' @export
#' @importFrom Biobase AnnotatedDataFrame
#'  
#' @examples
#'  # Check out the vignette metagenomeSeq for more details.
#'  # vignette("metagenomeSeq")
#'  data(soilrep)
#'  phyloseq_to_metagenomeSeq(soilrep)
phyloseq_to_metagenomeSeq = function(physeq, ...){
  # Enforce orientation. Samples are columns
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq)}
  # Coerce count data to vanilla matrix of integers
  countData = round(as(otu_table(physeq), "matrix"), digits=0)
  # Create sample annotation if possible
  if(!is.null(sample_data(physeq,FALSE))){
    ADF = AnnotatedDataFrame(data.frame(sample_data(physeq)))  
  } else { 
    ADF = NULL 
  }
  # Create taxa annotation if possible
  if(!is.null(tax_table(physeq,FALSE))){
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq),
                                        data.frame(tax_table(physeq)),row.names = taxa_names(physeq)))
  } else {
    TDF = AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq),
                                        row.names = taxa_names(physeq)))
  }
  # Create MRexperiment
  if(requireNamespace("metagenomeSeq")){
    mrobj = metagenomeSeq::newMRexperiment(counts = countData, phenoData = ADF, featureData = TDF,...)
    # Calculate normalization factor
    if (sum(colSums(countData > 0) > 1) < ncol(countData)) {
      p = suppressMessages(metagenomeSeq::cumNormStat(mrobj))
    }
    else {
      p = suppressMessages(metagenomeSeq::cumNormStatFast(mrobj))
    }
    mrobj = metagenomeSeq::cumNorm(mrobj, p = p)
    return(mrobj)
  }
}



#-------------------------------------------------------------------------------------------------#
#---------------------------------------#
# AUXILIAR FUNCTION 5: myLassoPredictor
#---------------------------------------#

# INPUT PARAMETERS:

#           x: matrix with explanatory variables (columns) for each individual (rows).
#           y: vector with the corresponding values referring to the assigned group.
#      n.fold: number of fold for each Cross - Validation group.
#      n.iter: number fo Cross - Validation groups.
#       C.vec: vector with REGULARIZATION STRENGTH values to analyze.
#   eval.crit: accuracy("acc"), auprc ("auprc"), auroc ("auroc"),
#              sensitivity/specifity ("sens.spe"), lift chart ("lift.chart")
# min.nonzero: minimum coefficients != 0 for the model selection to extract the Opt.C
#       Max.Q: accuracy (or quality values) for Opt.C
#        seed: seed value to replicate the results

# OUTPUT OBJECTS:

# A list with :
#                 x: the imput variable.
#                 y: the imput variable.
#         CV.groups: samples identification for CV-groups. 
#            n.fold: number of folds.
#            n.iter: number of iterations (n.fold * n.iter = nrow(CV.groups)).
#             Opt.C: vector of length n.iter with optimal REGULARIZATION STRENGTH value 
#                    for each FOLD or CV.group.
#             Max.Q: maximum values obtained for quality meassure for each fold.
#       min.nonzero: the input parameter.
#         eval.crit: value required for quality classification

#------------------------------------------------------------------------------------------------#


myLassoPredictor<-function(x,y,n.fold,n.iter,C.vec=10^seq(-2,2,length=17),
                          min.nonzero = 1,eval.crit="auroc",seed){
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  
  #------------------------------------------------#
  # Step 0: load required packages and modify data
  #------------------------------------------------#    
  # Starting time
  start.time = proc.time()[1]
  
  # Required packages:
  suppressMessages(require("LiblineaR")) # For LASSO logistic regression
  suppressMessages(library("LiblineaR"))
  suppressMessages(require("CMA"))       # For CV groups
  suppressMessages(library("CMA"))
  
  
  #-------------------------------------#
  # Modifications for response vector y
  #-------------------------------------#
  # Response vector as factor
  if (is.logical(y)){ newy<-as.factor(as.numeric(y))
  } else {
    # Define newy
    newy<-as.character(y)
    newy2<-as.factor(y)
    lev<-levels(newy2)
    
    # Redefine y with "0" or "1" levels
    newy[newy==lev[1]]<-0; newy[newy==lev[2]]<-1
    # Give factor format to y
    newy<-as.factor(newy)
  }
  # Define y as newy
  y<-newy
  # Levels of y
  lev<-levels(y)
  # Stop if there are no 2 levels
  stopifnot(length(lev)==2)
  
  
  #-------------------------------------#
  # Some parameters that will be needed
  #-------------------------------------#
  
  liblinear.type = 6      # LASSO
  class.weights = c(5, 1) # Class.weights
  eps = 1e-8              # Value for convergence
  
  
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#   
  
  #---------------------------------------------------#
  # Step 1: build different groups and folds for CV
  #---------------------------------------------------#
  
  # Fix a seed
  set.seed(seed)
  # Matrix where rows represent the individuals taking part for the learningsets
  CV.groups<-GenerateLearningsets(y=y,fold=n.fold,niter=n.iter,method="CV",strat=T)@learnmatrix
  
  # Vector for saving C.opt value for each n.iter
  Opt.C<-vector("numeric",length=n.iter)
  # Vector for saving the quality higher values
  Max.Q<-vector("numeric",length=n.iter)
  
  
  #----------------------------------------------------#
  # Step 2. for each learningset, get the Opt.C value
  #----------------------------------------------------#
  
  # For each global FOLD (as many as n.iter)
  for (k in 1:n.iter){
    
    #-------------------------------------#
    # Necessary matrices for the function
    #-------------------------------------#
    
    # Define c.vec as C.vec
    c.vec<-C.vec
    
    # folds associated to a "complete FOLD"
    CV.group<-CV.groups[((k-1)*n.fold + 1):(k*n.fold),]
    
    # Quality matrix
    Qual<-matrix(0,nrow=n.iter,ncol=length(C.vec))
    
    # Non-zero coefficients for model (used for a certain c.vec value and without the
    # corresponding fold)
    nonzero.coeff<-matrix(0,nrow=n.fold,ncol=length(C.vec))
    
    
    
    # For each c.vec value
    for (j in 1:length(c.vec)){
      
      # Load library  
      suppressMessages(require("CMA"))       # For CV groups
      suppressMessages(library("CMA"))
      
      # Print the model which the computer is working with
      cat(paste(" Computing results for TRAIN.DATASET number ",k,"\n . . . and ",j,
                "-th C.vec value\n"))
      
      # For each TRAIN dataset
      
      # Define a pred vector to save the predictions
      pred<-vector("numeric",length(y))
      for (i in 1:nrow(CV.group)){
        
        
        # Define training data.set and test data.set
        # Define training and test index
        train.idx<-CV.group[i,]; test.idx <- setdiff(1:nrow(x),train.idx)
        # Training dataset (x and y)
        x.train<-x[train.idx,]; y.train<-y[train.idx]
        # Test dataset (x and y)
        x.test<-x[test.idx,]; y.test<-y[test.idx]
        
        # If c.vec[j]==0
        if (c.vec[j]==0){
          # Model  
          model = glm(y.train~x.train, family=binomial(link="logit"))
          # Number of non-zero coefficients for the model
          nonzero.coeff[i,j]<-sum(model$coefficients!=0)-1
          
          # If accuracy required, . . . (working with predictions)
          if (eval.crit=='acc') {
            # Predictions for testx
            p = predict(model, data.frame(x.test), type="response")
            # Assign predictions in the correspondent place of pred vector
            pred[test.idx] = as.vector(p)
          } else {
            # If not accuracy required, . . .  (working with probabilities)
            # Predictions for testx
            p = predict(model, data.frame(x.test), type="response")
            # Assign predictions probabilities   
            pred[test.idx] = as.vector(p)
            # Check the lengths of results
            stopifnot(length(test.idx) == length(p))
          }  
          
        }else{
          class.weights = c(5, 1)
          # Model build from training dataset
          model = LiblineaR(x.train, y.train, type=liblinear.type, cost=c.vec[j],
                            bias=TRUE, epsilon=eps)
          # Number of non-zero coefficients for the model
          nonzero.coeff[i,j] = sum(model$W[1:(ncol(x.train)-1)] != 0)
          
          # If accuracy required, . . . (working with predictions)
          if (eval.crit=='acc') {
            # Predictions for testx
            p = predict(model, x.test, proba=FALSE)
            # Assign predictions in the correspondent place of pred vector
            pred[test.idx] = as.vector(p$predictions)
          } else {
            # If not accuracy required, . . .  (working with probabilities)
            # Predictions for testx
            p = predict(model, x.test, proba=TRUE)
            # Assign predictions probabilities   
            pred[test.idx] = as.vector(p$probabilities[,2])
            # Check the lengths of results
            stopifnot(length(test.idx) == length(p$probabilities[,2]))
          }
        }
        
      } # i
      
      # After working with all folds of a FOLD, compute "accuracy", "auroc", . . . with 
      # using pred infromation (prediction) and y (known real value)
      
      # If accuracy is required, . . .
      if (eval.crit=='acc') {
        # Percentage of good predictions 
        acc = mean(pred == y)
        # Write the result in Quality matrix
        Qual[k,j] = acc
        
        # If accuracy is not required, . . .    
        # if aupcr , . . .    
      } else if (eval.crit != 'acc') {
        # Load library
        # Detach CMA package
        detach("package:CMA", unload=TRUE)
        # Load ROCR package
        suppressMessages(require(ROCR)); suppressMessages(library(ROCR))
        suppressMessages(require(pracma)); suppressMessages(library(pracma))
        # Prediction obtained from pred (predictions) and y (labels)
        p<-prediction(predictions = pred, labels = y)
        # If "auprc" is required, . . .
        if (eval.crit=="auprc"){perf<-performance(p,measure = "prec", x.measure = "rec")
        } else if (eval.crit=="auroc"){perf<-performance(p,measure = "tpr", x.measure = "fpr")
        } else if (eval.crit=="sen.spe"){perf<-performance(p,measure = "sens",
                                                           x.measure = "spec")
        } else if (eval.crit=="lift.chart"){perf<-performance(p,measure = "lift",
                                                              x.measure = "rpp")
        } else {cat(" \n eval.crit variable is not one of the possible options.\n")}
        # Values to compute the area
        x.val<-perf@x.values[[1]]
        y.val<-perf@y.values[[1]]
        
        # If there is any NA, . . .
        if (sum(is.na(union(x.val,y.val)))>0){
          # Delete NAs
          U<-union(which(is.na(perf@x.values[[1]])),which(is.na(perf@y.values[[1]])))
          # Redefine x.values and y.values
          x.val<-perf@x.values[[1]][-U]; y.val<-perf@y.values[[1]][-U]
        }
        
        # Quality
        Qual[k,j]<-trapz(x.val,y.val)
        
      } # Quality indices
      
    } # j
    
    # Once nonzero.coeff and Quality matrix (for a specific k) are filled, Opt.C has to be chosen
    # for each iteration
    
    # Check whether these models have sufficiently many nonzero coefficients 
    # Minimum of non-zero values for a certain c.vec (between all possible folds)
    suff.nonzero = apply(nonzero.coeff, 2, min) > min.nonzero
    
    # Stop if all coefficients are 0 for porposed c.vec
    stopifnot(sum(suff.nonzero) > 0) 
    
    # Select the models (and the corresponding coefficients or values) with at least more nonzero
    # coefficients than min.nonzero
    # Select only coefficients of C.var values which supose non-zero values
    nonzero.coeff = nonzero.coeff[,suff.nonzero]
    # Select C.var values which implie non-zero values
    c.vec = c.vec[suff.nonzero]
    # Performance for nonzero values
    Q = Qual[k,suff.nonzero]
    # Index of the C.var which fits better
    opt.idx = which.max(Q)
    # Print how mani non-zero coefficients for each fold for this C.vec value
    cat(paste("\n \n  Minimum non-zero coefficients for all folds by C.vec values: \n",
              nonzero.coeff[,opt.idx]))
    # Which is c.vec associated value?
    Opt.C[k] = c.vec[opt.idx]
    # max(Q) values
    Max.Q[k]<-max(Q)
  } # k
  
  # How many seconds have been needed to run the code?
  cat('\nSuccessfully built ', n.iter*n.fold, ' LASSO models in ', proc.time()[1] - start.time,
      ' seconds\n', sep='')
  # Return (CV.groups, n.fold, n.iter and Opt.C)
  A<-list(x = x, y = y, CV.groups = CV.groups, n.fold = n.fold, n.iter = n.iter, Opt.C = Opt.C,
          Max.Q = Max.Q, min.nonzero = min.nonzero,eval.crit=eval.crit)
  return(A)
  
}

#-------------------------------------------------------------------------------------------------#
#---------------------------------------#
# AUXILIAR FUNCTION 6: myLASSOModelPlot
#---------------------------------------#
# INPUT PARAMETERS:
#         lp: object from lasso_predictor function with, . . .:
#                         * x: matrix with explanatory variables (columns) for each 
#                              individual (rows).
#                         * y: vector with the corresponding values referring to the assigned group.
#                   *  n.fold: number of fold for each Cross - Validation group.
#                   *  n.iter: number fo Cross - Validation groups.
#                   *   Opt.C: vector with the best REGULARIZATION STRENGTH values for each
#                              CV group (FOLD).
#                 * eval.crit: accuracy("acc"), auprc ("auprc"), auroc ("auroc"),
#                              sensitivity/specifity ("sens.spe"), lift chart ("lift.chart").
#               * min.nonzero: minimum coefficients != 0 for the model selection to
#                              extract the Opt.C.
#  zero.frac: maximum fraction of zeros in the models for a considered coefficient.
#       cols: vector of three colours to plot the Z-score
# OUTPUT OBJECTS:
#------------------------------------------------------------------------------------------------#
myLASSOModelPlot<-function(lp,zero.frac,cols=c("salmon4","white","seagreen"), file="Plot.pdf",sep=F){
  
  #----------------------------#
  # Step 0: extract parameters
  #----------------------------#
  
  # Parameters from lp object  
  x<-lp$x
  y<-lp$y
  lev<-levels(y)
  n.fold<-lp$n.fold
  n.iter<-lp$n.iter
  Opt.C<-lp$Opt.C
  CV.groups<-lp$CV.groups
  eval.crit<-lp$eval.crit
  min.nonzero<-lp$min.nonzero
  # Number of individuals
  n<-nrow(x)
  # Z-score limits
  z.score.lim<-c(-3,3)
  # Liblinear.type
  liblinear.type = 6       # LASSO
  eps = 1e-8               # Value for convergence
  
  
  
  #----------------------#
  # Step 1: predictions
  #----------------------#    
  
  # Define some variables
  # Quality vector (one value for each n.iter)
  Qual<-rep(0,n.iter); 
  # pred.M (Matrix for saving the predictions for each FOLD)
  pred.M<-matrix(0,nrow=n,ncol=n.iter)
  # Non-zero coefficients for model (used for a certain C.vec value and without the
  # corresponding fold)
  Coeff<-matrix(0,nrow=n.fold*n.iter,ncol=ncol(x)); colnames(Coeff)<-colnames(x)
  
  # For each CV-group
  for (k in 1:n.iter){
    # Select information
    # folds associated to a "complete FOLD"
    CV.group<-CV.groups[((k-1)*n.fold + 1):(k*n.fold),]
    
    # Optimum C
    C<-Opt.C[k]
    # Prediction vector
    pred<-rep(0,length(y))
    
    
    
    # For each fold (TRAIN i)  
    for (i in 1:nrow(CV.group)){          
      # Define training data.set and test data.set
      # Define training and test index
      train.idx<-CV.group[i,]; test.idx <- setdiff(1:nrow(x),train.idx)
      # Training dataset (x and y)
      x.train<-x[train.idx,]; y.train<-y[train.idx]
      # Test dataset (x and y)
      x.test<-x[test.idx,]; y.test<-y[test.idx]
      
      class.weights = c(5, 1)
      names(class.weights) = as.character(lev)   
      
      # Model build from training dataset
      model = LiblineaR(x.train, y.train, type=liblinear.type, cost=C,
                        bias=TRUE, epsilon=eps)
      
      # Number of non-zero coefficients for the model
      Coeff[((k-1)*n.fold + i),] = as.vector(model$W[1:(ncol(x))])
      
      # If accuracy required, . . . (working with predictions)
      if (eval.crit=='acc') {
        # Predictions for testx
        p = predict(model, x.test, proba=FALSE)
        # Assign predictions in the correspondent place of pred vector
        pred[test.idx] = as.vector(p$predictions)
      } else {
        # If not accuracy required, . . .  (working with probabilities)
        # Predictions for testx
        p = predict(model, x.test, proba=TRUE)
        # Assign predictions probabilities   
        pred[test.idx] = as.vector(p$probabilities[,2])
        # Check the lengths of results
        stopifnot(length(test.idx) == length(p$probabilities[,2]))
      }
      
    } # i
    
    # Add predictions values for k-th CV-group (FOLD) to pred matrix
    pred.M[,k]<-pred
    
    # If accuracy is required, . . .
    if (eval.crit=='acc') {
      # Percentage of good predictions 
      acc = mean(pred == y)
      # Write the result in Quality matrix
      Qual[k] = acc
      
      # If accuracy is not required, . . .    
      # if aupcr , . . .    
    } else if (eval.crit != 'acc') {
      # Load library
      # Load ROCR package
      suppressMessages(require(ROCR)); suppressMessages(library(ROCR))
      suppressMessages(require(pracma)); suppressMessages(library(pracma))
      # Prediction obtained from pred (predictions) and y (labels)
      p<-prediction(predictions = pred, labels = y)
      # If "auprc" is required, . . .
      if (eval.crit=="auprc"){perf<-performance(p,measure = "prec", x.measure = "rec")
      } else if (eval.crit=="auroc"){perf<-performance(p,measure = "tpr", x.measure = "fpr")
      } else if (eval.crit=="sen.spe"){perf<-performance(p,measure = "sens",
                                                         x.measure = "spec")
      } else if (eval.crit=="lift.chart"){perf<-performance(p,measure = "lift",
                                                            x.measure = "rpp")
      } else {cat(" \n eval.crit variable is not one of the possible options.\n")}
      # Values to compute the area
      x.val<-perf@x.values[[1]]
      y.val<-perf@y.values[[1]]
      
      # If there is any NA, . . .
      if (sum(is.na(union(x.val,y.val)))>0){
        # Delete NAs
        U<-union(which(is.na(perf@x.values[[1]])),which(is.na(perf@y.values[[1]])))
        # Redefine x.values and y.values
        x.val<-perf@x.values[[1]][-U]; y.val<-perf@y.values[[1]][-U]
      }
      
      # Quality
      Qual[k]<-trapz(x.val,y.val)
      
      
    }
    
  } # k
  
  
  #------------------------------#
  # Step 2: Variables importance
  #------------------------------#
  #-----------#
  # 1) Filter:  Delete variables whose coefficients are 0 more than a "zero.frac" proportions  
  #-----------#
  
  a<-Coeff
  
  # Selectec varaibles
  sel.var<-which(apply(Coeff==0,2,mean)<zero.frac)
  cat(paste("Variables passing the filter:", length(sel.var)))
  # Redefine Coeff matrix
  Coeff<-Coeff[,sel.var]
  
  #---------------#
  # 2) Importance:  Compute the importance of each variable (column) in each model (rows)
  #---------------#
  
  # Sum of absolute values of coefficients for each model (rows)
  Abs.sum<-apply(abs(Coeff),1,sum)
  # Weights
  Weights<-t(t(Coeff)/(Abs.sum))
  # Median of Weights for each variable  
  Imp<-apply(Weights,2,mean)
  Imp.med<-apply(Weights,2,median)
  
  cat(paste("\n Valores de Imp\n :",Imp))
  
  # Order variables by its mean importance
  Ord.w<-order(Imp,decreasing = T)
  # Order Imp values
  Imp<-Imp[Ord.w]
  Imp.med<-Imp.med[Ord.w]
  
  cat(paste("\n Valores de Imp",Imp))
  # Reorder sel.var
  sel.var<-colnames(Coeff)[Ord.w]
  # Number of selected variables
  num.sel.f<-length(sel.var)
  
  cat(paste("\n Selected variables",sel.var))
  cat(paste("\n Selected variables",names(Imp)))
  
  
  #------------------------------------------------------------------#
  # Step 3: Ordenation of individuals according to their prediction
  #------------------------------------------------------------------#
  
  # Matrix with mean prediction value (mean by individual)
  Mean.pred<-apply(pred.M,1,mean)
  # Sort individuals by the sum of the real label (y) and the mean predictor (Mean.pred)
  if (sep==F){predicti<- Mean.pred}
  else {predicti<- Mean.pred + as.numeric(y)}
  
  srt.idx = sort(predicti, index.return=TRUE)$ix
  
  
  
  # Order Mean.pred, predicti and y by srt.idx
  ord.Mean.pred<-Mean.pred[srt.idx]
  ord.predicti<-predicti[srt.idx]
  ord.y<-y[srt.idx]
  
  # To plot the heatmap
  # Select from x, "sel.var" variables and "srt.idx" ordered individuals
  heat.x<-x[srt.idx,sel.var]; cat(paste("\n \n Dimensi??n heat.x:",dim(heat.x)))
  # Standardize each variable
  m<-apply(heat.x,2,mean); cat(paste("\n \n Longitud medias: ",length(m)))
  sd<-apply(heat.x,2,sd); cat(paste("\n \n Longitud desviaciones: ",length(sd)))
  # For each individual
  for (i in 1:nrow(heat.x)){heat.x[i,]<-(heat.x[i,]-m)/sd}
  
  # Redefine heat.x not to exceed the z.score.limits
  heat.x[heat.x<z.score.lim[1]]<-z.score.lim[1]
  heat.x[heat.x>z.score.lim[2]]<-z.score.lim[2]
  
  # Load library
  suppressMessages(require(graphics))      
  suppressMessages(library(graphics))
  suppressMessages(require(RColorBrewer))      
  suppressMessages(library(RColorBrewer))
  
  # Color palette for heatmap
  color.scheme = colorRampPalette(cols)(100)
  # Names of what it supooses to be OTUs
  names<-colnames(x)
  
  # Load required libraries
  suppressMessages(require(ggplot2)); suppressMessages(library(ggplot2))
  suppressMessages(require(reshape)); suppressMessages(library(reshape))
  suppressMessages(require(gridExtra)); suppressMessages(library(gridExtra))
  suppressMessages(require(gtable)); suppressMessages(library(gtable))
  
  
  # PDF where save the plot
  pdf(file,height = 10, width = 12.5)
  
  #--------------------------#
  # Plotting options
  #--------------------------#
  
  # Maximum of 0.3 and 0.8 - 0.01*(# Selected variables)
  sel.f.cex = max(0.3, 0.8 - 0.01*num.sel.f)
  # Plot windows format
  lmat = rbind(c(4, 0, 1), c(2, 3, 5), c(0, 6, 0))
  h_t = 0.9 - 0.03*2 # dim(meta.data)[2])
  h_m = 0.1 + 0.03*2 # dim(meta.data)[2])
  h_b = 0.7 - 0.01*2 # dim(meta.data)[2])
  
  # Layout
  layout(lmat, widths=c(0.1, 0.65, 0.25), heights=c(0.1, h_t, h_m, h_b))
  # Margins for the plot (bottom, left, top, right)  
  par(oma=c(3, 4, 3, 4))
  
  
  #-----------------------#
  # Place 1: model header
  #-----------------------#
  
  # Plotting options
  par(mar=c(0, 8.1, 3.1, 1.1))
  # Plot (Empty plot)
  plot(NULL, type='n', xlim=c(-0.1,0.1), xaxt='n', xlab='',ylim=c(-0.1,0.1), yaxt='n', 
       ylab='', bty='n')
  # Text
  mtext('Linear model', side=3, line=2, at=0.04, cex=1, adj=0.5)
  mtext(paste('(Taxas = ', num.sel.f, ')', sep=''), side=3, line=1, at=0.04, cex=0.7, adj=0.5)
  
  
  #---------------------------------------------------------------------------------------#  
  # Place 2: boxplot displaying the poportion of weight per model that is actually shown
  #---------------------------------------------------------------------------------------#
  
  # Plotting options
  par(mar=c(0.1, 1.1, 0, 1.1))
  # Boxplot for weights proportion (proportion of total weights relative to selected ones ??
  # against all)
  boxplot(rowSums(abs(Coeff)) / rowSums(abs(a)), ylim=c(0,1))
  
  # Write text under the boxplot    
  mtext('Proportion of', side=1, line=1, at=1, adj=0.5, cex=0.7)
  mtext('weight shown', side=1, line=2, at=1, adj=0.5, cex=0.7)
  
  cat('  Finished plotting proportion of model weight shown.\n')
  
  
  #------------------------------------------------------------#
  # Place 3: feature heatmap with feature names to the right
  #------------------------------------------------------------#  
  
  # Plotting options  
  par(mar=c(0.1, 1.1, 0, 5.1))
  # Matrix to plot: heat.x
  
  # Load library
  suppressMessages(require(graphics))      
  suppressMessages(library(graphics))
  suppressMessages(require(RColorBrewer))      
  suppressMessages(library(RColorBrewer))
  
  
  # Heat.x as matrix
  heat.x<-as.matrix(heat.x)
  
  # Heatmap
  image(heat.x, zlim=c(z.score.lim[1],z.score.lim[2]), col=color.scheme, xaxt='n',
        yaxt='n', xlab='', ylab='')
  # Add features names
  for (f in 1:nrow(heat.x)) {
    mtext(colnames(heat.x)[f], side=4, line=1, at=(f-1)/(ncol(heat.x)-1), cex=sel.f.cex, las=2,
          col=ifelse(Imp[f]>0, "blue", "red"))
  }
  # Message
  cat('  finished plotting feature heatmap.\n')
  
  
  
  #------------------------------------------------------------------------------------------------#  
  # Place 4: additonally add to header row a corresponding color bar
  #------------------------------------------------------------------------------------------------#
  
  # Plotting options  
  par(mar=c(3.1, 1.1, 1.1, 1.1))
  
  # Barplot
  barplot(as.matrix(rep(1,100)), col=color.scheme, horiz=TRUE, border=0, ylab='', axes=FALSE)
  # Axis
  axis(side=1, at=seq(0,100,length.out=7), labels=seq(z.score.lim[1],z.score.lim[2],length.out=7))
  # Text
  mtext('Feature z-score', side=3, line=0.5, at=50, cex=0.7, adj=0.5)
  
  
  #------------------------------------------------------------------------------------------------#
  # Place 5: barplot of effect size associated with each feature
  #------------------------------------------------------------------------------------------------#
  
  # Plotting options
  par(mar=c(0.1, 12.1, 0, 1.1))
  # Minimum and maximum
  mi = (min(-Imp.med))*1.1
  mx = (max(-Imp.med))*1.1
  # Barplot
  barplot(-Imp.med, horiz=TRUE, width=1, space=0, yaxs='i', xlim=c(mi, mx), ylim=c(0, num.sel.f),
          xlab='', ylab='', yaxt='n', col=ifelse(-Imp>0,"red","deepskyblue2"),
          border=ifelse(-Imp>0,"firebrick","blue"))
  # Text
  mtext('Log-odds ratio', side=1, line=2, at=0, cex=0.7, adj=0.5)
  
  # Robustness indicated as percentage of models including a given feature 
  # For each selected feature
  for (f in 1:num.sel.f) {
    # Percentage of models inclulding this feature
    t = paste(format(100*mean(Coeff[,f] != 0), digits=1, scientific=F),
              '%', sep='')
    # Add as text into the plot  
    mtext(t, side=4, line=1.5, at=(f-0.5), cex=sel.f.cex, las=2, adj=1)
  }
  
  # Text (title)
  mtext('Effect size', side=3, line=1, at=0, cex=0.7, adj=1)
  mtext('Robustness', side=3, line=1, at=mx, cex=0.7, adj=0)
  
  # Message
  cat('  Finished plotting feature weights.\n')
  
  
  #-----------------------------------------------------------------------------------------------#  
  # Place 6 (middle): "heatmap" showing predictions and metadata (if given)
  #-----------------------------------------------------------------------------------------------# 
  
  # Plotting options
  par(mar=c(1.1, 1.1, 0.3, 5.1))
  # Define matrix
  img.data = as.matrix(cbind(ord.predicti,as.numeric(as.character(ord.y)))); colnames(img.data)<-c("Predictions","Group")
  
  
  
  # Define grays colours  
  #  grays = rev(gray(seq(0, 1, length.out=length(color.scheme))))
  Colores<-colorRampPalette(c("white", "black"))(100)
  
  
  
  # Plot the information (heatmap?)  
  image(as.matrix(img.data), col=Colores, xaxt='n', yaxt='n', xlab='', ylab='')
  # Define cex  
  meta.data.cex = max(0.3, 0.7 - 0.01*dim(img.data)[2])
  
  # Column names  
  for (m in 1:dim(img.data)[1]) {
    t = colnames(img.data)[m]
    t = gsub('\\.|_', ' ', t)
    mtext(t, side=4, line=1, at=(m-1)/(dim(img.data)[2]-1), cex=meta.data.cex, las=2)
  }
  
  # Mark missing values
  for (m in 1:dim(img.data)[1]) {
    idx = which(is.na(img.data[m,]))
    for (i in idx) {
      x = (i-1) / (dim(img.data)[1]-1)
      y = (m-1) / (dim(img.data)[2]-1)
      text(x, y, 'NA', col='red', cex=0.4)
    }
  }
  
  
  
  
  
  #-----------------------------------------------------------------------------------------------#  
  # Place 7 (right): AUC, roc curve and confidence intervals
  #-----------------------------------------------------------------------------------------------# 
  
  # Default graphical parameters
  par(mfrow=c(1,1))
  
  # Load library
  suppressMessages(library(pROC))
  
  # Start a ROC plot
  rocobj <- plot.roc(as.numeric(y),predicti,legacy.axes=T,xlim=c(1,0),
                     main="Cross-validation ROC for Cluster classification",
                     xlab="False positive rate",ylab="True positive rate")
  
  # CI
  ci.se.obj <- ci(rocobj, of="se", boot.n=5000,legacy.axes=F)
  plot(ci.se.obj,type="shape",col="grey87",border="white")
  text(x = 0.35,y=0.05,paste("Mean prediction AUC = ",round(rocobj$auc,2)))
  abline(1,-1)
  
  
  # Message  
  cat('  Finished plotting classification result and additional metadata.\n')
  
  # Close PDF file
  dev.off()
  
  # Print message:
  cat(paste("\n Results saved as 'Plot.pdf' in ",getwd()))
  
  
  
  return(list(a,ord.y,ord.predicti))
  
}


#########################
#### LefSE
#### From: https://rdrr.io/github/ying14/yingtools2/src/R/microbiota2.R
##########################
#### Call as:
 # lefse(phy, class, subclass = NA, subject = NA, anova.alpha = 0.05,
 #      wilcoxon.alpha = 0.05, lda.cutoff = 2,
 #      wilcoxon.within.subclass = FALSE, one.against.one = FALSE,
 #      mult.test.correction = 0, make.lefse.plots = FALSE,
 #      by_otus = FALSE, levels = rank_names(phy))


#' LEfSe
#'
#' Run LEfSe (LDA Effect Size) analysis using a phyloseq object.
#'
#' This function performs the analysis using the following steps.
#' (1) Creates lefse.txt from phyloseq data, a tab-delimited file in the format input required by LEfSe.
#' (2) Executes format_input.py to further format lefse.txt into lefse.in
#' (3) Executes run_lefse.py, which does the actual analysis and produces lefse.res.
#' (4) Executes plot_res.py and plot_cladogram.py, which create the graphics for LEfSe.
#' Note that you must have command-line lefse.py installed... this function is just an R wrapper for the original Huttenhower scripts.
#' @param phy the phyloseq object containing data
#' @param class variable to be tested by LEfSe. This must be a variable in sample_data(phy)
#' @param subclass variable to perform subclass testing. This step is skipped if it is not specified.
#' @param subject variable referring to the subject level designation. This is only necessary if multiple samples per subject.
#' @param anova.alpha alpha level of the kruskal-wallis testing. Default is 0.05
#' @param wilcoxon.alpha alpha level at which to perform wilcoxon testing of subclass testing. Default is 0.05.
#' @param lda.cutoff Cutoff LDA to be reported. Default is 2.0.
#' @param wilcoxon.within.subclass Set whether to perform Wilcox test only among subclasses with the same name (Default FALSE)
#' @param mult.test.correction Can be {0,1,2}. Set the multiple testing correction options. 0 no correction (more strict, default), 1 correction for independent comparisons, 2 correction for independent comparison
#' @param one.against.one for multiclass tasks, sets whether testing is performed one-against-one (TRUE - more strict) or one-against-all (FALSE - less strict)
#' @param levels Taxonomic levels to be tested. Default is to test all levels: rank_names(phy)
#' @return Returns data
#' @examples
#' lefse.tbl <- lefse(ph,class="CDI",subclass="Sex")
#' @author Ying Taur
#' @export
myLefse <- function(phy,class,subclass=NA,subject=NA,
                  anova.alpha=0.05,wilcoxon.alpha=0.05,lda.cutoff=2.0,
                  wilcoxon.within.subclass=FALSE,one.against.one=FALSE,
                  mult.test.correction=0,
                  make.lefse.plots=FALSE,by_otus=FALSE,
                  levels=rank_names(phy)) {
  #phy=ph.lefse;class="CDI";subclass=NA;subject=NA;anova.alpha=0.05;wilcoxon.alpha=0.05;lda.cutoff=2.0;wilcoxon.within.subclass=FALSE;one.against.one=FALSE;levels=rank_names(phy)
  #phy=ph.lefse;class="CDI";subclass=NA;subject=NA;anova.alpha=0.05;wilcoxon.alpha=0.05;lda.cutoff=2.0;wilcoxon.within.subclass=FALSE;one.against.one=FALSE;levels=rank_names(phy)
  #phy=ph.lefse;class="CDI";subclass="SampleType";subject="MRN";anova.alpha=0.05;wilcoxon.alpha=0.05;lda.cutoff=2.0;wilcoxon.within.subclass=FALSE;one.against.one=FALSE;levels=rank_names(phy)
  requireNamespace("phyloseq")
  keepvars <- c(class,subclass,subject,"sample")
  keepvars <- unique(keepvars[!is.na(keepvars)])
  samp <- get.samp(phy)[,keepvars]
  # note that lefse taxa names cannot have spaces, will replace ; and = with _
  # gsub("[;=]","_",xx$otu)
  if (by_otus) { #perform by otu only
    
    otu <- get.otu.melt(phy,sample_data=FALSE)
    otu.levels <- otu %>% mutate(taxon=otu) %>%
      group_by(sample,taxon) %>% summarize(pctseqs=sum(pctseqs)) %>%
      mutate(taxon=gsub(" ","_",taxon))
  } else { #divide by taxonomy
    otu <- get.otu.melt(phy,sample_data=FALSE)
    otu.list <- lapply(1:length(levels),function(i) {
      lvls <- levels[1:i]
      lvl <- levels[i]
      otu.level <- otu
      otu.level$taxon <- do.call(paste,c(lapply(lvls,function(l) otu[[l]]),sep="|"))
      otu.level$rank <- lvl
      otu.level2 <- otu.level %>% group_by(sample,taxon,rank) %>% summarize(pctseqs=sum(pctseqs)) %>% ungroup()
      return(otu.level2)
    })
    otu.levels <- bind_rows(otu.list) %>%
      mutate(taxon=gsub(" ","_",taxon))
  }
  
  otu.tbl <- otu.levels %>%
    dcast(sample~taxon,value.var="pctseqs",fill=0) %>%
    left_join(samp,by="sample") %>%
    select_(.dots=c(keepvars,lazyeval::interp(~everything())))
  if (is.na(subject) | subject!="sample") {
    otu.tbl <- otu.tbl %>% select(-sample)
  }
  tbl <- otu.tbl %>% t()
  write.table(tbl,"lefse.txt",quote=FALSE,sep="\t",col.names=FALSE)
  
  opt.class <- paste("-c",which(keepvars %in% class))
  opt.subclass <- ifelse(is.na(subclass),"",paste("-s",which(keepvars %in% subclass)))
  opt.subject <-ifelse(is.na(subject),"",paste("-u",which(keepvars %in% subject)))
  format.command <- paste("format_input.py lefse.txt lefse.in",opt.class,opt.subclass,opt.subject,"-o 1000000")
  system(format.command)
  #   -m {f,s}              set the policy to adopt with missin values: f removes
  #   the features with missing values, s removes samples
  #   with missing values (default f)
  #   -n int                set the minimum cardinality of each subclass
  #   (subclasses with low cardinalities will be grouped
  #   together, if the cardinality is still low, no pairwise
  #   comparison will be performed with them)
  
  lefse.command <- paste("run_lefse.py lefse.in lefse.res",
                         "-a",anova.alpha,
                         "-w",wilcoxon.alpha,
                         "-l",lda.cutoff,
                         "-e",as.numeric(wilcoxon.within.subclass),
                         "-y",as.numeric(one.against.one),
                         "-s",mult.test.correction)
  system(lefse.command)
  print("Wrote lefse.res")
  lefse.out <- read.table("lefse.res",header=FALSE,sep="\t") %>% rename(taxon=V1,log.max.pct=V2,direction=V3,lda=V4,p.value=V5)
  #   -a float        set the alpha value for the Anova test (default 0.05)
  #   -w float        set the alpha value for the Wilcoxon test (default 0.05)
  #   -l float        set the threshold on the absolute value of the logarithmic
  #   LDA score (default 2.0)
  #   --nlogs int     max log ingluence of LDA coeff
  #   --verbose int   verbose execution (default 0)
  #   --wilc int      wheter to perform the Wicoxon step (default 1)
  #   -r str          select LDA or SVM for effect size (default LDA)
  #   --svm_norm int  whether to normalize the data in [0,1] for SVM feature
  #   waiting (default 1 strongly suggested)
  #   -b int          set the number of bootstrap iteration for LDA (default 30)
  #   -e int          set whether perform the wilcoxon test only among the
  #   subclasses with the same name (default 0)
  #   -c int          set whether perform the wilcoxon test ing the Curtis's
  #                   approach [BETA VERSION] (default 0)
  #   -f float        set the subsampling fraction value for each bootstrap
  #                   iteration (default 0.66666)
  #   -s {0,1,2}      set the multiple testing correction options. 0 no correction
  #                   (more strict, default), 1 correction for independent
  #                   comparisons, 2 correction for independent comparison
  #   --min_c int     minimum number of samples per subclass for performing
  #                   wilcoxon test (default 10)
  #   -t str          set the title of the analysis (default input file without
  #                   extension)
  #   -y {0,1}        (for multiclass tasks) set whether the test is performed in
  #                   a one-against-one ( 1 - more strict!) or in a one-against-
  #                   all setting ( 0 - less strict) (default 0)
  
  if (make.lefse.plots) {
    system("plot_res.py lefse.res lefse_lda.png")
    print("Wrote lefse_lda.png")
    system("plot_cladogram.py lefse.res lefse_clado.pdf --format pdf")
    print("Wrote lefse_clado.pdf")
  }
  return(lefse.out)
}


#' Extract Phyloseq sample_data, companion function for myLefse
#'
#' Returns \code{sample_data} component from phyloseq object, as a data frame.
#'
#' This basically is similar to the function \code{phyloseq::sample_data}, but does a few extra things.
#' (1) Converts to a data frame
#' (2) The sample name is stored in a column called \code{sample}. (\code{phyloseq} normally stores as a row name)
#' (3) Calculates number of sequences and alpha diversity metrics, if desired.
#' This function is the opposite of \code{set.samp}, which converts the data frame back into a \code{sample_data}.
#' Note that if the \code{phyloseq} object does not contain \code{sample_data}, a data frame containing a single column, \code{sample}, is returned.
#' @param phy phyloseq object containing sample_data
#' @param stats logical, whether or not to include summary statistics of samples. Stores \code{nseqs}, and diversity metrics.
#' @param measures diversity measures to calculate, if stats is TRUE. Default: c("Observed","InvSimpson","Shannon")
#' @return Data frame containing \code{sample_data} data.
#' @export
get.samp <- function(phy,stats=FALSE,measures=c("Observed","InvSimpson","Shannon")) {
  requireNamespace(c("phyloseq","tibble"))
  
  if (is.null(sample_data(phy,FALSE))) {
    #if no sample_data, return single data frame with sample column
    sdata <- data.frame(sample=sample_names(phy),stringsAsFactors=FALSE)
  } else {
    if ("sample" %in% phyloseq::sample_variables(phy)) {stop("YTError: phyloseq sample_data already contains the reserved variable name \"sample\"")}
    sdata <- sample_data(phy) %>% data.frame(stringsAsFactors=FALSE) %>% tibble::rownames_to_column("sample")
  }
  if (stats) {
    dup.names <- intersect(c("nseqs",measures),names(sdata))
    if (length(dup.names)>0) {
      sdata <- sdata[,setdiff(names(sdata),dup.names)]
      warning("YTWarning: Following variables are duplicated. Deleting old values from phyloseq: ",paste(dup.names,collapse=", "))
    }
    sdata$nseqs <- phyloseq::sample_sums(phy)
    sdata <- cbind(sdata,estimate_richness(phy,measures=measures))
  }
  return(sdata)
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
myMultiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#### Multiomic data list, subset sample function
mySubsetSamplesFromList<-function(myList,sampleNames){
  return(lapply(myList, function(x) return(x[c(sampleNames),])))
}

#### Returns a color vector mapped to a numeric vector x
myMap2Color <- function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}


SparCC.count <- function(x, imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  # Dimension for w (latent variables)
  p <- ncol(x);
  n <- nrow(x);
  # Posterior distribution (alpha)
  x <- x + 1;
  # Store generate data
  y <- matrix(0, n, p);
  # Store covariance/correlation matrix
  cov.w <- cor.w <- matrix(0, p, p);
  indLow <- lower.tri(cov.w, diag = T);
  # Store covariance/correlation for several posterior samples
  covs <- cors <- matrix(0, p * (p + 1) / 2, imax);
  for(i in 1:imax) {
    # Generate fractions from posterior distribution
    y <- t(apply(x, 1, function(x) 
      gtools::rdirichlet(n = 1, alpha = x)));
    # Estimate covariance/correlation
    cov_cor <- SparCC.frac(x = y, kmax = kmax, alpha = alpha, Vmin = Vmin);
    # Store variance/correlation only low triangle 
    covs[, i] <- cov_cor$cov.w[indLow];
    cors[, i] <- cov_cor$cor.w[indLow];
  }
  # Calculate median for several posterior samples
  cov.w[indLow] <- apply(covs, 1, median); 
  cor.w[indLow] <- apply(cors, 1, median);
  #
  cov.w <- cov.w + t(cov.w);
  diag(cov.w) <- diag(cov.w) / 2;
  cor.w <- cor.w + t(cor.w);
  diag(cor.w) <- 1;
  #
  return(list(cov.w = cov.w, cor.w = cor.w));
}
#------------------------------------------------------------------------------#
# SparCC for fractions known
#   FUNCTION: SparCC.frac
#   INPUT:
#          x: nxp fraction data matrix, row is sample, col is variable
#       kmax: max iteration steps for SparCC. default is 10
#      alpha: the threshold for strong correlation. default is 0.1
#       Vmin: minimal variance if negative variance appears. default is 1e-4
#   OUTPUT: a list structure
#      cov.w: covariance estimation
#      cor.w: correlation estimation
#------------------------------------------------------------------------------#


SparCC.frac <- function(x, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  # Log transformation
  x <- log(x);
  # Number of variables  
  p <- ncol(x);
  # Variation matrix ( T0 = var(log(xi/xj)) ) function from stats library
  TT <- stats::var(x);
  T0 <- diag(TT) + rep(diag(TT), each = p) - 2 * TT;
  # Variance and correlation coefficients for Basic SparCC  
  rowT0 <- rowSums(T0);
  var.w <- (rowT0 - sum(rowT0) / (2 * p - 2))/(p - 2);
  var.w[var.w < Vmin] <- Vmin;
  
  #cor.w <- (outer(var.w, var.w, "+") - T0 ) / 
  #  sqrt(outer(var.w, var.w, "*")) / 2;
  
  Is <- sqrt(1/var.w);
  cor.w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5;
  # Truncated correlation in [-1, 1]
  cor.w[cor.w <= - 1] <- - 1; 
  cor.w[cor.w >= 1] <- 1;
  # Left matrix of estimation equation
  Lmat <- diag(rep(p - 2, p)) + 1; 
  # Remove pairs
  rp <- NULL;
  # Left components
  cp <- rep(TRUE, p);
  # Do loops until max iteration or only 3 components left
  k <- 0;  
  while(k < kmax && sum(cp) > 3) {
    # Left T0 = var(log(xi/xj)) after removing pairs
    T02 <- T0;
    # Store current correlation to find the strongest pair
    curr_cor.w <- cor.w;
    # Remove diagonal
    diag(curr_cor.w) <- 0;
    # Remove removed pairs
    if(!is.null(rp)) {
      curr_cor.w[rp] <- 0;
    }
    # Find the strongest pair in vector form
    n_rp <- which.max(abs(curr_cor.w));
    # Remove the pair if geater than alpha
    if(abs(curr_cor.w[n_rp]) >= alpha) {
      # Which pair in matrix form
      t_id <- c(arrayInd(n_rp, .dim = c(p, p)));
      Lmat[t_id, t_id] <- Lmat[t_id, t_id] - 1;
      # Update remove pairs
      n_rp <- c(n_rp, (p + 1) * sum(t_id) - 2 * p - n_rp);
      rp <- c(rp, n_rp);
      # Update T02
      T02[rp] <- 0;
      # Which component left
      cp <- (diag(Lmat) > 0);
      # Update variance and truncated lower by Vmin
      var.w[cp] <- solve(Lmat[cp, cp], rowSums(T02[cp, cp]));
      var.w[var.w <= Vmin] <- Vmin;
      # Update correlation matrix and truncated by [-1, 1]
      #cor.w <- (outer(var.w, var.w, "+") - T0 ) / 
      #  sqrt(outer(var.w, var.w, "*")) / 2;    
      Is <- sqrt(1/var.w);
      cor.w <- (var.w + rep(var.w, each = p) - T0) * 
        Is * rep(Is, each = p) * 0.5;
      # Truncated correlation in [-1, 1]
      cor.w[cor.w <= - 1] <- - 1;
      cor.w[cor.w >= 1] <- 1;
    }
    else {
      break;
    }
    # 
    k <- k + 1;
  }
  # Covariance
  Is <- sqrt(var.w);
  cov.w <- cor.w * Is * rep(Is, each = p);
  # Return the values
  return(list(cov.w = cov.w, cor.w = cor.w));
}
#------------------------------------------------------------------------------#
#SparCC.count(otu_table(myDADA2Phyloseq))




taxa_barplot<-function(Abundance_table,Taxonomy_table=NULL, metadata, tax_level=NULL, Nspecies,idvar, facet_vars = NULL){
  
  #' Creates a facetted composition barplot from abundance data, at the specified level.
  #' @Abundance_table Numerical matrix. For taxonomy abundance, relative abundance is recommended, although it teorically work anyway
  #' @Taxonomy_table Taxonomic assignment matrix just as the ones from Phyloseq. It can be set as NULL to just get the names from Abundance-table
  #' @Metadata The metadata matrix or Sample_data from phyloseq
  #' @tax_level if a Taxonomy_table is provided, which column (corresponds to taxonomic level) to take. One of the values from col.names(Taxonomy_table)
  #' @Nspecies The top N species to include in the plot
  #' @idvar Name of the variable which identifies the individuals (not the samples, but the patient itself)
  #' @facet_vars name of the variables from which the facetting formula will be constructed. If NULL, no facetting will be done. Otherwise it will facet following the variable levels. It can take up to two variables
  
  
  set.seed(1234)
  # Merge id and faceting variables to include in the plot
  VarsToChoose<-c(idvar,facet_vars)
  
  # Get number of species 
  
  TopNspecies<-as.vector(rownames(data.frame(sort(colSums(Abundance_table),decreasing=T)[0:Nspecies])))
  
  # Construct the inner data frame, cureate and melt
  
  if (!is.null(Taxonomy_table)){
    myInnerDataFrame<- merge( metadata[,VarsToChoose, drop=F],Abundance_table[,TopNspecies], by = "row.names",all = T) %>%
      
      transform(.,row.names=Row.names, Row.names=NULL) %>%
      
      rename_at(TopNspecies, ~ Taxonomy_table[TopNspecies,tax_level]) %>% 
      
      melt(id.vars = c(VarsToChoose)) %>%
      
      dplyr::rename(Species = variable, 
                    Abundance = value,
                    id_var = paste(idvar)) %>%
      
      mutate(id_var = as.factor(id_var))
  } else{
    myInnerDataFrame<- merge( metadata[,VarsToChoose, drop=F],Abundance_table[,TopNspecies], by = "row.names",all = T) %>%
      
      transform(.,row.names=Row.names, Row.names=NULL) %>%
      
      melt(id.vars = c(VarsToChoose)) %>%
      
      dplyr::rename(Species = variable, 
                    Abundance = value,
                    id_var = paste(idvar)) %>%
      
      mutate(id_var = as.factor(id_var))
  } 
  
  # # Taxa Level stands for taxonomical number (1 = Kingdom to 7 = Species)
  # tax_name<- tax_level
  
  
  # Extract and mix palettes from ColorBrewer
  
  qpalettes<- brewer.pal.info[which(brewer.pal.info$category=="qual"),]
  color_vector <- unlist(mapply(brewer.pal, qpalettes$maxcolors, rownames(qpalettes)))
  
  # Create palette according to the number of species
  n<-length(unique(TopNspecies))
  
  if (n<length(color_vector)){
    palette<-sample(color_vector, n, replace = F)
  }else{
    palette<-sample(color_vector, n, replace = T)
  }
  
  
  
  # Build the faceting formula in case it was inputted
  if (!is.null(facet_vars)){
    facet_formula<-formula(paste(facet_vars[2],"~",facet_vars[1]))
    
  }else{facet_formula<-NULL}
  
  
  barplot<- ggplot(myInnerDataFrame, aes(x=id_var,y=Abundance, fill=Species)) + 
    geom_bar(aes(color=Species, fill=Species), stat="identity", position="stack")+
    facet_grid(facets =  facet_formula, scales="free_x", space= "free",drop = TRUE)+
    theme_bw()+
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.text = element_text(size=10),
          legend.key.size = unit(0.7,"line"),
          axis.text.x = element_text(size=8,angle = 90 )) +
    scale_discrete_manual(aesthetics=c("color","fill"),values=palette) + 
    guides(fill=guide_legend(nrow=80, byrow=TRUE), col=guide_legend(nrow=40,byrow=TRUE)) + 
    labs(title = paste ("Abundance plot for top", Nspecies,tax_level, sep=" "), x= paste(idvar))
  
  
  return(barplot)
}


###### Wilcoxon test with hability to handle missing data in paired tests (need more work)
DA_wilcoxon<-function(AbundanceDF, TaxonomyDF, metadata, LinkVariable, SubsetVariable, TestVariable, VarSubset, TestSubset, threshold  = 1, is.longitudinal = F){
  # SubsetVariable<-CategoricalVariable
  # TestVariable<-LongitudinalVariable
  
  if (!is.null(TestSubset)){
    
    myVarDF<-metadata[c(LinkVariable,SubsetVariable, TestVariable)] %>%
      
      rename_all(., ~ c("LinkVariable","SubsetVariable","TestVariable")) %>%
      
      subset(TestVariable %in% TestSubset) %>%
      
      subset(SubsetVariable %in% VarSubset)
    
  } else {
    
    myVarDF<-metadata[c(LinkVariable,SubsetVariable, TestVariable)] %>% 
      
      rename_all(., ~ c("LinkVariable","SubsetVariable","TestVariable")) %>%
      
      subset(SubsetVariable %in% VarSubset)
    
  }
  
  
  # Remove samples without both TPs
  
  if (is.longitudinal == T){
    print(paste("is.longitudinal set as ", is.longitudinal,". This function will perform paired tests. Please, check if this is what you wanted", sep=))
    
    SelectedIDs<-c()
    for (id in levels(as.factor(myVarDF$LinkVariable))){
      if(nrow(myVarDF[which(myVarDF$LinkVariable==id),]) == 2){
        
        SelectedIDs<-c(SelectedIDs,id)
      }
      else{ print(paste("Individual with ID:",id, "didn't have all timepoints, discarded"))}
    } 
    
  } else{
    print(paste("is.longitudinal set as ", is.longitudinal,". This function will perform unpaired tests. Please, check if this is what you wanted", sep=))
    SelectedIDs<-myVarDF$LinkVariable
  }
  
  if (!is.null(TaxonomyDF)){
    myInnerDataFrame<-as.data.frame(myAbundanceDF) %>% rename_all(., ~ myTaxonomyDF$Species) %>%
      merge(., myVarDF,  by="row.names") %>%
      column_to_rownames("Row.names") %>%
      subset(LinkVariable %in% SelectedIDs)
  } else {
    myInnerDataFrame<-as.data.frame(myAbundanceDF) %>%
      merge(., myVarDF,  by="row.names") %>%
      column_to_rownames("Row.names") %>%
      subset(LinkVariable %in% SelectedIDs)
  }
  
  
  
  
  resultsDF<-data.frame()
  for (i in colnames(myInnerDataFrame)[1:ncol(myAbundanceDF)]){
    
    
    if (is.longitudinal ==T){
      
      testDF <- myInnerDataFrame[,c(i, "TestVariable", "LinkVariable")] %>% dplyr::rename("Abundance"=i) %>% mutate(TestVariable = as.factor(TestVariable)) %>% arrange(LinkVariable)
      myWilcox<- pairwise.wilcox.test( testDF$Abundance, g=testDF$TestVariable,p.adjust.method = "none", paired=T)
      
    } else{
      testDF <- myInnerDataFrame[,c(i, "TestVariable", "LinkVariable")] %>% dplyr::rename("Abundance"=i) %>% mutate(TestVariable = as.factor(TestVariable)) %>% arrange(LinkVariable)
      myWilcox<- pairwise.wilcox.test( testDF$Abundance, g=testDF$TestVariable,p.adjust.method = "none", paired=F)
      
    }
    
    resultsDF<-rbind(resultsDF,   melt(myWilcox$p.value) %>% mutate(name=i))
    # print(i)
    # print(myWilcox)
    
  }
  
  resultsDF<-resultsDF[,c(4,1,2,3)] %>% rename_all(., ~ c("Species","Group1","Group2","p.val")) %>%
    subset(.,p.val != "NaN") %>%
    mutate(p.val = as.numeric(p.val))%>%
    mutate(p.adj = p.adjust(p.val, method="fdr")) %>%
    # mutate(p.adj.method = myWilcox$p.adjust.method) %>%
    subset(., p.val <= threshold)
  
  return(resultsDF)
  
}




