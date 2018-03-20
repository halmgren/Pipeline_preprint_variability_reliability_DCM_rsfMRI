Work_dir='/media/hannes/Almgren_Disk3';

setwd(paste(Work_dir,'/Figures_paper_variability',sep=''));
library('R.matlab');
library(igraph);

A=readMat('Figure2_matrices.mat');

#vectorize matrices
S1_vec=as.vector(A$S1);
S2_vec=as.vector(A$S2);
S3_vec=as.vector(A$S3);
S4_vec=as.vector(A$S4);

for (sub in c(1,2,3,4)){
  print(sub)
  if(sub==1){
    Edges <- data.frame(
      from = rep(c('Prec','lIPC','mPFC','rIPC'),each=4),
      to = rep(c('Prec','lIPC','mPFC','rIPC'),times=4),
      thickness = abs(S1_vec),
      EC_value=S1_vec);}
  
  if(sub==2){
    Edges <- data.frame(
      from = rep(c('Prec','lIPC','mPFC','rIPC'),each=4),
      to = rep(c('Prec','lIPC','mPFC','rIPC'),times=4),
      thickness = abs(S2_vec),
    EC_value=S2_vec);}
  
  if(sub==3){
    Edges <- data.frame(
      from = rep(c('Prec','lIPC','mPFC','rIPC'),each=4),
      to = rep(c('Prec','lIPC','mPFC','rIPC'),times=4),
      thickness = abs(S3_vec),
      EC_value=S3_vec);}
  
  if(sub==4){
    Edges <- data.frame(
      from = rep(c('Prec','lIPC','mPFC','rIPC'),each=4),
      to = rep(c('Prec','lIPC','mPFC','rIPC'),times=4),
      thickness = abs(S4_vec),
      EC_value=S4_vec);}
  
  g <- graph.edgelist(as.matrix(Edges[,c(-3,-4)]));
  
  l <- layout.fruchterman.reingold(g); #location of plotted region
  
  l[1,1]=2;
  l[1,2]=3.8;
  l[2,1]=1.5;
  l[2,2]=4.2;
  l[3,1]=2;
  l[3,2]=4.8;
  l[4,1]=2.5;
  l[4,2]=4.2;
  
  c_scale <- colorRamp(c('blue','violetred','red'));
  E(g)$weight <- Edges$thickness * 0.75;
  E(g)$width <- Edges$thickness * 12;
  E(g)$loop.angle <- 0;
  E(g)$loop.angle[6]=3.14;
  E(g)$curved <- 0.3;
  E(g)$color = apply(c_scale(E(g)$weight), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) );
  
  
  
  V(g)$label.cex=2;
  V(g)$label.font=2;
  
  dev.new();
  g1=g;
  g1<-delete.edges(g1,E(g1)[S1_vec<0]);
  
  g2=g;
  g2<-delete.edges(g2,E(g2)[S1_vec>0]);
  
  c_scale_red<-colorRamp(c('red','darkred'));
  c_scale_blue<-colorRamp(c('light blue','dark blue'));
  
  E(g1)$color=apply(c_scale_red(E(g1)$weight),1,function(x) rgb(x[1]/255,x[2]/255,x[3]/255) );
  E(g2)$color=apply(c_scale_blue(E(g2)$weight),1,function(x) rgb(x[1]/255,x[2]/255,x[3]/255) );
  
  for (edge in 1:length(E(g1))){
    g1_tmp=delete.edges(g1,E(g1)[(1:length(E(g1)))[-edge]]);
    
    if(edge==1){
      E(g1_tmp)$arrow.width<-E(g1)$weight[edge]*4;
      plot.igraph(g1_tmp,vertex.size=60,layout=l,edge.arrow.size=E(g1)$thickness[edge]);
    }
    
    if(edge!=1){
      E(g1_tmp)$arrow.width<-E(g1)$weight[edge]*4;
      plot.igraph(g1_tmp,vertex.size=60,layout=l,add=TRUE);
    }
    
  }
  
  for (edge in 1:length(E(g2))){
    g2_tmp=delete.edges(g2,E(g2)[(1:length(E(g2)))[-edge]]);
    
      
        E(g2_tmp)$arrow.width<-E(g2)$weight[edge]*4;
        plot.igraph(g2_tmp,vertex.size=60,layout=l,add=TRUE);
      
    print(E(g2)$weight[edge]*4)
  }
  dev.print(bmp,filename=paste('DMN_figure_sub_',sub,'.bmp',sep=''),width = 1000, height = 500);
  dev.off;
  
}
