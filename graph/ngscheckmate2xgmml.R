
# randomly generate n distinct colors
distinct.colors<-function(n=40){

 n.pre=n*10 # number of colors to pre-generate before filtering to n colors.

 # 1. randomly generate n.pre colors.
 r=runif(n.pre)
 g=runif(n.pre)
 b=runif(n.pre)
 color=as.matrix(data.frame(r,g,b))

 # 2. iteratively remove a color that has the smaller minimum (Euclidean) distance in rgb space to other colors
 for(T in 1:(n.pre-n)){
  d=dist(color)
  to.be.removed = which.min(apply(as.matrix(d),1,function(xx)min(xx[order(xx)[2:length(xx)]])))
  color=color[setdiff(1:nrow(color),to.be.removed),]
 }

 # 3. convert to rgb.
 rgb=apply(color,1,function(xx)rgb(xx[1],xx[2],xx[3]))
 return(rgb)
}

## given hex code for a node color, determine a node label color as hex code (either black or white)
get.label.colors<-function(node.colors){
    node.colors.rgb=col2rgb(node.colors)
    apply(node.colors.rgb,2, function(node.rgb) {
       dist.to.black = dist(rbind(node.rgb,c(0,0,0)))
       dist.to.white = dist(rbind(node.rgb,c(255,255,255)))
       if(dist.to.black >= dist.to.white) { return("#000000") } else { return("#ffffff") }
    })
}

## This function takes nodes and edges (lists) and output xgmml file name as arguments and generates an xgmml file.
## The function requires four xgmml.template.xxx files in the current working directory.
## currently it requires two functions: distinct.colors and get.label.colors.
generate.xgmml<-function(nodes,edges,outfile){
   
   nNodes=length(nodes[[1]])
   nEdges=length(edges[[1]])
   
   # read node template 
   node.lines = paste(readLines("xgmml.template.node"),collapse="\n")
   
   # create an xgmml code for every node
   xgmml.lines.all.nodes = rep(node.lines,nNodes)
   for(i in 1:nNodes){
     for(att in names(nodes)){
       xgmml.lines.all.nodes[i] = gsub(paste("###",att,"###",sep=""),nodes[[att]][i],xgmml.lines.all.nodes[i],fixed=T)
     }
   }
   
   # read edge template 
   edge.lines = paste(readLines("xgmml.template.edge"),collapse="\n")
   
   # create an xgmml code for every edge
   xgmml.lines.all.edges = rep(edge.lines,nEdges)
   for(i in 1:nEdges){
     for(att in names(edges)){
       xgmml.lines.all.edges[i] = gsub(paste("###",att,"###",sep=""),edges[[att]][i],xgmml.lines.all.edges[i],fixed=T)
     }
   }
   
   # final xgml
   xgmml.1.lines = paste(readLines("xgmml.template.1"),collapse="\n")
   xgmml.2.lines = paste(readLines("xgmml.template.2"),collapse="\n")
   
   # date
   xgmml.1.lines = gsub("###DATE###",Sys.time(),xgmml.1.lines)
   
   
   write(c(xgmml.1.lines,xgmml.lines.all.nodes,xgmml.lines.all.edges,xgmml.2.lines),outfile)
   
}

## test example
#
#nodes=list()
#nodes$NODE_ID=1:3
#nodes$NODE_LABEL=paste("mynode",1:3,sep=":")
#nodes$NODE_BORDERCOLOR_R=c(0,0,0)
#nodes$NODE_BORDERCOLOR_G=c(0,0,0)
#nodes$NODE_BORDERCOLOR_B=c(0,0,0)
#nodes$NODE_BORDERCOLOR_HEX=c("#000000","#000000","#000000")
#nodes$NODE_COLOR_R=c(255,0,0)
#nodes$NODE_COLOR_G=c(0,255,0)
#nodes$NODE_COLOR_B=c(0,0,255)
#nodes$NODE_COLOR_HEX=c("#ff0000","#00ff00","#0000ff")
#nodes$NODE_LABELCOLOR_R=c(0,0,255)
#nodes$NODE_LABELCOLOR_G=c(0,0,255)
#nodes$NODE_LABELCOLOR_B=c(0,0,255)
#nodes$NODE_BORDER_WIDTH=c(1.0,1.0,5.0)
#nodes$NODE_X=c(-10,0,10)
#nodes$NODE_Y=c(-10,0,10)
#
#edges=list()
#edges$SOURCE_ID=1:2
#edges$TARGET_ID=c(1,1)
#edges$EDGE_LABEL=c("myedge1","myedge2")
#edges$EDGE_WIDTH=c(0.3,5.0)
#
#outfile="test.out.xgmml"
#
#generate.xgmml(nodes,edges,outfile)


## This function reads the NGSCheckMate output and turning it into nodes and edges.
## The function requires another function generate.xgmml
create.xgmml.from.ngscheckmateout<-function(label.file, ngscheckmateoutput.file, output.xgmml){
   labels = read.table(label.file,stringsAsFactors=F)
   if(ncol(labels)==2) labels[,3]=labels[,1]
   colnames(labels)=c("file","individual","tag")
   rownames(labels)=labels$file
   x = read.table(ngscheckmateoutput.file,stringsAsFactors=F)
   x=cbind(x,labels[x[,1],c("individual","tag")],labels[x[,3],c("individual","tag")])
   colnames(x)=c("file1","call","file2","r","depth","individual1","tag1","individual2","tag2")
   
   # all individuals that has ever been matched to a wrong individual
   individuals.containing.incorrect.matched = unique(c(x[which(x[,2]=="matched" & x[,6]!=x[,8]),6],x[which(x[,2]=="matched" & x[,6]!=x[,8]),8]))
   individuals.containing.incorrect.matched = unique(c(x[which(x$call=="matched" & x$individual1!=x$individual2),"individual1"],x[which(x$call=="matched" & x$individual1!=x$individual2),"individual2"]))
   # 21 individuals for luad usecase
   
   # all files associated with these individuals
   labels.filtered=labels[labels$individual%in%individuals.containing.incorrect.matched,]
   nNodes=nrow(labels.filtered)
   # 45 files for luad usecase (all of them will be included as nodes)
   
   # node IDs
   labels.filtered$nodeID=1:nNodes
   
   # all edges (matching pairs) associated with these individuals
   x.filtered = x[which((x$individual1%in%individuals.containing.incorrect.matched | x$individual2%in%individuals.containing.incorrect.matched) & x$call=="matched"),]
   nEdges=nrow(x.filtered)
   # 62 edges for luad usecase (all of them will be included as edges)
   
   # associate with node IDs
   x.filtered$nodeID1=labels.filtered[x.filtered$file1,"nodeID"]
   x.filtered$nodeID2=labels.filtered[x.filtered$file2,"nodeID"]
   
   
   # the nodes without edges will be represented as singletons.
   uniq.indiv= unique(labels.filtered$individual)
   nNodeColors = length(uniq.indiv)
   # 21 node colors for luad usecase
   source("~/Rscripts/Color.11.R")
   nodecolors.unq = distinct.colors(nNodeColors)
   names(nodecolors.unq)=uniq.indiv
   nodecolors = nodecolors.unq[labels.filtered$individual]
   labelcolors = get.label.colors(nodecolors)
   nodecolors.rgb=col2rgb(nodecolors)
   labelcolors.rgb=col2rgb(labelcolors)
   
   nodes=list()
   nodes$NODE_ID=1:nNodes
   nodes$NODE_LABEL=labels.filtered[,3]
   nodes$NODE_BORDERCOLOR_R=rep(0,nNodes)
   nodes$NODE_BORDERCOLOR_G=rep(0,nNodes)
   nodes$NODE_BORDERCOLOR_B=rep(0,nNodes)
   nodes$NODE_BORDERCOLOR_HEX=rep("#000000",nNodes)
   nodes$NODE_COLOR_R=nodecolors.rgb[1,]
   nodes$NODE_COLOR_G=nodecolors.rgb[2,]
   nodes$NODE_COLOR_B=nodecolors.rgb[3,]
   nodes$NODE_COLOR_HEX=nodecolors
   nodes$NODE_LABELCOLOR_R=labelcolors.rgb[1,]
   nodes$NODE_LABELCOLOR_G=labelcolors.rgb[1,]
   nodes$NODE_LABELCOLOR_B=labelcolors.rgb[1,]
   nodes$NODE_BORDER_WIDTH=rep(1,nNodes)
   nodes$NODE_X=runif(nNodes,-10,10)
   nodes$NODE_Y=runif(nNodes,-10,10)
   
   edges=list()
   edges$SOURCE_ID=x.filtered$nodeID1
   edges$TARGET_ID=x.filtered$nodeID2
   edges$EDGE_LABEL=x.filtered$r ## r
   edges$EDGE_WIDTH=(x.filtered$r-min(x.filtered$r))*4.5/(1-min(x.filtered$r))+0.5  ## e.g. (r-0.3)*5/0.7 (width ranges 0.5~5 for r=0.3-1 (if minimum r in the data is 0.3)
   
   outfile=output.xgmml
   
   generate.xgmml(nodes,edges,outfile)
   
   return(list(nodes=nodes,edges=edges))
   
}

## example for luad usecase
#label.file = "../usecase/LUAD_CGHUB_WGS_ALL.sample_label.txt"
#ngscheckmateoutput.file = "../usecase/LUAD_CGHUB_WGS_ALL_non_redandancy.txt"
#output.xgmml = "luad.out.xgmml"
#create.xgmml.from.ngscheckmateout(label.file,ngscheckmateoutput.file,output.xgmml)



