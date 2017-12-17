## TCR analysis- Aleksandar Obradovic

## read tables into R, then create cd4 and cd8 objects with columns (unstim,stim,Bx), and run

# data = read.table("fileName.tsv", header=T, sep="\t")       #can specify colClasses=c("factor", rep("numeric", numcolumns-1)) 
# colnames(data)
# cd4= data[,c(1,2,3)] #unstim,stim,sample
# cd8= data[,c(4,5,6)] #unstim,stim,sample
# run(cd4,cd8)


#remove donor/recipient ambiguous clones, either real shared clones, or ambiguous raised by CFSE-MLR sorting error
#Example: Pt15
#setwd('P:/CCTI_USERS/Jianing Fu/Adaptive samples and analysis')
#data=read.table('Pt15 05042017.tsv',header=T)

#setwd('P:/CCTI_USERS/Jianing Fu/Adaptive samples and analysis')
#data2=read.table('Pt15 05042017 table 2.tsv',header=T)

#p15=merge(data,data2,by="nucleotide", all=T)
#p15[is.na(p15)]=0

#write.table(p15,file ="P:/CCTI_USERS/Jianing Fu/Adaptive samples and analysis/Pt15 combined.tsv",quote=F,row.names=F,col.names=F, sep="\t")


# setwd('P:/CCTI_USERS/Jianing Fu/Adaptive samples and analysis')
# data=read.table('Pt15 combined.tsv',header=T)
# source('P:/CCTI_USERS/Jianing Fu/Adaptive samples and analysis/09-12-2017 JSD ileum colon native colon/alloreactivity_printtofile 09132017.R')
# rownames(data)=data[,1]
# x=data
# names(data)
#[1] "nucleotide"                        "total.x"                          
#[3] "Pt15MVTx.GVH.D4L"                  "Pt15MVTx.GVH.D8L"                 
#[5] "Pt15MVTx.GVH.R4L"                  "Pt15MVTx.GVH.R8L"                 
#[7] "Pt15MVTx.PBMCs.POD143.Donor_T"     "Pt15MVTx.PBMCs.POD143.Recipient_T"
#[9] "Pt15MVTx.pre.Tx.donor_MLN"         "Pt15MVTx.unstim.D4"               
#[11] "Pt15MVTx.unstim.D8"                "Pt15MVTx.unstim.R4"               
#[13] "Pt15MVTx.unstim.R8"                "Pt15_MVTx_PBMC_POD83"             
#[15] "Pt15_MVTx_PBMC_POD83_MLR_D3L"      "total.y"                          
#[17] "Pt15MVTx.Bx.POD17"                 "Pt15MVTx.Bx.POD27"                
#[19] "Pt15MVTx.Bx.POD55"                 "Pt15MVTx.PBMCs.POD11"             
#[21] "Pt15MVTx.PBMCs.POD19"              "Pt15MVTx.PBMCs.POD26"             
#[23] "Pt15MVTx_Bx_POD237"                "Pt15MVTx_PBMC_POD255"             


#ambiguous=which((x[,13]>0 | x[,12]>0)&(x[,11]>0 | x[,10]>0))
#ambiguous=union(ambiguous,which((x[,13]>0 | x[,12]>0)&(x[,3]>0 | x[,4]>0))) 
#ambiguous=union(ambiguous,which((x[,11]>0 | x[,10]>0)&(x[,5]>0 | x[,6]>0))) 
#x=x[setdiff(1:nrow(x),ambiguous),] 


# 5 fold default, stim vs unstim, run "fold=2" for alloreactive clones when needed.


run <- function(cd4, cd8, fold=5, freq=0.00001, clonality=FALSE, turnover=FALSE, absoluteexpansion=FALSE, ambiguityRatio=5, filename="Default") {
  write.table(paste("Number of unique clones in sample (from CD4 table):", length(which(cd4[,3]>0)),  sep="\t"),  quote=F, col.names = F, row.names=F, file=filename, append=T)
  write.table(paste("Number of unique clones in sample (from CD8 table):", length(which(cd8[,3]>0)),  sep="\t"),  quote=F, col.names = F, row.names=F, file=filename, append=T)
  
  cd4 = normalize(cd4)
  cd8 = normalize(cd8)
  
  rows1 = cleanup(cd4[,1], cd8[,1], ratio = ambiguityRatio)
  rows2 = cleanup(cd4[,2], cd8[,2], ratio = ambiguityRatio)
  
  cd4 = exclude(cd4, rows1[,1] | rows2[,1])
  cd8 = exclude(cd8, rows1[,2] | rows2[,2])
  
  cat(paste("\n",colnames(cd4)[3], ": cd4\n"), file=filename, append=T)
  reg(cd4, fold=fold, clonality=clonality, turnover=turnover,absoluteexpansion=absoluteexpansion, freq1=freq, filename=filename)
  cat(paste("\n",colnames(cd8)[3], ": cd8\n"), file=filename, append=T)
  reg(cd8, fold=fold, clonality=clonality, turnover=turnover, absoluteexpansion=absoluteexpansion, freq1=freq,filename=filename)  
}


reactiveClones <- function(data, fold=5, freq=1e-5) {
  rclones = data[data[,2] > freq & data[,2] > data[,1] * fold, ]
  return(rclones)
}


reg <- function(data, freq1 = 0.00001, freq2 = 0, fold = 5, clonality=FALSE, turnover=FALSE,  absoluteexpansion=FALSE, filename="Default") {
  ## columns in data:
  # unstim, stim, Bx
  
  data = normalize(data)
  
  allPre = data[data[,1] > freq2,]
  allHinBx = data[(data[,1] > freq2 | data[,2] > freq2 ) & data[,3] > freq2 , ]
  
  # all HvG clones
  allHvG = reactiveClones(data, fold=fold, freq= freq1)
  HvGinBx = reactiveClones(allHinBx, fold = fold, freq = freq1)
  
  npre = dim(allPre)[1]
  nHBx = dim(allHinBx)[1]
  
  nhvgpre = dim(allHvG[allHvG[,1]>0, ])[1]
  nhvgpost = dim(HvGinBx)[1]
  
  
  freqpre = sum(allPre[,1])
  freqBx = sum(allHinBx[,3])
  
  hfreqpre = sum(allHvG[,1])
  hfreqpost = sum(allHvG[,3])
  
  ratio1 = hfreqpre / freqpre
  ratio2 = hfreqpost / freqBx
  # Relative expansion analysis (default): Is a mappable clone (one that can be called alloreactive or nonalloreactive based on preTx MLR), more likely to be alloreactive in postTx than preTx?
  #(i.e.) is the alloreactivity rate higher in postTx than preTx? 
  # allopre, Npre-allopre
  # allopost, Npost-allopost
  write.table("Relative Expansion", quote=F, col.names=F, row.names=F, file=filename, append=T)
  write.table(paste("Samples","preTx_mappable_clones","alloreactive_clones","freq_mapped","freq_alloreactive","allo_cum_freq","allo_clone_fraction", sep="\t"),quote=F,col.names=F,row.names=F, file=filename, append=T)
  write.table(paste("preTx", npre, nhvgpre, round(freqpre, 5), round(hfreqpre, 5), round(ratio1, 5),round(nhvgpre/npre,5),  sep="\t"),  quote=F, col.names = F, row.names=F, file=filename, append=T)
  write.table(paste("postTx",  nHBx, nhvgpost, round(freqBx, 5) , round(hfreqpost,5), round(ratio2, 5),round(nhvgpost/nHBx,5),  sep="\t"),  quote=F, col.names = F, row.names=F, file=filename, append=T)
  
  x = c(nhvgpost, nHBx-nhvgpost)
  y = c(nhvgpre, npre-nhvgpre)
  ft = fisher.test(cbind(x,y))
  write.table(paste("Fisher'sExactTest:\tp-value=", signif(ft$p.value,3), "\tOR=",  signif(ft$estimate, 3),  " (95%CI:",  signif(ft$conf.int[1], 3), "-",signif(ft$conf.int[2], 3), ")", sep="" ),  sep="\t" , quote=F, col.names = F, row.names=F, file=filename, append=T)
  
  # repertoire turnover analysis: are preTx donor-reactive clones more likely to be persist than total preTx clones
  # N1,N0  (detected preTx clones, undetected preTx clones)
  # D1,D0 (donor-reactive clones, undetected donor-reactive clones)
  if(turnover){
    npre=dim(data[data[,1]>0,])[1]
    ndetected=dim(data[data[,1]>0 & data[,3]>0, ])[1]
    nhvgpre=dim(allHvG[allHvG[,1]>0,])[1]
    nhvgdetected=dim(allHvG[allHvG[,1]>0 & allHvG[,3]>0,])[1]
    write.table("Turnover", quote=F, col.names=F, row.names=F, file=filename, append=T)
    write.table(paste("clone sets","detected_in_sample","undetected_in_sample","persistent_clone_fraction", sep="\t"),quote=F,col.names=F,row.names=F, file=filename, append=T)
    write.table(paste("all preTx unstim clones", ndetected, npre-ndetected, ndetected/npre, sep="\t"),  quote=F, col.names = F, row.names=F, file=filename, append=T)
    write.table(paste("HvG clones",  nhvgdetected, nhvgpre-nhvgdetected, nhvgdetected/nhvgpre, sep="\t"),  quote=F, col.names = F, row.names=F, file=filename, append=T)
    x = c(ndetected, npre-ndetected)
    y = c(nhvgdetected,nhvgpre-nhvgdetected) 
    ft = fisher.test(cbind(y,x))
    write.table(paste("Fisher'sExactTest:\tp-value=", signif(ft$p.value,3), "\tOR=",  signif(ft$estimate, 3),  " (95%CI:",  signif(ft$conf.int[1], 3), "-",signif(ft$conf.int[2], 3), ")", sep="" ),  sep="\t" , quote=F, col.names = F, row.names=F, file=filename, append=T)
  }
  
  #absolute expansion analysis: are alloreative clones more likely to be detected in pre or post-transplant sample?
  # pre, N-pre
  # post, N-post   where N is the total number of clones termed alloreactive
  if(absoluteexpansion){
    nhvg=nrow(allHvG)
    write.table("Absolute Expansion", quote=F, col.names=F, row.names=F, file=filename, append=T)
    write.table(paste("clone sets","alloreactive_clones_detected","alloreactive_clones_undetected","detection_rate", sep="\t"),quote=F,col.names=F,row.names=F, file=filename, append=T)
    write.table(paste("preTx", nhvgpre, nhvg-nhvgpre, nhvgpre/nhvg, sep="\t"),  quote=F, col.names = F, row.names=F, file=filename, append=T)
    write.table(paste("postTx",  nhvgpost, nhvg-nhvgpost, nhvgpost/nhvg, sep="\t"),  quote=F, col.names = F, row.names=F, file=filename, append=T)
    x = c(nhvgpre, nhvg-nhvgpre)
    y = c(nhvgpost, nhvg-nhvgpost) 
    ft = fisher.test(cbind(y,x))
    write.table(paste("Fisher'sExactTest:\tp-value=", signif(ft$p.value,3), "\tOR=",  signif(ft$estimate, 3),  " (95%CI:",  signif(ft$conf.int[1], 3), "-",signif(ft$conf.int[2], 3), ")", sep="" ),  sep="\t" , quote=F, col.names = F, row.names=F, file=filename, append=T)
  }
  
  if(clonality==TRUE){
    write.table(paste("clonality\tPre:", cloneCal(allPre[,1]), "\tHvGinPre:", cloneCal(allHvG[allHvG[,1]>0,1]), "\tBx:", cloneCal(allHinBx[,3]), "\tHvGinBx:", cloneCal(HvGinBx[,3])), sep="\t", quote=F, col.names = F, row.names=F, file=filename, append=T)
    write.table(paste("entropy\tPre:", entropyCal(allPre[,1]), "\tHvGinPre:", entropyCal(allHvG[allHvG[,1]>0,1]), "\tBx:", entropyCal(allHinBx[,3]), "\tHvGinBx:", entropyCal(HvGinBx[,3])), sep="\t", quote=F, col.names = F, row.names=F, file=filename, append=T)
    write.table(paste("simpson's index\tPre:", simpsonCal(allPre[,1]), "\tHvGinPre:", simpsonCal(allHvG[allHvG[,1]>0,1]), "\tBx:", simpsonCal(allHinBx[,3]), "\tHvGinBx:", simpsonCal(HvGinBx[,3])), sep="\t", quote=F, col.names = F, row.names=F, file=filename, append=T)
  }
}


normalize <- function(data) {
  nc = ncol(data)
  for (i in 1:nc) {
    data[,i] = data[,i ] / sum(data[,i])
  }
  return(data)
}


cloneCal <- function(array) {
  x = array[array >0 ] / sum(array)
  #  x = sort(array, decreasing=T)
  l = length(x)
  entropy = sum(x * -1 * log2(x))
  maxentropy = -log2(1/l)
  return(signif(1 - entropy / maxentropy, 2))
}


entropyCal<-function(array) {
  x = array[array >0 ] / sum(array)
  entropy = sum(x * -1 * log2(x))
  return(signif(entropy,2))
}


simpsonCal<-function(array) {
  x = array[array >0 ] / sum(array)
  simpson = sum(x * x)
  return(signif(simpson,2))
}


cleanup <- function(cd4, cd8, ratio = 2) {
  ## remove contaminated clones  
  ambi = (cd4 > 0 & cd8 > 0 & cd4 / cd8 > 1/ratio & cd4 / cd8 < ratio) 
  cd4exclude = (cd4 > 0 & cd8 > 0 & cd4 / cd8 <= 1/ratio ) 
  cd8exclude = (cd4 > 0 & cd8 > 0 & cd4 / cd8 >= ratio )  
  return(cbind(ambi | cd4exclude, ambi | cd8exclude))
}


exclude <- function(data, excluderows) {
  data = data[excluderows == F, ]
  return(data)
}


#before running, must do rownames(data)=data[,1]
#allo[[1]] is the cd4 alloreactives
#allo[[2]] is the cd8 alloreactives
listAlloreactive<-function(cd4,cd8, fold=5, freq1=.00001, ambiguityRatio=5) # rows must be indexed by clone ID
{
  # need only column 1 unstim and column 2 stim
  cd4 = normalize(cd4)
  cd8 = normalize(cd8)
  
  rows1 = cleanup(cd4[,1], cd8[,1], ratio = ambiguityRatio)
  rows2 = cleanup(cd4[,2], cd8[,2], ratio = ambiguityRatio)
  
  cd4 = exclude(cd4, rows1[,1] | rows2[,1])
  cd8 = exclude(cd8, rows1[,2] | rows2[,2])
  
  allHvGcd4 = reactiveClones(normalize(cd4), fold=fold, freq=freq1)
  allHvGcd8 = reactiveClones(normalize(cd8), fold=fold, freq=freq1)
  
  return(list(rownames(allHvGcd4), rownames(allHvGcd8)))
}








#ckbmt_p4<-merge(ckbmt_p4_part1,ckbmt_p4_part2, by="nucleotide", all=T)
#ckmbt_p4[is.na(ckbmt_p4)]<-0


abundancePlot<-function(data, threshold=0, name="abundance plot"){
  par(mfrow=c(1,1)) #prepares space for plots, fill by row
  plot(0,type='n', main=name, xlab="log freq", ylab="log abundance", xlim=c(log10(min(data[data>0])),log10(max(data))), ylim=c(0,5))
  color=1
  if(is.null(ncol(data))){
    out=as.data.frame(table(data[data>threshold]))
    points(log10(as.numeric(as.character(out[,1]))), log10(out[,2]),pch=16, col=color)
  }
  else{
    for(index in 1:length(data))
    {
      out=as.data.frame(table(data[data[,index]>threshold,index]))
      points(log10(as.numeric(as.character(out[,1]))), log10(out[,2]),pch=16, col=color)
      color=color+1
    }
    legend("topright", legend=colnames(data), pch=16,col=1:length(data))
  }
}


#rownames(data)=data[,1]
#allo=listAlloreactive(data[,c(cd4unstim,cd4stim)],data[,c(cd8unstim,cd8stim)])
#allo=union(allo[[1]],allo[[2]])
#topclones(data[,sample],rownames(data),allo)
topclones<-function(sample, rownames,listAlloreactives, n=100)
{
  a<-order(sample,decreasing=TRUE)[1:n]
  color<-rep("blue", n)
  color[which(rownames[a] %in% listAlloreactives)]<-"red"
  #barplot(sample[a,3]/sum(sample[,3]), col=color)
  barplot(as.matrix((sample[a]/sum(sample))), col=color, beside=F)
  #barplot(p4[(order(p4[,3],decreasing=T)[1:100]),3]/sum(p4[,3]))
  #return(rownames[a])
}

resolveambiguous<-function(data,c1,c2,ratio=2){
  normalize(data[,c(c1,c2)])
  c1indices=which(data[,c1]>0 & data[,c2]>0 & data[,c1]>ratio*data[,c2])
  c2indices=which(data[,c1]>0 & data[,c2]>0 & data[,c2]>ratio*data[,c1])
  ambiindices=which(data[,c1]>0 & data[,c2]>0)
  ambiindices=setdiff(setdiff(ambiindices,c1indices),c2indices)
  data[c1indices,c2]=0
  data[c2indices,c1]=0
  data[ambiindices,c(c1,c2)]=0
  return(data)
}


listTreg<-function(samples, fold=5, freq1=.00001, ambiguityRatio=5) # rows must be indexed by clone ID
{
  # column 1 Treg and column 2 cd4nonTreg and column 3 cd8
  sample=normalize(samples)
  
  rows1 = cleanup(samples[,1], samples[,2], ratio = ambiguityRatio)
  rows2 = cleanup(samples[,1], samples[,3], ratio = ambiguityRatio)
  
  samples = exclude(samples, rows1[,1] | rows2[,1])
  samples=samples[samples[,1]>0,]
  
  return(rownames(samples))
}

getR20<-function(column,rval=0.2){
  column=column/sum(column)
  column=column[order(column,decreasing = T)]
  num=0
  total=0
  while(total<rval){
    total=total+column[num+1]
    num=num+1
  }
  return(num/length(column))
}



# Example of defining clone sets 
# a=rownames(data[data[,3]>0,]) # the clones present in column 3
# b=rownames(data[data[,4]>0,]) # the clones present in column 4
# c=rownames(data[data[,5]>0,]) # the clones present in column 5
#
# you can use allo=listAlloreactive(cd4,cd8) to define a list of alloreactive clones, where allo[[1]] is cd4allo and allo[[2]] is cd8 allo
#
# you can use the functions intersect(), union(), and setdiff() to build compound sets. e.g. using the above
# setdiff(a,b) # the clones present in column 3 but not column 4
# intersect(allo[[1]], c) # cd8 alloreactive clones present in column 5
# union(a,union(b,c)) # the clones present in column 3 or column 4 or column 5

venn_diagram<-function(a,b,c){
  print(c("a:", length(setdiff(setdiff(a,b),c))))
  print(c("b:" , length(setdiff(setdiff(b,a),c))))
  print(c("c:" , length(setdiff(setdiff(c,b),a))))
  print(c("a,b:" , length(setdiff(intersect(a,b),c))))
  print(c("a,c:" , length(setdiff(intersect(a,c),b))))
  print(c("b,c:" , length(setdiff(intersect(c,b),a))))
  print(c("a,b,c:" , length(intersect(intersect(a,b),c))))
}

venn_diagram_plot<-function(a,b,c){
  #dev.off()
  require("VennDiagram")
  draw.triple.venn(length(a),length(b),length(c), length(intersect(a,b)),length(intersect(b,c)),length(intersect(a,c)), length(intersect(intersect(a,b),c)),cex=2,fill = c('blue','red','green'), category = c("","",""))
  #draw.pairwise.venn
}






# data = read.table("fileName.tsv", header=T, sep="\t")
# colnames(data)
# jsdThresholded(data[,c(column1index, column2index, column3index)])


jsdReport<-function(data, topN=-1)
{
  out<-matrix(nrow=ncol(data),ncol=ncol(data), dimnames=list(colnames(data), colnames(data)))
  for(i in 1:ncol(data)){
    for(j in 1:ncol(data)){
      if(topN==-1){
        out[i,j]<-jensen_shannon(data[,i], data[,j])
      }
      else{
        a<-order(data[,i],decreasing=TRUE)[1:topN]
        b<-order(data[,j],decreasing=TRUE)[1:topN]
        z<-data[union(a,b),]
        out[i,j]<-jensen_shannon(z[,i],z[,j])
      }
    }
  }
  return(out)
}

shannon.entropy <- function(p)
{
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}

cloneCal <- function(p) {
  x = p[p>0] / sum(p)
  l = length(x)
  entropy = shannon.entropy(p)
  maxentropy = -log2(1/l)
  return(signif(1 - entropy / maxentropy, 2)) }

jensen_shannon <- function(p, q){
  ## JSD = H(0.5 *(p +q)) - 0.5H(p) - 0.5H(q)
  # H(X) = \sum x_i * log2(x_i)
  #  p = p[p >0 & q >0]
  #  q = q[p>0 & q>0]
  p = p / sum(p)
  q = q / sum(q)
  Hj = shannon.entropy(0.5 *(p+q)) 
  Hp = shannon.entropy(p) 
  Hq = shannon.entropy(q)
  
  jsd = Hj - 0.5*(Hp+Hq)
  #	cat(Hj, Hp, Hq, jsd, "\n")
  return(jsd)
}


jsdThresholded<-function(data)
{
  write.table(jsdReport(data),sep="\t")
  cat("\n")
  write.table(jsdReport(data,topN=5000),sep="\t")
  cat("\n")
  write.table(jsdReport(data,topN=1000),sep="\t")
  cat("\n")
  write.table(jsdReport(data,topN=500),sep="\t")
  cat("\n")
  write.table(jsdReport(data,topN=100), sep="\t")
}
