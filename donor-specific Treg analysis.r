## TCR analysis for donor-specific Tregs

## read table into R
## then use functions below as desired

## to get list of total donor-reactive sequences from CFSE-MLR, use listAlloreactive()
## to get list of Treg or non-Treg sequences after sorting error, use listTreg()
## to adjust CD4/CD8 samples for sorting error, use cleanup()
## for clonality and R20, use cloneCal() and getR20()
## for JSD calculations, use jsdThresholded()

reactiveClones <- function(data, fold=5, freq=1e-4) {
  rclones = data[data[,2] > freq & data[,2] > data[,1] * fold, ]
  return(rclones)
}


normalize <- function(data) {
  nc = ncol(data)
  for (i in 1:nc) {
    data[,i] = data[,i] / sum(data[,i])
  }
  return(data)
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


listAlloreactive<-function(cd4,cd8, fold=5, freq1=.0001, ambiguityRatio=5)
{
  # CD4 and CD8 objects with column 1 unstim and column 2 CFSE-low
  # rows must be indexed by clone ID 
  cd4 = normalize(cd4)
  cd8 = normalize(cd8)
  
  rows1 = cleanup(cd4[,1], cd8[,1], ratio = ambiguityRatio)
  rows2 = cleanup(cd4[,2], cd8[,2], ratio = ambiguityRatio)
  
  cd4 = exclude(cd4, rows1[,1] | rows2[,1])
  cd8 = exclude(cd8, rows1[,2] | rows2[,2])
  
  allHvGcd4 = reactiveClones(normalize(cd4), fold=fold, freq= freq1)
  allHvGcd8 = reactiveClones(normalize(cd8), fold=fold, freq=freq1)
  
  return(list(rownames(allHvGcd4), rownames(allHvGcd8)))
}


listTreg<-function(samples, ambiguityRatio=2)
{
  # make object "samples" with unstimulated Tregs, non-Tregs and unstimulated CD8s
  # for Tregs sequences: column 1 unstimulated Tregs, column 2 non-Tregs, column 3 unstimulated CD8s
  # for non-Tregs: column 1 unstimulated non-Tregs, column 2 Tregs, column 3 unstimulated CD8s
  # rows must be indexed by clone ID 
  
  samples=normalize(samples)
  
  rows1 = cleanup(samples[,1], samples[,2], ratio = ambiguityRatio)
  rows2 = cleanup(samples[,1], samples[,3], ratio = ambiguityRatio)
  samples = exclude(samples, rows1[,1] | rows2[,1])
  samples = samples[samples[,1]>0,]
  
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
  return(num/length(column[column>0]))
}


cloneCal <- function(array) {
  x = array[array >0 ] / sum(array)
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


jensen_shannon <- function(p, q){
  p = p / sum(p)
  q = q / sum(q)
  Hj = shannon.entropy(0.5 *(p+q)) 
  Hp = shannon.entropy(p) 
  Hq = shannon.entropy(q)
  
  jsd = Hj - 0.5*(Hp+Hq)
  return(jsd)
}

jsdThresholded<-function(data,topN=100)
{
  write.table(jsdReport(data,topN), sep="\t")
}