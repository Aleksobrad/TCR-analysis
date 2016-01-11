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
