## TCR analysis- Aleksandar Obradovic

## read tables into R, then create cd4 and cd8 objects with columns (unstim,stim,Bx), and run

# data = read.table("fileName.tsv", header=T, sep="\t")       #can specify colClasses=c("factor", rep("numeric", numcolumns-1)) 
#colnames(data)
# cd4= data[,c(1,2,3)]
# cd8= data[,c(4,5,6)]
# run(cd4,cd8)

run <- function(cd4, cd8, fold=5, freq=0.00001, clonality=FALSE, ambiguityRatio=5) {
  cd4 = normalize(cd4)
  cd8 = normalize(cd8)
  
  rows1 = cleanup(cd4[,1], cd8[,1], ratio = ambiguityRatio)
  rows2 = cleanup(cd4[,2], cd8[,2], ratio = ambiguityRatio)
  
  cd4 = exclude(cd4, rows1[,1] | rows2[,1])
  cd8 = exclude(cd8, rows1[,2] | rows2[,2])
  
  cat(paste(colnames(cd4)[3], ": cd4\n"))
  reg(cd4, fold=fold, clonality=clonality, freq1=freq)
  cat(paste(colnames(cd8)[3], ": cd8\n"))
  reg(cd8, fold=fold, clonality=clonality, freq1=freq)  
}


reactiveClones <- function(data, fold=5, freq=1e-5) {
  rclones = data[data[,2] > freq & data[,2] > data[,1] * fold, ]
  return(rclones)
}


reg <- function(data, freq1 = 0.00001, freq2 = 0, fold = 5, clonality=FALSE, expansion=FALSE,turnover=FALSE) {
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
  write.table(paste("Pre", npre, nhvgpre, round(freqpre, 5), round(hfreqpre, 5), round(ratio1, 5),  sep="\t"),  quote=F, col.names = F, row.names=F)
  write.table(paste("Bx",  nHBx, nhvgpost, round(freqBx, 5) , round(hfreqpost,5), round(ratio2, 5),  sep="\t"),  quote=F, col.names = F, row.names=F)
  
  x = c(nhvgpost, nHBx-nhvgpost)
  y = c(nhvgpre, npre-nhvgpre)
  ft = fisher.test(cbind(x,y))
  write.table(paste("Fisher'sExactTest:\tp-value=", signif(ft$p.value,3), "\tOR=",  signif(ft$estimate, 3),  " (95%CI:",  signif(ft$conf.int[1], 3), "-",signif(ft$conf.int[2], 3), ")", sep="" ),  sep="\t" , quote=F, col.names = F, row.names=F)
  
  # testing expansion and deletion: are donor-reactive clone more likely to be detectable in postTx than preTx
  # Dpre;N-Dpre
  # Dpost;N-Dpost
  if(expansion){
  write.table(paste("Pre", nhvgpre, dim(allHvG)[1]-nhvgpre, sep="\t"),  quote=F, col.names = F, row.names=F)
  write.table(paste("Post",  nhvgpost, dim(allHvG)[1]-nhvgpost,  sep="\t"),  quote=F, col.names = F, row.names=F)
  x = c(nhvgpre, dim(allHvG)[1]-nhvgpre)
  y = c(nhvgpost, dim(allHvG)[1]-nhvgpost)
  ft = fisher.test(cbind(y,x))
  write.table(paste("Fisher'sExactTest:\tp-value=", signif(ft$p.value,3), "\tOR=",  signif(ft$estimate, 3),  " (95%CI:",  signif(ft$conf.int[1], 3), "-",signif(ft$conf.int[2], 3), ")", sep="" ),  sep="\t" , quote=F, col.names = F, row.names=F)
  }
  # repertoire turnover analysis: are preTx donor-reactive clones more likely to be persist than total preTx clones
  # N1,N0  (detected preTx clones, undetected preTx clone)
  # D1,D0 (detected donor-reactive clone, undetecgted donor-reactive clones)
  if(turnover){
  npre=dim(data[data[,1]>0,])[1]
  ndetected=dim(data[data[,1]>0 & data[,3]>0, ])[1]
  nhvgpre=dim(allHvG[allHvG[,1]>0,])[1]
  nhvgdetected=dim(allHvG[allHvG[,1]>0 & allHvG[,3]>0,])[1]
  write.table(paste("Total", ndetected, npre-ndetected,  sep="\t"),  quote=F, col.names = F, row.names=F)
  write.table(paste("HvG",  nhvgdetected, nhvgpre-nhvgdetected,  sep="\t"),  quote=F, col.names = F, row.names=F)
  x = c(ndetected, npre-ndetected)
  y = c(nhvgdetected,nhvgpre-nhvgdetected) 
  ft = fisher.test(cbind(y,x))
  write.table(paste("Fisher'sExactTest:\tp-value=", signif(ft$p.value,3), "\tOR=",  signif(ft$estimate, 3),  " (95%CI:",  signif(ft$conf.int[1], 3), "-",signif(ft$conf.int[2], 3), ")", sep="" ),  sep="\t" , quote=F, col.names = F, row.names=F)
  }
  
  if(clonality==TRUE){
    write.table(paste("clonality\tPre:", cloneCal(allPre[,1]), "\tHvGinPre:", cloneCal(allHvG[allHvG[,1]>0,1]), "\tBx:", cloneCal(allHinBx[,3]), "\tHvGinBx:", cloneCal(HvGinBx[,3])), sep="\t", quote=F, col.names = F, row.names=F)
    write.table(paste("entropy\tPre:", entropyCal(allPre[,1]), "\tHvGinPre:", entropyCal(allHvG[allHvG[,1]>0,1]), "\tBx:", entropyCal(allHinBx[,3]), "\tHvGinBx:", entropyCal(HvGinBx[,3])), sep="\t", quote=F, col.names = F, row.names=F)
    write.table(paste("simpson's index\tPre:", simpsonCal(allPre[,1]), "\tHvGinPre:", simpsonCal(allHvG[allHvG[,1]>0,1]), "\tBx:", simpsonCal(allHinBx[,3]), "\tHvGinBx:", simpsonCal(HvGinBx[,3])), sep="\t", quote=F, col.names = F, row.names=F)
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




#ckbmt_p4<-merge(ckbmt_p4_part1,ckbmt_p4_part2, by="nucleotide", all=T)
#ckmbt_p4[is.na(ckbmt_p4)]<-0