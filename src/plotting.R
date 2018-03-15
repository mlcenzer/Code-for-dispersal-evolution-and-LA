## *************** nice little pdf function ****************
pdf.f <- function(f, file, ...) {
  cat(sprintf("Writing %s\n", file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}


## plot a landscape
plot.landscape <- function(prms.ls, label) {

  ## generalize this code so that it can stack two different shading
  ## hues
  ## I just cannot make it look good.
prms.ls<-prms.ls$landscape[,,2]
K<-10
frac<-0.3
#total number of individuals the landscape supports
tot<-K*frac*128*128

#normalize to tot
prms.ls<-prms.ls/(sum(prms.ls)/tot)
l<-8

#make breaks consistent between landscapes
breaks<-seq(0,max(prms.ls),length=(l+1))

#varying saturation
sat<-seq(1,0.1,length=l)

#blue palette varying in saturation
colsB<-sapply(1:l, FUN=function(x)
	hsv(0.66,sat[x],1))

#for a passable red
colsR<-sapply(1:l, FUN=function(x)
	hsv(0.03,sat[x],1))


  image(prms.ls, col=colsB,
        xlab='', ylab='', xaxt='n', yaxt='n', breaks=breaks)
        
      text(x=0.039,y=0.949, label, cex=1.75)
}

