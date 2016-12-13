source('fish4mich.R')
source('import4mich.R')
print(filein)
fileout<-paste(filein,".BP.out",sep="")
filename<-paste(filein,".BP.names",sep="")
print(fileout)
print(filename)
importfisherdata(fileout)->fisherdatatable
out<-fish(fisherdatatable)
read.table(filename)->IronBP.names
as.data.frame(IronBP.names)->IronBP.names.df
as.data.frame(out)->out.p
df1 = data.frame(IronBP.names.df$V1, out.p[,1])
write.table(df1,file=paste("A",".BP.output",sep=""),sep="\t",quote=F)

fileoutMF<-paste(filein,".MF.out",sep="")
filenameMF<-paste(filein,".MF.names",sep="")
importfisherdata(fileoutMF)->fisherdatatableMF
outMF<-fish(fisherdatatableMF)
read.table(filenameMF)->IronMF.names
as.data.frame(IronMF.names)->IronMF.names.df
as.data.frame(outMF)->out.pMF
df1MF = data.frame(IronMF.names.df$V1, out.pMF[,1])
write.table(df1MF,file=paste(filein,".MF.output",sep=""),sep="\t",quote=F)

