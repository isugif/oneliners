fish<-function(datamatrix){
#	dmat<-dim(inputmatrix)
#	pvalues<-matrix(999,dmat[1],sum(1:(dmat[2]-1)))
#	index<-0
dataout<-NULL
	for(i in 1:dim(datamatrix)[1]){
		fi<-fisher.test(matrix(c(datamatrix[i,1],datamatrix[i,2],datamatrix[i,3],datamatrix[i,4]),nrow=2),"two.sided")$p.value
#final<-append(fi)
#print(final)
#return(fi)
dataout[i]<-fi
#	write.table(fi,file="fisher.txt",sep="\t",quote=F)
	}
	return(dataout)
}



