
#This function will give the index number in a matrix given the rowname.	
IndexFromGeneCall<-function(inputmatrix,rowname){	
	match(rowname,rownames(inputmatrix),nomatch=0)
}

#########################################################	

#function for exporting the genes in each bin
#appendtofilename,significantIntervals,glymacoords,glymaIndex,Anorm
exportGeneLists<-function(appendtofilename,significantIntervals,geneCoordinates,geneIndexValues,NormalizedExpressionValues){		
		#print(significantIntervals)
		significantIntervals<-significantIntervals[which(significantIntervals[,3]>1),]
		if (is.vector(significantIntervals)){
		sigbins<-significantIntervals[1:2]
		}else{
			if (is.matrix(significantIntervals)){
				significantIntervals<-significantIntervals[sort(significantIntervals[,1],index.return=T)$ix,]
				sigbins<-significantIntervals[,1:2]
			}else{
			return	
			}
		}
		
		if (is.vector(sigbins)){
			sigbins<-t(as.matrix(sigbins,rows=1,col=1))
		}
		#print(geneIndexValues)
		for (k in 1:dim(sigbins)[1]){
			#print(which(((geneCoordinates[geneIndexValues,3]+geneCoordinates[geneIndexValues,4])/2)>sigbins[j,1] & ((geneCoordinates[geneIndexValues,3]+geneCoordinates[geneIndexValues,4])/2)<sigbins[j,2]))
			#print(dim(geneCoordinates))
			#print(length(geneIndexValues))
			genePositions<-as.matrix(geneCoordinates[geneIndexValues[which(((geneCoordinates[geneIndexValues,3]+geneCoordinates[geneIndexValues,4])/2)>sigbins[k,1] & ((geneCoordinates[geneIndexValues,3]+geneCoordinates[geneIndexValues,4])/2)<sigbins[k,2])],])	
			#print(genePositions)
			
			geneClusters<-rownames(genePositions)
			geneExpressionProfile<-NormalizedExpressionValues[which(as.numeric(rownames(NormalizedExpressionValues) %in% geneClusters)==1),]
			GbrowseAddress<-GbrowseAddressFun(genePositions[1,2],sigbins[k,1],sigbins[k,2])

			#print(soybaseLocation)
			binStart<-sigbins[k,1]
			binEnd<-sigbins[k,2]
			clusterdata<-(cbind(genePositions,binStart,binEnd,GbrowseAddress,geneExpressionProfile))
			
			#write out only one set of column variables
			if (k == 1){
				write.table(clusterdata,file=paste(i,"/clusterTable",i,appendtofilename,".txt",sep=""),append=T,quote=F,col.names=T)
			#print(clusterdata)
			}else{
				write.table(clusterdata,file=paste(i,"/clusterTable",i,appendtofilename,".txt",sep=""),append=T,quote=F,col.names=F)
			#print(clusterdata)
			}
			
		}
}

#########################################################	
#this function provides a text sequence that can be placed into GBrowse to view a particular region
GbrowseAddressFun<-function(chromosome,startposition,endposition){
			GbrowseAddress<-paste(chromosome,startposition,sep=":")
			GbrowseAddress<-paste(GbrowseAddress,endposition,sep="..")
			return(GbrowseAddress)
}
#########################################################	
#plotting function
#this function will generate a histogram like file at the scaled position on the chromosome.

ChromosomePlot<-function(chromosomelength,intervalstart,intervalend,intervalheight,maxNumGenesInCluster,maxchromosomelength,binScales,currentBinScale,geneCoordinatesOfInterest){
		
		#The axis labels require resizing depending on the width of the plot  
		#the following is a crude estimate of the resizing required
			if (maxchromosomelength<750000){
				XcexVar<-1
			}
			if (maxchromosomelength>750000 & maxchromosomelength<12000000){
				XcexVar<-(-0.3125)*log(maxchromosomelength/1000000)+0.91
			}
			if (maxchromosomelength>12000000){
				XcexVar<-.1
			}
		
		a=(1/maxchromosomelength)*intervalstart
		b1=0
		width=(intervalend-intervalstart)/maxchromosomelength
		averagepos<-(intervalend+intervalstart)/(2*maxchromosomelength)
		
		average<-(intervalend+intervalstart)/(2)
		cv=a+width
		d=intervalheight/maxNumGenesInCluster
		rect(0,0,chromosomelength/maxchromosomelength,-0.04)
		rect(geneCoordinatesOfInterest/maxchromosomelength,0,geneCoordinatesOfInterest/maxchromosomelength ,-0.04,lwd=.1)
		rect(a,b1,cv,d,border="black",col=rainbow(length(binscales))[which(binscales==currentBinScale)])
		text(chromosomelength/(2*maxchromosomelength),-.02,labels=paste("chromosome",i),cex=.5)
		text(geneCoordinatesOfInterest/maxchromosomelength-2000/maxchromosomelength,-.02,labels=rownames(geneCoordinatesOfInterest),cex=XcexVar,srt=90)
		axis(1,tick=T,at=intervalstart/maxchromosomelength,labels=intervalstart,cex.axis=XcexVar,las = 2,lwd=.5)
		axis(2,tick=T,at=seq(0,1,1/maxNumGenesInCluster),labels=0:maxNumGenesInCluster,cex.axis=.8)
		ablineMulti<-function(i){abline(i,0,col="darkgrey",lwd=.1)}
		sapply(seq(0,1,2/maxNumGenesInCluster),ablineMulti)
}


identifyGenesOnChromosomeForSoybean<-function(GeneList){
		#this loop identifies all gene model names (glymas) on a specified chromosome
		if (i<10){
		glymas<-GeneList[grep(paste("0",i,"g",sep=""),GeneList)]
		#print(glymas)
		}else{
		glymas<-GeneList[grep(paste(i,"g",sep=""),GeneList)]	
		#print(glymas)
		}
		return(glymas)
}
			
	
	
				#genesAll<-AllGeneCalls[grep(paste(i,"g",sep=""),AllGeneCalls)]
	#**********
	#this is repeated every binscale though doesn't need to be repeated.
	#keep in mind the scope will erase it unless I assign and input it everytime.
	#new function?	
	
#this section will take the same number of genes and simulate how the genes will fall into the bins based on the 1000(or numofsims) random collections of genes
	simulateData<-function(numofsims,AllGeneCalls,geneCoordinates,genes){		
							

			#matrix for storing simulations
			#chromBinsSample<-matrix(0,numofsims,length(breaks1))
			
				genesAll<-identifyGenesOnChromosomeForSoybean(AllGeneCalls)
				genesSample<-sample(genesAll,length(genes)*numofsims,replace=T)
				AllindexCoords<-function(i){IndexFromGeneCall(geneCoordinates,genesSample[i])}
				sampleIndex<-sapply(1:length(genesSample),AllindexCoords)
				genesSample<-matrix(genesSample,ncol=length(genes))
				sampleIndex<-matrix(sampleIndex,ncol=length(genes))
				#print(genesSample)
				#print(sampleIndex)
				
				return(list(genesSample=genesSample,sampleIndex=sampleIndex))
	}
				
	#************		
			
####################
#Boostrap function for clustering on a chromosome

#generation of the bin sizes across the chromosome

clusterByBoostrap<-function(chromosomelength,binsize,geneIndexValues,geneCoordinatesAve,AllGeneCalls,numofsims,bootData){
	
	numBins<-floor(chromosomelength/binsize)
	
#this section will calcluate how many of of the genes we are interested in fall into each bin
	breaks1<-seq(0,chromosomelength,binsize) 
	chromBinsFind<-findInterval(geneCoordinatesAve[geneIndexValues],breaks1, rightmost.closed=T)

  	chromBins<-hist(chromBinsFind,breaks=seq(0,length(breaks1),1),plot=F)$counts
		
	#genesSample<-bootData$genesSample
	sampleIndex<-bootData$sampleIndex
			
			
			chromBinsSample<-matrix(0,numofsims,length(breaks1))
			#genesAll<-identifyGenesOnChromosomeForSoybean(AllGeneCalls)	
				
			for (j in 1:numofsims){	
			#generation of the bin sizes across the chromosome for the simulated data
				breaks1<-seq(0,chromosomelength,binsize) 
				chromBinsSam<-findInterval(geneCoordinatesAve[sampleIndex[j,]],breaks1, rightmost.closed=T)
				chromBinsSample[j,]<-hist(chromBinsSam,breaks=seq(0,length(breaks1),1),plot=F)$counts	
			}
			#Average and standard deviation of the simulated data
			chrombinsAve<-colMeans(chromBinsSample)
			chrombinsSD<-sd(chromBinsSample)
	
			
#this section is the determination of the bins that are significant
			over3stdev<-which((chromBins-(chrombinsAve+3*chrombinsSD))>0)
			over3stdevBy<-(chromBins-(chrombinsAve+3*chrombinsSD))[over3stdev]
			over3stdevZscore<-((chromBins-(chrombinsAve))[over3stdev])/chrombinsSD[over3stdev]
			
			#this if statement is required in case no intervals are found to be significant

			if (length(over3stdev)==0){
				print("No significant intervals found")
				return(0)		
			}else{
			significantIntervals<-matrix(c(breaks1[over3stdev],breaks1[over3stdev]+binsize),length(over3stdev),2)
			}
			
			#append number of genes in the bin and the zscore to significantIntervals
			significantIntervals<-matrix(cbind(significantIntervals,chromBins[over3stdev],over3stdevZscore),ncol=4)
			significantIntervals<-matrix(significantIntervals[sort(significantIntervals[,1],index.return=T)$ix,],ncol=4)

			#this if statement resolves the case where it is the first time through
#			if (length(significantIntervalsOrig)==1){
#				significantIntervalsOrig<-significantIntervals
#			}else{
#				#removing those that are not clusters ie contain only 1 gene
#				significantIntervals<-matrix(significantIntervals[which(significantIntervals[,3]>1),],ncol=4)
#				significantIntervalsOrig<-matrix(significantIntervalsOrig[which(significantIntervalsOrig[,3]>1),],ncol=4)
#				
#				#determine which smaller significant bins overlap with larger significant bins
#				Intervalsig<-Intervals(significantIntervals[,1:2],closed = c( FALSE, FALSE ))
#				IntervalsigOrig<-Intervals(significantIntervalsOrig[,1:2],closed = c( FALSE, FALSE ))
#				IntervalOverlap<-distance_to_nearest(IntervalsigOrig,Intervalsig)
#				print(IntervalOverlap)
#				
#				#combine the list of intervals that didn't overlap with the smaller intervals that are significant
#				significantIntervals<-matrix(rbind(significantIntervalsOrig[which(IntervalOverlap>0),], significantIntervals), ncol=4)
#				
#				#sorting of the intervals so they go from first to last interval
#				significantIntervals<-matrix(significantIntervals[sort(significantIntervals[,1],index.return=T)$ix,],ncol=4)
#				significantIntervalsOrig<-significantIntervals
#			}

			significantIntervals<-matrix(significantIntervals[which(significantIntervals[,3]>1),],ncol=4)

			print(significantIntervals)
			return(significantIntervals)
}










