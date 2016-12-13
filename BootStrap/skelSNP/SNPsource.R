#This file contains small scripts required to run SNPscript.R  I have found many of these functions are also useful for other analyses.
#Created on 04/06/10 by Andrew Severin andrewseverin@gmail.com 
#Iowa State University


#########################################################
#This function will give the index number in a matrix given the rowname.	
IndexFromGeneCall<-function(inputmatrix,rowname){	
	match(rowname,rownames(inputmatrix),nomatch=0)
}


#########################################################	
#plotting function
#this function will generate intervals of significant clustering of SNPs on soybean chromosomes scaled to the longest chromosome.

ChromosomePlot<-function(chromosomelength,intervalstart,intervalend,intervalheight,maxNumGenesInCluster,maxchromosomelength,binScales,currentBinScale,SNPcoords){
		
		#The axis labels require resizing depending on the width of the plot  
		#the following is an estimate of the resizing required
			if (maxchromosomelength<750000){
				XcexVar<-1
			}
			if (maxchromosomelength>750000 & maxchromosomelength<12000000){
				XcexVar<-(-0.3125)*log(maxchromosomelength/1000000)+0.91
			}
			if (maxchromosomelength>12000000){
				XcexVar<-.1
			}

		
		#I make use of the rect function that has input as (xleft, ybottom, xright, ytop)
		#keep in mind everything is scaled to Gm18, the largest chromosome
		xleft=(1/maxchromosomelength)*intervalstart
		ybottom=0
		width=(intervalend-intervalstart)/maxchromosomelength
		averagepos<-(intervalend+intervalstart)/(2*maxchromosomelength)
		
		average<-(intervalend+intervalstart)/(2)
		xright=xleft+width
		ytop=intervalheight/maxNumGenesInCluster
		rect(0,0,chromosomelength/maxchromosomelength,-0.04)																			#this draws the chromosome on the bottom
		rect(SNPcoords/maxchromosomelength,0,SNPcoords/maxchromosomelength ,-0.04,lwd=.1)												#this draws the location of each SNP
		rect(xleft,ybottom,xright,ytop,border="black",col=rainbow(length(binscales))[which(binscales==currentBinScale)])				#Significant Intervals of SNPs
		text(chromosomelength/(2*maxchromosomelength),-.02,labels=paste("chromosome",i),cex=.5)											#chromosome name
		text(SNPcoords/maxchromosomelength-2000/maxchromosomelength,-.02,labels=rownames(SNPcoords),cex=XcexVar,srt=90)					#location of each Interval
		axis(1,tick=T,at=intervalstart/maxchromosomelength,labels=intervalstart,cex.axis=XcexVar,las = 2,lwd=.5)						#axis 1
		axis(2,tick=T,at=seq(0,1,1/maxNumGenesInCluster),labels=0:maxNumGenesInCluster,cex.axis=.8)										#axis 2
		ablineMulti<-function(i){abline(i,0,col="darkgrey",lwd=.1)}																		#Creates a grid at 2 SNP intervals
		sapply(seq(0,1,2/maxNumGenesInCluster),ablineMulti)																				#sapply to create the grid	
		}

#########################################################	

#this function is not used in the SNPscript but is a handy little function.
#Used in a similar script for clustering genes.
identifyGenesOnChromosomeForSoybean<-function(GeneList){
		#this loop identifies all gene model names (glymas) on a specified chromosome
		if (i<10){
		glymas<-GeneList[grep(paste("0",i,"g",sep=""),GeneList)]
		}else{
		glymas<-GeneList[grep(paste(i,"g",sep=""),GeneList)]	
		}
		return(glymas)
}

identifySNPSOnChromosomeForSoybean<-function(snpList,chromosomeNumber){
			#this loop identifies all SNPs on a specified chromosome
			if (chromosomeNumber<10){
			snps<-snpList[grep(paste("Gm0",chromosomeNumber,sep=""),rownames(snpList)),]
			}else{
			snps<-snpList[grep(paste("Gm",chromosomeNumber,sep=""),rownames(snpList)),]	
			}
			return(snps) #this will return the snp matrix that corresponds only to the chromosome of interest
	}
			
	
#########################################################	

	
#this section will take the same number of genes and simulate how the genes will fall into the bins based on the 1000(or numofsims) random collections of genes
	simulateData<-function(numofsims,AllGeneCalls,geneCoordinates,SNPS,chromosomeNumber){		
			#matrix for storing simulations
				genesAll<-identifyGenesOnChromosomeForSoybean(AllGeneCalls)
				genesSample<-sample(genesAll,dim(SNPS)[1]*numofsims,replace=T)
				AllindexCoords<-function(i){IndexFromGeneCall(geneCoordinates,genesSample[i])}
				sampleIndex<-sapply(1:length(genesSample),AllindexCoords)
				genesSample<-matrix(genesSample,ncol=dim(SNPS)[1])
				sampleIndex<-matrix(sampleIndex,ncol=dim(SNPS)[1])
				return(list(genesSample=genesSample,sampleIndex=sampleIndex))
	}

#########################################################	
#Boostrap function for clustering on a chromosome

#generation of the bin sizes across the chromosome
clusterByBoostrap<-function(chromosomelength,binsize,geneIndexValues,SNPCoordinates,AllGeneCalls,numofsims,bootData){
	print(SNPCoordinates)
	numBins<-floor(chromosomelength/binsize)
	
#this section will calcluate how many of of the SNPs we are interested in fall into each bin
	breaks1<-seq(0,chromosomelength,binsize) 
	chromBinsFind<-findInterval(SNPCoordinates,breaks1, rightmost.closed=T)

  	chromBins<-hist(chromBinsFind,breaks=seq(0,length(breaks1),1),plot=F)$counts
	print(chromBins)	
	sampleIndex<-bootData$sampleIndex												#actual data
	
	chromBinsSample<-matrix(0,numofsims,length(breaks1))							#simulated data
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
			over3stdevZscore<-round(((chromBins-(chrombinsAve))[over3stdev])/chrombinsSD[over3stdev],2)
			
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
			significantIntervals<-matrix(significantIntervals[which(significantIntervals[,3]>1),],ncol=4)
			
			print(significantIntervals)
			if(dim(significantIntervals)[1]==0){
			return(0)}else{
			return(significantIntervals)
			}
}










