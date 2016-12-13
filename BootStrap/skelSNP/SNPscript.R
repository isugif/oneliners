#This script is used to cluster SNPs onto Soybean chromosomes.  Requires SNPsource.R and .RDataSNP
#Created on 04/06/10 by Andrew Severin andrewseverin@gmail.com 
#Iowa State University

#required libraries
	library(gplots)
	source('SNPsource.R')
	load(".RDataSNP")
	
#starting parameters 
	dir<-"./"  
	numofsims<-3
	SNPsofInterest<-read.table('./exampleSNPsFile.txt') 		#list with SNPs of interest
	numberOfChromosome<-20  									#this variable will allow you to loop through the first X chromosomes (see for loop below)

#this section is optional if you would like to have multiple bin sizes uncomment	
	#StartingBinsize<-6000000									#important the the vector in this forloop results in binsizes that include the binsizes before it
	#binscales<-c(1,2,6,12,60,120)								#for binsize 6M 3M 1M 500K 100K 50k
#For multiple bin sizes comment out this block of code
	StartingBinsize<-500000									#Here I chose just one binsize
	binscales<-c(1)	
	
	
#variables calculated from the input parameters
	geneCoordinatesAve<-matrix(round((geneCoordinates[,3]+geneCoordinates[,4])/2),ncol=1)
	rownames(geneCoordinatesAve)<-rownames(geneCoordinates)
	chromosomelengthAll<-chrom[,4]
	maxchromosomelength<-max(chrom[,4])
	significantIntervalsOrig<-0


#this for loop will cycle through each chromosomes.
	for (i in 1:numberOfChromosome){
	#for (i in numberOfChromosome:numberOfChromosome){			#This line can be uncommented if you want to run it on a specific chromosome
		dir.create(paste("./",i,sep=""))						#create directory to export outfiles
		chromosomelength<-chromosomelengthAll[i]
	
	maxNumGenesInCluster<-0										#initiate a variable that will be needed later for plotting
	SNPs<-identifySNPSOnChromosomeForSoybean(SNPsofInterest,i)	#this function will identify the SNPs on each chromosome as it goes through the loop
	
	if (dim(SNPs)[1]<3){										#No need to look at chromosomes that do not have at least 3 SNPs
	print(SNPs)
	next()
	}

	bootData<-simulateData(numofsims,AllGeneCalls,geneCoordinates,SNPs,i)			#generate the simulated data  (See SNPsource for code)



#binSize (for loop) will cycle through the binsizes determined above
	for (b in binscales){
		
		binsize<-StartingBinsize/b
		print(binsize)
		appendtofilename<-paste("_",binsize/binscales,sep="")						#this variable is used for the outputfiles to distinguish between bins
		SNPcoords<-SNPs[,1]															#for retrieval of the coordinates of the SNPs of interest

#function to do bootstrap method
		significantIntervals<-clusterByBoostrap(chromosomelength,binsize,geneCoordinatesAve,SNPcoords,AllGeneCalls,numofsims,bootData)
		print(significantIntervals)
		
		if (b==binscales[1]){
		#open a pdf file
		pdf(file=paste(i,"/chrom",i,"ALL",".pdf",sep=""),paper="special",height=7,width=100)
		plot(0:1, 0:1, type="n", axes=FALSE, ann=FALSE)
		}
#if there are no significant Intervals go to the next binsize in the loop
		
		if(significantIntervals==0){
			next
		}
#This block estimates the required Y dimension for plotting and works for most cases
			if (maxNumGenesInCluster==0){												
				maxNumGenesInCluster<-max(significantIntervals[,3])+1
			}

			
#required input variables for the plotting function.  ChromosomePlot can be found in SNPsource.
			intervalstart<-significantIntervals[,1]
			intervalend<-significantIntervals[,2]
			intervalheight<-significantIntervals[,3]
			currentBinScale<-b
			ChromosomePlot(chromosomelength,intervalstart,intervalend,intervalheight,maxNumGenesInCluster,maxchromosomelength, binScales, currentBinScale,SNPcoords)
			
			if(b==binscales[length(binscales)]){
			#now that the plotting is finished, close the pdf file
			dev.off()
			}
			
			#write to file gene lists with intervals that are significant
			colnames(significantIntervals)<-c('intervalstart','intervalend','numberInInterval','ZscoreaboveBootstrap')
			write.table(significantIntervals,file=paste(i,"/clusterTable",i,appendtofilename,".txt",sep=""),append=T,quote=F,col.names=T)

	}
	
	#this commands save the R sessions for each chromosome into each chromosome folder respectively.
	save.image(file = paste(i,"/.RData",i,sep=""))

}




