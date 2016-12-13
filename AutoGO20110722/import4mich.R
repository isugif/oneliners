importfisherdata<-function (filename){
#filename is the name of a file containing selexa data in the following format
#filename='/Volumes/I\ AM\ Jenna/Itpk1.prn'
#           leaf flower cm_pod ...
#Glyma09g33   54     39     46 ...    
#Glyma16g13   67     52    110 ...     
#Glyma17g09   13     17     12 ...   
#...

#scan in header and other data separately
columndata<-scan(filename,list(0,0,0,0))
#headers<-scan(filename,nlines=1,list(""))

#convert data to a usable matrix and assign row and column names
datamatrix<-sapply(columndata[1:length(columndata)],c)
#rownames(datamatrix)<-sapply(columndata[1],c)
#colnames(datamatrix)<-t(sapply(headers,c))
print('datamatrix= ')
return(datamatrix)
}

