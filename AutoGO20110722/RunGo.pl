#!/bin/csh
#to run this script execute it in the following way
#ls *.list | xargs -I xx ./RunGo.pl xx 
perl lookupGOinfoBPstrict.pl $1 Glyma2GODB_TAIR10 > $1.BP
perl lookupGOinfoMFstrict.pl $1 Glyma2GODB_TAIR10> $1.MF

perl countGO.pl $1.BP > $1.BP.count
perl countGO.pl $1.MF > $1.MF.count

perl compareGOcount.pl AllGlymaIDs.BPcount $1.BP.count > $1.BP.compare
perl compareGOcount.pl AllGlymaIDs.MFcount $1.MF.count > $1.MF.compare

perl createfisher.pl $1.BP.compare > $1.BP.fisher
perl createfisher.pl $1.MF.compare > $1.MF.fisher
#regenerate the files with the new list
awk '{print $2,$3,$4,$5}' $1.BP.fisher > $1.BP.out #IronBP.out
awk '{print $1}' $1.BP.fisher > $1.BP.names
awk '{print $2,$3,$4,$5}' $1.MF.fisher > $1.MF.out #IronMF.out
awk '{print $1}' $1.MF.fisher > $1.MF.names

R CMD BATCH '--args filein="'$1'"' fishertest.R
perl combinefilesbyGO.pl $1.BP.fisher ATH_GO_GOSLIM.txt $1.BP.output > $1.BP.final
perl combinefilesbyGO.pl $1.MF.fisher ATH_GO_GOSLIM.txt $1.MF.output > $1.MF.final

mkdir $1.folder
mv $1* $1.folder
