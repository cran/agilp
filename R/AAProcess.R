#AAProcess. A function to extract and optionaly normalise expression data for Agilent arrays.

# Copyright 2010 Benny Chain
  #This program is free software: you can redistribute it and/or modify
  #  it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation, either version 3 of the License, or
  # (at your option) any later version.
#
  #  This program is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  #GNU General Public License for more details.

  #  You should have received a copy of the GNU General Public License
  # along with this program.  If not, see <http://www.gnu.org/licenses/>.

AAProcess<-function(input,output,baseline,s=9,norm = TRUE) {
#processing all files in current sdirectory; extracts green and red median signals, 
#replaces duplicates by means, log2 transforms and stores in new outputfiles 
#(rawrFILENAME.txt and rawgFILENAME.txt)). 

#reads the directory of file names into directory
directory<-unlist(dir(input))[]

#read in first data set skipping first "s" rows
inputpath<-c(rep(1,length(directory)))
inputpath[1]<-paste(input,directory[1],"" ,sep = "")
data<-read.table(inputpath[1],skip = s, header = FALSE, quote = "",comment.char = "", sep = "\t", fill = TRUE,  stringsAsFactors=FALSE )

#assign column names to row 1, and then remove it
colnames(data) <- data[1,]
data<-data[-1,]

#extract probenames, and average of green and red for each replicate probe 
probenames <- data[,"ProbeName"]
green<-tapply(as.numeric(data[,"gMedianSignal"]),probenames,mean)
red<-tapply(as.numeric(data[,"rMedianSignal"]),probenames,mean)
probes<-levels(factor(probenames))

rownames(green)<-probes
rownames(red)<-probes

#set up output files preprocesg and preprocessr
preprocessg<-green
preprocessr<- preprocessg

preprocessg<-merge(preprocessg, green, by = "row.names")
preprocessr<-merge(preprocessr, red, by = "row.names")

row.names(preprocessg)<-preprocessg[,1]
row.names(preprocessr)<-preprocessr[,1]

preprocessg<-preprocessg[,-(1:2)]
preprocessr<-preprocessr[,-(1:2)]


#if there are more than one file, loop through all files in directory
if (length(directory) > 1) {
for (i in 2:length(directory)){

inputpath[i]<-paste(input,directory[i],"" ,sep = "")
data<-read.table(inputpath[i],skip = s, header = FALSE, quote = "",comment.char = "", sep = "\t", fill = TRUE,  stringsAsFactors=FALSE )

colnames(data) <- data[1,]
data<-data[-1,]
probenames <- data[,"ProbeName"]
green<-tapply(as.numeric(data[,"gMedianSignal"]),probenames,mean)
red<-tapply(as.numeric(data[,"rMedianSignal"]),probenames,mean)
probes<-levels(factor(probenames))

rownames(green)<-probes
rownames(red)<-probes

preprocessg<-merge(preprocessg, green, by = "row.names")
preprocessr<-merge(preprocessr, red, by = "row.names")

row.names(preprocessg)<-preprocessg[,1]
row.names(preprocessr)<-preprocessr[,1]

preprocessg<-preprocessg[,-(1)]
preprocessr<-preprocessr[,-(1)]
}
colnames(preprocessg) <- unlist(dir(input))
colnames(preprocessr) <- unlist(dir(input))
} else print("Only one file to be processed")

#convert to log2

datag<<-log2(preprocessg)
datar<<-log2(preprocessr)

#check for non-numeric entries
print ("All entries in preprocessg are numeric")
print(is.numeric(as.matrix(preprocessg)))
print ("All entries in preprocessr are numeric")
print(is.numeric(as.matrix(preprocessr)))

#number of rows and columns of data
N<-nrow(datag)
M<-length(directory)


#print raw log data
if (M == 1){
dg<-c(rep(1,length(directory)))
dr<-c(rep(1,length(directory)))

dg[1]<-paste(output,"rawg",directory[1],"" ,sep = "")
dr[1]<-paste(output,"rawr",directory[1],"" ,sep = "")

write.table(datag,dg[1],sep="\t", col.names = directory[1], row.names = rownames(datag))
write.table(datar,dr[1],sep="\t", col.names = directory[1], row.names = rownames(datag))
}else {

for (i in 1:M){
dg<-c(rep(1,length(directory)))
dr<-c(rep(1,length(directory)))

dg[i]<-paste(output,"rawg",directory[i],"" ,sep = "")
dr[i]<-paste(output,"rawr",directory[i],"" ,sep = "")

write.table(datag[,i],dg[i],sep="\t", col.names = directory[i], row.names = rownames(datag))
write.table(datar[,i],dr[i],sep="\t", col.names = directory[i], row.names = rownames(datag))
}
	}
#stop run if norm = FALSE
if (norm == FALSE) {
print ("No normalisation carried out")
}
if (norm == TRUE){
#read in baseline file
rowmeans<-read.table(file = baseline,header=TRUE, row.names = 1, fill = TRUE, sep = "\t")

	#Normalisation by loess against baseline
	

	datag<-merge(rowmeans, datag, by = "row.names")
	datar<-merge(rowmeans, datar, by = "row.names")

	row.names(datag)<-datag[,1]
	row.names(datar)<-datar[,1]

	datag<-datag[,-(1:2)]
	datar<-datar[,-(1:2)]
	
	normg<-datag
	normr<-datar
	
#ensure that rowmeans and the data have the same number of probes

	rowm<-merge(rowmeans, normg, by = "row.names")
	rownew<-rowm[,2]
	dim(rownew)<-c(length(rownew),1)
	row.names(rownew)<-rowm[,1]
	rowmeans<-rownew

		if (M > 1){
		for (i in 1:M) {
		#fits loess model to ith column of data versus means
		loessg<-loess(as.matrix(datag[,i])~ as.matrix(rowmeans), span = 0.1)
		loessr<-loess(as.matrix(datar[,i])~ as.matrix(rowmeans), span = 0.1)

		#calculate the new normalised data; z is difference between predicted and actual

		zg<-rowmeans-predict(loessg,rowmeans)
		zr<-rowmeans-predict(loessr,rowmeans)

		normg[,i]<-datag[,i]+zg[,1]
		normr[,i]<-datar[,i]+zr[,1]


		#write file with normalised data

		dg[i]<-paste(output,"ng",directory[i],"" ,sep = "")
		dr[i]<-paste(output,"nr",directory[i],"" ,sep = "")

		write.table(normg[,i],dg[i],sep="\t", col.names = directory[i], row.names = rownames(rowmeans))
		write.table(normr[,i],dr[i],sep="\t", col.names = directory[i], row.names = rownames(rowmeans))

		}
	normg<<-normg
	normr<<-normr
	} else {
		loessg<-loess(as.matrix(datag)~ as.matrix(rowmeans), span = 0.1)
		loessr<-loess(as.matrix(datar)~ as.matrix(rowmeans), span = 0.1)

		#calculate the new normalised data; z is difference between predicted and actual

		zg<-rowmeans-predict(loessg,rowmeans)
		zr<-rowmeans-predict(loessr,rowmeans)

		normg<-datag + as.vector(zg)
		normr<-datar + as.vector(zr)
		
		normg<<-normg
		normr<<-normr

		#write file with normalised data

		dg[1]<-paste(output,"ng",directory[1],"" ,sep = "")
		dr[1]<-paste(output,"nr",directory[1],"" ,sep = "")

		write.table(normg,dg[1],sep="\t", col.names = directory[1], row.names = rownames(rowmeans))
		write.table(normr,dr[1],sep="\t", col.names = directory[1], row.names = rownames(rowmeans))

		}
	}
rm(data,rowmeans,N,M,zg,zr,s,input,output)
#end of program
}
