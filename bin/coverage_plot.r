#; USAGE: coverage_plot.r [input dir] [output file] [target size] [sample1 name [...]]

options(warn=-1)
stdin <- commandArgs(TRUE) 

if(length(stdin) < 4){
	stop("ERROR! Incorrect number of arguments. \nUSAGE: coverage_plot.r [input dir] [output file] [target size] [sample1 name [...]]")
}

input_dir <- stdin[1]
output_file <- stdin[2]
target_size <- as.double(stdin[3])
sample_size <- length(stdin)-3
samplenames <- stdin[4:length(stdin)]
extension <- ".coverage.grp"

quantile_coverage <- as.double(c())
cov_col <- 0
# calculate 95 percentile depth for each sample
for(i in 1:sample_size){
	coveragefile <- paste(input_dir,"/",samplenames[i],extension, sep="", collapse=NULL)
	samplematrix <- as.matrix(read.table(file=coveragefile,header=T))
	if(cov_col == 0){
		cov_col <- ncol(samplematrix)
	}
	percent <- as.double(c())
	for(j in 1:nrow(samplematrix)){
		percent[j] <- as.double(((as.double(sum(as.numeric(samplematrix[samplematrix[,1]>=j,cov_col]))))/target_size)*100)
	}
	quantile_coverage <- cbind(quantile_coverage,length(percent[percent>=5]))
}
depth_cutoff <- round(mean(quantile_coverage))

# create plot
png(file=output_file, units="in", width=11, height=8.5, res=300)
if(sample_size > 1){
	limSize = ifelse(sample_size<=1024, sample_size, 1024)
	palette(rainbow(limSize))
}

target_size_units <- paste(target_size," bp ",sep="")
if(target_size >= 1000){
	target_size_units <- sprintf("%.2f kbp ",target_size/1000)
}
if(target_size >= 1000000){
	target_size_units <- sprintf("%.2f Mbp ",target_size/1000000)
}
if(target_size >= 1000000000){
	target_size_units <- sprintf("%.2f Gbp ",target_size/1000000000)
}
for(i in 1:sample_size){
	coveragefile <- paste(input_dir,"/",samplenames[i],extension, sep="", collapse=NULL)
	samplematrix <- as.matrix(read.table(file=coveragefile,header=T))
	percent <- as.double(c())
	for(j in 1:depth_cutoff){
		percent[j] <- as.double(((as.double(sum(as.numeric(samplematrix[samplematrix[,1]>=j,cov_col]))))/target_size)*100)
	}
	linecolor <- i
	sym <- i
	if(i == 1){
		plot(percent,xlim=c(1,depth_cutoff),ylim=c(0,100),main=paste("Coverage Across ",target_size_units,"Target Region",sep=""),xlab="Depth of coverage",ylab="Percent coverage", pch=sym,col=linecolor,type="o",cex=0.6,lwd=1)
	}else{
		lines(percent,col=linecolor,pch=sym,type="o",cex=0.6,lwd=1)
	}
}

# add legend
linecolors <- seq(1:sample_size)
symbol <- seq(1:sample_size)
legend(x="topright",samplenames,inset=0.02,col=linecolors,pch=symbol,lty=1,lwd=1,cex=0.6)

dev.off()


