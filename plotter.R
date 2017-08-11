rm(list=ls())
setwd(getwd())

breaker <- function(pcsc, bin_num) {

        cut(pcsc, breaks=unique(quantile(pcsc, seq(0,1,by=1/bin_num), include.lowest = TRUE, na.rm=T)))

}

binner <- function(score, breaks) {

        tapply(score, breaks, mean)

}


dict <- read.table("dummy_dictionary.txt", sep="\t", row.names=1, as.is=T)

for(label in unique(dict[,2])) {
	
	print(paste("First working with", label))
	d <- dict[which(dict[,2]== label),]
	
	pdf_file_name <- paste(label,"_phenotypes_4pops.pdf", sep="")
	table_file_name <- paste(label,"_table_4pops.txt", sep="")

	pdf(pdf_file_name, height=12, width=8)
	m <- data.frame()
	par(mfrow=c(4,3))
	n_breaks <- 50
	for(element in row.names(d)) {
		
		filename <- paste("merged.4pops.", element, sep="")
		if(filename %in% list.files()) {
			
	        x <- read.table(filename, header=T)
	        x <- x[complete.cases(x),]
	        x <- x[which(is.finite(x$chi.sq)),]
	
	        breaks <- breaker(x[,"bs.ld.score"], n_breaks)
	        bs.ld.score <- binner(x[,"bs.ld.score"], breaks)
	        chi.sq <- binner(x[,"chi.sq"], breaks)
	        plot(bs.ld.score, chi.sq, pch=16, col="#4169e180", main=d[element,1], cex.main=1, xlab="Bulik-Sullivan LD score", ylab="Chi-squared")
	        abline(lm(chi.sq ~ bs.ld.score), col="red")
	        m[d[element,1],"BS_intercept_50_bins"] <- lm(chi.sq ~ bs.ld.score)$coefficients[1]
	        m[d[element,1],"BS_intercept_unbinned"] <- lm(x$chi.sq ~ x$bs.ld.score)$coefficients[1]
	
	        breaks <- breaker(x[,"new.ld.score"], n_breaks)
	        new.ld.score <- binner(x[,"new.ld.score"], breaks)
	        chi.sq <- binner(x[,"chi.sq"], breaks)
	        plot(new.ld.score, chi.sq, pch=16, col="#4169e180", main=d[element,1], cex.main=1, xlab="New LD score", ylab="Chi-squared")
	        abline(lm(chi.sq ~ new.ld.score), col="red")
	        m[d[element,1],"OUR_intercept_50_bins"] <- lm(chi.sq ~ new.ld.score)$coefficients[1]
	        m[d[element,1],"OUR_intercept_unbinned"] <- lm(x$chi.sq ~ x$new.ld.score)$coefficients[1]
	
	
	        breaks <- breaker(x[,"pca.score"], n_breaks)
	        pca.score <- binner(x[,"pca.score"], breaks)
	        chi.sq <- binner(x[,"chi.sq"], breaks)
	        plot(pca.score, chi.sq, pch=16, col="#4169e180", main=d[element,1], cex.main=1, xlab="PCA score", ylab="Chi-squared")
	        abline(lm(chi.sq ~ pca.score), col="red")
	        lm(chi.sq ~ pca.score)$coefficients[1]
	        lm(x$chi.sq ~ x$pca.score)$coefficients[1]
	
			
		}
		
	}
	dev.off()
	
	write.table(m, file=table_file_name, sep="\t", quote=F)
	
}
