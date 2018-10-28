#######
#voy a correr para un dado perc para varios stdd 

setwd("C:\\Users\\Elizabeth\\Desktop\\Tesis-Redes\\Practica-R\\PracticaRTesis\\VarianceNeighbourhood\\rev3\\conShuffle\\MatrixLowDim\\buenShuffle\\ahoraSi\\pruebas2")

library("igraph")

t = seq(0,80,1)    
klim = 11  
genesTotal = 500 
perc = 0.2 
std = c(0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3)
PercFilt = vector()
PercInter = vector()
for(stdd in std){
	matrizExpr = dget("matrixExp.R")
	adjConectado = dget("adjConectado.R")
	#
	matrixx = matrizExpr(genesTotal, t, stdd, perc)
	matrix = matrixx[,1:length(t)]
	genesShuffled = matrixx[,length(t)+1]
	grupos = 5
	adjConectado = adjConectado(matrix,t,klim,stdd, perc)
	
	source("filtradoFinal.R")

	filt = filtradoFinal(adjConectado,matrix,stdd, perc)

	matrixFiltr = filt[,1:length(t)]
	genesFilt = filt[,length(t)+1]
	
	percFilt = length(genesFilt[ genesFilt != 0 ])/genesTotal 
	percInter = (length(intersect(genesShuffled,genesFilt))-1)/(length(genesShuffled)*perc)

	PercFilt = c(PercFilt,percFilt)
	PercInter = c(PercInter,percInter)
}

Corner_text = function(text, location="topright"){
legend(location,legend=text, bty ="n", pch=NA) 
}
percFiltEsperado = rep(perc,length(std))

infile = paste("FiltradosVsShufflePerc0.2", ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(std,PercFilt, pch=19, cex=1.5, col="black", xlab="Std", ylab="% Genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,1.2), xlim = c(0,3.2))
lines(std,percFiltEsperado , pch=19, cex=1.5, col="black", xlab="Std", ylab="% Genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim=c(0,1.2), xlim = c(0,3.2))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
#legend("bottomleft", legend= ("std = 0.5, 500 genes, 5 grupos funcionales, 4 vecinos"),font=2, cex=1.5)
Corner_text(text="Perc = 0.2, 500 genes, 5 grupos funcionales",location= "topright")

dev.off()
		
infile = paste("CoincVsShufflePerc0.2", ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(std,PercInter, pch=19, cex=1.5, col="black", xlab="Std", ylab="% Overlap genes aleatorios y filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim=c(0,1.2), xlim = c(0,3.2))
lines(std,rep(1,length(std)), pch=19, cex=1.5, col="black", xlab="Std", ylab="% Overlap genes aleatorios y filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim=c(0,1.2), xlim = c(0,3.2))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
text(0.5, 1.1 , ("Perc = 0.2, 500 genes, 5 grupos funcionales"),font=2, cex=1.5)

dev.off()