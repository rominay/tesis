
#voy a correr para un dado std para varios perc 

setwd("C:\\Users\\Elizabeth\\Desktop\\Tesis-Redes\\Practica-R\\PracticaRTesis\\VarianceNeighbourhood\\ahoraSi\\EstudioOrdenTemporal")
####llamar funciones
#matrixExp = dget("matrixExp.R")

#matrixx = matrixExp(genesTotal, t, stdd, perc)
###llamar funciones 

library("igraph")

t = seq(0,80,1)    
klim = 11  
genesTotal = 500 
stdd = 1
percen = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
PercFilt = vector()
PercInter = vector()
for(perc in percen){
	matrizExpr = dget("matrixExp.R")
	adjConectado = dget("adjConectado.R")
	#
	matrixx = matrizExpr(genesTotal, t, stdd, perc)
	matrix = matrixx[,1:length(t)]
	genesShuffled = matrixx[,length(t)+1]
	grupos = 5
	adjConectado = adjConectado(matrix,t,klim,stdd, perc)
	adj = adjConectado + t(adjConectado) 

	source("histoCaminoTemporal.R")
	histoCaminoTemporal(adj)
	
	source("filtradoFinal.R")

	filt = filtradoFinal(adjConectado,matrix,stdd, perc)

	matrixFiltr = filt[,1:length(t)]
	genesFilt = filt[,length(t)+1]
	
	percFilt = length(genesFilt[ genesFilt != 0 ])/genesTotal 
	percInter = (length(intersect(genesShuffled,genesFilt))-1)/(length(genesShuffled)*perc)

	PercFilt = c(PercFilt,percFilt)
	PercInter = c(PercInter,percInter)

	adjClean = adjConectado(matrixFiltr,t,klim,stdd, perc)
	adjCleanComplt = adjClean + t(adjClean)
}

if (length(PercFilt) != 10){
	PercFilt = c(PercFilt, 1)
	PercInter = c(PercInter, 1)
}


Corner_text = function(text, location="topright"){
legend(location,legend=text, bty ="n", pch=NA) 
}

infile = paste("FiltradosVsShufflestd", stdd, ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(percen,PercFilt, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="% Genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,1.2), xlim = c(0,1.2))
lines(percen,percen, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="% Genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,1.2), xlim = c(0,1.2))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
#legend("bottomleft", legend= ("std = 0.5, 500 genes, 5 grupos funcionales, 4 vecinos"),font=2, cex=1.5)
Corner_text(text="std = 1, 500 genes, 5 grupos funcionales",location= "topright")

dev.off()
		
infile = paste("CoincVsShufflestd",stdd, ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(percen,PercInter, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="% Overlap genes aleatorios y filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0.3,1.2), xlim = c(0,1.2))
lines(percen,rep(1,10), pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="% Overlap genes aleatorios y filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,1.2), xlim = c(0,1.2))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
#text(0.5, 1.1 , ("std = 2, 500 genes, 5 grupos funcionales"),font=2, cex=1.5)
Corner_text(text="std = 1, 500 genes, 5 grupos funcionales",location= "topright")

dev.off()











#######
setwd("C:\\Users\\Elizabeth\\Desktop\\Tesis-Redes\\Practica-R\\PracticaRTesis\\VarianceNeighbourhood\\rev3\\conShuffle\\MatrixLowDim\\buenShuffle\\ahoraSi")

t = seq(0,80,1)      
genesTotal = 500 
stdd = 3 
perc = 0.1
matrixx = matrizExpr(genesTotal, t, stdd, perc)
grupos = 5
adjDadosVecinos = adjDadosVecinos(matrixx,t,10,stdd)
filtrado = filtrado(adjDadosVecinos,matrixx,stdd, perc)

length(filtrado[,1])