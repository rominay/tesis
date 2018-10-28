#######
#voy a correr para un dado perc y std para varios valores de k 
#vi que para std=1 y perc=0.2 funciona bien... 
#uso una funcion que no se queda con el k que hace el grafo conectado sino que lo toma de input 


setwd("C:\\Users\\Elizabeth\\Desktop\\Tesis-Redes\\Practica-R\\PracticaRTesis\\VarianceNeighbourhood\\rev3\\conShuffle\\MatrixLowDim\\buenShuffle\\ahoraSi\\pruebas2")

library("igraph")

t = seq(0,80,1)      
genesTotal = 500 
perc = 0.5 
std = 2
ks = 2:11
PercFilt = vector()
PercInter = vector()
for(k in ks){
	matrizExpr = dget("matrixExp.R")
	adjDadosVecinos = dget("adjDadosVecinos.R")
	#
	matrixx = matrizExpr(genesTotal, t, stdd, perc)
	matrix = matrixx[,1:length(t)]
	genesShuffled = matrixx[,length(t)+1]
	grupos = 5
	adjDadosVecinos = adjDadosVecinos(matrix,t,k,stdd, perc)
	
	source("filtradoFinal.R")

	filt = filtradoFinal(adjDadosVecinos ,matrix,stdd, perc)

	matrixFiltr = filt[,1:length(t)]
	genesFilt = filt[,length(t)+1]
	
	percFilt = length(genesFilt[ genesFilt != 0 ])/genesTotal 
	percInter = (length(intersect(genesShuffled,genesFilt))-1)/(length(genesShuffled)*perc)

	PercFilt = c(PercFilt,percFilt)
	PercInter = c(PercInter,percInter)
}

###para los graficos
Corner_text = function(text, location="topright"){
legend(location,legend=text, bty ="n", pch=NA) 
}
###para los graficos

percFiltEsperado = rep(perc,length(ks))

infile = paste("FiltradosVsShufflePerc0.2Std1Ks", ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(ks,PercFilt, pch=19, cex=1.5, col="black", xlab="k", ylab="% Genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,1.2), xlim = c(0,11))
lines(ks,percFiltEsperado , pch=19, cex=1.5, col="black", xlab="k", ylab="% Genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim=c(0,1.2), xlim = c(0,11))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
#legend("bottomleft", legend= ("std = 0.5, 500 genes, 5 grupos funcionales, 4 vecinos"),font=2, cex=1.5)
Corner_text(text="Perc = 0.5, std = 2, 500 genes, 5 grupos funcionales, kc =4",location= "topright")

dev.off()
		
infile = paste("CoincVsShufflePerc0.2Std1Ks", ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(ks,PercInter, pch=19, cex=1.5, col="black", xlab="k", ylab="% Overlap genes aleatorios y filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim=c(0,1.2), xlim = c(0,11))
lines(ks,rep(1,length(ks)), pch=19, cex=1.5, col="black", xlab="k", ylab="% Overlap genes aleatorios y filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim=c(0,1.2), xlim = c(0,11))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
Corner_text(text=("Perc = 0.5, std = 2, 500 genes, 5 grupos funcionales, kc = 4"),location= "topright")

dev.off()