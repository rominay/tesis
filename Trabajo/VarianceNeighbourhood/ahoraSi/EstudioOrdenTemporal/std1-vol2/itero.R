
#voy a correr para un dado std para varios perc 

setwd("C:\\Users\\Elizabeth\\Desktop\\Tesis-Redes\\Practica-R\\PracticaRTesis\\VarianceNeighbourhood\\ahoraSi\\EstudioOrdenTemporal\\iteraciones")

###llamar funciones 

library("igraph")

t = seq(0,80,1)    
klim = 11  
genesTotal = 500 
stdd = 1
percen = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
SInt = vector()
SAft = vector()
PercFilt = vector()
PercInter = vector()
it = seq(1,20,1)

sIntIt = vector()
SdInt = vector()

percFiltIt = vector()
percInterIt = vector()
SdpercFilt = vector()
SdpercInter = vector()

SAftIt = vector()
SdAft = vector()
for(perc in percen){
	for (i in it){
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
		#histoCaminoTemporal(adj,stdd, perc, 'init')


		GrafoConectado = dget("GrafoConectado.R")
		GrafoInt = GrafoConectado(matrix,t,klim,stdd, perc)


		sInt = shortest.paths(GrafoInt, 1,81, algorithm = "dijkstra")
		sIntIt = c(sIntIt, sInt) 
	
		source("filtradoFinal.R")

		filt = filtradoFinal(adjConectado,matrix,stdd, perc)

		matrixFiltr = filt[,1:length(t)]
		genesFilt = filt[,length(t)+1]
	
		percFilt = length(genesFilt[ genesFilt != 0 ])/genesTotal 
		percInter = (length(intersect(genesShuffled,genesFilt))-1)/(length(genesShuffled)*perc)

		percFiltIt = c(percFiltIt, percFilt)
		percInterIt = c(percInterIt, percInter)

		
		adjConectado = dget("adjConectado.R")
		adjClean = adjConectado(matrixFiltr,t,klim,stdd, perc)
		adjCleanComplt = adjClean + t(adjClean)
	
		#histoCaminoTemporal(adjCleanComplt,stdd, perc, 'filt')

		GrafoAft = GrafoConectado(matrixFiltr,t,klim,stdd, perc)

		sAft = shortest.paths(GrafoAft, 1,81, algorithm = "dijkstra")
		SAftIt = c(SAftIt, sAft)
	}
	sIntProm = mean(sIntIt)
	sdInt = sd(sIntIt)
	SInt = c(SInt, sIntProm)
	SdInt = c(SdInt, sdInt)

	percFilt = mean(percFiltIt)
	sdpercFilt = sd(percFiltIt)
	SdpercFilt = c(SdpercFilt,sdpercFilt)
	PercFilt = c(PercFilt,percFilt)

	percInter = mean(percInterIt)
	sdpercInter = sd(percFiltIt)
	SdpercInter = c(SdpercFilt,sdpercFilt)
	PercInter = c(PercInter,percInter)

	sAftProm = mean(SAftIt)
	sdAft = sd(SAftIt)
	SAft = c(SAft, sAftProm)
	SdAft = c(SdAft, sdAft)
}




#if (length(PercFilt) != 10){
#	PercFilt = c(PercFilt, 1)
#	PercInter = c(PercInter, 1)
#}


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
plot(percen,PercFilt, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="% Genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,1.2), xlim = c(0,1))
lines(percen,percen, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="% Genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,1.2), xlim = c(0,1))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
#legend("bottomleft", legend= ("std = 0.5, 500 genes, 5 grupos funcionales, 4 vecinos"),font=2, cex=1.5)
Corner_text(text="std = 1, 20 iteraciones, 500 genes, 5 grupos funcionales",location= "topright")

dev.off()

infile = paste("STPercFiltstd", stdd, ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(percen,SdpercFilt, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="ST % genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)#, ylim = c(0,1.2), xlim = c(0,1))
#lines(percen,percen, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="% Genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,1.2), xlim = c(0,1))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
#legend("bottomleft", legend= ("std = 0.5, 500 genes, 5 grupos funcionales, 4 vecinos"),font=2, cex=1.5)
Corner_text(text="std = 1, 20 iteraciones, 500 genes, 5 grupos funcionales",location= "topright")

dev.off()
		
infile = paste("CoincVsShufflestd",stdd, ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(percen,PercInter, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="% Overlap genes aleatorios y filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0.3,1.2), xlim = c(0,1))
lines(percen,rep(1,9), pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="% Overlap genes aleatorios y filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,1.2), xlim = c(0,1))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
#text(0.5, 1.1 , ("std = 2, 500 genes, 5 grupos funcionales"),font=2, cex=1.5)
Corner_text(text="std = 1, 20 iteraciones, 500 genes, 5 grupos funcionales",location= "topright")

dev.off()

SdpercInter = SdpercInter[1:9]
infile = paste("STCoincstd", stdd, ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(percen,SdpercInter, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="ST overlap",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)#, ylim = c(0,1.2), xlim = c(0,1))
#lines(percen,percen, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="% Genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,1.2), xlim = c(0,1))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
#legend("bottomleft", legend= ("std = 0.5, 500 genes, 5 grupos funcionales, 4 vecinos"),font=2, cex=1.5)
Corner_text(text="std = 1, 20 iteraciones, 500 genes, 5 grupos funcionales",location= "topright")

dev.off()

infile = paste("ShortestPath",stdd, ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(percen,SInt, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="Shortest path global",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(5,85))#, xlim = c(0,1.2))
lines(percen,rep(80,9), pch=19, cex=1.5, col="red", xlab="% Genes aleatorios", ylab="Shortest path global",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(5,85))#, xlim = c(0,1.2))
points(percen,SAft, pch=19, cex=1.5, col="red", xlab="% Genes aleatorios", ylab="Shortest path global",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(5,85))#, xlim = c(0,1.2))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
#text(0.5, 1.1 , ("std = 2, 500 genes, 5 grupos funcionales"),font=2, cex=1.5)
legend(0.1,75,legend=c("Sin filtrar","Filtrado"), col=c("black","red"),pch=19, cex=1.5,box.lty = 2)
Corner_text(text="std = 1, 20 iteraciones, 500 genes, 5 grupos funcionales",location= "topright")

dev.off()

infile = paste("STShortPathInt", stdd, ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(percen,SdInt, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="ST short path sin filtro",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)#, ylim = c(0,1.2), xlim = c(0,1))
#lines(percen,percen, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="% Genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,1.2), xlim = c(0,1))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
#legend("bottomleft", legend= ("std = 0.5, 500 genes, 5 grupos funcionales, 4 vecinos"),font=2, cex=1.5)
Corner_text(text="std = 1, 20 iteraciones, 500 genes, 5 grupos funcionales",location= "topright")

dev.off()

infile = paste("STShortPathAft", stdd, ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(percen,SdAft, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="ST short path filtrado",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)#, ylim = c(0,1.2), xlim = c(0,1))
#lines(percen,percen, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="% Genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,1.2), xlim = c(0,1))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
#legend("bottomleft", legend= ("std = 0.5, 500 genes, 5 grupos funcionales, 4 vecinos"),font=2, cex=1.5)
Corner_text(text="std = 1, 20 iteraciones, 500 genes, 5 grupos funcionales",location= "topright")

dev.off()
