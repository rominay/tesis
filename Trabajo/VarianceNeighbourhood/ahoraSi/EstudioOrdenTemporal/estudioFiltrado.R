#para un dado valor de std y de perc calculo el grado de cada nodo y hago la distribucion
setwd("C:\\Users\\Elizabeth\\Desktop\\Tesis-Redes\\Practica-R\\PracticaRTesis\\VarianceNeighbourhood\\ahoraSi\\EstudioOrdenTemporal")

library("igraph")

t = seq(0,80,1)    
klim = 11  
genesTotal = 500 
stdd = 1
perc = 0.8

matrizExpr = dget("matrixExp.R")
adjConectado = dget("adjConectado.R")
GrafoConectado = dget("GrafoConectado.R")

matrixx = matrizExpr(genesTotal, t, stdd, perc)
matrix = matrixx[,1:length(t)]
genesShuffled = matrixx[,length(t)+1]

#hago un chequeo de que está haciendo las cosas bien 
#para un gen que no fue shuffleado 
15 %in% genesShuffled 

Corner_text = function(text, location="topright"){
legend(location,legend=text, bty ="n", pch=NA) 
}

infile = paste("ExpreGenes", stdd, "perc", perc, ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(t, matrix[15,], pch="*", cex=2, col="black", xlab="Genes", ylab="",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(2,14))#, xlim = c(0,500)
points(t, matrix[2,],col="red", pch="o", cex=2)			
legend(1,14,legend=c("Varianza vecinos","std"), col=c("black","red"),pch=c("*","o"), cex=1.5,box.lty = 2)
#text(40,14,paste(" eliminados=",  percRechazado, "k=", kc), font=2, cex=1.5)
Corner_text(text="std = 0.2, 500 genes, 5 grupos funcionales, perc = 0.2",location= "topright")

dev.off()


grupos = 5
adjConectado = adjConectado(matrix,t,klim,stdd, perc)
adj = adjConectado + t(adjConectado)
GrafoConectado = GrafoConectado(matrix,t,klim,stdd, perc)

s = shortest.paths(GrafoConectado, 1,81 algorithm = "dijkstra")



kc = sum(adjConectado[1,])
#veo la cantidad maxima de vecinos 
adj[adj==2]=1

ks = sum(adj[1,])
for (i in seq(2,length(t),1)){
	ks = c(ks, sum(adj[i,]) )
}
kmax = max(ks)

distTiempo = function(Vector, fila){
  distanciasTiempo = vector()
  for (j in seq(1,length(Vector),1)){
    if( Vector[j] != 0){
      dist = j - fila
      distanciasTiempo = c(distanciasTiempo,dist)
    }
  }
  distanciasTiempo
}	


DistanciasTiempo = matrix(0L, nrow=length(t), ncol=kmax)

for (i in seq(1,length(t),1)){
    Vector = adj[i,]
    fila = i 
    dist = distTiempo(Vector, fila)
    DistanciasTiempo[i,0:length(dist)] = dist
}


#####PONERLO COMO HISTOGRAMA
DistanciasTiempoVector = vector()
for (i in seq(1,length(t),1)){
  Vector = adj[i,]
  fila = i 
  dist = distTiempo(Vector, fila)
  DistanciasTiempoVector = c(DistanciasTiempoVector,dist)
}

#caso particular para std 0.5 perc 0.6
#DistanciasTiempoVector = c(DistanciasTiempo[,1], DistanciasTiempo[,2])


infile = paste("DistGradoHistograma", stdd, "perc", perc, ".png",sep="")
png(filename=infile, 
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=72)
hist(DistanciasTiempoVector,
     breaks = seq(-12,11,1),
     main="Histograma de la densidad de pasos temporales",
     xlab="Pasos temporales",
     ylab = "Densidad nodos vecinos",
     ylim=c(0,0.13),
     #xlim=c(-7,7),
     col="darkmagenta",
     freq=FALSE,
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
)
Corner_text(text="std genes = 1, 500 genes, 5 grupos funcionales, shuffle = 0.8",location= "topright")

dev.off()
#ahora cuando filtro veo que pasa 
source("filtradoFinal.R")

filt = filtradoFinal(adjConectado,matrix,stdd, perc)

matrixFiltr = filt[,1:length(t)]
genesFilt = filt[,length(t)+1]

adjConectado = dget("adjConectado.R")
adjClean = adjConectado(matrixFiltr,t,klim,stdd, perc)
adjCleanComplt = adjClean + t(adjClean)


DistanciasTiempoVectorClean = vector()
for (i in seq(1,length(t),1)){
  Vector = adjCleanComplt[i,]
  fila = i 
  dist = distTiempo(Vector, fila)
  DistanciasTiempoVectorClean = c(DistanciasTiempoVectorClean,dist)
}

#####PONERLO COMO HISTOGRAMA

infile = paste("DistGradoHistogramaClean", stdd, "perc", perc, ".png",sep="")
png(filename=infile, 
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=72)
hist(DistanciasTiempoVectorClean,
     breaks = seq(-3,2,1),
     main="Histograma de la densidad de pasos temporales",
     xlab="Pasos temporales",
     ylab = "Densidad nodos vecinos",
     ylim=c(0,0.6),
     #xlim=c(-2,2),
     col="darkmagenta",
     freq=FALSE,
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
)
Corner_text(text="std genes = 1, 500 genes, 5 grupos funcionales, shuffle = 0.8",location= "topright")

dev.off()




####para ponerlo como puntos discretos... 

wmax = which.max(DistanciasTiempo)
distMax = DistanciasTiempo[row(DistanciasTiempo)[wmax], col(DistanciasTiempo)[wmax]]

wmin = which.min(DistanciasTiempo)
distMin = DistanciasTiempo[row(DistanciasTiempo)[wmin], col(DistanciasTiempo)[wmin]]

distPosibles = seq(distMin, distMax, 1)

Cantidades = vector()
for (dist in distPosibles){
  Cantidad = 0
  for (i in seq(1,length(DistanciasTiempo[,1]),1)){
    for (j in seq(1,kc,1)){
      if (DistanciasTiempo[i,j] == dist){
        Cantidad = Cantidad + 1 
      }
    }
  }
  Cantidades = c(Cantidades,Cantidad)
}

#divido por 2 para no considerar el doble y normalizo 
Cantidades = Cantidades/(2*length(t))

###prueba
m <- diag(4-abs(-4:4))  # test matrix
wm = which.max(m)
m[row(m)[wm], col(m)[wm]]   # c(5,5) 
###prueba


Corner_text = function(text, location="topright"){
legend(location,legend=text, bty ="n", pch=NA) 
}

infile = paste("DistGrado2STD", stdd, "perc", perc, ".png",sep="")
png(filename=infile, 
	units="in", 
	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
plot(distPosibles,Cantidades, pch=19, cex=1.5, col="black", xlab="Dist tiempo", ylab="Cantidad muestras",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)#, ylim = c(0,1.2), xlim = c(0,1.2))
#lines(percen,percen, pch=19, cex=1.5, col="black", xlab="% Genes aleatorios", ylab="% Genes filtrados",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,1.2), xlim = c(0,1.2))
#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
#legend("bottomleft", legend= ("std = 0.5, 500 genes, 5 grupos funcionales, 4 vecinos"),font=2, cex=1.5)
Corner_text(text="std = 0.5, 500 genes, 5 grupos funcionales, perc = 0.6",location= "topright")

dev.off()
