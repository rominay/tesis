#toma de argumento la adj total
histoCaminoTemporal = function(adj, stdd, perc, string){
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
	DistanciasTiempoVector = vector()
	for (i in seq(1,length(t),1)){
  		Vector = adj[i,]
  		fila = i 
  		dist = distTiempo(Vector, fila)
  		DistanciasTiempoVector = c(DistanciasTiempoVector,dist)
	}

	min = min(DistanciasTiempoVector)-1
	max = max(DistanciasTiempoVector)


	infile = paste("DistGradoHistograma", stdd, "perc", perc, string, ".png",sep="")
	png(filename=infile, 
   	units="in", 
    	width=10, 
    	height=8, 
    	pointsize=12, 
    	res=72)
	hist(DistanciasTiempoVector,
     	breaks = seq(min,max,1),
     	main="Histograma de la densidad de pasos temporales",
     	xlab="Pasos temporales",
     	ylab = "Densidad nodos vecinos",
     	#ylim=c(0,0.13),
     	#xlim=c(-7,7),
     	col="darkmagenta",
     	freq=FALSE,
     	cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
	)
	text(5, 1.2, paste0("std=", stdd," perc=", perc),font=2, cex=1.5)
	dev.off()
	
}
