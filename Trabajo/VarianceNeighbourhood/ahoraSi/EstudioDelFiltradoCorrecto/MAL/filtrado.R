filtrado = function(adjConectado,matrixx,stdd, perc){

			kc = sum(adjConectado[,1])

			#a la matriz de adyacencia debo hacerla simetrica 
			adjConectado = adjConectado + t(adjConectado)


			devSTD2 = matrix(nrow = length(matrixx[,1]), ncol = 2)

			for (g in seq(1,length(matrixx[,1]),1)){
				devSTD2[g,] = neighbourhoodVarSTD(matrixx, adjConectado,g)
			}

			genes = seq(1,length(matrixx[,1]),1)

			cantGenes = length(matrixx[,1])

			
			genesSi = 0
			for (i in genes){
				if (devSTD2[,1][i] < devSTD2[,2][i]){
					genesSi = genesSi + 1
				}
			}

			percRechazado = (cantGenes - genesSi)/cantGenes 


			png(filename="Neighbourhood.png", 
    				units="in", 
    				width=10, 
    				height=8, 
    				pointsize=12, 
    				res=72)
			plot(genes,devSTD2[,1], pch="*", cex=2, col="black", xlab="Genes", ylab="", ylim = c(0,70), xlim = c(0,100),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
			points(genes,devSTD2[,2],col="red", pch="o", cex=2)			
			legend(1,70,legend=c("Varianza vecinos","std"), col=c("black","red"),pch=c("*","o"), cex=1.5,box.lty = 2)
			text(80, 70, paste(" eliminados=",  percRechazado), font=2, cex=1.5)
			dev.off()



			#para ver cuales no fueron aceptados 

			genesRechazados = vector()
			for (i in genes){
				if (devSTD2[,1][i] > devSTD2[,2][i]){
					genesRechazados = c(genesRechazados,i)
				}
			}

			#ahora tengo que ver como quedarme solo con esos genes y tener esa matrix
			matrixFinal = matrix(nrow = 0, ncol= length(t))
			for (i in genes){
				if (devSTD2[,1][i] < devSTD2[,2][i]){
					matrixFinal = rbind(matrixFinal,matrixx[i,])
				}
			}
			
			matrixOut = matrix(0L, nrow=length(matrixFinal[,1]), ncol=length(matrixFinal[1,]) + 1)
			matrixOut[,1:length(matrixFinal[1,])]	= matrixFinal 
			matrixOut[1:length(genesRechazados),(length(matrixFinal[1,])+1)] = genesRechazados 
			matrixOut
}






