adjDadosVecinos = function(matrixExprShuffled,t,k,stdd){
			distEuclideana = function(vectA, vectB){
  				d = sqrt(sum( (vectA - vectB)^2 ))
  				d
			}

			matrizRangos = function(matrizExpr, t){
				samples = matrix(nrow = length(t), ncol = length(t)) #uso que es simetrica
				for (i in seq(1,length(t),1)){
  					for (j in seq(1,length(t),1)){
    						samples[i,j] = samples[j,i] = distEuclideana(matrizExpr[,i],matrizExpr[,j])
  					}
				}
				matrizRangos = t(apply(samples,1,rank))
				matrizRangos 
			}



			adyacencia = function(matrizRangos, k){
				matrizRangos[matrizRangos <= k] = 1
				matrizRangos[matrizRangos > k] = 0
				matrizRangos
			}

			#le saco los 1 de la diagonal en la adjMatrix 
			diag = matrix(nrow = length(t), ncol = length(t))
			for (i in seq(1,length(adjMatrix[1,]),1)){
				for (j in seq(1,length(adjMatrix[1,]),1)){	
					if (i==j){
						diag[i,i] = 1
					}
					if (i != j) {
						diag[i,j] = 0	
					}
				}
			}

			rangos = matrizRangos(matrixExprShuffled,t)

			#lo hago hasta que me diga TRUE 
			adjFinal = adyacencia(rangos,i)
			sinSelf = adjFinal - diag			
			png(filename="GrafoConectado.png", 
    			units="in", 
    			width=10, 
    			height=8, 
    			pointsize=12, 
    			res=72)
			plot(grafo, layout = layout_on_sphere(grafo), edge.color = "black", edge.arrow.size=2, edge.width = 2, edge.arrow.width = 2, edge.lty = 1, vertex.size = 10)
			text(1, 1, paste0("std=", stdd," k=", k),font=2, cex=1.5)
			dev.off()
			sinSelf 
}
