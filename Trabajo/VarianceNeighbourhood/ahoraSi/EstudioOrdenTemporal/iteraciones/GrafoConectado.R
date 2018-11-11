###dada la matriz de expresion genetica...

###hago la matriz de adyacencia tal que queda un grafo conexo 
GrafoConectado = function(matrixExprShuffled,t,klim,stdd, perc){
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
			for (i in seq(1,length(t),1)){
				for (j in seq(1,length(t),1)){	
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
			adjFinal = function(rangos,klim){ 
					for (i in seq_along(1:klim)){
						adjMatrix = adyacencia(rangos,i)
						grafo = graph_from_adjacency_matrix(adjMatrix, mode = c("undirected"), weighted = NULL, diag = FALSE,
 	 					add.colnames = NULL, add.rownames = NA)
						if (is_connected(grafo) == TRUE){
							adjFinal = adyacencia(rangos,i)
							break	
						}
					}
					sinSelf = adjFinal - diag
					sinSelf 			
			}

			adjFinal = adjFinal(rangos, klim)
			grafo = graph_from_adjacency_matrix(adjFinal, mode = c("undirected"), weighted = NULL, diag = FALSE,
 	 					add.colnames = NULL, add.rownames = NA)
			kc = sum(adjFinal[1,])
			infile = paste("GrafoConectado", perc, "k",kc, "std", stdd, ".png",sep="")
			png(filename=infile, 
    			units="in", 
    			width=10, 
    			height=8, 
    			pointsize=12, 
    			res=72)
			plot(grafo, layout = layout_on_sphere(grafo), edge.color = "black", edge.arrow.size=2, edge.width = 2, edge.arrow.width = 2, edge.lty = 1, vertex.size = 10)
			text(1, 1, paste0("std=", stdd," kc=", kc, " perc=", perc),font=2, cex=1.5)
			dev.off()
			grafo
			#return(list(adjFinal=adjFinal,g=grafo))
}

#adjConectado = adjConectado(matrixx,t,5,3)
