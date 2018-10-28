neighbourhoodVarSTD = function(expMatrix, adjMatrix,g){
					t = t = seq(1,81,1)
					expres = vector()
					js = vector()
					is = vector()
					S = 0 
					std = 0
					n = length(expMatrix[1,])
					#kc = 2
					kc = apply(adjMatrix,1,sum)[1]
					#i representa a la sample1,sample2,..
					for (i in seq(1 , length(expMatrix[1,]) , 1)){
						egi = expMatrix[g,i]
						#quiero sumar sobre los tiempos vecinos..
						for (j in seq(1 , length(expMatrix[1,]) , 1)){
							s  = adjMatrix[i,][j]
							if (s != 0){ # o sea que es mi vecino
								#ese j es el tiempo que me interesa
								egj = expMatrix[g,j]
								#print(j)
								js = c(js,j)
								is = c(is,i)
								expres = c(expres,egi) 
								S = S + ((egi - egj)^2)
							}
						}
					}
				SN = S/(n*kc-1)
				std = desvStandarCuadrado(expMatrix[g,])/length(expMatrix[1,])
				c(SN,std)
				#info = matrix(nrow = length(js), ncol = 3)
				#info[,1] = expres
				#info[,2] = is
				#info[,3] = js
				#info
}

