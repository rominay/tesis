f1 = function(t, c){
	(5 * cos(t/5) + 8) + c 
}

f2 = function(t, c){
  (5 * sin(t/5) + 8) + c 
}

f3 = function(t, c){
  (t)^0.5 + c 
}

f4 = function(t, c){
  ( (t/20)^2 ) + c 
}

f5 = function(t, c){
  (16 - (t/20))^2 + c 
}

t = seq(0,80,1)      #tiempo    

genesTotal = 100 

f = matrix(nrow = genesTotal, ncol = length(t)) #va a ser mi matriz de expresion genetica

#tomo 100 arbitrario asi tengo 5 grupitos de 20 que son de cada una de las trayectorias   

#armo una matriz donde cada fila es una de las funciones y cada columna el tiempo..entonces 
#si me quedo con una fila ya tengo la base de cada gen.
grupos = 5
curv = matrix(nrow = grupos, ncol = length(t))
curv [1,] = f1(t,0)
curv [2,] = f2(t,0)
curv [3,] = f3(t,0)
curv [4,] = f4(t,0)
curv [5,] = f5(t,0)

matrizExpGenetica = function(g, stdd){
  for (i in seq(1,genesTotal,1)){
    if (i <= 20){                                 #el primer grupito
      f[i,] = g[1,] + rnorm(length(t), mean = 0, sd = stdd)
    }
    if (20 < i & i <=40){                                 #el primer grupito
      f[i,] = g[2,] + rnorm(length(t), mean = 0, sd = stdd)
    }
    if (40 < i & i <= 60){                                 #el primer grupito
      f[i,] = g[3,] + rnorm(length(t), mean = 0, sd = stdd)
    }
    if (60 < i & i <= 80){                                 #el primer grupito
      f[i,] = g[4,] + rnorm(length(t), mean = 0, sd = stdd)
    }
    if (80 < i & i <= 100){                                 #el primer grupito
      f[i,] = g[5,] + rnorm(length(t), mean = 0, sd = stdd)
    }
  }
  f
}

matrizFSueve = matrizExpGenetica(curv,0.2)  

distEuclideana = function(vectA, vectB){
  d = sqrt(sum( (vectA - vectB)^2 ))
  d
}

samples = matrix(nrow = length(t), ncol = length(t)) #uso que es simetrica
for (i in seq(1,length(t),1)){
  for (j in seq(1,length(t),1)){
    samples[i,j] = samples[j,i] = distEuclideana(matrizFSueve[,i],matrizFSueve[,j])
  }
}

matrizRangos = t(apply(samples,1,rank))

adyacencia = function(matrizRangos, k){
			matrizRangos[matrizRangos <= k] = 1
			matrizRangos[matrizRangos > k] = 0
			matrizRangos
		}


adjMatrix = adyacencia(matrizRangos, 2)

library("igraph") 
grafo = graph_from_adjacency_matrix(adyacenciaMatrix, mode = c("undirected"), weighted = NULL, diag = FALSE,
  add.colnames = NULL, add.rownames = NA)
#simplify(grafo, remove.multiple = F, remove.loops = T)
#plot(grafo)
#l = layout_on_sphere(grafo) 
plot(grafo, layout = layout_on_sphere(grafo), edge.color = "black", edge.arrow.size=2, edge.width = 2, edge.arrow.width = 2, edge.lty = 1, vertex.size = 10)
text(1, 1, labels = "std = 0.2, k=1")

expMatrix  = matrizFSueve 

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


adjMatrixNoSelf = adjMatrix - diag
adjMatrix = adjMatrixNoSelf 

desvStandarCuadrado = function(e){
	desvStandr = 0
	eMean = mean(e)
	#desvStandr = vector()
	for (i in  seq(1,(length(e)),1)){
		#print(i)
		desvStandr = desvStandr + ( eMean - e[i] )^2
	}
	std = desvStandr
	std
}


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




devSTD2 = matrix(nrow = length(expMatrix[,1]), ncol = 2)

for (g in seq(1,length(expMatrix[,1]),1)){
	devSTD2[g,] = neighbourhoodVarSTD(expMatrix, adjMatrix,g)
}

genes = seq(1,length(expMatrix[,1]),1)



plot(genes,devSTD2[,1], pch="*", col="black", xlab="Genes", ylab="", ylim = c(0,24))#, xlim = c(40,60))
points(genes,devSTD2[,2],col="red", pch="o")
# Adding a legend inside box at the location (2,40) in graph coordinates.
# Note that the order of plots are maintained in the vectors of attributes.
legend(1,20,legend=c("Varianza vecinos","std"), col=c("black","red"),
       pch=c("o","*"),lty=c(1,2), ncol=1)
text(50, 20, labels = "std = 0.2, kc=1")




