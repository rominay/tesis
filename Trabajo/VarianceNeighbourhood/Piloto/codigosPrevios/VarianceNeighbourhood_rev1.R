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

#para armarme los genes lo que voy a hacer es pasarle a cada una un c que es una gaussiana random 

t = seq(0,80,1)      #tiempo    

genesTotal = 100 
f = matrix(nrow = genesTotal, ncol = length(t)) #va a ser mi matriz de expresion genetica

#tomo 100 arbitrario asi tengo 5 grupitos de 20 que son de cada una de las trayectorias   

#armo una matriz donde cada fila es una de las funciones y cada columna el tiempo..entonces 
#si me quedo con una fila ya tengo la base de cada gen.
grupos = 5
g = matrix(nrow = grupos, ncol = length(t))
g[1,] = f1(t,0)
g[2,] = f2(t,0)
g[3,] = f3(t,0)
g[4,] = f4(t,0)
g[5,] = f5(t,0)

#cada fila es un gen entonces..voy llenando 
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

matrizFSueve = matrizExpGenetica(g,6)   


#ya tengo mi matriz de expresion genetica
plot(t,matrizFSueve[1,], pch="o", col="black", xlab="samples (tiempo)", ylab="expresion", ylim = c(-10,40))#, xlim = c(0,8))
#plot(t,matrizFSueve[1,])#,ylim=c(-100, 100)) 
points(t,matrizFSueve[22,],col="red", pch="*")
points(t,matrizFSueve[42,],col="blue", pch="*")
points(t,matrizFSueve[62,],col="green", pch="*")
points(t,matrizFSueve[82,],col="pink", pch="*")
# Adding a legend inside box at the location (2,40) in graph coordinates.
# Note that the order of plots are maintained in the vectors of attributes.
legend(1,40,legend=c("Gen tipo 1","Gen tipo 2","Gen tipo 3","Gen tipo 4","Gen tipo 5"), col=c("black","red","blue", "green", "pink"),
       pch=c("o","*","*","*","*"),lty=c(1,2,3,4,5), ncol=1)
text(50, 40, labels = "std = 16")



#ahora me armo la matriz distancia euclideana entre samples 

 

#defino la funcion distancia euclideana entre vectores 


distEuclideana = function(vectA, vectB){
  d = sqrt(sum( (vectA - vectB)^2 ))
  d
}

#B = matrix( c(2, 4, 3, 1, 5, 7), nrow=2, ncol=3, byrow = TRUE)ra        # fill matrix by rows 
 


#armo la matriz E dado que E_ij = d(S_i, S_j)
samples = matrix(nrow = length(t), ncol = length(t)) #uso que es simetrica
for (i in seq(1,length(t),1)){
  for (j in seq(1,length(t),1)){
    samples[i,j] = samples[j,i] = distEuclideana(matrizFSueve[,i],matrizFSueve[,j])
  }
}

###chequeo que mi funcion distEuclideana funca bien
samples2 = matrix(nrow = length(t), ncol = length(t)) #uso que es simetrica 

for (i in seq(1,length(t),1)){
  for (j in seq(1,length(t),1)){
    samples2[i,j] = samples2[j,i] = sqrt(sum((matrizFSueve[,i]-matrizFSueve[,j])*(matrizFSueve[,i]-matrizFSueve[,j])))
  }
}
###fin del chequeo 

###chequeo que tiene traza nula 
sum(diag(samples))
###fin del chequeo 

#####armar la matriz de adyacencia 
matrizRangos = t(apply(samples,1,rank))

#ahora dependiendo de la cantidad de vecinos que quiero permetir que tenga cada nodo(samples) voy a armar
#distintas matrices de adjacencia.
#hagamoslo para k 

adyacencia = function(matrizRangos, k){
			matrizRangos[matrizRangos <= k] = 1
			matrizRangos[matrizRangos > k] = 0
			matrizRangos
		}

#por ejemplo aca tengo la matriz de adyacencia si consideo 
adyacenciaMatrizSampleos = adyacencia(matrizRangos, 2)

###chequeo que cada muestra tiene 5 vecinos 
apply(adyacenciaMatrizSampleos,1,sum)
###fin del chequeo 

##########################################ahora hacer grafo para distintos valores de k y ver que es conexo 

###chequeo tener instalado igraph

"xtable" %in% rownames(installed.packages())
###fin chequeo 

library("igraph") 
grafo = graph_from_adjacency_matrix(adyacenciaMatrizSampleos, mode = c("undirected"), weighted = NULL, diag = FALSE,
  add.colnames = NULL, add.rownames = NA)
#simplify(grafo, remove.multiple = F, remove.loops = T)
#plot(grafo)
#l = layout_on_sphere(grafo) 
plot(grafo, layout = layout_on_sphere(grafo), edge.color = "black", edge.arrow.size=2, edge.width = 2, edge.arrow.width = 2, edge.lty = 1, vertex.size = 10)
text(1, 1, labels = "std = 16, k=1")

###tengo self loop porque cuando hace el rank cuenta al cero de las diagonales.. 
length(adyacenciaMatrizSampleos[,1])
adyacenciaMatrizSampleos[60,63]
###fin de eso 

###para un std = 3 el k para que quede conectado es 3 vecinos. 

###entoncesssss.... calculemos la neighbourhood variance. 

adjMatrix = adyacencia(matrizRangos, 2)
#expMatrix = matrizExpGenetica(g,3)   

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

#para un dado gen, digamos el primero 
#g = 1

##############
#ahora calculo la desviacion estandar....es distinto a std... 

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
					S = 0 
					std = 0
					n = length(expMatrix[1,])
					#kc = 2
					kc = apply(adjMatrix,1,sum)[1]
					#i representa a la sample1,sample2,..
					for (i in seq(1 , length(expMatrix[1,]) , 1)){
						egi = expMatrix[g,i]
						#quiero sumar sobre los tiempos vecinos..
						for (j in seq(1,length(adjMatrix[g,]),1)){
							s  = adjMatrix[g,][j]
							if (s != 0){ # o sea que es mi vecino
								#ese j es el tiempo que me interesa
								egj = expMatrix[g,j]
								S = S + ((egi - egj)^2)
							}
						}
					}
				SN = S/(n*kc-1)
				std = desvStandarCuadrado(expMatrix[g,])/length(expMatrix[1,])
				c(SN,std) 						
}



SN = neighbourhoodVarSTD(expMatrix, adjMatrix,3)
SN2  = neighbourhoodVarSTD(expMatrix, adjMatrix,3)



#std = sd(expMatrix[g,])
#para todos los genes, voy armando un vector

dev = matrix(nrow = length(expMatrix[1,]), ncol = 2)

for (g in seq(1,length(expMatrix[1,]),1)){
	dev[g,] = neighbourhoodVarSTD(expMatrix, adjMatrix,g)
}


#ahora hago un plot de ambas columnas, la neighvariance y la std en funcion de los genes 
genes = seq(1,length(expMatrix[1,]),1)
length(dev[,1])

plot(genes,dev[,1], pch="*", col="black", xlab="Genes", ylab="", ylim = c(-6,130))#, xlim = c(40,60))
points(genes,dev[,2],col="red", pch="o")
# Adding a legend inside box at the location (2,40) in graph coordinates.
# Note that the order of plots are maintained in the vectors of attributes.
legend(1,120,legend=c("Varianza vecinos","std"), col=c("black","red"),
       pch=c("o","*"),lty=c(1,2), ncol=1)
text(50, 120, labels = "std = 3, kc=1")


########################
std = sum(desvStandarNoNormal(expMatrix[g,]))
##############
####################################pruebas
egi = expMatrix[g,1]
j = 2
s  = adjMatrix[g,][j]

for (s in adjMatrix[1,]){
	if (s != 0){
		print(s)
	}	
}


for (s in c(1,2,3)){
	print(s)
}
####################################fin pruebas



################################################################################################

############################################################################################3

#para ver la trayectoria que busco en el tiempo... me quedo con una columna 
#(dado valor de c, dado un gen)

plot(t,f[,1])

#cada columna es un gen y cada fila un tiempo (muestra) 
#where e_ij is the expression of gene j in sample i

#NEIGHBOURHOOD VARIANCE 
#tengo que hallar primero los primeros vecinos de un dado gen (digamos aquellos 2 de menos dist euclideana)
#hagamos para el gen 1

#pienso en un plano que es expresion vs muestra para todos los genes. Veamoslo para los 3 primeros 
expsGen1= f[,1]
expsGen2= f[,2]
expsGen3= f[,3]

plot(t,expsGen1,ylim=c(-100, 100)) 
points(t,expsGen2,col="red", pch="*")
points(t,expsGen3,col="blue", pch="*")
# Adding a legend inside box at the location (2,40) in graph coordinates.
# Note that the order of plots are maintained in the vectors of attributes.
legend(1,100,legend=c("gen1","gen2","gen3"), col=c("black","red","blue"),
                                   pch=c("o","*","+"),lty=c(1,2,3), ncol=1)

plot(t,expsGen1,xlim=c(18, 22), ylim=c(-100, 100)) 
points(t,expsGen2,col="red", pch="*")
points(t,expsGen3,col="blue", pch="*")

###busco los primeros vecinos para el gen 2 a tiempo 20 considerando el entorno de 
#tiempos entre 18 y 22. y 3 primeros genes. 

#################como tengo que calcular distancia euclideana me puede definir esa funcion 
#dado un gen y un tiempo me calcula los vecinos y termino teniendo las samples neighbour.


distEuclideanaPuntos = function(t, e, tGenTarget, eGenTarget){
	( (t - tGenTarget)^2 + (e - eGenTarget)^2 ) ^ (0.5)
}

#################busco los primeros vecinos para el gen 2 a tiempo 20 considerando el entorno de 
#tiempos entre 18 y 22. y 3 primeros genes. 

entornoGen2 = f[t[19:23], 1:3] #armo submatriz de 5 x 3.. un criterio... 

euclideana = matrix(nrow = length(entornoGen2[,1]), ncol = length(entornoGen2[1,]))
for (i in seq(1,length(entornoGen2[,1]),1)){
	for (j in seq(1,length(entornoGen2[1,]),1)){	
		#print(i)
		euclideana[i,j] = distEuclideana(t[18+i],entornoGen2[i,j], t[21], entornoGen2[3,2])
	}
}

#notese que el elemento target que es en el medio de la matriz, el euclideana[3,2] es cero

#ahora puedo decir que me quedo con los 3 que estan a menor distancia euclideana como mis vecinos mas cercanos

euclideanaMenorAMayor = sort(euclideana, decreasing = FALSE)

vecinoPrimero = which(euclideana==euclideanaMenorAMayor[2], arr.ind=TRUE)
vecinoSegundo = which(euclideana==euclideanaMenorAMayor[3], arr.ind=TRUE)
vecinoTercero = which(euclideana==euclideanaMenorAMayor[4], arr.ind=TRUE)
vecinoCuarto = which(euclideana==euclideanaMenorAMayor[5], arr.ind=TRUE)

#entonces ya tengo las coordenadas de los vecinos, lo que me interessa es el tiempo para la formula 
#armo un vector con los indices de esos tiempos. 

indexesTiemposVecinos = c(18+vecinoPrimero[1]-1, 18+vecinoSegundo[1]-1, 18+vecinoTercero[1]-1, 18+vecinoCuarto[1]-1)


###dado el tiempo 20 calculo la primera sumatoria.. 

#sum_j (e_ig - e_N(i,j),g)^2 

g = 2
i = t[21]

primeraSuma = 0 

for (k in seq(1,length(indexesTiemposVecinos),1)){
	primeraSuma = primeraSuma + (f[i , g] - f[indexesTiemposVecinos[k] , g])^2 
	#print(primeraSuma)
} 

####################################################################################################
	
##los vecinos dependen del tiempo(dado un gen)..hago una funcion generica ahora 
## #armo submatriz de 5 x 3.. un criterio...  entorno de mi punto target (tTarget,gTarget) el medio de la matriz es mi target


entornoGenTarget = function(indexTGenTarget, indexEGenTarget){
	#agarro un gen adelante y otro detras, y dos tiempos delantes y dos tiempos detras.. 
	f[(indexTGenTarget-2):(indexTGenTarget+2), (indexEGenTarget-1):(indexEGenTarget+1)] 
}

indexesTiemposVecinos = function(indexTGenTarget, indexEGenTarget){
	entornoGen = entornoGenTarget(indexTGenTarget, indexEGenTarget)
	euclideana = matrix(nrow = length(entornoGen[,1]), ncol = length(entornoGen[1,]))

	for (i in seq(1,length(entornoGen[,1]),1)){       #las samples 
		for (j in seq(1,length(entornoGen[1,]),1)){	
			#print(i)
			indexTiempoBase = indexTGenTarget-2
			indexTiempoTargetMedio = 3  
			indexETargetMedio = 2
			euclideana[i,j] = distEuclideana(t[indexTiempoBase+i],entornoGen[i,j], t[indexTGenTarget+1], entornoGen[indexTiempoTargetMedio ,indexETargetMedio])
		}
	}

	#notese que el elemento target que es en el medio de la matriz, el euclideana[3,2] es cero

	#ahora puedo decir que me quedo con los 3 que estan a menor distancia euclideana como mis vecinos mas cercanos

	euclideanaMenorAMayor = sort(euclideana, decreasing = FALSE)
	#tomo el criterio de quedarme con los primeros 4 cuatro vecinos mas cercanos en dist euclidea
	vecinoPrimero = which(euclideana==euclideanaMenorAMayor[2], arr.ind=TRUE)
	vecinoSegundo = which(euclideana==euclideanaMenorAMayor[3], arr.ind=TRUE)
	vecinoTercero = which(euclideana==euclideanaMenorAMayor[4], arr.ind=TRUE)
	vecinoCuarto = which(euclideana==euclideanaMenorAMayor[5], arr.ind=TRUE)

	#entonces ya tengo las coordenadas de los vecinos, lo que me interessa es el tiempo para la formula 
	#armo un vector con los indices de esos tiempos. 

	c(18+vecinoPrimero[1]-1, 18+vecinoSegundo[1]-1, 18+vecinoTercero[1]-1, 18+vecinoCuarto[1]-1)
}

#######

distEuclideanaConVecinos = function(indexTGenTarget, indexEGenTarget){
	entornoGen = entornoGenTarget(indexTGenTarget, indexEGenTarget)
	euclideana = matrix(nrow = length(entornoGen[,1]), ncol = length(entornoGen[1,]))

	for (i in seq(1,length(entornoGen[,1]),1)){       #las samples 
		for (j in seq(1,length(entornoGen[1,]),1)){	
			#print(i)
			indexTiempoBase = indexTGenTarget-2
			indexTiempoTargetMedio = 3  
			indexETargetMedio = 2
			euclideana[i,j] = distEuclideana(t[indexTiempoBase+i],entornoGen[i,j], t[indexTGenTarget+1], entornoGen[indexTiempoTargetMedio ,indexETargetMedio])
		}
	}
	euclideana
}

distEuclideanaConVecinos(3,2)
distEuclideanaConVecinos(4,2)
distEuclideanaConVecinos(5,2)

plot(t,f[,1], type="b", pch=19, col="red", xlab="samples", ylab="expresion", ylim = c(0,10), xlim = c(0,8))
points(t,f[,2],pch=18, col="blue", type="b", lty=2)
points(t,f[,3],pch=18, col="black", type="b", lty=2)
legend(0,100, legend=c("Gen 1", "Gen 2 target", "Gen 3"),
       col=c("red", "blue", "black"),lty=c(1,2,3), ncol=1)


###dado el tiempo 20 calculo la primera sumatoria.. 
#sum_j (e_ig - e_N(i,j),g)^2 

indexesTiemposVecinos2 = indexesTiemposVecinos(20,2)
g = 2
i = t[21]

primeraSuma = 0 

for (k in seq(1,length(indexesTiemposVecinos2),1)){
	primeraSuma = primeraSuma + (f[i , g] - f[indexesTiemposVecinos2[k] , g])^2 
	#print(primeraSuma)
} 

###lo hago para todos los tiempos para el gen g

g = 2 
neighbourhoodVariance = vector()
samples = vector()
indexesTiemposVecinosDelTargetVector = list()
 
#tengo que ver como hacer con los 3 primeros tiempos... deberian tomar como vecinos en tiempo anteriores de una forma especial

for (i in seq(3,length(t)-2,1)){
	indexesTiemposVecinosDelTarget = indexesTiemposVecinos(i,g) 
	indexesTiemposVecinosDelTargetVector[[i-2]] = indexesTiemposVecinosDelTarget 
	primeraSuma = 0 
	#print(i)
	for (k in seq(1,length(indexesTiemposVecinosDelTarget),1)){
		primeraSuma = primeraSuma + (f[i , g] - f[indexesTiemposVecinosDelTarget[k] , g])^2 
	}
	samples[i-2] = i-2
	neighbourhoodVariance[i-2] = primeraSuma
} 

###prueba 

i=3
indexesTiemposVecinosDelTarget = indexesTiemposVecinos(i,g) 
indexesTiemposVecinosDelTargetVector[[i-2]] = indexesTiemposVecinosDelTarget 
primeraSuma = 0 
#print(i)
for (k in seq(1,length(indexesTiemposVecinosDelTarget),1)){
	primeraSuma = primeraSuma + (f[i , g] - f[indexesTiemposVecinosDelTarget[k] , g])^2 
}
samples[i] = i
neighbourhoodVariance[i] = primeraSuma
 

###fin de prueba


length(neighbourhoodVariance)
length(samples)
length(t)
length(indexesTiemposVecinosDelTargetVector)
#esta lista indexesTiemposVecinosDelTargetVector tiene en cada elemento los indices de los primeros vecinos 
#para el gen g en la muestra i esima. 
#de las 81 samples, hacemos sobre 78 

plot(samples,neighbourhoodVariance)

#ahora calculo la desviacion estandar. 

desvStandarNoNormal = function(e){
	eMean = mean(e)
	desvStandr = vector()
	for (i in  seq(1,(length(e)),1)){
		print(i)
		desvStandr[i] = ( eMean - e[i] )^2
	}
	desvStandr
}

eGen2 = f[,2]
eGen2Recortado = eGen2[3:(length(t)-2)]
length(eGen2Recortado)
desvStandrGen2 = desvStandarNoNormal(eGen2Recortado)
length(desvStandrGen2)		

plot(samples,neighbourhoodVariance,type="b", pch=19, col="red", xlab="samples", ylab="",ylim=c(-100, 37000))
points(samples,desvStandrGen2,pch=18, col="blue", type="b", lty=2)
legend(1,37000, legend=c("Vecinos", "Estandar"),title = "Gen 2",
       col=c("red", "blue"), lty=1:2, cex=0.8)

#ahora lo ultimo para ver si nos quedamos con este gen o no lo que nos importa es la suma total 

standard = sum(desvStandrGen2)
neighbour = sum(neighbourhoodVariance)

if (neighbour > standard){
	print("Rechazado gen")
}
########################################################################
S = 0 
for (i in c){
  for (j in vec){  
    S = (f[i,g] - f[j,g])^2 + S 
  }
}



 