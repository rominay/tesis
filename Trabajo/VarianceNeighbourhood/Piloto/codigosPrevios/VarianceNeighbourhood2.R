f1 = function(t, c){
	5 * c * cos(t/5)
}

t = seq(0,80,1)      #tiempo    
c = seq(1,600,10)    #genes
f = matrix(nrow = length(t), ncol = length(c)) 

for (i in seq(1,length(t),1)){
  for (j in seq(1,length(c),1)){  
    f[i,j] = f1(t[i],c[j])
  }
}


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

plot(t,f[,1], type="b", pch=19, col="red", xlab="samples", ylab="expresion")#, ylim = c(0,10), xlim = c(0,8))
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



 