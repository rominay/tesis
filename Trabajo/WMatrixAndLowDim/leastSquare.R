setwd("C:\\Users\\Elizabeth\\Desktop\\Tesis-Redes\\Practica-R\\PracticaRTesis\\VarianceNeighbourhood\\rev3\\conShuffle\\MatrixLowDim")

library("limSolve") 

%%%%%%%%%%%%pruebaaaaa
b = c(1,2,3)

A = matrix(nrow = 3, ncol = 3)
A[,1] = c(0.25,0.50,3/4)
A[,2] = c(0.1,2/10,3/10)
A[,3] = c(1,2,3)

E = rep(1, length(b))
f = 1


lsei(A, b, E, f)


%%%%%%%%%%%%%%%%%%

%la onda es que dada la matriz ya mas acotada en genes de expresion calculo para muestra el vector de los W


%tengo la matriz de expresion ya con los genes sacados...matrixFinal
%la que me dice los vecinos es la....adjMatrix 
#la condicion es que sume 1 cada columna de W sume 1
adjMatrixFull = t(adjMatrix) + adjMatrix 

kc = 3                                                   #es clave que tengo que tomar ese valor porque necesito que todo tenga mismo tamano
f = 1

samples = length(matrixFinal[1,])

W = matrix(0L, nrow = samples, ncol= samples)

for (i in seq(1 , length(matrixFinal[1,]) , 1)){
		A = matrix(nrow = length(matrixFinal[,1]), ncol= 0)
		b = matrixFinal[,i]
		vecinos = vector()
		for (j in seq(1 , length(expMatrix[1,]) , 1)){
				s  = adjMatrixFull[i,][j]
				if (s != 0){ # o sea que es mi vecino
					A = cbind(A,matrixFinal[,j])
					vecinos = c(vecinos,j)
				}
		}
		E = rep(1,length(vecinos))
		ww = lsei(A, b, E, f)
		w = array(as.numeric(unlist(ww[1])), dim=length(vecinos))
		#ahora tengo que pones esos coeficientes donde corresponde 
		#y en los otros lugares donde los vecinos no aportan poner 0
		#vector vecinos me dice cuales son los que aportan ese coeficiente 
		for (k in seq(1, length(vecinos),1)){
			W[i,vecinos[k]] = w[k]
		}
}

#entoncs W es una matrix de 81*81 que es la cantidad de samples con 
#construimos M 

###antes prueba de como multiplica matrices esto 
A = matrix(nrow = 2, ncol = 2)
A[1,] = c(1,2)
A[2,] = c(3,4)

B = matrix(nrow = 2, ncol = 2)
B[1,] = c(1,2)
B[2,] = c(3,4)
###fin
I = diag(samples)
M = t(I - W)%*%(I -W)

#ahora tengo que diagonalizarla 

eigen = eigen(M)
eigenValues = eigen$values
eigenVectors = eigen$vectors


#armo L la matriz de baja dimensionalidad 
#W es la matriz que en cada fila describe 
#los coeficientes de sus respectivos vecinos. A cada fila es de los coeficientes 
#para la sample de esa fila i esima. 

#ahora L es de 2 * samples y sus filas son los autovectores de menor autovalor 

#primero eliminamos el 0 

##ojoque vi que tengo un autovalor que va a la -16 entonces lo saco a manopla y 
##tomo los otros dos que le siguen. 

##se que me losordena de mayor a menor entonces... el menor es la columna 
##81
L = matrix(nrow = 2, ncol=samples)
L[1,] = eigenVectors[,80]
L[2,] = eigenVectors[,79]

#ya tengo la matriz de baja dimensionalidad!!!!

##calculo devuelva la distancia euclideana entre vecinos de muestras..

distEuclideana = function(vectA, vectB){
  d = sqrt(sum( (vectA - vectB)^2 ))
  d
}

samples = matrix(nrow = samples, ncol = samples) #uso que es simetrica
for (i in seq(1,length(t),1)){
  for (j in seq(1,length(t),1)){
    samples[i,j] = samples[j,i] = distEuclideana(L[,i],L[,j])
  }
}

matrizRangos = t(apply(samples,1,rank))

adyacencia = function(matrizRangos, k){
			matrizRangos[matrizRangos <= k] = 1
			matrizRangos[matrizRangos > k] = 0
			matrizRangos
		}

adjMatrix = adyacencia(matrizRangos,4)

library("igraph") 
grafo = graph_from_adjacency_matrix(adjMatrix, mode = c("undirected"), weighted = NULL, diag = FALSE,
  add.colnames = NULL, add.rownames = NA)
#simplify(grafo, remove.multiple = F, remove.loops = T)
#plot(grafo)
#l = layout_on_sphere(grafo) 
png(filename="KNearNeighBajaDim.png", 
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=72)
plot(grafo, layout =  layout_with_mds(grafo), edge.color = "black", edge.arrow.size=2, edge.width = 2, edge.arrow.width = 6, edge.lty = 1, vertex.size = 20)
text(1, 1, labels = "std = 3, k=3, low dim d=2",font=2, cex=1.5)
dev.off()

is_connected(grafo)

####################
rbPal = colorRampPalette(c('red','blue'))
datCol = rbPal(10)[as.numeric(cut(t,breaks = 10))]

plot(L[1,],L[2,],pch = 20,col = datCol)
