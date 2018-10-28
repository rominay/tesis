setwd("C:\\Users\\Elizabeth\\Desktop\\Tesis-Redes\\Practica-R\\PracticaRTesis\\VarianceNeighbourhood\\rev3\\conShuffle\\MatrixLowDim\\buenShuffle")


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

genesTotal = 500 

f = matrix(nrow = genesTotal, ncol = length(t)) #va a ser mi matriz de expresion genetica


#armo una matriz donde cada fila es una de las funciones y cada columna el tiempo..entonces 
#si me quedo con una fila ya tengo la base de cada gen.
grupos = 5
curv = matrix(nrow = grupos, ncol = length(t))
curv [1,] = f1(t,0)
curv [2,] = f2(t,0)
curv [3,] = f3(t,0)
curv [4,] = f4(t,0)
curv [5,] = f5(t,0)

matrizExpGenetica = function(g, stdd,genesGrupo){
  for (i in seq(1,genesTotal,1)){
    if (i <= genesGrupo){                                 #el primer grupito
      f[i,] = g[1,] + rnorm(length(t), mean = 0, sd = stdd)
    }
    if (genesGrupo < i & i <= 2*genesGrupo){                                 #el primer grupito
      f[i,] = g[2,] + rnorm(length(t), mean = 0, sd = stdd)
    }
    if (2*genesGrupo < i & i <= 3*genesGrupo){                                 #el primer grupito
      f[i,] = g[3,] + rnorm(length(t), mean = 0, sd = stdd)
    }
    if (3*genesGrupo < i & i <= 4*genesGrupo){                                 #el primer grupito
      f[i,] = g[4,] + rnorm(length(t), mean = 0, sd = stdd)
    }
    if (4*genesGrupo < i & i <= 5*genesGrupo){                                 #el primer grupito
      f[i,] = g[5,] + rnorm(length(t), mean = 0, sd = stdd)
    }
  }
  f
}

matrizFSuave = matrizExpGenetica(curv,3,100)

#tengo que hacer shuffle...pruebo la funcion sample apra eso 
data = matrix(1:100,nrow=10) #mtrix of 5 rows and 20 columns
x = data[sample(1:10),] #permute rows
y =data[,sample(1:10)]# the columns

#el tema de la probabilidad
N = 100
n = length(50:90)

set.seed(1)
random.sample1 = sample(50:90, N, replace=TRUE, prob=rep(1/n, times=n))

set.seed(1)
random.sample2 = sample(50:90, N, replace=TRUE)

plot(density(random.sample1),col="blue")
lines(density(random.sample2),col="red")
#hace una probabilidad equiprobable si el replace es TRUE. 

#la onda es ir agarrando grupitos de 5, una de cada serie, hacer un reshuffle y devolverlas a su lugar. 
matrixExpression = matrizFSueve
###ya no lo uso 
for (i in seq(1,length(genes)/grupos-1,1)){
	matrixForReshuffle = matrix(nrow = grupos, ncol=length(t))
	matrixReshuffled = matrix(nrow = grupos, ncol=length(t))
	matrixForReshuffle[1,] = matrixExpression[i,]
	matrixForReshuffle[2,]= matrixExpression[i+20,]
	matrixForReshuffle[3,]= matrixExpression[i+40,]
	matrixForReshuffle[4,]= matrixExpression[i+60,]
	matrixForReshuffle[5,]= matrixExpression[i+80,]
	matrixReshuffled = matrixForReshuffle[sample(1:grupos),]
	matrixExpression[i,] = matrixReshuffled[1,]
	matrixExpression[i+20,] = matrixReshuffled[2,]
	matrixExpression[i+40,] = matrixReshuffled[3,]
	matrixExpression[i+60,] = matrixReshuffled[4,]
	matrixExpression[i+80,] = matrixReshuffled[5,]
}


matrizFSueve = matrixExpression #ahora ya con el reshuffle
### 

##ahora uso esto 

random = function(grupos, genesTotal, matrix, perc){
  for(i in 1:grupos){
    genesGrupo = genesTotal/grupos
    FilasRnd = sample((1:genesGrupo) + (i-1)*genesGrupo,perc*genesGrupo ) #(i-1)*genesGrupo es un offset al ir pasando
    for(j in seq_along(FilasRnd)){
      matrix[FilasRnd[j],] = sample(matrix[FilasRnd[j],], ncol(matrix)) #hago un shuffle de esas filas rnd que agarre. 
    }
  }
  matrix
}

matrixExprShuffled = random(5, 500, matrizFSuave, 0.1)

##a ver elshuffleado 

png(filename="Shuffle1Std3.png", 
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=72)
plot(matrizFSuave,matrixExprShuffled, pch="*", col="black", xlab="No shuffle", ylab="Shuffle",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,270), xlim = c(0,270))
#legend(100,200,legend=c("Gen1","Gen2","Gen3"), col=c("black","red","blue"),
#       pch=c("*","o","o"),lty=c(1,2,3), ncol=1,cex=1.5,box.lty = 2)
text(100, 265, labels = "500 genes, std = 3, shuffle 100%",font=2, cex=1.5)
dev.off()


png(filename="ExpresionShuffle0.1Std3.png", 
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=72)
plot(t,matrixExprShuffled[1,], pch="*", col="black", xlab="Tiempo", ylab="Expres",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,20))#, xlim = c(0,100))
points(t,matrixExprShuffled[5,],col="red", pch="o")
points(t,matrixExprShuffled[8,],col="blue", pch="o")
# Adding a legend inside box at the location (2,40) in graph coordinates.
# Note that the order of plots are maintained in the vectors of attributes.
legend(65,20,legend=c("Gen1","Gen2","Gen3"), col=c("black","red","blue"),
       pch=c("*","o","o"),lty=c(1,2,3), ncol=1,cex=1.5,box.lty = 2)
text(20, 20, labels = "Función coseno,std = 3, sin shuffle.",font=2, cex=1.5)
dev.off()


distEuclideana = function(vectA, vectB){
  d = sqrt(sum( (vectA - vectB)^2 ))
  d
}

matrizRangos = function(matrizExpr, t){
	matrix(nrow = length(t), ncol = length(t)) #uso que es simetrica
	for (i in seq(1,length(t),1)){
  		for (j in seq(1,length(t),1)){
    			samples[i,j] = samples[j,i] = distEuclideana(matrizFSueve[,i],matrizFSueve[,j])
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

###prueba
adjMatrix = adyacencia(matrizRangos,4)

library("igraph") 
grafo = graph_from_adjacency_matrix(adjMatrix, mode = c("undirected"), weighted = NULL, diag = FALSE,
  add.colnames = NULL, add.rownames = NA)
#simplify(grafo, remove.multiple = F, remove.loops = T)
#plot(grafo)
#l = layout_on_sphere(grafo) 
png(filename="GrafoConShuffleStd3k3.png", 
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=72)
plot(grafo, layout = layout_on_sphere(grafo), edge.color = "black", edge.arrow.size=2, edge.width = 2, edge.arrow.width = 2, edge.lty = 1, vertex.size = 10)
text(1, 1, labels = "std = 3, k=3",font=2, cex=1.5)
dev.off()

###fin prueba

library("igraph") 


#lo hago hasta que me diga TRUE 

for (i in seq_along(1:10)){
	adjMatrix = adyacencia(matrizRangos,i)
	grafo = graph_from_adjacency_matrix(adjMatrix, mode = c("undirected"), weighted = NULL, diag = FALSE,
 	 add.colnames = NULL, add.rownames = NA)
	if ()
}




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

#a la matriz de adyacencia debo hacerla simetrica 
adjMatrix = adjMatrix + t(adjMatrix)


devSTD2 = matrix(nrow = length(expMatrix[,1]), ncol = 2)

for (g in seq(1,length(expMatrix[,1]),1)){
	devSTD2[g,] = neighbourhoodVarSTD(expMatrix, adjMatrix,g)
}

genes = seq(1,length(expMatrix[,1]),1)

#setwd("C:\\Users\\Elizabeth\\Desktop\\Tesis-Redes\\Practica-R\\PracticaRTesis\\VarianceNeighbourhood\\rev3\\conShuffle")

png(filename="NeighbourhoodConShuffleAdjMatrixCompletaSISI.png", 
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=72)
#png(filename="pruebar4.png", res=300)
plot(genes,devSTD2[,1], pch="*", cex=2, col="black", xlab="Genes", ylab="", ylim = c(0,70), xlim = c(0,100),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
points(genes,devSTD2[,2],col="red", pch="o", cex=2)
# Adding a legend inside box at the location (2,40) in graph coordinates.
# Note that the order of plots are maintained in the vectors of attributes.
legend(1,70,legend=c("Varianza vecinos","std"), col=c("black","red"),pch=c("*","o"), cex=1.5,box.lty = 2)
text(80, 70, labels = "std = 3, kc=3, 80/100 genes", font=2, cex=1.5)
dev.off()



genesSi = 0
for (i in genes){
	if (devSTD2[,1][i] < devSTD2[,2][i]){
		genesSi = genesSi + 1
	}
}

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
		matrixFinal = rbind(matrixFinal,expMatrix[i,])
	}
}

#tengo que calcular la matriz de adyacencia devuelta..

samples = matrix(nrow = length(t), ncol = length(t)) #uso que es simetrica
for (i in seq(1,length(t),1)){
  for (j in seq(1,length(t),1)){
    samples[i,j] = samples[j,i] = distEuclideana(matrixFinal[,i],matrixFinal[,j])
  }
}

matrizRangos = t(apply(samples,1,rank))

adjMatrix = adyacencia(matrizRangos,4)

grafo = graph_from_adjacency_matrix(adjMatrix, mode = c("undirected"), weighted = NULL, diag = FALSE,
  add.colnames = NULL, add.rownames = NA)
#simplify(grafo, remove.multiple = F, remove.loops = T)
#plot(grafo)
#l = layout_on_sphere(grafo) 
png(filename="GrafoConShuffleStd3k3Quitado.png",
    units="in", 
    width=10, 
    height=8, 
    pointsize=12, 
    res=72)
plot(grafo, layout = layout_on_sphere(grafo), edge.color = "black", edge.arrow.size=2, edge.width = 2, edge.arrow.width = 2, edge.lty = 1, vertex.size = 10)
text(1, 1, labels = "std = 3, k=3",font=2, cex=1.5)
dev.off()

is_connected(grafo)


#le saco la diagonal 

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

#djMatrix = adjMatrix + t(adjMatrix)
#no hago eso para que sea un numerofijo para lo que sige la cantidadde vecinos

write.table(matrixFinal, file="M.txt", row.names=FALSE, col.names=FALSE)

abrir = read.table("M.txt",header=TRUE,sep=" ")
#prueba
MatrixA = matrix(nrow=2,ncol=2)
MatrixA[,1] = c(1,2)
MatrixA[,2] = c(2,3)

MatrixB =cbind(MatrixA, c(5,6))
###fin prueba


