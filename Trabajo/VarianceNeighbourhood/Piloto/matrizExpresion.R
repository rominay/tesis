
matrizExpr = function(genesTotal, t, stdd, perc){

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

		f = matrix(nrow = genesTotal, ncol = length(t)) #va a ser mi matriz de expresion genetica

		grupos = 5
		genesGrupo = genesTotal / grupos
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

		matrizFSuave = matrizExpGenetica(curv,stdd,genesGrupo)



		##el shuffle

		random = function(grupos, genesTotal, matrix, perc){
			matrixOut = matrix(0L, nrow=genesTotal,ncol=(length(matrix[1,])+1))
			genesShuffled = vector()
  			for(i in 1:grupos){
				FilasRnd = sample((1:genesGrupo) + (i-1)*genesGrupo,perc*genesGrupo ) #(i-1)*genesGrupo es un offset al ir pasando
    				for(j in seq_along(FilasRnd)){
      					matrix[FilasRnd[j],1:length(matrix[1,])] = sample(matrix[FilasRnd[j],], ncol(matrix)) #hago un shuffle de esas filas rnd que agarre. 
    				}
				genesShuffled = c(genesShuffled, FilasRnd)
  			}
			matrixOut[,1:length(matrix[1,])] = matrix
  			matrixOut[1:length(genesShuffled),(length(matrix[1,])+1)] = genesShuffled 
			matrixOut
		}

		matrixOut = random(grupos, genesTotal, matrizFSuave, perc)
		matrixExprShuffled = matrixOut[,1:length(matrizFSuave[1,])]

		png(filename="ShuffleVsNo.png", 
    		units="in", 
    		width=10, 
    		height=8, 
    		pointsize=12, 
    		res=72)
		plot(matrizFSuave,matrixExprShuffled, pch="*", col="black", xlab="No shuffle", ylab="Shuffle",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, ylim = c(0,270), xlim = c(0,270))
		#text(100, 265, labels = "genesString genes, std = 3, shuffle 100%",font=2, cex=1.5)
		text(125, 265, paste0("genes=", genesTotal," grupos=", grupos, " std=", stdd, " shuffle=", perc),font=2, cex=1.5)
		dev.off()
		
		matrixOut 
}

####################################
#t = seq(0,80,1)      
#genesTotal = 500 
#matrixx = matrizExpr(genesTotal, t, 3, 0.1)

##################################


