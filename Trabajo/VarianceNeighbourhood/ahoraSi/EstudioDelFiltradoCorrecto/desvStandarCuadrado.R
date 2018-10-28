#dada la matriz de adyacencia adjConectado 
# y la matriz de expresion matrixx


###defino para el filtrado las funciones correspondientes 
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