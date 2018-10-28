####################IGRAPH##############

Graf1 = graph(edges = c(1,2, 2,3, 1,3), n = 3, directed = T)
plot(Graf1)

###########para una red dirigida importan el orden en que escribo en el vector de las edges
###########1,3 implica flecha de 1 a 3 

#########grafo con nombre en los nodos 

Graf2 = graph( c("Romi","Mai", "Mai","Solcha", "Solcha","Yam", "Yam","Romi", "Yam", "Mai", "Solcha", "Romi"))
plot(Graf2)

Graf3 = graph( c("Romi","Mai", "Mai","Solcha", "Solcha","Yam", "Mai","Yam"),
               isolates = c("Juan", "Nico", "Fede", "Nacho") )

##########forma copada de plotear##############

plot(Graf3, edge.arrow.size = 0.5, vertex.color = 'gold', vertex.size = 15, vertex.frame.color = 'gray',vertex.label.color = 'black',vertex.label.cex = 0.8,vertex.label.dist = 1,edge.curved = 0.2)

#########ANALISIS#############################

#enlaces 
E(Graf3)                                  #hay 3 enlaces y te muestra de que nodo a que nodo es

#nodos 
V(Graf3)

#la matriz
Graf3[]

#le puedo pedir los nombres de los nodos 
V(Graf3)$name

#le puedo agregar un campo escalar como un atributo a mi grafo 
#sobre los nodos
V(Graf3)$gender = c("female", "female", "female", "male", "male", "male", "male", "male")
vertex_attr(Graf3)

#sobre los enlaces 
E(Graf3)$type = "email"
E(Graf3)$weight = c(2,3,4,5)
edge_attr(Graf3)

############como ploteaar con el campo impuesto

plot(Graf3, edge.arrow.size = 5, vertex.label.color = "black", vertex.label.dist = 1.5, vertex.color = c("pink", "skyblue")[1+(V(Graf3)$gender=="male")])