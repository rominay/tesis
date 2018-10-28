#############LISTAS#################

v1 = c(1,2,3)
v2 = 1:10
v3 = 1 
v4 = "holas"
l1 = list(boo = v1, foo = v2, moo = v3, hola = v4)
###dos formas de acceder a sus elementos
typeof(l1["boo"])                                     #devuelve como lista 
l1[["boo"]]                                           #solo devuelve elementos como vector 
##########DATA FRAMES###############


DataFr1 = data.frame( ID = 1:4, 
                      FirstName = c("John", "Jim", "Jane", "Jill"),
			    Female = c(F,F,T,T),
                      Age = c(22,33,44,55) )

###puedo acceder a cada una, que es un factor 

DataFr1$FirstName                                   # y me dice el level 
#lo puedo pedir como un vector 
as.vector(DataFr1$FirstName)

#algo super util es poder filtrar de acuerdo a alguna condicion sobre alguna variable 

DataFr1[DataFr1$Age > 30,2]   #el 2 me dice que quiero los nombres (el 2do factor del data frame)
