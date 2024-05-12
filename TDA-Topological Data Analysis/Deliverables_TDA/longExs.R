library('TDA')

#### DELIVERY 4


X1 <- matrix(c(1.05, 2.00, 2.16, 1.59, 1.42, 2.66,
               2.50, 2.74, 2.07, 1.92, 1.06, 1.64, 
               1.19, 1.83, 1.55, 2.05, 1.76, 2.98, 
               2.89, 1.84, 2.16, 1.05, 2.05, 1.59, 
               1.76, 1.09, 2.21, 0.99, 1.90, 1.84, 
               2.35, 2.64, 2.69, 2.49, 2.15, 2.89, 
               1.67, 1.15, 1.95, 2.96, 1.84, 2.09, 
               2.27, 1.59, 1.16, 2.55, 1.59, 1.30,
               2.29, 1.97, 1.01, 1.53, 2.80, 2.24, 
               1.41, 1.96, 1.29, 2.16, 1.63, 1.08, 
               2.60, 1.91, 2.82, 0.99, 1.66, 2.29, 
               2.88, 2.63, 2.05, 2.52, 2.94, 2.19, 
               1.37, 2.46, 1.35, 2.45, 2.27, 1.14, 
               2.49, 2.77, 1.87, 1.69, 1.24, 2.60, 
               1.80, 1.15, 1.47, 2.25, 2.94, 2.29,
               1.52, 2.82, 1.81, 1.70, 2.51, 2.81, 
               1.62, 1.27, 2.63, 1.21, 1.68, 2.64, 
               2.57, 1.95, 1.17, 2.13, 1.96, 1.08,
               1.68, 1.25, 2.55, 2.21, 1.05, 2.40, 
               1.99, 2.42, 2.82, 2.22, 2.76, 2.53,
               1.31, 2.40, 1.44, 2.38, 2.57, 2.71,
               1.66, 2.58, 2.66, 1.79, 1.74, 2.88,
               2.91, 2.27, 2.29, 1.67, 3.00, 1.80,
               1.31, 1.45, 2.63, 2.62, 2.76, 2.47,
               1.03, 2.11, 1.74, 1.20, 1.73, 2.36
               ), ncol = 3, byrow = TRUE)
X1

## Compute the Vietoris-Rips persistence module:
Diag <- ripsDiag(X1, maxdimension = 3, maxscale = 5, library = 'GUDHI')

par(mfrow = c(1,2))
plot(X1)
plot(Diag[['diagram']], barcode = FALSE)
plot(Diag[['diagram']], barcode = TRUE)

Diag
#persistence silhouette
DiagLim = 5
tseq <- seq(0, DiagLim, length = 1000)
Sil <- silhouette(Diag[['diagram']], p = 1,  dimension = 3, tseq)
#Sil <- silhouette(Diag[['diagram']], p = 1,  dimension = 3)
Land <- landscape(Diag[['diagram']], dimension = 3, KK=2,tseq)
plot(Sil)
plot(Land, type = 'l')

par(mfrow = c(1,1))
plot(Diag[['diagram']], barcode = FALSE)

par(mfrow = c(1,1))
plot(silhouette(Diag[['diagram']], p = 1,  dimension = 3, tseq), type = 'l')




###Landscape:
## dimension <-  the dimension of the topological features under consideration. The default value is 1 (loops).
## KK <- a vector: the order of the landscape function. The default value is 1. (First Landscape function).

Diag

maxscale <- 1.5
tseq <- seq(0, maxscale, length = 1000)   #domain
Land <- landscape(Diag[["diagram"]], dimension = 1, KK = 1, tseq)
Sil <- silhouette(Diag[["diagram"]], p = 1, dimension = 1, tseq)
par(mfrow = c(1,2))
plot(tseq, Land, type = "l")
plot(tseq, Sil, type = "l")


par(mfrow = c(2,2))
for (i in 1:4){
  Land1 = landscape(Diag[['diagram']], dimension = 1, KK=i,tseq)
  plot(tseq, Land1, type= 'l')
}
par(mfrow = c(1,1))
par(mfrow = c(1,3))
for (i in 1:3){
  Land1 = landscape(Diag[['diagram']], dimension = 1, KK=i,tseq)
  plot(tseq, Land1, type= 'l')
}

maxscale <- 3
tseq <- seq(0, maxscale, length = 1000)   #domain
Land <- landscape(Diag[["diagram"]], dimension = 1, KK = 1, tseq)
Sil <- silhouette(Diag[["diagram"]], p = 1, dimension = 1, tseq)
par(mfrow = c(1,2))
plot(tseq, Land, type = "l")
plot(tseq, Sil, type = "l")

par(mfrow = c(2,3))
for (i in 1:5){
  Land1 = landscape(Diag[['diagram']], dimension = 1, KK=i,tseq)
  plot(tseq, Land1, type= 'l')
}

par(mfrow = c(2,2))
for (i in 1:4){
  Land1 = landscape(Diag[['diagram']], dimension = 1, KK=i,tseq)
  plot(tseq, Land1, type= 'l')
}


maxscale <- 5
tseq <- seq(0, maxscale, length = 1000)   #domain
Land <- landscape(Diag[["diagram"]], dimension = 2, KK = 1, tseq)
Sil <- silhouette(Diag[["diagram"]], p = 1, dimension = 2, tseq)
par(mfrow = c(1,2))
plot(tseq, Land, type = "l")
plot(tseq, Sil, type = "l")

par(mfrow = c(2,3))
for (i in 1:5){
  Land1 = landscape(Diag[['diagram']], dimension = 2, KK=i,tseq)
  plot(tseq, Land1, type= 'l')
}


par(mfrow = c(1,1))
plot(tseq, Sil, type = 'l')


##### DELIVERY 5

X1 <- matrix(
  c(0.83089090, 0.92139106, 1.05135118, 0.69658666,
  1.40570229, 1.45601583, 0.05738789, 1.58270914,
  1.83601599, 2.31975042, 0.47738787, 0.41686889,
  0.58929427, 1.24252248, 1.65529594, 0.98534007,
  0.82368247, 1.00514813, 1.49864810, 0.92041849,
  0.14943596, 0.59539077, 1.40450974, 1.63377801,
  0.81189552, 1.27512527, 2.08298134, 0.87369037,
  2.53595820, 1.29246603, 1.01228192, 0.15565445, 
  2.05421121, 2.05534683, 0.75427762, 0.75383665,
  -0.14872117, 0.02096184, 1.47630053, 0.87510234,
  -0.42530458, 0.49917853, 1.20266640, 0.79526767,
  1.23633497, 1.81306232, 0.42736194, 1.51483777, 
  0.86969203, 0.99283119, 2.11973316, 1.44747695, 
  1.37592899, 0.96259255, 0.99710074, 1.67116645, 
  0.61740303, 1.30389238, 1.81213878, 1.91489774), 
  ncol = 2, byrow = TRUE
)
X1

X2 <- matrix(c(
  1.40754745, 2.647683, 1.24523782, 2.820778, 
  1.55517672, 3.303610, 1.97133972, 2.321545, 
  1.39970730, 2.477310, 1.15300714, 3.251407,
  2.04172416, 4.218545, 2.00200034, 2.743398,
  1.82156736, 3.471864, 1.75719936, 3.097291,
  1.90579751, 3.350476, 2.33118508, 3.963780, 
  0.60439609, 2.430479, 1.98517772, 2.459581,
  2.81785877, 3.541370, -0.02755078, 1.812042,
  2.06339041, 3.167594, 0.75285330, 3.240254, 
  2.15114208, 2.660630, 2.28021529, 2.584600,
  1.86389493, 3.231760, 2.58430981, 3.650785,
  2.20839977, 3.496688, 3.53630008, 3.469574,
  2.13849006, 3.351805, 2.78733827, 3.396214,
  1.92692528, 2.975792, 3.52340370, 3.839808, 
  1.95290166, 3.071435, 0.99651700, 2.791795),
  ncol = 2, byrow = TRUE
)

X2 

X1_2 <- matrix(
  c(0.83089090, 0.92139106, 1.05135118, 0.69658666,
    1.40570229, 1.45601583, 0.05738789, 1.58270914,
    1.83601599, 2.31975042, 0.47738787, 0.41686889,
    0.58929427, 1.24252248, 1.65529594, 0.98534007,
    0.82368247, 1.00514813, 1.49864810, 0.92041849,
    0.14943596, 0.59539077, 1.40450974, 1.63377801,
    0.81189552, 1.27512527, 2.08298134, 0.87369037,
    2.53595820, 1.29246603, 1.01228192, 0.15565445, 
    2.05421121, 2.05534683, 0.75427762, 0.75383665,
    -0.14872117, 0.02096184, 1.47630053, 0.87510234,
    -0.42530458, 0.49917853, 1.20266640, 0.79526767,
    1.23633497, 1.81306232, 0.42736194, 1.51483777, 
    0.86969203, 0.99283119, 2.11973316, 1.44747695, 
    1.37592899, 0.96259255, 0.99710074, 1.67116645, 
    0.61740303, 1.30389238, 1.81213878, 1.91489774,
    2.29024248, 2.526838, 2.13586188, 2.792384, 
    2.64627954, 2.317721, 2.01157856, 2.221607,
    2.94471778, 3.495130, 1.59348125, 3.347204,
    2.50426440, 3.819559, 1.75890949, 2.685264, 
    0.93981079, 2.512890, 2.37635825, 2.883139),
  ncol = 2, byrow = TRUE
)
X1_2

X2_2 <- matrix(c(
  1.40754745, 2.647683, 1.24523782, 2.820778, 
  1.55517672, 3.303610, 1.97133972, 2.321545, 
  1.39970730, 2.477310, 1.15300714, 3.251407,
  2.04172416, 4.218545, 2.00200034, 2.743398,
  1.82156736, 3.471864, 1.75719936, 3.097291,
  1.90579751, 3.350476, 2.33118508, 3.963780, 
  0.60439609, 2.430479, 1.98517772, 2.459581,
  2.81785877, 3.541370, -0.02755078, 1.812042,
  2.06339041, 3.167594, 0.75285330, 3.240254, 
  2.15114208, 2.660630, 2.28021529, 2.584600,
  1.86389493, 3.231760, 2.58430981, 3.650785,
  2.20839977, 3.496688, 3.53630008, 3.469574,
  2.13849006, 3.351805, 2.78733827, 3.396214,
  1.92692528, 2.975792, 3.52340370, 3.839808, 
  1.95290166, 3.071435, 0.99651700, 2.791795, 
  2.29024248, 2.526838, 2.13586188, 2.792384, 
  2.64627954, 2.317721, 2.01157856, 2.221607,
  2.94471778, 3.495130, 1.59348125, 3.347204,
  2.50426440, 3.819559, 1.75890949, 2.685264, 
  0.93981079, 2.512890, 2.37635825, 2.883139),
  ncol = 2, byrow = TRUE
)
X2_2

## Compute the Vietoris-Rips persistence module:
#Diag <- ripsDiag(X1, maxdimension = 3, maxscale = 5, library = 'GUDHI')

par(mfrow = c(1,3))
plot(X1)
Diag <- ripsDiag(X1, maxdimension = 2, maxscale = 5, library = 'GUDHI')
plot(Diag[['diagram']], barcode = FALSE)
plot(Diag[['diagram']], barcode = TRUE)

par(mfrow = c(1,3))
plot(X1_2)
Diag <- ripsDiag(X1_2, maxdimension = 2, maxscale = 5, library = 'GUDHI')
plot(Diag[['diagram']], barcode = FALSE)
plot(Diag[['diagram']], barcode = TRUE)


par(mfrow = c(1,3))
plot(X2)
Diag <- ripsDiag(X2, maxdimension = 2, maxscale = 5, library = 'GUDHI')
plot(Diag[['diagram']], barcode = FALSE)
plot(Diag[['diagram']], barcode = TRUE)


par(mfrow = c(1,3))
plot(X2_2)
Diag <- ripsDiag(X2_2, maxdimension = 2, maxscale = 5, library = 'GUDHI')
plot(Diag[['diagram']], barcode = FALSE)
plot(Diag[['diagram']], barcode = TRUE)

distFct(X1, X1_2)
max(distFct(X1, X1_2))

distFct(X2, X2_2)
max(distFct(X2, X2_2))


