##  Software compatible con ISO/TS 28037:2010(E)

## Introducción

# Este documento ejecuta un ejemplo numérico para ilustrar el algoritmo de
# regresión de distancia generalizada (GDR) descrito en la *cláusula
# 8 (Modelo para  incertidumbres asociadas a las x_i y las y_i y
# covarianzas asociadas a los pares (x_i, y_i))*.
#


# Para cerrar todas las ventanas de gráficos abiertas
dev.off()
# Para limpiar todas las variables en el entorno de trabajo
rm(list = ls())
# Para establecer el formato de salida numérica
options(digits = 3)

## Asignar datos de medición

##

x <- c(1.2, 1.9, 2.9, 4.0, 4.7, 5.9)
m <- length(x)

##
# Asignar valores y.
y <- c(3.4, 4.4, 7.2, 8.5, 10.8, 13.5)

##
# Asignar incertidumbres asociadas a los valores x.
ux <- c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2)

##
# Asignar incertidumbres asociadas a los valores y.
uy <- c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4)

##
# Asignar covarianzas asociadas a pares de valores (x, y).
uxy <- c(0.036, 0.036, 0.036, 0.036, 0.036, 0.036)

# Obtener estimaciones de los parámetros de la función de calibración de línea
# recta y sus incertidumbres estándar asociadas y covarianza.
# Resolver el problema de regresión de distancia generalizada para
# obtener los mejores parámetros de ajuste de línea recta.

# Paso 1. Aproximación inicial utilizando mínimos cuadrados ponderados:
# Llamada a la función algm_wls_steps_1_to_8 con los argumentos x, y, uy

GDR2<- function(x,ux,y,uy,uxy){

resultados <- algm_wls_steps_1_to_8(x, y, uy)

# Extracción de los resultados individuales
ai <- resultados$ai
bi <- resultados$bi
u2ai <- resultados$u2ai
u2bi <- resultados$u2bi
uabi <- resultados$uabi
wi <- resultados$wi
g0i <- resultados$g0i
h0i <- resultados$h0i
gi <- resultados$gi
hi <- resultados$hi
G2i <- resultados$G2i
ri <- resultados$ri
Ri <- resultados$Ri

at <- list() # Inicializar lista para almacenar valores de 'ai'
bt <- list() # Inicializar lista para almacenar valores de 'bi'

ai <- round(10000 * ai) / 10000
at <- list(ai)
bi <- round(10000 * bi) / 10000
bt <- list(bi)

##
# Repita los pasos 2 a 5 hasta que se cumplan los criterios de convergencia.
# (En este ejemplo, se considera que se ha alcanzado la convergencia cuando
# las magnitudes de los incrementos deltaa y deltab en a y b no son
# superiores a 0,00005.) Para los datos generales del usuario, se sugiere
# que se considere que se ha alcanzado la convergencia cuando las magnitudes
# de todos los incrementos *relativos a las aproximaciones iniciales de los
# parámetros de la recta* no sean superiores a 1e-7. En este caso, se considera
# que se ha alcanzado la convergencia cuando las magnitudes de los incrementos
# deltaa y deltab en a y b no sean superiores a 0,00005. En este caso, la
# tolerancia puede asignarse utilizando el comando tol <- 1e-7 * norm(c(ai, bi))

# Paso 2: Inicialización de variables y asignar tolerancias e inicializar variables
  # Asignar tolerancias e inicializar variables.

tol=0.00005

# Crear las listas vacías da, db, t, xs, z, f, g, h, F2, g0, h0, gt, ht, Gt2 y r
da <- list()
db <- list()
t <- list()
xs <- list()
z <- list()
f <- list()
g <- list()
h <- list()
F2 <- list()
g0 <- list()
h0 <- list()
gt <- list()
ht <- list()
Gt2 <- list()
r <- list()

# Paso 2 a 5: Ejecución iterativa

ind<-1

# Llamada a la función 'algm_gdr1_steps_2_to_5' con los argumentos correspondientes
algm_gdr1_steps_2_to_5(x, ux, y, uy, at, bt, da, db, t, xs, z, f, g, h, F2, g0, h0, gt, ht, Gt2, r, ind)

# Bucle while para iterar hasta que la condición sea falsa
while ((abs(da[[ind]]) > tol) || (abs(db[[ind]]) > tol)) {
  # Actualizar el número de iteración
  ind <- ind + 1
  # Paso 6: Bucle while para iterar hasta que se alcance la convergencia
  # Llamada a la función algm_gdr1_steps_2_to_5 con los argumentos correspondientes
  algm_gdr1_steps_2_to_5(x, ux, y, uy, at, bt, da, db, t, xs, z, f, g, h, F2, g0, h0, gt, ht, Gt2, r, ind)
}

# Obtener valores finales de a y b
a <- at[[ind + 1]]
b <- bt[[ind + 1]]

# Paso 7: Evaluar las incertidumbres
u2a <- 1 /F2[[ind]] + g0[[ind]]^2/Gt2[[ind]]
u2b <- 1 /Gt2[[ind]]
uab <- -g0[[ind]]/Gt2[[ind]]

# Paso 8: Calcular el valor observado del chi-cuadrado y los grados de libertad
chi_sq_obs <- sum(r[[ind]] * r[[ind]])
nu <- (m - 2)

##
# Paso 9: Calcular el percentil del 95% de la distribución chi-cuadrado
# Ver: calc_chi_sq_95_percent_quantile.R.
chi_sq <- calc_chi_sq_95_percent_quantile(nu);

## Mostrar información en pantalla y generar cifras
# Modelo de medición.
cat("\nMODELO PARA LAS INCERTIDUMBRES ASOCIADAS A LAS XI Y LAS YI Y LAS COVARIANZAS ASOCIADAS A LOS PARES (XI, YI)\n\n")

##
# Datos de medición.
cat("AJUSTE \n")
cat(paste("Datos que representan", m, "puntos de medición\n"))
for (i in 1:m) {
  cat(sprintf("%8.1f %8.1f %8.1f %8.1f %8.3f\n", x[i], ux[i], y[i], uy[i], uxy[i]))
}
cat("\n")

##
# Tabla de cálculo para el problema WLS (mínimos cuadrados ponderados).
cat("Aproximaciones iniciales a los parámetros\n")
write_wls_tableau(x, y, wi, g0i, h0i, gi, hi, ai, bi, ri, '%9.4f ')

##
# Aproximación inicial a la intercepción
cat("Aproximación inicial a la intercepción\n")
cat(sprintf("%9.4f \n\n", ai))

# Aproximación inicial a la pendiente
cat("Aproximación inicial a la pendiente\n")
cat(sprintf("%9.4f \n\n", bi))



##
# Tabla de cálculo para cada iteración: ver write_gdr2_tableaux.R
# El tablero es similar al del problema de regresión por distancia
# generalizada (GDR) descrito en la *cláusula 7*, pero incluye una
# columna adicional (columna 5) en el primer tablero de cada iteración
# que contiene los valores de las covarianzas uxy).
niter<-1
# Bucle para cada iteración
for (niter in 1:ind) {
  write_gdr2_tableaux(x, ux, y, uy, uxy, t, xs, z, f, g, h, g0, h0, gt, ht,
                      at, bt, da, db, r, niter, '%9.4f')
}

##
# Solucioes estimadas:
cat('Estimación de la intersección:\n')
cat(sprintf('%9.4f\n', a), '\n')

cat('Estimación de la pendiente:\n')
cat(sprintf('%9.4f\n', b), '\n')

# Incertidumbres estándar asociadas a las estimaciones de las soluciones.
cat('Incertidumbre estándar asociada con la estimación del intercepto:\n')
cat(sprintf('%9.4f\n', sqrt(u2a)), '\n')

cat('Incertidumbre estándar asociada con la estimación de la pendiente:\n')
cat(sprintf('%9.4f\n', sqrt(u2b)), '\n')

cat('Covarianza asociada con las estimaciones de la intersección y la pendiente:\n')
cat(sprintf('%9.4f\n', uab), '\n')

##
cat('Grados de libertad:\n')
cat(sprintf('%9.4f\n', nu), '\n')

cat('Valor observado del chi-cuadrado:\n')
cat(sprintf('%9.3f\n', chi_sq_obs), '\n')

cat('Cuantil del 95% de la distribución del chi-cuadrado con', nu,'grados de libertad\n')
cat(sprintf('%9.3f\n', chi_sq), '\n')

# Prueba de chi-cuadrado
if (chi_sq_obs > chi_sq) {
  cat('FALLA LA PRUEBA CHI-CUADRADO: SE RECHAZA EL MODELO RECTILÍNEO\n\n')
} else {
  cat('PRUEBA DE CHI-CUADRADO SUPERADA: SE ACEPTA EL MODELO LINEAL\n\n')
}

##
# Cargar librerías necesarias (si es necesario)
# install.packages("ggplot2") #  si ggplot2 no está instalado
library(ggplot2)

# Configurar parámetros de visualización
options(scipen = 999) # Evitar notación científica en los ejes

# Figura 1

xlabel <- expression(italic(x)) # Etiqueta del eje x en itálica
ylabel <- expression(italic(y)) # Etiqueta del eje y en itálic

plot_errorbar2(x, ux, y, uy,  marker=1,  linestyle=1,main="Función de calibración con incertidumbre asociada",
               xlab=xlabel,ylab=ylabel) # Puntos con barras de error
lines(xs[[ind]], at[[ind]] + bt[[ind]] * xs[[ind]], col = "red", pch = "o") # Línea de ajuste
lines(xs[[ind]], at[[ind]] + bt[[ind]] * xs[[ind]], col = "blue", lwd = 2) # Línea de ajuste en azul

# Figura 2
axis1 <- range(x, na.rm = TRUE)
xlabel <- expression(italic(x)) # Etiqueta del eje x en itálica
ylabel <- expression(italic(r)) # Etiqueta del eje y en itálica

plot(NULL, xlim = c(axis1[1], axis1[2]), type = "n",ylim = c(min(r[[ind]]), max(r[[ind]])),
xlab = xlabel, ylab = ylabel, main = "Residuales del modelo") # Crear el lienzo de la figura
for (i in 1:m) {
  lines(c(x[i], x[i]), c(0, r[[ind]][i]), col = "black", lty = "solid") # Líneas verticales
  points(x[i], r[[ind]][i], pch = "o", col = "white", bg = "white",cex = 0.6) # Puntos con círculos blancos
}
abline(h = 0, lty = 2, col = "blue") # Línea horizontal en azul
}

## Agradecimientos
# El trabajo aquí descrito ha contado con el apoyo de la Oficina Nacional de Medición
# del Departamento de Empresa, Innovación y Cualificaciones del Reino Unido como parte
# de su programa NMS Software Support for Metrology.
