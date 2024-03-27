## Software compatible con ISO/TS 28037:2010(E)

# Introducción
# Este documento ejecuta el ejemplo numérico de regresión de Gauss-Markov (GMR)
# descrita en la *Cláusula 9 (Modelo para incertidumbres y covarianzas
# asociadas con y_i)*.

# Los usuarios de MATLAB pueden ejecutar directamente el código en el archivo M
# correspondiente para  obtener los resultados indicados en la norma ISO/TS 28037:2010(E)
# y también pueden modificar los datos y ejecutar el código con los nuevos datos.

# Para los usuarios que no tengan acceso a MATLAB, el software puede utilizarse como
# base para preparar implementaciones en otros lenguajes de programación.

# El software se suministra con un <NPL_TS_28037(2010)_MSC_L_10_001.pdf
# acuerdo de licencia de software> (REF: MSC/L/10/001) y el uso del
# software está sujeto a los términos establecidos en dicho acuerdo.
# Al ejecutar el código M, el usuario acepta los términos del acuerdo.

# Limpiamos entorno de trabajo
# rm(list=ls())

# Asignación de datos de medición

# Asignar valores de x
# x <- c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0)
#
#
# # Asignar valores de y
# y <- c(1.3, 4.1, 6.9, 7.5, 10.2, 12.0, 14.5, 17.1, 19.5, 21.0)
#
# # Asignar matriz de covarianza asociada a los valores de y
# Uy <- matrix(c(2,1,1,1,1,0,0,0,0,0,
#                1,2,1,1,1,0,0,0,0,0,
#                1,1,2,1,1,0,0,0,0,0,
#                1,1,1,2,1,0,0,0,0,0,
#                1,1,1,1,2,0,0,0,0,0,
#                0,0,0,0,0,5,4,4,4,4,
#                0,0,0,0,0,4,5,4,4,4,
#                0,0,0,0,0,4,4,5,4,4,
#                0,0,0,0,0,4,4,4,5,4,
#                0,0,0,0,0,4,4,4,4,5), nrow = 10, byrow = TRUE)


# Obtener estimaciones de los parámetros de la función de calibración
# en línea recta y las incertidumbres y covarianzas estándar asociadas.
# Resolver el problema de regresión de Gauss-Markov para obtener los
# parámetros de la recta de mejor ajuste.



  m <- length(x)
##
# Paso 1.
Ly <- chol(Uy, lower = TRUE)


##
# Paso 2.
e <- rep(1, m)
f <- solve(Ly, e)
g <- solve(Ly, x)
h <- solve(Ly, y)

##
# Paso 3.
F2 <- sum(f * f)

##
# Paso 4.
g0 <- sum(f*g)/F2
h0 <- sum(f*h)/F2

##
# Paso 5.
gt <- g - g0 * f
ht <- h - h0 * f

## Paso 6.
Gt2 <- sum(gt*gt)

##
# Paso 7.
b <- sum(gt*ht)/Gt2
a <- h0 - b*g0

##
# Paso 8.
u2a <- 1/F2 + g0^2/Gt2
u2b <- 1/Gt2
uab <- -g0/Gt2


##
# Paso 9.
r <- ht - b*gt

##
# Paso 10.
chi_sq_obs <- sum(r*r)
nu <- m - 2

##
# Paso 11.
chi_sq <- calc_chi_sq_95_percent_quantile(nu);

## Visualizar información en pantalla y generar cifras

## Modelo de medición.
cat('\nMODELO PARA LAS INCERTIDUMBRES Y COVARIANCIAS ASOCIADAS AL YI\n\n')
cat('ISO/TS 28037:2010(E) CLÁUSULA 9 \n')
cat('EJEMPLO \n\n')

## Datos de medición.
cat('AJUSTE\n')
cat(paste0('Datos que representan ', m,' puntos de medida  \n'))
for (i in 1:m){
  cat(sprintf('%8.1f %8.1f \n', x[i], y[i]))
}
cat('\n')

## Matriz de covarianza Uy
cat("Matriz de covarianza Uy\n")
print(Uy)

## Factor de Cholesky Ly de Uy
cat("Factor de Cholesky Ly de Uy\n")
print(Ly)


##
# Tabla de cálculo: véase write_gmr_tableaux.R
write_gmr_tableaux(f, g, h, g0, h0, gt, ht, a, b, r, '%9.4f')

##
# Estimaciones de la solución.
cat("Estimación del intercepto\n", format(a, digits=4), "\n\n")
cat("Estimación de la pendiente\n", format(b, digits=4), "\n\n")

##
# Incertidumbres estándar asociadas con la estimación de la solución.
cat("Incertidumbre estándar asociada con la estimación del intercepto\n")
cat(sprintf("%9.4f \n\n", sqrt(u2a)))
cat("Incertidumbre estándar asociada a la estimación de la pendiente \n")
cat(sprintf("%9.4f \n\n", sqrt(u2b)))
cat("Covarianza asociada con las estimaciones del intercepto y la pendiente\n")
cat(sprintf("%9.4f \n\n", uab))

##
# Validación del modelo.
cat("VALIDACIÓN\n")
cat("Grados de libertad\n")
cat(nu, "\n\n")
cat("Valor chi-cuadrado observado\n")
cat(sprintf("%9.3f", chi_sq_obs), "\n\n")
cat("Cuantil del 95% de la distribución chi-cuadrado con", nu,"grados de libertad", "\n")
  cat(sprintf("%9.3f", chi_sq), "\n\n")
if (chi_sq_obs > chi_sq) {
  cat("FALLÓ LA PRUEBA DE CHI CUADRADO - SE RECHAZA EL MODELO DE LÍNEA RECTA\n\n")
} else {
  cat("PRUEBA CUADRADA CHI APROBADA - MODELO DE LÍNEA RECTA ACEPTADO\n\n")
}

##
# Figuras
plot(x, y, pch='o', col='black', cex=1.5, lwd=1, xlab='x', ylab='y')
error <- sqrt(diag(Uy))
arrows(x, y-error, x, y+error, length=0.05, angle=90, code=3)
abline(a=a, b=b, col='blue', lwd=2)

par(lwd=2, cex.lab=1.2, font.lab=2)
axis1 <- range(x, na.rm = TRUE)
plot(axis1, c(0, 0), type = "n", xlab = xlabel, ylab = ylabel, main="Residuales del modelo",
     xlim=c(axis1[1], axis1[2]), ylim=c(min(r), max(r)))
for (i in 1:m) {
  segments(x[i], 0, x[i], r[i], lwd=1)
  points(x[i], r[i], pch='o', col='black', bg='white', cex=1.5)
}
abline(h=0, lty=2, col='blue')
segments(axis1[1], 0, axis1[2], 0, lty=2, col='blue')
points(x, r, pch='o', col='black', bg='white', cex=1.5)

# par(lwd=2, cex.lab=1.2, font.lab=2)
# plot(x, y, pch='o', col='black', cex=1.5, lwd=1, xlab='x', ylab='y')
# points(x, y, pch='o', col='black', cex=1.5)
# arrows(x, y-r, x, y+r, length=0.05, angle=90, code=3)
# abline(a=a, b=b, col='blue', lwd=2)


# Agradecimientos El trabajo aquí descrito ha contado con el apoyo de la Oficina Nacional
# de Medición del Departamento de Empresa, Innovación y Competencias del Reino Unido,
# como parte de su programa NMS Software Support for Metrology.





