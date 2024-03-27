## Software compatible con ISO/TS 28037:2010(E)

## Introducción
# Este documento ejecuta el ejemplo numérico de mínimos cuadrados ponderados (WLS)
# con *pesos iguales desconocidos* descrito en el *Anexo E (Incertidumbres conocidas
#  hasta un factor de escala)*.

# Los usuarios de MATLAB pueden ejecutar directamente el código en el archivo
# M correspondiente para obtener los resultados indicados en la norma ISO/TS 28037:2010(E)
# y también pueden modificar los datos y ejecutar el código con los nuevos datos.

# Para los usuarios que no tengan acceso a MATLAB, el software puede utilizarse
# como base para preparar implementaciones en otros lenguajes de programación.

# El software se suministra con un <NPL_TS_28037(2010)_MSC_L_10_001.pdf
# acuerdo de licencia de software> (REF: MSC/L/10/001) y el uso del
# software está sujeto a los términos establecidos en dicho acuerdo.
# Al ejecutar el código M, el usuario acepta los términos del acuerdo.


# Borrar todas las variables previas en R
# rm(list = ls())
#
# ## Asignar datos de medición
#
# ##
# # Asignar valores x.
# x <- c(1.000, 2.000, 3.000, 4.000, 5.000, 6.000)
#
#
# ## Asignar valores y.
# #
# y <- c(3.014, 5.225, 7.004, 9.061, 11.201, 12.762)
#
# ##
# # Asignar incertidumbres asociadas a los valores y.
# uy <- c(1, 1, 1, 1, 1, 1)



WLS3 <- function(x, y, uy=rep(1,length(x))){
  m <- length(x)
  #INSERTAR WARNING
  if(length(x)<=4)print("Estimación no válida para m menor a 4")

## Obtener estimaciones de los parámetros de la función de calibración en
# línea recta y las incertidumbres estándar y covarianza asociadas.
# Resolver el problema de mínimos cuadrados ponderados para obtener
# los parámetros de la recta de mejor ajuste.

##
# Paso 1.
w <- 1/uy
F2 <- sum(w*w)

##
# Paso 2.
g0 <- sum(w*w*x)/F2
h0 <- sum(w*w*y)/F2

##
# Paso 3.
g <- w*(x - g0)
h <- w*(y - h0)

##
# Paso 4.
G2 <- sum(g*g)

##
# Paso 5.
b <- sum(g*h)/G2
a <- h0 - b*g0

##
# Paso 6.
u2a <- 1/F2 + g0^2/G2
u2b <- 1/G2
uab <- -g0/G2

##
# Paso 7.
r <- w*(y - a - b*x)

##
# Paso 8.
chi_sq_obs <- sum(r*r)


## Visualizar información en pantalla y generar cifras
# Modelo de medición.
cat('\nMODELO PARA LAS INCERTIDUMBRES ASOCIADAS CON LOS YI \n\n')
cat('ISO/TS 28037:2010(E) ANEXO E \n')
cat('EJEMPLO (PESOS DESCONOCIDOS) \n\n')

##
# Datos de medición
cat("Ajuste\n")
cat(paste("Datos que representan", m, "puntos de medición, pesos iguales desconocidos \n"))
for (i in 1:m) {
  cat(sprintf("%.1f %.1f %.1f \n", x[i], y[i], uy[i]))
}
cat("\n")

##
# Tabla de cálculo: véase 'write_wls_tableau.R'
write_wls_tableau(x, y, w, g0, h0, g, h, a, b, r, '%.3f ');

##
# Solución de estimaciones.
cat('Estimación del intercepto \n', a, '\n\n')
cat('Estimación de la pendiente \n', b, '\n\n')

##
# Incertidumbres estándar asociadas con las estimaciones de la solución

# Incertidumbre estándar asociada a la estimación del intercepto
cat("Incertidumbre estándar asociada a la estimación del intercepto\n")
cat(sprintf("%8.3f \n\n", sqrt(u2a)))
# Incertidumbre estándar asociada a la estimación de la pendiente
cat("Incertidumbre estándar asociada a la estimación de la pendiente \n")
cat(sprintf("%8.3f \n\n", sqrt(u2b)))
# Covarianza asociada a las estimaciones del intercepto y pendiente
cat("Covarianza asociada a las estimaciones del intercepto y pendiente \n")
cat(sprintf("%8.3f \n\n", uab))

##
# Figuras
par(lwd = 2, cex.lab = 1.2, font.lab = 2, font.main = 2)
plot(x, y, pch = 20, col = "black", ylim = c(0, max(y+uy)), xlim = c(0, max(x)+1),
     xlab = expression(x), ylab = expression(y))
segments(x, y - uy, x, y + uy, col = "black")
abline(a, b, col = "blue")
legend("topright", c("Data", "Fit"), col = c("black", "blue"), pch = c(20, NA), lty = c(NA, 1))

plot(x, r, type = "n", ylim = c(-max(abs(r)), max(abs(r))), xlim = c(0, max(x)+1),
     xlab = expression(x), ylab = expression(r))
for (i in 1:m){
  segments(x[i], 0, x[i], r[i], col = "black")
  points(x[i], r[i], pch = 20, col = "white", cex = 1.2)
}
abline(h = 0, lty = 2, col = "blue")

##
# Estimaciones posteriores de las incertidumbres estándar
sigmah2 <- chi_sq_obs/(m - 2)
sigmah <- sqrt(sigmah2)
u2ah <- sigmah2*u2a
u2bh <- sigmah2*u2b
uabh <- sigmah2*uab
cat("Estimación posterior de las incertidumbres estándar asociadas a los valores y \n")
cat(sprintf("%8.3f \n\n", sigmah))
cat("Incertidumbre estándar escalada asociada con la estimación del intercepto \n")
cat(sprintf("%8.3f \n\n", sqrt(u2ah)))
cat("Incertidumbre estándar escalada asociada a la estimación de la pendiente \n")
cat(sprintf("%8.3f \n\n", sqrt(u2bh)))
cat("Covarianza escalada asociada con las estimaciones del intercepto y la pendiente \n")
cat(sprintf("%8.3f \n\n", uabh))

## Mejores estimaciones posteriores de la incertidumbre
sigmat2 <- chi_sq_obs/(m - 4)
sigmat <- sqrt(sigmat2)
u2at <- sigmat2*u2a
u2bt <- sigmat2*u2b
uabt <- sigmat2*uab
cat("Mejor estimación posterior de las incertidumbres estándar asociadas con los valores y\n", sigmat, "\n\n")
cat("Mejor incertidumbre estándar asociada con la estimación de la intercepción\n", sqrt(u2at), "\n\n")
cat("Mejor incertidumbre estándar asociada con la estimación de la pendiente\n", sqrt(u2bt), "\n\n")
cat("Mejor covarianza asociada con las estimaciones del intercepto y la pendiente\n", uabt, "\n\n")
}
## Agradecimientos
# El trabajo aquí descrito ha contado con el apoyo de la Oficina Nacional
# de Medición del Departamento de Empresa, Innovación y Cualificaciones
# del Reino Unido como parte de su programa  NMS Software Support for Metrology.


