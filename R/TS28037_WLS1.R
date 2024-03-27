## Basado en la ISO/TS 28037:2010

# Borrar todas las variables previas en R
rm(list = ls())

cat("Este documento ejecuta el ejemplo numérico de mínimos cuadrados ponderados (WLS)
 con *pesos iguales conocidos* descrito en la *Cláusula 6 (Modelo para
 incertidumbres asociadas a y_i)*, y realiza la predicción descrita en
 *11.1 EJEMPLO 1* y la evaluación hacia delante descrita en *11.2 EJEMPLO*.")


####################### DEFINICION DE FUNCIONES ###################################

write_wls_tableau <- function(x, y, w, g0, h0, g, h, a, b, r, fstr){
  m <- length(x)
  C1 <- c(0, 0, 0, 0, g0, h0, 0, 0, a, 0)
  C2 <- cbind(w, w*w, w*w*x, w*w*y, g, h, g*g, g*h, r, r*r)
  C3 <- colSums(C2)
  C3[1] <- 0
  C3[5] <- 0
  C3[6] <- 0
  C3[9] <- b
  Ca <- rbind(C1, C2, C3)
  cat("Cuadro de cálculo\n")
  print(round(Ca,4))

}


# Generar tabla de consulta.
quantile = c(3.841,5.991,7.815,9.488,11.070,12.592,14.067,15.507,16.919,18.307,
             19.675,21.026,22.362,23.685,24.996,26.296,27.587,28.869,30.144,31.410
             ,32.671,33.924,35.172,36.415,37.652,38.885,40.113,41.337,42.557,43.773
             ,44.985,46.194,47.400,48.602,49.802,50.998,52.192,53.384,54.572,55.758
             ,56.942,58.124,59.304,60.481,61.656,62.830,64.001,65.171,66.339,
             67.505,68.669,69.832,70.993,72.153,73.311,74.468,75.624,76.778,77.931
             ,79.082,80.232,81.381,82.529,83.675,84.821,85.965,87.108,88.250,89.391
             ,90.531,91.670,92.808,93.945,95.081,96.217,97.351,98.484,99.617,100.749
             ,101.879,103.010,104.139,105.267,106.395,107.522,108.648,109.773,110.898
             ,112.022,113.145,114.268,115.390,116.511,117.632,118.752,119.871,120.990
             ,122.108,123.225,124.342)
# length(quantile)

table <- as.data.frame(cbind(nu = c(1:100), quantile))

# Busqua el valor de nu en la primera columna de la tabla de búsqueda. Si 'nu' aparece en
# la tabla, establezca 'chi_sq' para que sea el cuantil correspondiente, de lo contrario use el
# aproximación de Wilson y Hilferty.

calc_chi_sq_95_percent_quantile <- function(nu) {
  table <- data.frame(nu = c(1, 2), quantile = c(3.841, 5.991))
  ind <- match(nu, table$nu)
  if (!is.na(ind)) {
    chi_sq <- table$quantile[ind]
  } else {
    chi_sq <- nu * ((1 - 2/(9*nu) + 1.644854*sqrt(2/(9*nu)))^3)
  }
  return(chi_sq)
}





################################################################################
################## MODELOS DE REGRESION PONDERADA #############################
###############################################################################





#### Asignar datos de medición ####
# # Asignar valores de x
 x <- c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
#
#
# # Asignar valores de y
 y <- c(3.3, 5.6, 7.1, 9.3, 10.7, 12.1)
#
# # Asignar incertidumbres asociadas a los valores de y
 uy <- c(0.5, 0.5, 0.5, 1, 1, 1)


cat(" El objetivo es obtener estimaciones de los parámetros de la función de calibración
 en línea recta y las incertidumbres y covarianzas estándar asociadas.
 Para esto se resuelve el problema de mínimos cuadrados ponderados para obtener
 los parámetros de la recta de mejor ajuste.")

# tamaño de los datos #
m <- length(x)


# Paso 1.
w <- 1/uy
F2 <- sum(w*w)

# Paso 2.
g0 <- sum(w * w * x) / F2
h0 <- sum(w * w * y) / F2

# Paso 3.
g <- w * (x - g0)
h <- w * (y - h0)

# Paso 4.
G2 <- sum(g * g)

# Paso 5.
b <- sum(g * h) / G2
a <- h0 - b * g0

# Paso 6.
u2a <- 1/F2 + g0^2/G2
u2b <- 1/G2
uab <- -g0/G2

# Paso 7.
r <- w * (y - a - b * x)

# Paso 8.
chi_sq_obs <- sum(r * r)
nu <- (m - 2)

# Paso 9.
chi_sq <- calc_chi_sq_95_percent_quantile(nu)

## Mostrar información en pantalla y generar figuras ####



cat("MODELO PARA LAS INCERTIDUMBRES ASOCIADAS AL Y
    ISO/TS 28037:2010(E) CLÁUSULA 6
    EJEMPLO (PESOS IGUALES) \n\n")

##
# Measurement data.

cat("Los datos que representan", m, "puntos de medida, con pesos iguales\n")
print(cbind(x,y,uy))

##
# Tabla de cálculo:
# Ver: write_wls_tableau.R
write_wls_tableau(x, y, w, g0, h0, g, h, a, b, r)

##
# Estimaciones de solución.

cat(sprintf("Estimación del intercepto : %8.3f\n", a))
cat(sprintf("Estimación de la pendiente: %8.3f\n", b))
cat(sprintf("Incertidumbre estándar asociada a la estimación del intercepto: %.3f\n", sqrt(u2a)))
cat(sprintf("Incertidumbre estándar asociada con la estimación de la pendiente: %.3f\n", sqrt(u2b)))
cat(sprintf("Covarianza asociada con las estimaciones del intercepto y la pendiente: %.3f\n", uab))

##
# Validación del modelo.
cat("VALIDACIÓN \n
    Grados de libertad \n")
cat(sprintf("%4f \n\n", nu))
cat("Valor chi-cuadrado observado \n")
cat(sprintf("%9.3f \n\n", chi_sq_obs))
cat("Cuantil del 95% de la distribución chi-cuadrado con", nu ," grados de libertad")
cat(sprintf("\n%9.3f \n\n", chi_sq))
if (chi_sq_obs > chi_sq) {
  cat("FALLÓ LA PRUEBA DE CHI CUADRADO - SE RECHAZA EL MODELO DE LÍNEA RECTA \n\n")
} else {
  cat("PRUEBA DE CHI CUADRADO APROBADA - MODELO ACEPTADO \n\n")
}


# Configuración de fuente y tamaño de letra
par(lwd=2, cex=1.2, font.lab=2)

# Gráfico de líneas con barras de error
plot(x, y, type="n", xlab=expression(italic(x)), ylab=expression(italic(y)),main="Modelo WLS")
arrows(x, y-uy, x, y+uy, length=0.05, angle=90, code=3, col="black")
lines(x, a+b*x, col="blue")

# Gráfico de residuales
plot(x, r, type="n", xlab=expression(italic(x)), ylab=expression(italic(r)),main="Residuales WLS")
for (i in 1:m) {
  segments(x[i], 0, x[i], r[i], col="black")
  points(x[i], r[i], pch=21, bg="white", cex=1.2)
}
segments(par("usr")[1], 0, par("usr")[2], 0, lty=2, col="blue")


### Estimación de nuevas observaciones

# Prediccíon
# Un valor y se desea  predecir el x.
yp <- 10.5
uyp <- 0.5
xp <- (yp - a)/b
ca <- -1/b
cb <- -(yp - a)/(b^2)
cy <- 1/b
u2x <- ca^2*u2a + cb^2*u2b + 2*ca*cb*uab + cy^2*uyp^2
cat("ISO/TS 28037:2010(E) 11.1
    EJEMPLO 1
    PREDICCIÓN
    Valor 'y' medido", y)
cat("Incertidumbre estándar asociada al valor 'y' medido", uy)
cat("Estimación de x", x)
cat("Incertidumbre estándar asociada al valor 'x' estimado", sqrt(u2x))
cat("Incertidumbre asociada con el estimador de a es:", ca)
cat("Incertidumbre asociada con el estimador de b es:", cb)



## Evaluación previa
# con un valor x se desea  predecir el y.
xp <- 3.5
uxp <- 0.2
yp <- a + b*xp
ca <- 1
cb <- xp
cx <- b
u2y <- ca^2*u2a + cb^2*u2b + 2*ca*cb*uab + cx^2*uxp^2
cat("ISO/TS 28037:2010(E) 11.2 \n")
cat("EJEMPLO
    EVALUACIÓN ADELANTADA \n")
cat("Valor x medido \n")
cat(sprintf("%8.3f \n\n", xp))
cat("Incertidumbre estándar asociada al valor x medido \n")
cat(sprintf("%8.3f \n\n", uxp))
cat("Estimación de y \n")
cat(sprintf("%8.3f \n\n", yp))
cat("Incertidumbre estándar asociada a la estimación de y \n")
cat(sprintf("%8.3f \n\n", sqrt(u2y)))
cat("Incertidumbre asociada con el estimador de a es: \n")
cat(sprintf("%8.3f \n\n", ca))
cat("Incertidumbre asociada con el estimador de b es: \n")
cat(sprintf("%8.3f \n\n", cb))






wls1 <- data.frame(Método = c("WLS"),
                   Intercept = 0, U_Intrcpt = 0, Slope = 0, U_Slope = 0, Cor.est=0)

wls1$Intercept[1] <-a # Intercept
wls1$U_Intrcpt[1] <- sqrt(u2a)   # U-Intercept
wls1$Slope[1]     <- b # Slope
wls1$U_Slope[1]   <- sqrt(u2b) # U-Slope
wls1$Cor.est[1]  <-(uab) #Correlación entra a y b


print(wls1)


