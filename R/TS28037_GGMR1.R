## Software compatible con ISO/TS 28037:2010(E)
#
# Introducción
# Este documento ejecuta el ejemplo numérico de regresión de Gauss Markov
# generalizada (GGMR) descrito en la *Cláusula 10 (Modelo para incertidumbres
# y covarianzas asociadas a las x_i y las y_i)*.
#
# Los usuarios de MATLAB pueden ejecutar directamente el código en el
# archivo M correspondiente para obtener los resultados indicados en
# la norma ISO/TS 28037:2010(E) y también pueden modificar los datos
# y ejecutar el código con los nuevos datos.
#
# Para los usuarios que no tengan acceso a MATLAB, el software puede
# utilizarse como base para preparar implementaciones en otros
# lenguajes de programación.
#
# El software se suministra con 'NPL_TS_28037(2010)_MSC_L_10_001.pdf'
# acuerdo de licencia de software> (REF: MSC/L/10/001) y el uso del
# software está sujeto a los términos establecidos en dicho acuerdo.
# Al ejecutar  el código M, el usuario acepta los términos del acuerdo.


# Cerrar todas las ventanas de gráficos abiertas
graphics.off()

# Limpiar el entorno de trabajo
# rm(list = ls())

# Establecer el formato de salida numérica a corto
options(digits = 4)


## Asignar datos de medición
##
# Asignar valores de x
# x <- c(50.4, 99.0, 149.9, 200.4, 248.5, 299.7, 349.1)
#
#
# ##
# # Asignar valores de y
# y <- c(52.3, 97.8, 149.7, 200.1, 250.4, 300.9, 349.2)
#
# ##
# # Asignar matriz de covarianza asociada a los valores de x
# Ux <- matrix(c(0.5, 0, 0.25, 0, 0.25, 0, 0.25,
#                0, 1.25, 1, 0, 0, 1, 1,
#                0.25, 1, 1.5, 0, 0.25, 1, 1.25,
#                0, 0, 0, 1.25, 1, 1, 1,
#                0.25, 0, 0.25, 1, 1.5, 1, 1.25,
#                0, 1, 1, 1, 1, 2.25, 2,
#                0.25, 1, 1.25, 1, 1.25, 2, 2.5), nrow = 7, byrow = TRUE)
#
#
#
# # Asignar matriz de covarianza asociada a los valores de y
# Uy <- matrix(c(5, 1, 1, 1, 1, 1, 1,
#                1, 5, 1, 1, 1, 1, 1,
#                1, 1, 5, 1, 1, 1, 1,
#                1, 1, 1, 5, 1, 1, 1,
#                1, 1, 1, 1, 5, 1, 1,
#                1, 1, 1, 1, 1, 5, 1,
#                1, 1, 1, 1, 1, 1, 5), nrow = 7, byrow = TRUE)

##
# Construya una matriz de covarianza completa para los valores x e y
# (suponiendo que la correlación asociada a cada par (x e y) es cero).


GGMR<- function(x,Ux,y,Uy){

m <- length(x)

U <- diagU <- diag(2*m)
U[1:m, 1:m] <- Ux
U[(m+1):(2*m), (m+1):(2*m)] <- Uy


# Obtener estimaciones de los parámetros de la función de calibración
# en línea recta y las incertidumbres y covarianzas estándar asociadas.
#  Resolver el problema de regresión de Gauss-Markov generalizado
# para obtener los parámetros de la recta de mejor ajuste.


##
# Paso 1. Aproximación inicial usando mínimos cuadrados ponderados:
# ver algm_wls_steps_1_to_8.R


# Llamar a la función algm_wls_steps_1_to_8 para obtener los valores ai y bi
result <- algm_wls_steps_1_to_8(x, y, sqrt(diag(Uy)))
ai <- result$ai
bi <- result$bi




# Redondear las aproximaciones iniciales de los parámetros a cuatro decimales
# (Este paso se incluye para producir los resultados dados en ISO/TS 28037:2010:
# normalmente no se realiza este paso.)
ai <- round(ai, 4)
at <- list(ai)
bi <- round(bi, 4)
bt <- list(bi)

##
# Aproximación inicial
tt <- list()
tt[[1]] <- c(x, at[[1]], bt[[1]])

# Bucle a través de los pasos 2 a 9 hasta que se cumplan los criterios de convergencia
# (En este ejemplo, se considera que se ha alcanzado la convergencia cuando las
# magnitudes de todos los incrementos no son mayores que 1e-7.
# Para datos de usuario general, se sugiere que se considere que se ha alcanzado
# la convergencia cuando las magnitudes de todos los incrementos
# *relativos a las aproximaciones iniciales de los parámetros de la recta* no son mayores que 1e-7.
# En este caso, la tolerancia se puede asignar usando el comando tol = 1e-7*norm(c(ai, bi));)

# Asignar tolerancias e inicializar variables
tol <- 1e-7
dt <- list()
f <- list()
J <- list()
ft <- list()
Jt <- list()
g <- list()
H <- list()
M <- list()
q <- list()




##
# Pasos 2 a 9 ver :algm_ggmr_cholesky_steps_2_to_9.m
#

ind <- 1

algm_ggmr_cholesky_steps_2_to_9(x, y, U, tt, dt, f, J, ft, Jt, g, H, M, q, ind)


while (any(abs(dt[[ind]]) > tol)) {
    # Actualiza el número de iteración.
    ind <- ind + 1
  algm_ggmr_cholesky_steps_2_to_9(x, y, U, tt, dt, f, J, ft, Jt, g, H, M, q, ind)
}



##
# Paso 10. Repita los pasos 2 a 9 hasta alcanzar la convergencia: véase
# 'algm_ggmr_cholesky_steps_2_to_9.R'

# Calcular los valores de a y b en la última iteración
a <- tt[[ind + 1]][m + 1]
b <- tt[[ind + 1]][m + 2]


##
# Paso 11. Particionar el factor de Cholesky y evaluar las incertidumbres
M22 <- M[[ind]][(m + 1):(m + 2), (m + 1):(m + 2)]
m11 <- M22[1, 1]
m21 <- M22[2, 1]
m22 <- M22[2, 2]
u2a <- (m22^2 + m21^2)/(m11^2 * m22^2)
u2b <- 1/(m22^2)
uab <- -m21/(m11 * m22^2)

##
# Paso 12. Formar el valor observado del chi-cuadrado y los grados de libertad
chi_sq_obs <- sum(ft[[ind]] * ft[[ind]])
nu <- (m - 2)

##
# Paso 13. Calcular el percentil del 95% de la distribución chi-cuadrado
# Ver: 'calc_chi_sq_95_percent_quantile.R
chi_sq <- calc_chi_sq_95_percent_quantile(nu)

## Mostrar información en pantalla

##
# Modelo de medición
cat("\nMODELO DE INCERTIDUMBRES Y COVARIANZAS ASOCIADAS A LAS XI Y LAS YI\n\n")
cat("ISO/TS 28037:2010(E) CLÁUSULA 10\n")
cat("EJEMPLO\n\n")

##
# Datos de medición.
cat("AJUSTE\n")
cat("Datos que representan",num2str(m), "puntos de medición\n")
for (i in 1:m) {
  cat(sprintf("%8.1f %8.1f\n", x[i], y[i]))
}
cat("\n")

##
# Matrix de varianzas-covarianzas: 'Ux'
cat("Covariance matrix Ux\n")
print(Ux)

##
# Factor Cholesky 'Lx' de 'Ux'
cat("Factor Cholesky Lx de Ux\n")
print(chol(Ux, lower = TRUE))


##
# Matrix de varianzas-covarianzas: 'Uy'
cat("Covariance matrix Uy\n")
print(Uy)

##
# Factor Cholesky 'Ly' de 'Uy'
cat("Factor Cholesky Ly de Uy\n")
print(chol(Uy, lower = TRUE))

##
# Aproximaciones iniciales a los parámetros
cat("Aproximación inicial a la intercepción\n")
print(ai)
cat("\n")
cat("Aproximación inicial a la pendiente\n")
print(bi)
cat("\n")

##
# Resumen del procedimiento iterativo: ver
# 'write_ggmr_tableau.m'
# Crear una matriz con los datos de tt y dt
write_ggmr_tableau(tt, dt, ind, '%14.4e')



# Escribir el tableau en la consola:
cat("Summary of iterative procedure:\n")
cat(sprintf("%14.4e ", tabla), "\n")


##
# Matriz M en iteración final
cat("Matriz M en iteración final:\n")
print(M[[ind]])

# Matriz M22 en iteración final
cat("Matriz M22 en iteración final:\n")
print(M[[ind]][(m+1):(m+2), (m+1):(m+2)])



##
# Soluciones estimadas.
# Estimación del intercepto
cat("Estimación del intercepto:\n")
cat(sprintf("%9.4f\n", a), "\n")

# Estimación de la pendiente
cat("Estimación de la pendiente:\n")
cat(sprintf("%9.4f\n", b), "\n")

##
# Incertidumbres estándar asociadas a las estimaciones de las soluciones.
# Incertidumbre estándar asociada a la estimación del intercepto
cat("Incertidumbre estándar asociada a la estimación del intercepto:\n")
cat(sprintf("%9.4f\n", sqrt(u2a)), "\n")

# Incertidumbre estándar asociada a la estimación de la pendiente
cat("Incertidumbre estándar asociada a la estimación de la pendiente:\n")
cat(sprintf("%9.4f\n", sqrt(u2b)), "\n")

# Covarianza asociada a las estimaciones del intercepto y la pendiente
cat("Covarianza asociada a las estimaciones del intercepto y la pendiente:\n")
cat(sprintf("%9.4f\n", uab), "\n")

##
# Validación del modelo.
cat("VALIDACION:\n")
# Grados de libertad
cat("Grados de libertad:\n")
cat(sprintf("%9.1f", nu), "\n")

# Valor observado del estadístico chi-cuadrado
cat("Valor chi-cuadrado observado:\n")
cat(sprintf("%9.3f", chi_sq_obs),"\n")

# Cuantil del 95% de la distribución chi-cuadrado con grados de libertad nu
cat("Cuantil del 95% de la distribución chi-cuadrado con", nu, "grados de libertad:\n")
cat(sprintf("%9.3f", chi_sq), "\n")

# Resultado de la prueba chi-cuadrado
if (chi_sq_obs > chi_sq) {
  cat("FALLÓ LA PRUEBA DE CHI CUADRADO - SE RECHAZA EL MODELO DE LÍNEA RECTA")
} else {
  cat("LA PRUEBA DE CHI CUADRADO PASÓ - SE ACEPTA EL MODELO DE LÍNEA RECTA")
}

# Agradecimientos
# El trabajo aquí descrito ha contado con el apoyo de la Oficina Nacional de
# Medición del Departamento de Empresa, Innovación y Cualificaciones del Reino
# Unido como parte de su programa  NMS Software Support for Metrology.






xlabel <- expression(italic(x)) # Etiqueta del eje x en itálica
ylabel <- expression(italic(y)) # Etiqueta del eje y en itálic

plot_errorbar2(x, diag(Ux), y, diag(Uy),  marker=1,  linestyle=1,main="Función de calibración con incertidumbre asociada",
               xlab=xlabel,ylab=ylabel) # Puntos con barras de error
lines(x, a+b*x, col="blue")


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
