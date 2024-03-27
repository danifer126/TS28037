## Basado en la ISO/TS 28037:2010

# Borrar todas las variables previas en R
rm(list = ls())

cat("Regresión de distancia generalizada (GDR) descrito en la *Cláusula 7
 Software compatible con ISO/TS 28037:2010(E)
 Este documento ejecuta el ejemplo numérico de regresión de distancia
 generalizada (GDR) descrito en la *Cláusula 7 (Modelo para incertidumbres
 asociadas con la x_i y la y_i)*")


####################### DEFINICION DE FUNCIONES ###################################


### primera función es el algoritmo de regresion ponderda.
algm_wls_steps_1_to_8 <- function(x, y, uy){
  # Step 1
  m <- length(x)
  w <- 1/uy
  F2 <- sum(w * w)

  # Step 2
  g0 <- sum(w * w * x) / F2
  h0 <- sum(w * w * y) / F2

  # Step 3
  g <- w * (x - g0)
  h <- w * (y - h0)

  # Step 4
  G2 <- sum(g * g)

  # Step 5
  b <- sum(g * h) / G2
  a <- h0 - b * g0

  # Step 6
  u2a <- 1/F2 + (g0^2) / G2
  u2b <- 1 / G2
  uab <- -g0 / G2

  # Step 7
  r <- w * (y - a - b * x)

  # Step 8
  R <- sum(r * r)

  return(list(ai = a, bi = b, u2ai = u2a, u2bi = u2b, uabi = uab, wi= w,
              g0i = g0, h0i = h0, gi = g, hi = h, G2i = G2, ri = r, Ri = R))
}



### Segunda funcion para llevar a cabo los pasos 2 a 5 del algoritmo GDR (Cláusula 7)
algm_gdr1_steps_2_to_5 <- function(x, ux, y, uy, at, bt, da, db,t, xs, z, f,
                                   g, h, F2, g0, h0, gt, ht, Gt2, r, ind){



  # Paso 2.
  m <- length(x)
  t[[ind]] <- 1 / (uy^2 + bt[[ind]]^2 * ux^2)
  xs[[ind]] <- (x * uy^2 + bt[[ind]] * ux^2 * (y - at[[ind]])) * t[[ind]]
  z[[ind]] <- y - at[[ind]] - bt[[ind]] * x


  # Paso 3.
  f[[ind]] <- sqrt(t[[ind]])
  g[[ind]] <- f[[ind]] * xs[[ind]]
  h[[ind]] <- f[[ind]] * z[[ind]]

  # Paso 4.
  # Substep (i).
  F2[[ind]] <- sum(f[[ind]] * f[[ind]])

  # Substep (ii).
  g0[[ind]] <- sum(f[[ind]] * g[[ind]]) / F2[[ind]]
  h0[[ind]] <- sum(f[[ind]] * h[[ind]]) / F2[[ind]]

  # Substep (iii).
  gt[[ind]] <- g[[ind]] - g0[[ind]] * f[[ind]]
  ht[[ind]] <- h[[ind]] - h0[[ind]] * f[[ind]]

  # Substep (iv).
  Gt2[[ind]] <- sum(gt[[ind]] * gt[[ind]])

  # Substep (v).
  db[[ind]] <- sum(gt[[ind]] * ht[[ind]]) / Gt2[[ind]]
  da[[ind]] <- h0[[ind]] - db[[ind]] * g0[[ind]]

  # Paso 5.
  at[[ind+1]] <- at[[ind]] + da[[ind]]
  bt[[ind+1]] <- bt[[ind]] + db[[ind]]
  r[[ind]] <- ht[[ind]] - db[[ind]] * gt[[ind]]


  # Devuelve los resultados actualizados
  return(list(at <<- at, bt <<- bt, da <<- da, db <<- db,
              t <<- t, xs <<- xs, z <<- z, f <<- f, g <<- g, h <<- h,
              F2 <<- F2, g0 <<- g0, h0 <<- h0, gt <<- gt, ht <<- ht,
              Gt2 <<- Gt2, r <<- r))

}


### Segunta funcion tabla de cálculo de la regresión generalizada de distancias (GDR)

write_gdr1_tableaux <- function(x, ux, y, uy, t, xs, z, f, g, h, g0, h0, gt, ht,
                                at, bt, da, db, r, niter) {

  m <- length(x)
  ci1<-c(0, 0, 0, 0, at[[niter]], bt[[niter]], 0, 0, 0, 0)
  ci2<-cbind(x, ux, y, uy, t[[niter]], xs[[niter]], z[[niter]], f[[niter]], g[[niter]], h[[niter]])
  Ci<-rbind(ci1,ci2)
  cat("Iteración ", niter, " para determinar fi, gi, y hi")
  print(round(Ci,4))

  C1 <- c(0, 0, 0, g0[[niter]], h0[[niter]], 0, 0, da[[niter]], 0)
  C2 <- cbind(f[[niter]]^2, f[[niter]] * g[[niter]], f[[niter]] * h[[niter]],
              gt[[niter]], ht[[niter]], gt[[niter]]^2, gt[[niter]] * ht[[niter]],
              r[[niter]], r[[niter]]^2)
  C3 <- colSums(C2)
  C3[4] <- 0
  C3[5] <- 0
  C3[8] <- db[[niter]]
  Ca <- rbind(C1, C2, C3)
  cat("Iteración ", niter, " para determinar los incrementos delta_a y delta_b \n")
  print(round(Ca,4))
}



# Tercera Función  trazar datos con barras de error en x e y
plot_errorbar2 <- function(x, ux, y, uy, marker, linestyle, marker_color = "black",
                           main = "", xlab = "", ylab = "") {
  teex <- (max(x) - min(x)) / 100
  teey <- (max(y) - min(y)) / 100
  plot(x,y,pch = marker, col = marker_color, main = main, xlab = xlab, ylab = ylab) # Modificar para agregar el argumento col para el color de marcador
  for (i in 1:length(x)) {
    lines(c(x[i] - ux[i], x[i] + ux[i]), c(y[i], y[i]), lty = linestyle)
    lines(c(x[i] - ux[i], x[i] - ux[i]), c(y[i] - teey, y[i] + teey), lty = linestyle)
    lines(c(x[i] + ux[i], x[i] + ux[i]), c(y[i] - teey, y[i] + teey), lty = linestyle)
    lines(c(x[i], x[i]), c(y[i] - uy[i], y[i] + uy[i]), lty = linestyle)
    lines(c(x[i] - teex, x[i] + teex), c(y[i] - uy[i], y[i] - uy[i]), lty = linestyle)
    lines(c(x[i] - teex, x[i] + teex), c(y[i] + uy[i], y[i] + uy[i]), lty = linestyle)
  }
}
# Final de plot_errorbar2.R


# Cuarte Función de visualización del cuadro de cálculo de los mínimos cuadrados ponderados (WLS)

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

# Fin de write_wls_tableau.R


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
#### REGRESION DE DISTANCIA GENERALIZADA (Regresion ortogonal ponderada) ######
##############################################################################




# mediciones
  x <- c(1.2, 1.9, 2.9, 4.0, 4.7, 5.9)
# señales
  y <- c(3.4, 4.4, 7.2, 8.5, 10.8, 13.5)
# Incertidumbres de las mediciones
  ux <- c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
# Incertidumbres de las señales
  uy <- c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4)


  # tamaño de los datos #
 m <- length(x)

cat(" Obtener estimaciones de los parámetros de la función de calibración de línea recta
 y las incertidumbres estándar asociadas y la covarianza.
 Resolver el problema de regresión de distancia generalizada para obtener los
 mejores parámetros de ajuste de la línea recta.")

# Paso 1. Aproximación inicial usando mínimos cuadrados ponderados: ver
#  algm_wls_steps_1_to_8(x, y, uy)
  # Código de la función algm_wls_steps_1_to_8 en R
  # Reemplazar con la lógica y los cálculos necesarios



# Llamada a la función algm_wls_steps_1_to_8 con los argumentos x, y, uy
resultados <- algm_wls_steps_1_to_8(x, y, uy)

# Extracción de los resultados individuales
ai <- resultados$ai
bi <- resultados$bi


# Redondee las aproximaciones a los parámetros en el paso 1 a cuatro lugares decimales.
# (Este paso se incluye para producir los resultados proporcionados en ISO/TS 28037:2010:
#el paso generalmente no se realiza).

ai <- round(10000 * ai) / 10000
at <- list(ai)
bi <- round(10000 * bi) / 10000
bt <- list(bi)

# Repita los pasos 2 a 5 hasta que se cumplan los criterios de convergencia.
# (En este ejemplo, se considera que se ha alcanzado la convergencia cuando
# las magnitudes de los incrementos deltaa y deltab en a y b no son superiores a 0,00005.)
# Para los datos generales del usuario, se sugiere que  se considere que se ha alcanzado
# la convergencia cuando las magnitudes de todos los incrementos *relativos a las
# aproximaciones iniciales a los parámetros de la recta* no sean superiores a 1e-7.
#  En este caso, la tolerancia puede asignarse utilizando el comando tol = 1e-7*norm(c(ai, bi))
# Calcular la norma euclidiana del vector que contiene ai y bi: norm(c(ai, bi))

##

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



# Pasos 2 a 5 del wls ####
# Llamada a la función 'algm_gdr1_steps_2_to_5' con los argumentos correspondientes
# para la primera iteracion
ind <- 1
algm_gdr1_steps_2_to_5(x, ux, y, uy, at, bt, da, db, t, xs, z, f, g, h, F2,
                       g0, h0, gt, ht, Gt2, r, ind)



# Paso 6: Bucle while para iterar hasta que se alcance la convergencia ####
# Asignar tolerancias e inicializar variables.
tol <- 0.00005
# Bucle while para iterar hasta que la condición sea falsa
while ((abs(da[[ind]]) > tol) || (abs(db[[ind]]) > tol)) {
  # Actualizar el número de iteración
  ind <- ind + 1

  # Llamada a la función algm_gdr1_steps_2_to_5 con los argumentos correspondientes
  algm_gdr1_steps_2_to_5(x, ux, y, uy, at, bt, da, db, t, xs, z, f, g, h, F2,
                         g0, h0, gt, ht, Gt2, r, ind)
}


##


# Obtener los resultados finales de a y b
a <- at[[ind + 1]]
b <- bt[[ind + 1]]

##
# Paso 7: Evaluar las incertidumbres u2a, u2b, uab
u2a <- 1/F2[[ind]] + g0[[ind]]^2/Gt2[[ind]]
u2b <- 1/Gt2[[ind]]
uab <- -g0[[ind]]/Gt2[[ind]]



##
# Paso 8: Calcular el valor observado de la chi-cuadrado y los grados de libertad
chi_sq_obs <- sum(r[[ind]] * r[[ind]])
nu <- (m-2)

##
# Paso 9. Calcular el 95% del percentil de la distribución chi-cuadrado
chi_sq <- calc_chi_sq_95_percent_quantile(nu);

## Visualizar información en pantalla y generar cifras

##
# Modelo de medición.
cat("\n MODELO PARA LAS INCERTIDUMBRES ASOCIADAS AL XI Y AL Y_i
    ISO/TS 28037:2010(E) CLAUSULA 7
    EJEMPLO Tabla 10")

##
# Datos de medición.

cat("AJUSTE
    Datos correspondientes a ", m, "puntos de medición\n")
for (i in 1:m) {
  cat(sprintf("%8.1f %8.1f %8.1f %8.1f \n", x[i], ux[i], y[i], uy[i]))
}
cat("\n")

##
# Tabla de cálculo para el problema WLS: ver
#  <write_wls_tableau.html write_wls_tableau.m>.
cat("Aproximaciones iniciales a los parámetros  \n")
write_wls_tableau(x, y, resultados$wi, resultados$g0i, resultados$h0i, resultados$gi,
                  resultados$hi, resultados$ai, resultados$bi, resultados$ri)


##
# Aproximaciones iniciales a los parámetros.
cat("Aproximación inicial a la intercepción ",sprintf("%9.4f ", ai))
cat("Aproximación inicial a la pendiente ",sprintf("%9.4f ", bi))


# Cuadro de cálculo para cada iteración:
# <write_gdr1_tableaux.html write_gdr1_tableaux.m>.

for(niter in 1:ind){
  cat("Iteration ", niter, "\n")
  write_gdr1_tableaux(x, ux, y, uy, t, xs, z, f, g, h, g0, h0, gt, ht, at, bt, da, db, r, niter)
}

##
# Estimaciones de la solución
cat("Estimación del intercepto", a)

cat("Estimación de la pendiente", b)

# Incertidumbres estándar asociadas a las estimaciones de las soluciones
cat("Incertidumbre estándar asociada a la estimación del intercepto", sqrt(u2a))
cat("Incertidumbre estándar asociada a la estimación de la pendiente", sqrt(u2b))
cat("Covarianza asociada a las estimaciones del intercepto y la pendiente",uab)



##
# Validación del modelo
cat("VALIDACIÓN
    Grados de libertad", nu)

cat("Valor chi-cuadrado observado\n",chi_sq_obs)

cat("Cuantil del 95% de la distribución chi-cuadrado con",nu,"grados de libertad", chi_sq)

## qchisq(.95,nu)

if (chi_sq_obs > chi_sq) {
  cat("LA PRUEBA DE CHI-CUADRADO HA FALLADO - SE RECHAZA EL MODELO")
} else {
  cat("SE SUPERA LA PRUEBA DE CHI-CUADRADO - SE ACEPTA EL MODELO")
}



# Figuras
par(lwd = 2, cex = 0.8, font = 2)

plot_errorbar2(x, ux, y, uy,  marker=1,  linestyle=1,main="Función código implementado")
points(xs[[ind]], at[[ind]] + bt[[ind]] * xs[[ind]], col = "red", pch = "o")
lines(xs[[ind]], at[[ind]] + bt[[ind]] * xs[[ind]], col = "blue")

lab<-paste0(round(at[[ind]],4),'+',paste0(round(bt[[ind]],4),'x'))

text(10,max(y) ,  lab )

xlabel <- expression(italic(x))
ylabel <- expression(italic(y))
axis1 <- range(x, na.rm = TRUE)
  plot(axis1, c(0, 0), type = "n", xlab = xlabel, ylab = ylabel,
       ylim=c(min(r[[1]]),max(r[[1]])),  main="Residuales del modelo")

for (i in 1:m) {
  segments(x[i], 0, x[i], r[[ind]][i], col = "black", lty = "solid")
  points(x[i], r[[ind]][i], pch = "o", col = "black", bg = "white", cex = 0.6)
}
abline(h = 0, col = "blue", lty = "dashed")


Gdr1 <- data.frame(Método = c("GDR"),
                  Intercept = 0, U_Intrcpt = 0, Slope = 0, U_Slope = 0, Cor.est=0)

Gdr1$Intercept[1] <- at[[ind]] # Intercept
Gdr1$U_Intrcpt[1] <- sqrt(u2a)   # U-Intercept
Gdr1$Slope[1]     <- bt[[ind]] # Slope
Gdr1$U_Slope[1]   <- sqrt(u2b) # U-Slope
Gdr1$Cor.est[1]  <-(uab)

print(Gdr1)



################################################################################
###################  PREDICCIONES DEL MODELO  ##################################
################################################################################



# Evaluación inversa - Prediccíon
yp <- 6
uyp <- 0.3

a<-rslt$coefficients[1]
b<-rslt$coefficients[2]


xp <- (yp - a)/b
ca <- -1/b
cb <- -((yp - a)/(b^2))
cy <- 1/b
u2xp <- ca^2*u2a + cb^2*u2b + 2*ca*cb*uab + cy^2*uyp^2
uxp<-sqrt(u2xp)
int<-x+c(-1,1)*pt(0.025,m-2)*sqrt(u2xp/m)



plot_errorbar2(x,  ux, y,  uy,  marker=1,  linestyle=1,main="Nuevas observaciones",
               xlab=xlabel,ylab=ylabel) # Puntos con barras de error
lines(x, a+b*x, col="blue")

lines(c(xp - uxp, xp + uxp), c(yp, yp), lty = 1,col="red")
lines(c(xp - uxp, xp - uxp), c(yp - teey, yp + teey), lty = 1,col="red")
lines(c(xp + uxp, xp + uxp), c(yp - teey, yp + teey), lty = 1,col="red")
lines(c(xp, xp), c(yp - uyp, yp + uyp), lty = 1,col="red")
lines(c(xp - teex, xp + teex), c(yp - uyp, yp - uyp), lty = 1,col="red")
lines(c(xp - teex, xp + teex), c(yp + uyp, yp + uyp), lty = 1,col="red")

lines(c(xp + uxp, xp + uxp), c(0, yp + teey), lty = 3,col="red")
lines(c(xp - uxp, xp - uxp), c(0, yp + teey), lty = 3,col="red")

arrows(0,yp,xp - 1.2*uxp,yp,col="pink",length = 0.1)

lab<-paste0(round(a,4),'+',paste0(round(b,4),'x'))
text(mean(x),max(y) ,  lab )
#puntual
teex <- (max(x) - min(x)) / 100
teey <- (max(y) - min(y)) / 100
lines(c(xp - sqrt(u2xp), xp + sqrt(u2xp)), c(yp, yp), lty = 1)
lines(c(xp - sqrt(u2xp), xp - sqrt(u2xp)), c(yp - teey, yp + teey), lty = 1)
lines(c(xp + sqrt(u2xp), xp + sqrt(u2xp)), c(yp - teey, yp + teey), lty = 1)
lines(c(xp, xp), c(yp - uyp, yp + uyp), lty = 1)
lines(c(xp - teex, xp + teex), c(yp - uyp, yp - uyp), lty = 1)
lines(c(xp - teex, xp + teex), c(yp + uyp,yp + uyp), lty = 1)




## Evaluación directa

xf <- 3.5
uxf <- 0.2
a<-rslt$coefficients[1]
b<-rslt$coefficients[2]

yf <- a + b*xf
ca <- 1
cb <- xf
cx <- b
u2yf <- ca^2*u2a + cb^2*u2b + 2*ca*cb*uab + cx^2*uxf^2


yf +c(-1,1)*u2yf


rnumbers<-rnorm(3,0.1,0.02)
shapiro.test(rnumbers)
plot(rnumbers)
abline(h=0.1)

