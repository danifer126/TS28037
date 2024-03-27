





# Función para llevar a cabo los pasos 2 a 5 del algoritmo GDR (Cláusula 7)
algm_gdr1_steps_2_to_5 <- function(x, ux, y, uy, at, bt, da, db,t, xs, z, f, g, h, F2, g0, h0, gt, ht, Gt2, r, ind){

  # Step 2.
  m <<- length(x)
  t[[ind]] <<- 1 / (uy^2 + bt[[ind]]^2 * ux^2)
  xs[[ind]] <<- (x * uy^2 + bt[[ind]] * ux^2 * (y - at[[ind]])) * t[[ind]]
  z[[ind]] <<- y - at[[ind]] - bt[[ind]] * x


  # Step 3.
  f[[ind]] <<- sqrt(t[[ind]])
  g[[ind]] <<- f[[ind]] * xs[[ind]]
  h[[ind]] <<- f[[ind]] * z[[ind]]

  # Step 4.
  # Substep (i).
  F2[[ind]] <<- sum(f[[ind]] * f[[ind]])

  # Substep (ii).
  g0[[ind]] <<- sum(f[[ind]] * g[[ind]]) / F2[[ind]]
  h0[[ind]] <<- sum(f[[ind]] * h[[ind]]) / F2[[ind]]

  # Substep (iii).
  gt[[ind]] <<- g[[ind]] - g0[[ind]] * f[[ind]]
  ht[[ind]] <<- h[[ind]] - h0[[ind]] * f[[ind]]

  # Substep (iv).
  Gt2[[ind]] <<- sum(gt[[ind]] * gt[[ind]])

  # Substep (v).
  db[[ind]] <<- sum(gt[[ind]] * ht[[ind]]) / Gt2[[ind]]
  da[[ind]] <<- h0[[ind]] - db[[ind]] * g0[[ind]]

  # Step 5.
  at[[ind+1]] <<- at[[ind]] + da[[ind]]
  bt[[ind+1]] <<- bt[[ind]] + db[[ind]]
  r[[ind]] <<- ht[[ind]] - db[[ind]] * gt[[ind]]


  # Devuelve los resultados actualizados
  return(list(at <- at, bt<- bt, da <- da, db <- db,
              t <- t, xs <- xs, z <- z, f <- f, g <- g, h <- h,
              F2 <- F2, g0 <- g0, h0 <- h0, gt <- gt, ht <- ht,
              Gt2 <- Gt2, r <- r))

}

algm_ggmr_cholesky_steps_2_to_9 <- function(x, y, U, tt, dt, f, J, ft, Jt, g, H, M, q, ind) {

  # Step 2
  m <- length(x)
  f[[ind]] <- c(x - tt[[ind]][1:m], y - (tt[[ind]][m+1] + tt[[ind]][m+2]*tt[[ind]][1:m]))
  J[[ind]] <- rbind( cbind(-diag(m), matrix(0, nrow = m, ncol = 2)),cbind(
    -tt[[ind]][m+2] * diag(m), -matrix(1, nrow = m, ncol = 1),-tt[[ind]][1:m]))

  # Step 3 (Lower)
  L1 <- chol(U[1:m,1:m])
  L2 <- chol(U[(m+1):(2*m),(m+1):(2*m)])
  L<-diag(2*m)
  L[1:m,1:m] <- t(L1)
  L[(m+1):(2*m),(m+1):(2*m)] <- t(L2)


  # Step 4
  ft[[ind]] <- solve(L,f[[ind]])
  Jt[[ind]] <- solve(L, J[[ind]])

  # Step 5
  g[[ind]] <- t(Jt[[ind]]) %*% ft[[ind]]
  H[[ind]] <- t(Jt[[ind]]) %*% Jt[[ind]]

  # Step 6
  M[[ind]] <- t(chol(H[[ind]]))

  # Step 7
  q[[ind]] <- solve(M[[ind]], -g[[ind]])

  # Step 8
  dt[[ind]] <- solve(t(M[[ind]]), q[[ind]])

  # Step 9
  tt[[ind + 1]] <- tt[[ind]] + dt[[ind]]

  # Return the updated variables
  return(list(tt <- tt, dt <- dt, f<- f, J <- J, ft<- ft, Jt <- Jt, g <- g, H <- H, M <- M, q <- q))
}



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

  return(list(ai = a, bi = b, u2ai = u2a, u2bi = u2b, uabi = uab, wi= w, g0i = g0, h0i = h0, gi = g, hi = h, G2i = G2, ri = r, Ri = R))
}



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

table <- cbind(nu = c(1:100), quantile)
table<- as.data.frame(table)

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



plot_errorbar2 <- function(x, ux, y, uy, marker, linestyle, marker_color = "black", main = "", xlab = "", ylab = "") {
  teex <- (max(x) - min(x)) / 100
  teey <- (max(y) - min(y)) / 100
  plot(x,y,pch = marker, col = marker_color, main = main, xlab = xlab, ylab = ylab,xlim=c(0,max(x)),ylim=c(0,max(y))) # Modificar para agregar el argumento col para el color de marcador
  for (i in 1:length(x)) {
    lines(c(x[i] - ux[i], x[i] + ux[i]), c(y[i], y[i]), lty = linestyle)
    lines(c(x[i] - ux[i], x[i] - ux[i]), c(y[i] - teey, y[i] + teey), lty = linestyle)
    lines(c(x[i] + ux[i], x[i] + ux[i]), c(y[i] - teey, y[i] + teey), lty = linestyle)
    lines(c(x[i], x[i]), c(y[i] - uy[i], y[i] + uy[i]), lty = linestyle)
    lines(c(x[i] - teex, x[i] + teex), c(y[i] - uy[i], y[i] - uy[i]), lty = linestyle)
    lines(c(x[i] - teex, x[i] + teex), c(y[i] + uy[i], y[i] + uy[i]), lty = linestyle)
  }
}



write_gdr1_tableaux <- function(x, ux, y, uy, t, xs, z, f, g, h, g0, h0, gt, ht, at, bt, da, db, r, niter, fstr) {

  m <- length(x)
  ci1<-c(0, 0, 0, 0, at[[niter]], bt[[niter]], 0, 0, 0, 0)
  ci2<-cbind(x, ux, y, uy, t[[niter]], xs[[niter]], z[[niter]], f[[niter]], g[[niter]], h[[niter]])

  Ci<-rbind(ci1,ci2)

  cat("Iteración ", niter, " para determinar fi, gi, y hi \n")
  cat("\n")
  print(round(Ci,4))



  C1 <- c(0, 0, 0, g0[[niter]], h0[[niter]], 0, 0, da[[niter]], 0)
  C2 <- cbind(f[[niter]]^2, f[[niter]] * g[[niter]], f[[niter]] * h[[niter]], gt[[niter]], ht[[niter]], gt[[niter]]^2, gt[[niter]] * ht[[niter]], r[[niter]], r[[niter]]^2)
  C3 <- colSums(C2)
  C3[4] <- 0
  C3[5] <- 0
  C3[8] <- db[[niter]]
  Ca <- rbind(C1, C2, C3)
  cat("Iteration ", niter, " para determinar los incrementos deltaa y deltab \n")
  print(round(Ca,4))


}

write_gdr2_tableaux <- function(x, ux, y, uy, uxy, t, xs, z, f, g, h, g0, h0, gt, ht,
                                at, bt, da, db, r, niter, fstr) {
  m <- length(x)
  ci1 <- c(0, 0, 0, 0, 0, at[[niter]], bt[[niter]], 0, 0, 0, 0)
  ci2 <-cbind(  x, ux, y, uy, uxy, t[[niter]], xs[[niter]], z[[niter]], f[[niter]],
                g[[niter]], h[[niter]])
  ci<-rbind(ci1,ci2)

  cat(paste0('Iteración', niter, 'para determinar fi, gi, y hi\n'))
  for (i in 1:(m+1)) {
    cat(sprintf(fstr, ci[i,]), '\n')
  }
  cat('\n')
  C1 <- c(0, 0, 0, g0[[niter]], h0[[niter]], 0, 0, da[[niter]], 0)
  C2 <- cbind(f[[niter]] * f[[niter]], f[[niter]] * g[[niter]], f[[niter]] * h[[niter]],
              gt[[niter]], ht[[niter]], gt[[niter]] * gt[[niter]], gt[[niter]] * ht[[niter]],
              r[[niter]], r[[niter]] * r[[niter]])
  C3 <- colSums(C2)
  C3[4] <- 0
  C3[5] <- 0
  C3[8] <- db[[niter]]
  Ca <- rbind(C1, C2, C3)
  cat(paste0('Iteración ', niter, 'para determinar los incrementos deltaa y deltab \n'))
  for (i in 1:(m+2)) {
    cat(sprintf(fstr, Ca[i,]), '\n')
  }
  cat('\n')
}


write_ggmr_tableau <- function(tt, dt, ind, fstr) {
  m <- length(tt[[1]]) - 2
  TT <- matrix(0, nrow = m + 2, ncol = ind + 2)
  TT[, 1] <- tt[[1]]

  for (niter in 1:(ind-1)) {
    TT[, niter + 1] <- dt[[niter]]
  }
  TT[, ind + 2] <- tt[[ind + 1]]

  cat('Resumen del procedimiento iterativo (aproximaciones iniciales, correcciones y estimaciones finales)\n')
  print(TT)
  cat('\n')
}

write_gmr_tableaux <- function(f, g, h, g0, h0, gt, ht, a, b, r, fstr) {

  # Cuadro de cálculo inicial.
  cat("Cuadro de cálculo inicial.\n")
  m <- length(f)
  for (i in 1:m) {
    cat( f[i], g[i], h[i],'\n')
  }
  cat('\n')

  # Cuadro de cálculo principal.
  C1 <- c(0, 0, 0, g0, h0, 0, 0, a, 0)
  C2 <- cbind(f*f, f*g, f*h, gt, ht, gt*gt, gt*ht, r, r*r)
  C3 <- colSums(C2)
  C3[4] <- 0
  C3[5] <- 0
  C3[8] <- b
  Ca <- rbind(C1, C2, C3)
  cat("Cuadro de cálculo principal\n")
  for (i in 1:(m+2)) {
    cat(sprintf( '%9.4f', Ca[i, ]),'\n')
  }
  cat('\n')
}

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


### Métodos



GDR1<- function(x,ux,y,uy){


  # Asignar datos de medición
  #  x <- c(1.2, 1.9, 2.9, 4.0, 4.7, 5.9)
  # # #
  # # #
  #  y <- c(3.4, 4.4, 7.2, 8.5, 10.8, 13.5)
  # # #
  #  ux <- c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
  # # #
  #  uy <- c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4)



  #  x<-c(5.032,15.068,24.905,35.074,44.781,55.230,65.366,75.021,84.387,94.862)
  #  y<-c(0.02193, 0.072, 0.12834, 0.17902, 0.21606, 0.28298, 0.32895, 0.35127, 0.40098, 0.47454)
  #
  #
  #
  #
  # #
  # ux<-c(0.004,0.005,0.008,0.01,0.013,0.016,0.019,0.021,0.024,0.027)
  # uy<-c(0.00008,0.00014,0.00022,0.00029,0.00023,0.00045,0.00045,0.00084,0.00056,0.00055)
  # m <- length(x)

  #, Obtener estimaciones de los parámetros de la función de calibración de línea recta
  # y las incertidumbres estándar asociadas y la covarianza.
  # Resolver el problema de regresión de distancia generalizada para obtener los
  # mejores parámetros de ajuste de la línea recta.

  # Paso 1. Aproximación inicial usando mínimos cuadrados ponderados: ver
  # <algm_wls_steps_1_to_8.html algm_wls_steps_1_to_8.m>.

  #  algm_wls_steps_1_to_8(x, y, uy)
  # Código de la función algm_wls_steps_1_to_8 en R
  # Reemplazar con la lógica y los cálculos necesarios




  # Llamada a la función algm_wls_steps_1_to_8 con los argumentos x, y, uy
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
  # Asignar tolerancias e inicializar variables.
  tol <- 0.000000000001;
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


  ##
  # Pasos 2 a 5: véase

  ind <- 1


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
  cat("\n MODELO PARA LAS INCERTIDUMBRES ASOCIADAS AL XI Y AL YI\n\n")
  cat("ISO/TS 28037:2010(E) CLAUSULA 7 \n")
  cat("EJEMPLO \n\n")

  ##
  # Datos de medición.

  cat("AJUSTE \n")
  cat("Datos correspondientes a ", m, "puntos de medición\n")
  for (i in 1:m) {
    cat(sprintf("%8.1f %8.4f %8.1f %8.4f \n", x[i], ux[i], y[i], uy[i]))
  }
  cat("\n")

  ##
  # Tabla de cálculo para el problema WLS: ver
  #  <write_wls_tableau.html write_wls_tableau.m>.
  cat("Aproximaciones iniciales a los parámetros  \n")
  write_wls_tableau(x, y, wi, g0i, h0i, gi, hi, ai, bi, ri, fstr = "%9.4f")


  ##
  # Aproximaciones iniciales a los parámetros.
  cat("Aproximación inicial a la intercepción  \n")
  cat(sprintf("%9.4f \n\n", ai))

  cat("Aproximación inicial a la pendiente \n")
  cat(sprintf("%9.4f \n\n", bi))

  ##
  # Cuadro de cálculo para cada iteración: ver
  # <write_gdr1_tableaux.html write_gdr1_tableaux.m>.

  for(niter in 1:ind){
    cat("Iteration ", niter, "\n")
    write_gdr1_tableaux(x, ux, y, uy, t, xs, z, f, g, h, g0, h0, gt, ht, at, bt, da, db, r, 2, '%9.4f')
  }

  ##
  # Estimaciones de la solución
  cat("Estimación del intercepto\n")
  cat(sprintf("%9.4f\n", a))
  cat("\n")

  cat("Estimación de la pendiente\n")
  cat(sprintf("%9.4f\n", b))
  cat("\n")


  ##
  # Incertidumbres estándar asociadas a las estimaciones de las soluciones
  cat("Incertidumbre estándar asociada a la estimación del intercepto\n")
  cat(sprintf("%9.4f\n", sqrt(u2a)))
  cat("\n")

  cat("Incertidumbre estándar asociada a la estimación de la pendiente\n")
  cat(sprintf("%9.4f\n", sqrt(u2b)))
  cat("\n")

  cat("Covarianza asociada a las estimaciones del intercepto y la pendiente\n")
  cat(sprintf("%9.4f\n", uab))
  cat("\n")

  ##
  # Validación del modelo
  cat("VALIDACIÓN\n")
  cat("Grados de libertad\n")
  cat(sprintf("%9.3f", nu))
  cat("\n")
  cat("Valor chi-cuadrado observado\n")
  cat(sprintf("%9.3f", chi_sq_obs))
  cat("\n")
  cat("Cuantil del 95% de la distribución chi-cuadrado con",nu,"grados de libertad\n")
  cat(sprintf("%9.3f", chi_sq))
  cat("\n")

  if (chi_sq_obs > chi_sq) {
    cat("LA PRUEBA DE CHI-CUADRADO HA FALLADO - SE RECHAZA EL MODELO RECTILÍNEO\n\n")
  } else {
    cat("SE SUPERA LA PRUEBA DE CHI-CUADRADO - SE ACEPTA EL MODELO RECTILÍNEO\n\n")
  }


  ##
  # Figuras
  par(lwd = 2, cex = 0.8, font = 2)

  plot_errorbar2(x, ux, y, uy,  marker=1,  linestyle=1,main="Función de calibración con incertidumbre asociada")
  points(xs[[ind]], at[[ind]] + bt[[ind]] * xs[[ind]], col = "red", pch = "o")
  lines(xs[[ind]], at[[ind]] + bt[[ind]] * xs[[ind]], col = "blue")
  lab<-paste0(round(at[[ind]],4),'+',paste0(round(bt[[ind]],4),'x'))

  text(10,max(y) ,  lab )

  xlabel <- expression(italic(x))
  ylabel <- expression(italic(y))
  axis1 <- range(x, na.rm = TRUE)
  plot(axis1, c(0, 0), type = "n", xlab = xlabel, ylab = ylabel, main="Residuales del modelo")
  for (i in 1:m) {
    segments(x[i], 0, x[i], r[[ind]][i], col = "black", lty = "solid")
    points(x[i], r[[ind]][i], pch = "o", col = "black", bg = "white", cex = 0.6)
  }
  abline(h = 0, col = "blue", lty = "dashed")
}

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


GMR<- function(x,y,Uy){
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
}

### REGRESIONES PODERADAS ###

WLS1 <- function(x, y, uy){
  m <- length(x)
  # Paso 1.
  w <- 1/uy
  F2 <- sum(w * w)

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

  ## Mostrar información en pantalla y generar figuras


  # Measurement model.
  cat("\nMODELO PARA LAS INCERTIDUMBRES ASOCIADAS AL YI  \n\n")
  cat("ISO/TS 28037:2010(E) CLÁUSULA 6 \n")
  cat("EJEMPLO (PESOS IGUALES) \n\n")

  ##
  # Measurement data.
  cat("Ajuste\n")
  cat("Datos que representan", m, "puntos de medida, con pesos iguales\n")
  # cat(" x y uy\n")
  for (i in 1:m) {
    cat(sprintf("%8.1f %8.1f %8.1f\n", x[i], y[i], uy[i]))
  }
  cat("\n")


  ##
  # Tabla de cálculo:
  # Ver: write_wls_tableau.R
  write_wls_tableau(x, y, w, g0, h0, g, h, a, b, r, '%.3f ');

  ##
  # Estimaciones de solución.

  cat(sprintf("Estimación del intercepto : %8.3f\n", a))
  cat(sprintf("Estimación de la pendiente: %8.3f\n", b))
  cat(sprintf("Incertidumbre estándar asociada a la estimación del intercepto: %.3f\n", sqrt(u2a)))
  cat(sprintf("Incertidumbre estándar asociada con la estimación de la pendiente: %.3f\n", sqrt(u2b)))
  cat(sprintf("Covarianza asociada con las estimaciones del intercepto y la pendiente: %.3f\n", uab))

  ##
  # Validación del modelo.
  cat("VALIDACIÓN \n")
  cat("Grados de libertad \n")
  cat(sprintf("%4f \n\n", nu))
  cat("Valor chi-cuadrado observado \n")
  cat(sprintf("%9.3f \n\n", chi_sq_obs))
  cat("Cuantil del 95% de la distribución chi-cuadrado con", nu ," grados de libertad")
  cat(sprintf("\n%9.3f \n\n", chi_sq))
  if (chi_sq_obs > chi_sq) {
    cat("FALLÓ LA PRUEBA DE CHI CUADRADO - SE RECHAZA EL MODELO DE LÍNEA RECTA \n\n")
  } else {
    cat("PRUEBA DE CHI CUADRADO APROBADA - MODELO DE LÍNEA RECTA ACEPTADO \n\n")
  }


  # Configuración de fuente y tamaño de letra
  par(lwd=2, cex=1.2, font.lab=2)

  # Gráfico de líneas con barras de error
  plot(x, y, type="n", xlab=expression(italic(x)), ylab=expression(italic(y)))
  arrows(x, y-uy, x, y+uy, length=0.05, angle=90, code=3, col="black")
  lines(x, a+b*x, col="blue")

  # Segundo gráfico
  plot(x, r, type="n", xlab=expression(italic(x)), ylab=expression(italic(r)))
  for (i in 1:m) {
    segments(x[i], 0, x[i], r[i], col="black")
    points(x[i], r[i], pch=21, bg="white", cex=1.2)
  }
  segments(par("usr")[1], 0, par("usr")[2], 0, lty=2, col="blue")

  # Prediccíon
  yp <- 10.5
  uyp <- 0.5
  xp <- (yp - a)/b
  ca <- -1/b
  cb <- -(yp - a)/(b^2)
  cy <- 1/b
  u2xp <- ca^2*u2a + cb^2*u2b + 2*ca*cb*uab + cy^2*uyp^2
  cat("ISO/TS 28037:2010(E) 11.1\n")
  cat("EJEMPLO 1\n\n")
  cat("PREDICCIÓN\n")
  cat("Valor y medido\n")
  cat(sprintf("%8.3f\n\n", yp))
  cat("Incertidumbre estándar asociada al valor y medido\n")
  cat(sprintf("%8.3f\n\n", uyp))
  cat("Estimación de x\n")
  cat(sprintf("%8.3f\n\n", xp))
  cat("Coeficiente de sensibilidad wrt a\n")
  cat(sprintf("%8.3f\n\n", ca))
  cat("Coeficiente de sensibilidad para b\n")
  cat(sprintf("%8.3f\n\n", cb))
  cat("Coeficiente de sensibilidad para y\n")
  cat(sprintf("%8.3f\n\n", cy))
  cat("Incertidumbre estándar asociada a la estimación de x\n")
  cat(sprintf("%8.3f\n\n", sqrt(u2xp)))


  ## Evaluación previa
  xf <- 3.5
  uxf <- 0.2
  yf <- a + b*xf
  ca <- 1
  cb <- xf
  cx <- b
  u2yf <- ca^2*u2a + cb^2*u2b + 2*ca*cb*uab + cx^2*uxf^2
  cat("ISO/TS 28037:2010(E) 11.2 \n")
  cat("EJEMPLO \n\n")
  cat("EVALUACIÓN ADELANTADA \n")
  cat("Valor x medido \n", cat(sprintf("%8.3f \n\n", xf)))
  cat("Incertidumbre estándar asociada al valor x medido \n", cat(sprintf("%8.3f \n\n", uxf)))
  cat("Estimación de y \n", cat(sprintf("%8.3f \n\n", yf)))
  cat("Coeficiente de sensibilidad respecto a a \n", cat(sprintf("%8.3f \n\n", ca)))
  cat("Coeficiente de sensibilidad wrt b \n", cat(sprintf("%8.3f \n\n", cb)))
  cat("Coeficiente de sensibilidad wrt x \n", cat(sprintf("%8.3f \n\n", cx)))
  cat("Incertidumbre estándar asociada a la estimación de y \n", cat(sprintf("%8.3f \n\n", sqrt(u2yf))))
}

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


