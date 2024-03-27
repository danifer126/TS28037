# Función de visualización del cuadro de cálculo de los mínimos cuadrados ponderados (WLS)

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





