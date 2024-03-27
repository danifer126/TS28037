## Función para visualizar la tabla de cálculo de la regresión de Gauss Markov (GMR)


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



## Fin de write_gmr_tableaux.R
