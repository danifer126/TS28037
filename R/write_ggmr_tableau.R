## Función para mostrar la tabla de cálculo para la regresión generalizada
# de Gauss Markov (GMR) 'function write_ggmr_tableau(tt, dt, ind, fstr)'

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

##
# End of write_ggmr_tableau.R
