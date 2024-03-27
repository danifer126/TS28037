
#' @title Función para resolver el problema de mínimos cuadrados ponderados (WLS) (Cláusula 6)
#' @description Este documento ejecuta el ejemplo numérico de mínimos cuadrados ponderados (WLS)
#' con pesos iguales conocidos, descrito en la Cláusula 6 (Modelo para incertidumbres asociadas a los y_i),
#' y realiza la predicción descrita en el Ejemplo 1 de la Cláusula 11.1 y la evaluación hacia adelante
#' descrita en el Ejemplo 2 de la Cláusula 11.2.


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
