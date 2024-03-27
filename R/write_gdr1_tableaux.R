## Función de visualización de la tabla de cálculo de la regresión generalizada de distancias (GDR)

write_gdr1_tableaux <- function(x, ux, y, uy, t, xs, z, f, g, h, g0, h0, gt, ht, at, bt, da, db, r, niter) {

  m <- length(x)
  ci1<-c(0, 0, 0, 0, at[[niter]], bt[[niter]], 0, 0, 0, 0)
  ci2<-cbind(x, ux, y, uy, t[[niter]], xs[[niter]], z[[niter]], f[[niter]], g[[niter]], h[[niter]])
  Ci<-rbind(ci1,ci2)
  cat("Iteración ", niter, " para determinar fi, gi, y hi")
  print(round(Ci,4))

  C1 <- c(0, 0, 0, g0[[niter]], h0[[niter]], 0, 0, da[[niter]], 0)
  C2 <- cbind(f[[niter]]^2, f[[niter]] * g[[niter]], f[[niter]] * h[[niter]], gt[[niter]], ht[[niter]], gt[[niter]]^2, gt[[niter]] * ht[[niter]], r[[niter]], r[[niter]]^2)
  C3 <- colSums(C2)
  C3[4] <- 0
  C3[5] <- 0
  C3[8] <- db[[niter]]
  Ca <- rbind(C1, C2, C3)
  cat("Iteración ", niter, " para determinar los incrementos delta_a y delta_b \n")
  print(round(Ca,4))
}
## End of write_gdr1_tableaux.R


