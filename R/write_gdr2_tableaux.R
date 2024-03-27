## Función para mostrar la tabla de cálculo de la regresión de distancia generalizada (GDR)

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

## Fin de write_gdr2_tableaux.m
