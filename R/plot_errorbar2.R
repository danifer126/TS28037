# Función para trazar datos con barras de error en x e y
plot_errorbar2 <- function(x, ux, y, uy, marker, linestyle, marker_color = "black", main = "", xlab = "", ylab = "") {
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


