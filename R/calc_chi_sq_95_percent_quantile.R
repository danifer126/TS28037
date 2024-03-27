# Función para calcular el cuantil del 95 % de la distribución chi-cuadrado
## función chi_sq = calc_chi_sq_95_percent_quantile(nu)
##
# Para valores de $\nu \leq 100$, los cuantiles se proporcionan en una tabla de consulta.
# Para valores de $\nu > 100$, el resultado en La distribución de chi-cuadrado E B Wilson y M M Hilferty

# Actas de la Academia Nacional de Ciencias de los Estados Unidos de
# América, 17 (1931):684-688
# que $(\chi^2/\nu)^{1/3}$ tiene una distribución aproximadamente normal con media
# Se usa $1 - 2/(9 \nu)$ y la varianza $2/(9 \nu)$.

# En el rango de $\nu = 101$ a $\nu = 1000$, la aproximación de Wilson y Hilferty
# tiene un error absoluto máximo de 0,002 5 y un error porcentual relativo máximo de 0,002 %.

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

# Fin de 'calc_chi_sq_95_percent_quantile'

