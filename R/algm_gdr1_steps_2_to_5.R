# Función para llevar a cabo los pasos 2 a 5 del algoritmo GDR (Cláusula 7)
algm_gdr1_steps_2_to_5 <- function(x, ux, y, uy, at, bt, da, db,t, xs, z, f, g, h, F2, g0, h0, gt, ht, Gt2, r, ind){



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
