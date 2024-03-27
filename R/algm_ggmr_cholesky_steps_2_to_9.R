# Function to undertake steps 2 to 9 of the GGMR algorithm (Clause 10)
algm_ggmr_cholesky_steps_2_to_9 <- function(x, y, U, tt, dt, f, J, ft, Jt, g, H, M, q, ind) {

  # Step 2
  m <- length(x)
  f[[ind]] <- c(x - tt[[ind]][1:m], y - (tt[[ind]][m+1] + tt[[ind]][m+2]*tt[[ind]][1:m]))
  J[[ind]] <- rbind( cbind(-diag(m), matrix(0, nrow = m, ncol = 2)),cbind(
                    -tt[[ind]][m+2] * diag(m), -matrix(1, nrow = m, ncol = 1),
                    -tt[[ind]][1:m]))

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
  return(list(tt <<- tt, dt <<- dt, f<<- f, J <<- J, ft<<- ft, Jt <<- Jt, g <<- g, H <<- H, M <<- M, q <<- q))
 }


