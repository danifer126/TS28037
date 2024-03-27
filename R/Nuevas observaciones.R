

x<-c(5.032,15.068, 24.905,35.074, 44.781, 55.23, 65.366, 75.021, 84.387, 94.862)
y<-c(0.02193,  0.072, 0.12834,0.17902, 0.21606, 0.28298,0.32895,0.35127,0.40098,0.47454)
ux<-c(0.004,  0.005,0.008,0.01, 0.013,0.016,0.019, 0.021,0.024,  0.027)
uy<-c(0.00008,0.00014,0.00022,0.00029,0.00023,0.00045,0.00045,0.00084,0.00056,  0.00055)


GDR1(x,ux,y,uy)

library(ggmr)

rslt <- ggmr::ggmr(x,y,diag(ux^2),diag(uy^2) )
print(rslt)
sqrt(diag(rslt$cov))

?ggmr.pred()


plot(rslt$chisq.validation)


#rslt <- ggmr::ggmr(x,y, Ux= diag(0, length(x)),Uy = diag(1, length(x)))



GRD_t <- data.frame(Método = c("GDR PAQUETE GGMR"),
                    Intercept = 0, U_Intrcpt = 0, Slope = 0, U_Slope = 0, Cor.est=0)
#, R2_ggpVal = 0,Rslt_Intp = 0, U_c = 0, U_exp = 0, row.names = 1) # Data to estimate the uncertainty using the first function (funcionUP).

GRD_t$Intercept<- rslt$coefficients[1] # Intercept
GRD_t$U_Intrcpt<- rslt$cov[1,1]  # U-Intercept
GRD_t$Slope<- rslt$coefficients[2]      # Slope
GRD_t$U_Slope<- rslt$cov[2,2]    # U-Slope
Gdr1$Cor.est[1]  <-rslt$cov[1,2]


# Evaluación inversa - Prediccíon
yp <- 6
uyp <- 0.3

a<-rslt$coefficients[1]
b<-rslt$coefficients[2]


xp <- (yp - a)/b
ca <- -1/b
cb <- -((yp - a)/(b^2))
cy <- 1/b
u2xp <- ca^2*u2a + cb^2*u2b + 2*ca*cb*uab + cy^2*uyp^2
uxp<-sqrt(u2xp)
int<-x+c(-1,1)*pt(0.025,m-2)*sqrt(u2xp/m)



plot_errorbar2(x,  ux, y,  uy,  marker=1,  linestyle=1,main="Nuevas observaciones",
               xlab=xlabel,ylab=ylabel) # Puntos con barras de error
lines(x, a+b*x, col="blue")

lines(c(xp - uxp, xp + uxp), c(yp, yp), lty = 1,col="red")
lines(c(xp - uxp, xp - uxp), c(yp - teey, yp + teey), lty = 1,col="red")
lines(c(xp + uxp, xp + uxp), c(yp - teey, yp + teey), lty = 1,col="red")
lines(c(xp, xp), c(yp - uyp, yp + uyp), lty = 1,col="red")
lines(c(xp - teex, xp + teex), c(yp - uyp, yp - uyp), lty = 1,col="red")
lines(c(xp - teex, xp + teex), c(yp + uyp, yp + uyp), lty = 1,col="red")

lines(c(xp + uxp, xp + uxp), c(0, yp + teey), lty = 3,col="red")
lines(c(xp - uxp, xp - uxp), c(0, yp + teey), lty = 3,col="red")

arrows(0,yp,xp - 1.2*uxp,yp,col="pink",length = 0.1)

lab<-paste0(round(a,4),'+',paste0(round(b,4),'x'))
text(mean(x),max(y) ,  lab )
#puntual
teex <- (max(x) - min(x)) / 100
teey <- (max(y) - min(y)) / 100
lines(c(xp - sqrt(u2xp), xp + sqrt(u2xp)), c(yp, yp), lty = 1)
lines(c(xp - sqrt(u2xp), xp - sqrt(u2xp)), c(yp - teey, yp + teey), lty = 1)
lines(c(xp + sqrt(u2xp), xp + sqrt(u2xp)), c(yp - teey, yp + teey), lty = 1)
lines(c(xp, xp), c(yp - uyp, yp + uyp), lty = 1)
lines(c(xp - teex, xp + teex), c(yp - uyp, yp - uyp), lty = 1)
lines(c(xp - teex, xp + teex), c(yp + uyp,yp + uyp), lty = 1)




## Evaluación directa

xf <- 3.5
uxf <- 0.2
a<-rslt$coefficients[1]
b<-rslt$coefficients[2]

yf <- a + b*xf
ca <- 1
cb <- xf
cx <- b
u2yf <- ca^2*u2a + cb^2*u2b + 2*ca*cb*uab + cx^2*uxf^2


yf +c(-1,1)*u2yf


rnumbers<-rnorm(3,0.1,0.02)
shapiro.test(rnumbers)
plot(rnumbers)
abline(h=0.1)

