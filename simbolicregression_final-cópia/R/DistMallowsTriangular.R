#Defini??o de constantes e par?metros
eps<-1e-6
#Fim - Defini??o de constantes e par?metros



###Defini??o de Fun??es

safesqrt <- function(x, prec=1e-6) {
  if (length(x) == 1) return (safesqrt0(x,prec))
  if (is.null(dim(x))) return (sapply(x,safesqrt0,prec=prec))
  apply(x,1:lenght(dim(x)),safesqrt0,prec=prec)
}


safesqrt0 <- function(x, prec=1e-6) {
  if (x >= 0.) return (sqrt(x))
  if (abs(x) < prec) return(0.)
  if (abs(x) >= prec) stop(cat("Negative argument",x,"for safesqrt0\n"))
}


safearcsin <- function(x, prec=1e-6) {
  if (abs(x) <= 1.) return(asin(x))
  if ((abs(x)-1.) > prec) stop(cat("Invalid argument for safearcsin:",x,"\n"))
  if (x > 1.) return(asin(1.))
  if (x < -1.) return(asin(-1.))
}


DistMallows2<-function(c1,r1,m1,c2,r2,m2)
  #  Fun??o que calcula o quadrado da dist?ncia de Mallows entre intervalos com distribui??o triangular geral

  #  Argumentos:
  #     c1 -- centro do intervalo 1
  #     r1 -- raio do intervalo 1
  #     m1 -- moda do intervalo 1
  #     c2 -- centro do intervalo 2
  #     r2 -- raio do intervalo 2
  #     m2 -- moda do intervalo 2

{ #c1 <- round(c1,12)
  #r1 <- round(r1,12)
  #m1 <- round(m1,12)
  #c2 <- round(c2,12)
  #r2 <- round(r2,12)
  #Ym2 <- round(m2,12)

  aux11<-0
  aux22<-0
  if (abs(r1)<eps) {
    if (abs(r2)<eps) distM2<-(c1-c2)^2
    else if (abs(r2)>=eps) distM2<-(c1-c2)^2-r2^2/3-4/3*(m2-c2)^2+2/3*(m2-c2)*(c2-c1)+(m2-c2+r2)^3/(4*r2)+(-m2+c2+r2)^3/(4*r2)
  }
  else if (abs(r1)>=eps) {
    if (abs(r2)<eps) distM2<-(c1-c2)^2-r1^2/3-4/3*(m1-c1)^2+2/3*(m1-c1)*(c1-c2)+(m1-c1+r1)^3/(4*r1)+(-m1+c1+r1)^3/(4*r1)
    else if (abs(r2)>=eps)  {
      aux11<-(m1-c1)/r1
      aux22<-(m2-c2)/r2
      z1<-safesqrt(m1-c1+r1)
      s1<-safesqrt(c1+r1-m1)
      z2<-safesqrt(m2-c2+r2)
      s2<-safesqrt(c2+r2-m2)
      sr<-safesqrt(r1*r2)
      if ((aux11+1)/2<=(aux22+1)/2) {
        distM2<-(c1-c2)^2+(r1-r2)^2/6+(m1-c1)^2/6+(m2-c2)^2/6-5/3*r1*r2+2/3*(m1-c1)*(c1-c2+r2)-2/3*(m2-c2)*(c1-c2+r1)+sr*z1*z2*(5-aux11)/6+sr*s1*s2*(5+aux22)/6+sr*s1*z2*(safearcsin(aux22)-safearcsin(aux11))/2
      }
      else if ((aux11+1)/2>(aux22+1)/2) {
        distM2<-(c1-c2)^2+(r1-r2)^2/6+(m1-c1)^2/6+(m2-c2)^2/6-5/3*r1*r2+2/3*(m1-c1)*(c1-c2-r2)-2/3*(m2-c2)*(c1-c2-r1)+sr*z1*z2*(5-aux22)/6+sr*s1*s2*(5+aux11)/6+sr*s2*z1*(safearcsin(aux11)-safearcsin(aux22))/2
      }
    }
  }
  return(distM2)
}
#Fim - DistMallows2


SumDist<-function(par,m,nvar,mcv,mrv,mmv,mlambdas)
  #  SumDist -- Soma dos Quadrados da Dist?ncia de Mallows
  #  Fun??o que calcula a soma dos quadrados da dist?ncia de Mallows entre as vari?veis e as combina??es lineares dos fatores

  #  Argumentos:
  #     par -- Vetor dos par?metros: li, dlim, dlsm
  #           li   -- Limites inferiores dos intervalos dos fatores
  #           dlim -- Dist?ncia dos limites inferiores dos intervalos dos fatores ?s modas
  #           dlsm -- Dist?ncia dos limites superiores dos intervalos dos fatores ?s modas
  #     m        -- N?mero de fatores
  #     nvar     -- N?mero de vari?veis
  #     mcv      -- Matriz dos centros das vari?veis
#     mrv      -- Matriz dos raios das vari?veis
#     mmv      -- Matriz das modas das vari?veis
#     mlambdas -- Matriz dos coeficientes dos factores

{
  li<-par[1:m]
  dlim<-par[(m+1):(2*m)]
  dlsm<-par[(2*m+1):(3*m)]
  cf<-c(rep(0,m))
  rf<-c(rep(0,m))
  mf<-c(rep(0,m))
  cclf<-c(rep(0,nvar))
  rclf<-c(rep(0,nvar))
  mcrf<-c(rep(0,nvar))
  cf<-li+(dlim+dlsm)/2        # Vetor dos centros dos factores
  rf<-(dlim+dlsm)/2           # Vetor dos raios dos factores
  mf<-li+dlim                 # Vetor das modas dos factores
  cclf<-mlambdas%*%cf         # Vetor dos centros das combinacoes lineares de factores
  rclf<-abs(mlambdas)%*%rf    # Vetor dos raios das combinacoes lineares de factores
  mclf<-mlambdas%*%mf         # Vetor das modas das combinacoes lineares de factores
  distM<-0
  if (op_met==1)
    for (j in 1:nvar) distM <- distM + DistMallows2(mcv[j],mrv[j],mmv[j],cclf[j],rclf[j],mclf[j])
  else
    for (j in 1:nvar)
      distM <- distM + DistMallows2(mcv[j],mrv[j],mmv[j],cclf[j],rclf[j],mclf[j])/varerros[j]
  return(distM)
}
#Fim - SumDist

###Fim - Defini??o de Fun??es
