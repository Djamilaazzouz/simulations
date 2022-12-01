#MINI PROJET 1 : simulations et copules
#On cherche à calculer une fonction en dimension 2 qu'on note phi.

#methode 1
library(rgl)
n=1000
N= 200 
phi= function(x,y)(x^2)*exp(-y) #la fonction qu'on a choisi d'éstimer par MC
X=seq(0,2,length=1100) #intervalle dans 0 et 2 avec 1100 de pixel 
Y=X
z=outer(X,Y,phi)
persp3d(X,Y,z, col = "red") #la version 3d de la fonction phi
library(pracma)
I=quad2d(phi,0,2,0,2) #calcul de la double integral 
print(I)
U= runif(n,0,2)
V= runif(n,0,2)
W= runif(n,0,4)
I_chap= (4*max(z))*mean(W<= phi(U,V))  #max= 2X4
evol1=(4*max(z))*cumsum(W<= phi(U,V))/(1:n) 
#Evolution de l'estimateur I1_chap:
plot(evol1,lwd='2',type = 'l', col='blue', main="evolution de l'estimateur I_chap")
abline(h=I, col='red') # I pour référence par rapport à la convergence

#calcul de l'estimateur
# de N échantillonnage
M_1= matrix(runif(N*n,0,2),ncol=N)   # taille n et echantillon N
M_2= matrix(runif(N*n,0,2),ncol=N)
M_3= matrix(runif(N*n,0,4),ncol=N)
M_12=phi(M_1,M_2)
I1_chap= (4*max(z))*M_3[,1]<= phi(M_1[,1],M_2[,1])
evol1=(4*max(z))*cumsum(M_3[,1]<= phi(M_1[,1],M_2[,1]))/(1:n)#convergence de I1_chap le premier estimateur
M=(M_3<= phi(M_1,M_2))
I1_vect= (4*max(z))*apply(M, MARGIN=2, FUN=mean)
hist(I1_vect, main="Histogramme N estimation de I", col='green')
#on vérifie la qualité de notre estimateur en calculant l'équart type
ecart_type1=sd(I1_vect)
print(ecart_type1)

####################################################################################

#methode 2
I2_chap= 4*phi(M_1[,1],M_2[,1])#avec la densité de la uniforme du couple (u,v)=1/4
evol2= 4*cumsum(phi(M_1[,1],M_2[,1]))/(1:n)#convergence du deuxième estimateur I2_chap
#Evolution de I2_chap:
plot(evol2,type='l',lwd=2,col="blue")
abline(h=I)
#Supperposition des deux estimateurs
plot(evol2,lwd=2,lty=1,type='l',col="blue")#on montre la convergence de I2_chap vers I 
lines(evol1,lwd=2,lty=1,type='l',col='purple')#on compare la convergence des deux estimateurs vers I
abline(h=I,lwd=2,col="red")
legend(x="topright", legend=c("I2chap","I1chap","I"),
       lty=c(1,1), col=c("blue","purple","red"), lwd=2)
#Histogramme de N estimation de I
I2_vect=4*apply(M_12, MARGIN=2, FUN=mean)
hist(I2_vect, main="Histogramme N estimation de I", col='blue')
#observation: on remarque que I2_chap converge plus rapidement vers I que I1_chap
#Calcul de l'écart type:
ecart_type2=sd(I2_vect)
print(ecart_type2)#en effet on voit que la valeur de l'écart-type de I2_chap est plus petite que celle de I1_chap


################################################################################

#methode 3
#simulation du couple de (Z1,Z2) 
#cherchons d'abord p(x)et p(y) tel que p(x)et p(y) sont independant et différent de 0
#le choix des p est basé sur le fait que p(x)*p(y) doit ressembler graphiquement à notre fonction phi
px=function(x)(1/(exp(2)-1))*exp(x) #densité de Z1
py=function(y)(1/2)*y #densité de Z2
#par la méthode d'inversion on va simuler Z1 et Z2
inv_px=function(n)log((exp(2)-1)* runif(n)+1)#la fonction inverse de Z1
inv_py=function(n)2*sqrt(runif(n))#la fonction inverse de Z2
dp=function(x,y)px(x)*py(y)
#on vérifie que p(x)*p(y) ressemble à phi
X=seq(0,2, length=1100) 
z=outer(X,Y,dp)
persp3d(X, Y, z,col='yellow')

#création matrice pour simuler les Z1_i et Z2_i pour tout i=[1:n]
M3= matrix(inv_px(N*n),ncol=N)
M4= matrix(inv_py(N*n),ncol=N)
I3_chap= phi(M3[,1],M4[,1])/(px(M3[,1])*py(M4[,1]))
evol3= cumsum(phi(M3[,1],M4[,1])/(px(M3[,1])*py(M4[,1])))/(1:n)
#Evolution de I3_chap:
#plot(evol3,type='l',lwd=2,col="purple")
#abline(h=I)
#Construction de l'histogramme de N estimation de Ichap3
calc3=function(Z1,Z2) mean(phi(Z1,Z2)/(px(Z1)*py(Z2)))
I3_vect=apply(M3,M4, MARGIN=2, FUN=calc3)
#Superposition des 3 estimateurs:
plot(evol3,lwd=2,lty=1,type='l',col="blue")
lines(evol2,lwd=2,lty=1,type='l',col='red')
lines(evol1,lwd=2,lty=1,type='l',col='purple')
abline(h=I,lwd=2,col="black")
legend(x="topright", legend=c("I3chap","I2chap", "I1chap","I"),
lty=c(1,1,1), col=c("blue","red","purple","black"), lwd=2)
#comparaison des 3 histogrammes 
par(mfrow=c(2,2))
hist(I1_vect, col = "pink")
hist(I2_vect, col = "green" )
hist(I3_vect, col = "aquamarine")
#observation: on remarque que l'histogramme associé à l'éstimateur I3_chap est plus centré que les autres.
#En effet, d'apres le TCL, l'estimateur I3_chap converge vers la loi normale centrée réduite
#plus que les deux premiers estimateurs

#comparaison des écart-types pour les trois méthodes
print(ecart_type1)
print(ecart_type2)
ecart_type3= sd(I3_vect)
print(ecart_type3)
#le but de la methode de Monte Carlo est d'approcher une intégrale, en utilisant les trois méthodes présentées au-dessus
#nous cherchons à trouver le meilleur estimateur de cette intégrale en minimisant l'écart-type.
#on conclut que la methode 3 est la plus efficace car son écart-type est plus petit que les deux autres.

#####################################################################################
#Méthode 4:
theta=1/n
F_clay=function(U,Z,theta)
  ((Z^(-theta/(theta+1))-1)*(U^-theta)+1)^(-1/theta)
rcop_clay=function(n,theta){
  U=runif(n)
  Z=runif(n)
  V=F_clay(U,Z,theta)
  return(data.frame(U,V))
}
#Tout d'abord nous allons faire en sorte d'avoir plus de points en (2,0)
#On pose donc  (U1,V1)=(1-U1,V2).
#On note (1-U1,V2)=(U2,V2)
#on simule les marginales  selon la copule de clayton
rcop=rcop_clay(n,theta)
marge1=1-rcop$U #marginale de 1-U
marge2=rcop$V#marginale de V
plot(marge1,marge2, col="brown")
#transformer les marginales  simulées sur [0,2]
U2=qunif(marge1,0,2) # U2=F^(-1)(1-U)
V2=qunif(marge2,0,2) # U2=G^(-1)(V)
plot(U2,V2,col="purple")
#La densité de Clayton suivant le couple  (1-U,V) est la suivante:
dcop_clayton=function(U,V)
  (1+theta)*((1-U)*V)^((-theta)-1)*(((1-U)^(-theta)+V^(-theta)-1)^((-1/theta)-2))
#Présentation en 3d de dcop_clayton:
x=seq(0,1,length=1100)
y=x
z=outer(x,y,dcop_clayton)
persp3d(x,y,z, col="green")
#Calcul de l'intégrale de cette couple:
quad2d(dcop_clayton,0,1,0,1)
#Observation :
#l'intégrale de cette copule vaut bien 1 donc la copule définie au dessus est 
#bien une denisté
##################################################################################
#Ecrivons la fonction de répartition du couple (U2,V2) selon clayton
#on sait que d'après le théorème de Sklar on a :
#H(x,y)=C(F(x),G(y)) or F(x)=1-U et G(y)=V
#La densité du couple (U2,V2) est donc la dérivée de H 
#h(x,y)=c(F(x),G(y))*f(x)*g(y) où f(x) et g(y) sont les densités respectives
#de U2 et V2.
f=function(x)dunif(x,0,2)
K=function(y)punif(y,0,2)
g=function(y)dunif(y,0,2)
G=function(x)punif(x,0,2)
#densité du couple:
h=function(x,y) f(x)*g(y)*dcop_clayton(K(y),G(x))
#Visualisation de la densité h en 3D 
a=seq(0,2,length=50) 
b=a
c=outer(a,b,h)
persp3d(a,b,c, col = "red")
#calcul de l'intégrale de h
quad2d(h,0,2,0,2)
#Obervation: 
#h est bien une densité.
#Constrution de l'estimateur de I
Ichap4=mean(phi(U2,V2)/h(U2,V2))
evol4= cumsum(phi(U2,V2)/h(U2,V2))/(1:n)
#Présentation graphique de l'évolution de l'estimateur Ichap4.
plot(evol4,type='l',lwd=2,col="red")
abline(h=I)
# Constrution de l'histogramme de N estimation De ichap4
rcop=rcop_clay(n*N,theta)
a1=1-rcop$U
b1=rcop$V
d=qunif(a1,0,2)
s=qunif(b1,0,2)
M10= matrix(phi(d,s)/h(d,s),nrow=n)
I4_vect= apply(M10, MARGIN=2, FUN=mean)
hist(I4_vect,col="red")
#Calcul de l'écart type de cette méthode
ecart_type4=sd(I4_vect)
print(ecart_type4)
#Superposition des 4 estimateurs 
plot(evol4,lwd=2,lty=1,type='l',col="blue")
lines(evol3,lwd=2,lty=1,type='l',col='red')
lines(evol2,lwd=2,lty=1,type='l',col='purple')
lines(evol1,lwd=2,lty=1,type='l',col="#009999")
abline(h=I,lwd=2,col="black")
legend(x="topright", legend=c("I4chap","I3chap", "I2chap","I1chap","I"),
lty=c(1,1,1,1), col=c("blue","red","purple","#009999","black"), lwd=2)
#Projection des histogrammes avec les 4 méthodes 
par(mfrow=c(2,2))
hist(I1_vect, col = "pink",main="méthode 1")
hist(I2_vect, col = "green",main="méthode 2" )
hist(I3_vect, col = "aquamarine",main="méthode 3")
hist(I4_vect, col= "orange",main="méthode 4")
#Superposition des 4 histogrammes 
hist(I1_vect,ylim=c(0,50), col = "pink",main="Superpositions des 3 histogrammes")
hist(I2_vect,add=T, col = "green" )
hist(I3_vect,add=T, col = "aquamarine")
hist(I4_vect,add=T ,col= "white")



