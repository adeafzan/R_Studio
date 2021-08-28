#Matriks/Matrix/Matrices on R

#Assalamu'alaikum Warahmatullahi Wabarakaatuh

#Perkalian Matriks
x <- matrix(c(2,-1,1,-2,3,-2,-4,4,-3),nrow = 3, ncol = 3)
x
xtrc <- sum(diag(x))
xtrc
xx <- x%*%x
xx

#Invers Matriks
xinv <- solve(x) #Tidak punya invesrs krn merupakan matriks singular
xinv
xdet <- det(x) #singular matriks --> det == 0
xdet

y <- matrix(c(2,4,0,-1),nrow = 2, ncol = 2)
y
print(det(y))
print(solve(y))
ytran <- t(y)
ytraninv <- solve(ytran)
ytraninv
yinvtran <- t(solve(y))
yinvtran == ytraninv  ### Trans(inv(y)) = Inv (trans(y))

#Eigen Value & Eigen Vector dari Matriks
A <- matrix(c(1,-2,1,4),nrow = 2,ncol = 2)
A
eiA <- eigen(A, only.values = TRUE)
eiA
eiA <- eigen(A)   #menghasilkan eigen value dan eigen vector
eiA

#Singular Value Decompotition (SVD)
A <- matrix(c(1,2,2,2,-2,2), nrow = 3, ncol = 2)
A
As <- svd(A) #Hasil --> Matriks D (Eigen Value), U, V(Eigen Vector)
As
Ab <- As$u %*% diag(As$d) %*% t(As$v) # A=U*Diag(D)*Vt
Ab
Ab == A #Tidak sama persis krn pembulatan

#Cukup disini dulu, sampai jumpa di catatan berikutnya
#Selamat mencatat dan membaca
#Salam sehat selalu
#Wassalamu"alaikum warahmatullahi wabarakaatuh




