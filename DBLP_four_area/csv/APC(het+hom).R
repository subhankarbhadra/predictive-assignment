a = read.csv("~/csv/author_label.csv", header = F, col.names = c("a.id","a.lab","a.name"))
a <- a[order(a$a.id),]
a$a.id = as.factor(a$a.id); a$a.lab = as.factor(a$a.lab);
a.list = levels(a$a.id)

ap = read.csv("~/csv/paper_author.csv", header = F, col.names = c("p.id","a.id"))
pc = read.csv("~/csv/paper_conf.csv", header = F, col.names = c("p.id","c.id"))
conf = rep(0, dim(ap)[1])
for (i in 1:dim(ap)[1]){
id = ap[i,]$p.id
conf[i] = subset(pc, p.id == id)$c.id
}

pac = cbind(ap, conf)
pac1 = subset(pac, a.id %in% a.list, select = c(a.id, conf))
pac1$a.id = as.factor(pac1$a.id); pac1$conf = as.factor(pac1$conf)
n.a = length(levels(pac1$a.id)); n.c = length(levels(pac1$conf))

pac2 = subset(pac, a.id %in% a.list, select = c(a.id, p.id))
pac2$a.id = as.factor(pac2$a.id); pac2$p.id = as.factor(pac2$p.id)
n.p = length(levels(pac2$p.id))

pac3 = subset(pac, a.id %in% a.list, select = c(p.id, conf))
pac3 = unique(pac3)
pac3$p.id = as.factor(pac3$p.id); pac3$conf = as.factor(pac3$conf)

#Adjacency matrices
write.csv(levels(pac1$a.id), "~/csv/a-levels.csv")
write.csv(levels(pac1$conf), "~/csv/c-levels.csv")
write.csv(levels(pac2$p.id), "~/csv/p-levels.csv")
levels(pac1$a.id) <- 1:n.a
levels(pac1$conf) <- 1:n.c
levels(pac2$a.id) <- 1:n.a
levels(pac2$p.id) <- 1:n.p
levels(pac3$conf) <- 1:n.c
levels(pac3$p.id) <- 1:n.p

AP = matrix(0, nrow = n.a, ncol = n.p)
for (i in 1:n.p){
AP[subset(pac2, p.id == i)$a.id,i] = 1
}
APA = matrix(0, nrow = n.a, ncol = n.a)
for (i in 1:(n.a-1)){
for (j in (i+1):n.a){
if (t(AP[i,])%*%AP[j,] > 0) {APA[i,j] = 1; APA[j,i] = 1}
}}
AC = matrix(0, nrow = n.a, ncol = n.c)
for (i in 1:n.a){
id = unique(subset(pac1, a.id == i)$conf)
AC[i, id] = 1
}
ACA = matrix(0, nrow = n.a, ncol = n.a)
for (i in 1:(n.a-1)){
for (j in (i+1):n.a){
if (t(AC[i,])%*%AC[j,] > 0) {ACA[i,j] = 1; ACA[j,i] = 1}
}}
PC = matrix(0, nrow = n.p, ncol = n.c)
for (i in 1:n.p){
PC[i,subset(pac3, p.id == i)$conf] = 1
}

A = matrix(0, nrow = n.a+ n.p+ n.c, ncol = n.a+ n.p+ n.c)
A[1:n.a,(n.a+1):(n.a+n.p)] = AP
A[(n.a+1):(n.a+n.p),1:n.a] = t(AP)
A[(n.a+1):(n.a+n.p),(n.a+n.p+1):(n.a+n.p+n.c)] = PC
A[(n.a+n.p+1):(n.a+n.p+n.c),(n.a+1):(n.a+n.p)] = t(PC)

library("clue", lib ="~/R")
pMatrix.min <- function(A, B) { 
        n <- nrow(A) 
        D <- matrix(NA, n, n) 
        for (i in 1:n) { 
        for (j in 1:n) { 
        D[j, i] <- (sum((B[j, ] - A[i, ])^2)) 
        } } 
vec <- c(solve_LSAP(D)) 
list(A=A[vec,], pvec=vec) 
} 

require(clue)  # need this package to solve the LSAP 

n = ncol(A)
d = colSums(A)
L = matrix(0, nrow = n, ncol = n)
for (i in 1:(n-1)){
for (j in (i+1):n){
L[i,j] = A[i,j]/sqrt(d[i]*d[j]) 
L[j,i] = L[i,j]
}}
k = 4

# we need to pick the eigenvectors for top eigenvalues of L by abs value
val = abs(eigen(L, symmetric = TRUE)$values)
id = order(-val)
id_vec = id[1:(2*k)]
X = eigen(L)$vectors[,id_vec]
N1 = 4057; #CHANGE this
X1 = X[1:N1,]
K1 = kmeans(X1, centers = k, nstart = 10);K1$size
Err1 = K1$cluster

### Homogeneous version(s) ####
A1 = APA; d1 = colSums(A1); 
null1 = which(d1 == 0); 

# do the next patch only when there are 0 values in d1
if (length(null1) > 0 ) {
A1 = A1[-null1,-null1]; 
d1 = d1[-null1]; 
}
N11 = N1 - length(null1)
L1 = matrix(0, nrow = N11, ncol = N11)
for (i in 1:(N11-1)){
for (j in (i+1):N11){
L1[i,j] = A1[i,j]/sqrt(d1[i]*d1[j]) 
L1[j,i] = L1[i,j]
}}
val1 = abs(eigen(L1, symmetric = TRUE)$values)
id1 = order(-val1)
id_vec1 = id1[1:k]
X3 = eigen(L1)$vectors[,id_vec1]	
K3 = kmeans(X3, centers = k, nstart = 10);K3$size
Err3 = K3$cluster


A2 = ACA; d2 = colSums(A2); 
null2 = which(d2 == 0); 

# do the next patch only when there are 0 values in d1 and d2
if (length(null2) > 0 ) {
A2 = A2[-null2,-null2]; 
d2 = d2[-null2]; 
}
N12 = N1 - length(null2)
L2 = matrix(0, nrow = N12, ncol = N12)
for (i in 1:(N12-1)){
for (j in (i+1):N12){
L2[i,j] = A2[i,j]/sqrt(d2[i]*d2[j]) 
L2[j,i] = L2[i,j]
}}
val2 = abs(eigen(L2, symmetric = TRUE)$values)
id2 = order(-val2)
id_vec2 = id2[1:k]
X4 = eigen(L2)$vectors[,id_vec2]	
K4 = kmeans(X4, centers = k, nstart = 10);K4$size
Err4 = K4$cluster

### OUTPUT ANALYSIS #######
C <- a$a.lab
levels(C) <- 1:k

N_err = 0
E = factor(x = Err1, levels = 1:k)
A = table(E,C)
X <- pMatrix.min(A,diag(1,k)) 
A = X$A
N_err_het1 = sum(A) - sum(diag(A))
N_err = N_err_het1
P = N_err/N1 
write.csv(A, "~/all-het-output(APC).csv")

N_err = 0 
if (length(null1) > 0){C1 = C[-null1]
} else {C1 = C} 
E = factor(x = Err3, levels = 1:k)
A1 = table(E,C1)
X <- pMatrix.min(A1,diag(1,k)) 
A1 = X$A
N_err_hom1 = sum(A1) - sum(diag(A1))
N_err = N_err_hom1 + ((k-1)/k)*length(null1)
P1 = N_err/N1
write.csv(A1, "~/all-hom1-output(APA).csv")

N_err = 0 
if (length(null2) > 0){C2 = C[-null2]
} else {C2 = C} 
E = factor(x = Err4, levels = 1:k)
A2 = table(E,C2)
X <- pMatrix.min(A2,diag(1,k)) 
A2 = X$A
N_err_hom2 = sum(A2) - sum(diag(A2))
N_err = N_err_hom2 + ((k-1)/k)*length(null2)
P2 = N_err/N1
write.csv(A2, "~/all-hom2-output(ACA).csv")