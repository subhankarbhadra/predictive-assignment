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
write.csv(A, "~/all-het-output(AC).csv")

