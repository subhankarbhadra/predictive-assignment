# Code for cleaning data from the raw version available at: https://github.com/subhankarbhadra/predictive-assignment/tree/main/DBLP_four_area

library(Matrix)

a = read.csv("DBLP_four_area/csv/author_label.csv", header = F, col.names = c("a.id","a.lab","a.name"))
a <- a[order(a$a.id),]
a$a.id = as.factor(a$a.id); a$a.lab = as.factor(a$a.lab);
a.list = levels(a$a.id)

ap = read.csv("DBLP_four_area/csv/paper_author.csv", header = F, col.names = c("p.id","a.id"))
pc = read.csv("DBLP_four_area/csv/paper_conf.csv", header = F, col.names = c("p.id","c.id"))
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

#Adjacency matrices
write.csv(levels(pac1$a.id), "DBLP_four_area/csv/a-levels.csv")
write.csv(levels(pac1$conf), "DBLP_four_area/csv/c-levels.csv")
write.csv(levels(pac2$p.id), "DBLP_four_area/csv/p-levels.csv")
levels(pac1$a.id) <- 1:n.a
levels(pac1$conf) <- 1:n.c
levels(pac2$a.id) <- 1:n.a
levels(pac2$p.id) <- 1:n.p

AC = matrix(0, nrow = n.a, ncol = n.c)
for (i in 1:n.a){
  id = unique(subset(pac1, a.id == i)$conf)
  AC[i, id] = 1
}

ACA = matrix(0, nrow = n.a, ncol = n.a)
for (i in 1:(n.a-1)){
  for (j in (i+1):n.a){
    if (sum(AC[i,]*AC[j,]) > 0) {ACA[i,j] = 1; ACA[j,i] = 1}
  }}

n <- nrow(ACA)
classes <- as.numeric(a$a.lab)

saveRDS(list(adjacency = as(ACA, "dgCMatrix"), classes = classes), "dblp.RData")
