install.packages("qtl")
library(qtl)

#read genotype data
cross <- read.cross(format = "csvs", F.gen = 8,
                    genfile = "genotype_115markers_144lines.csv", 
                    phefile = "pheno_46lines.csv")
summary(cross)
cross$geno$`1`$map


#construct linkage map
map <- est.map(cross, map.function="kosambi")
plotMap(map, show.marker.names=F)


#add linkage map to the model
cross <- replace.map(cross, map)
cross$geno$`1`$map
cross$geno$`1`$data[1:40, ]


#format to RIL
cross <- convert2riself(cross)
summary(cross)
cross$geno$`1`$data[1:40, ]


#simulate pseudo markers
cross <- calc.genoprob(cross, step = 2, map.function = "kosambi")
cross <- sim.geno(cross, step = 2, n.draws = 100, map.function = "kosambi")
cross$geno$`1`$prob


#Simple Interval Mapping (Your data)
sim <- scanone(cross, pheno.col = 2, method = "hk")
sim.perm <- scanone(cross, pheno.col = 2, method = "hk", n.perm = 100)
plot(sim)
abline(h = summary(sim.perm, alpha = 0.05))


#Simple Interval Mapping (Previous data)
sim <- scanone(cross, pheno.col = 4, method = "hk")
sim.perm <- scanone(cross, pheno.col = 4, method = "hk", n.perm = 100)
plot(sim)
abline(h = summary(sim.perm, alpha = 0.05))


#Simple Interval Mapping (Complete data)
sim <- scanone(cross, pheno.col = 5, method = "hk")
sim.perm <- scanone(cross, pheno.col = 5, method = "hk", n.perm = 100)
plot(sim)
abline(h = summary(sim.perm, alpha = 0.05))
#Focus on chr.6
plot(sim, chr = 6)
abline(h = summary(sim.perm, alpha = 0.05))


#Composite interval mapping
cim <- cim(cross, pheno.col = 5, method = "hk", n.marcovar = 2, window = 20)
cim.perm <- cim(cross, pheno.col = 5, method = "hk", n.marcovar = 2,
                   window = 20, n.perm = 100)
plot(cim)
abline(h = summary(cim.perm, alpha = 0.05))
#Focus on chr.6
plot(cim, chr = 6)
abline(h = summary(cim.perm, alpha = 0.05))


#Calculate effects of QTLs
summary(cim)
cim[cim$lod > as.numeric(summary(cim.perm, alpha = 0.05)), ]
qtl <- makeqtl(cross, chr = 6, pos = 42.59, what = "prob")
res <- fitqtl(cross, pheno.col = 5, qtl = qtl, get.ests = T, method="hk")
summary(res)



###Additional figures

#Compare your data and previous data
data <- read.csv("pheno_46lines.csv")
plot(x = data$Absorption, y = data$Absorption_2)
plot(x = data$Absorption, y = data$Absorption_3)
cor(data[, c(2,4,5)], use = "complete.obs")

#Association of the trait and Alk genotype
geno <- read.csv("genotype_115markers_144lines.csv")
geno.alk <- geno$Alk[-1]
data <- cbind(data, geno=as.factor(geno.alk))
list <- list(data$Absorption[which(data$geno == "A")], data$Absorption[which(data$geno == "B")])
boxplot(list)
list <- list(data$Absorption_2[which(data$geno == "A")], data$Absorption_2[which(data$geno == "B")])
boxplot(list)
list <- list(data$Absorption_3[which(data$geno == "A")], data$Absorption_3[which(data$geno == "B")])
boxplot(list)
