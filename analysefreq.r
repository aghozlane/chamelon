#freq
data=as.matrix(read.table("a_Frequences.txt",header=T))
A=gl(20,1,20,labels=c("L","A","G","V","E","S","D","T","K","I","R","P","N","F","Q","Y","H","M","C","W"))
rownames(data) =levels(A)
barplot(data[,1],names.arg=levels(A))
png("PDB2008NRfreq.png", width = 800, height = 500) 
barplot(t(data),beside=T,ylim=c(0,15),legend.text=colnames(data),col=c("red","blue","green", "yellow"),ylab="Fréquence(%)", xlab="Acides Aminés", main="Fréquence des acides aminés et de leurs structures secondaires dans le jeu de données non redondant\n(PDB 2008)")
dev.off()

#propensions
data2=as.matrix(read.table("a_Propension.txt",header=T))
A=gl(20,1,20,labels=c("L","A","G","V","E","S","D","T","K","I","R","P","N","F","Q","Y","H","M","C","W"))
rownames(data2) =levels(A)
barplot((data2[,1]-1),names.arg=levels(A))
png("PDB2008NRprop.png", width = 800, height = 500)
barplot(t(data2),beside=T,ylim=c(-1,1),legend.text=colnames(data2),col=c("blue","green","yellow"),ylab="Propension P(aa s2DB)", xlab="Acides Aminés", main="Proprension des acides aminés en hélice, boucle et en feuillet dans le jeu de données non redondant\n(PDB 2008)")
dev.off()

#freqglobales
data3=as.matrix(read.table("a_FrequencesGlobales.txt",header=T))
A=gl(4,1,4,labels=c("1","2","3","4"))
png("PDB2008NRfreqglobal.png", width = 800, height = 500)
barplot(t(data3),ylim=c(0,60),beside=T,names.arg=levels(A),col=c("grey10","grey30","grey50", "grey80"),legend.text=colnames(data3), ylab="Fréquence (%)", main="Fréquence des structures dans le jeu de donnéess non redondant\n(PDB 2008)")
dev.off()

#distribution
data4=as.matrix(read.table("a_Distribution.txt",header=T))
png("PDB2008NRdistribution.png", width = 600, height = 500)
hist(data4, main="Distribution de la longueur des séquences dans le jeu de donnéess non redondant\n(PDB 2008)", xlab="Longueur des séquences (Nombre d'acides aminés)",breaks=c(hist4$breaks))
dev.off()

#CAMELEON :)
data4=as.matrix(read.table("Cameleon.txt",header=T))
A=gl(5,1,5,labels=c("4","5","6","7","8"))
png("Cameleon1.png", width = 800, height = 500)
barplot(t(data4),beside=T,names.arg=levels(A),col=c("blue","green"),legend.text=colnames(data4), ylab="Nombre d'occurence", main="Nombre de séquences caméléons identifiés (Assignation de structure secondaire réalisée par DSSP)")
dev.off()

data4=as.matrix(read.table("Cameleon.txt",header=T))
A=gl(5,1,5,labels=c("4","5","6","7","8"))
png("Cameleon2.png", width = 800, height = 500)
barplot(t(data4),beside=T,names.arg=levels(A),col=c("blue","green"),legend.text=colnames(data4),ylim=c(0,5000), ylab="Nombre d'occurence", main="Nombre de séquences caméléons identifiés (Assignation de structure secondaire réalisée par DSSP)")
dev.off()

#Probabilite
data5=as.matrix(read.table("a_SeqCameleon_4_lettres_Probabilite.txt",header=T))
png("PDB2008NRPROBCameleon.png", width = 12000, height = 500)
barplot(t(data5),beside=T,col=c("red","blue","yellow"),legend.text=colnames(data5),ylab="Probabilite de conformation", xlab="Séquences cameleon", main="Probabilité des séquences caméléon à adopter une conformation en hélice, feuillet et boucle dans le jeu de données non redondant (Assignation de structure secondaire réalisée par DSSP)")
dev.off()

data5=as.matrix(read.table("a_SeqCameleon_4_lettres_Probabilite.txt",header=T))
png("PDB2008NRPROBCameleonHelice.png", width = 5000, height = 500)
barplot(t(data5[,1]),col="red",legend.text=colnames(data5[,1]),ylab="Probabilite de conformation", xlab="Séquences cameleon", main="Probabilité des séquences caméléon à adopter une conformation en hélice, feuillet et boucle dans le jeu de données non redondant (Assignation de structure secondaire réalisée par DSSP)")
dev.off()

data5=as.matrix(read.table("a_SeqCameleon_4_lettres_Probabilite.txt",header=T))
png("PDB2008NRPROBCameleonFeuillet.png", width = 5000, height = 500)
barplot(t(data5[,2]),col="blue",legend.text=colnames(data5[,1]),ylab="Probabilite de conformation", xlab="Séquences cameleon", main="Probabilité des séquences caméléon à adopter une conformation en feuillet dans le jeu de données non redondant (Assignation de structure secondaire réalisée par DSSP)")
dev.off()

data5=as.matrix(read.table("a_SeqCameleon_4_lettres_Probabilite.txt",header=T))
png("PDB2008NRPROBCameleonCoil.png", width = 5000, height = 500)
barplot(t(data5[,2]),col="green",legend.text=colnames(data5[,1]),ylab="Probabilite de conformation", xlab="Séquences cameleon", main="Probabilité des séquences caméléon à adopter une conformation en boucle dans le jeu de données non redondant (Assignation de structure secondaire réalisée par DSSP)")
dev.off()
#V2
#4lettres
data5=as.matrix(read.table("a_SeqCameleon_4_lettres_Probabilite.txt",header=T))
png("PDB2008NRPROBCameleonHelice4lettres.png", width = 800, height = 500)
hist(t(data5[,1]),col="red",xlab="Probabilité de conformation en hélice",  ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 4) à adopter une conformation en hélice\ndans le jeu de données non redondant (PDB 2008 - Assignation de structure secondaire réalisée par DSSP)",breaks=seq(0.1,2.5,0.01), ylim=c(0,1000),xlim=c(0,2))
dev.off()

data5=as.matrix(read.table("a_SeqCameleon_4_lettres_Probabilite.txt",header=T))
png("PDB2008NRPROBCameleonFeuillet4lettres.png", width = 800, height = 500)
hist(t(data5[,2]),col="blue",xlab="Probabilité de conformation en feuillet", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 4) à adopter une conformation en feuillet\ndans le jeu de données non redondant (PDB 2008 - Assignation de structure secondaire réalisée par DSSP)",breaks=seq(0.5,1.5,0.01), ylim=c(0,1000),xlim=c(0.8,1.2))
dev.off()

data5=as.matrix(read.table("a_SeqCameleon_4_lettres_Probabilite.txt",header=T))
png("PDB2008NRPROBCameleonCoil4lettres.png", width = 800, height = 500)
hist(t(data5[,3]),col="green",xlab="Probabilité de conformation en boucle", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 4) à adopter une conformation en boucle\ndans le jeu de données non redondant (PDB 2008 - Assignation de structure secondaire réalisée par DSSP)",breaks=seq(0.5,1.5,0.01), ylim=c(0,1000),xlim=c(0.8,1.2))
dev.off()

data5=as.matrix(read.table("a_SeqCameleon_4_lettres_Probabilite.txt",header=T))
png("PDB2008NRPROBCameleonplot4lettres.png", width = 800, height = 500)
plot((data5[,1]),(data5[,2]), xlab="Probabilité de conformation en hélice", ylab="Probabilité de conformation en feuillet", main="Probabilité de conformation en hélice en fonction de la probabilité de conformation en feuillet\ndes séquences caméléon de longueur 4 (PDB non redondant 2008)",xlim=c(0.90,1.10),ylim=c(0.5,1.7))
dev.off()

#5lettres
data5=as.matrix(read.table("a_SeqCameleon_5_lettres_Probabilite.txt",header=T))
png("PDB2008NRPROBCameleonHelice5lettres.png", width = 800, height = 500)
hist(t(data5[,1]),col="red",xlab="Probabilité de conformation en hélice",  ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 5) à adopter une conformation en hélice\ndans le jeu de données non redondant (PDB 2008 - Assignation de structure secondaire réalisée par DSSP)",breaks=seq(0.5,1.5,0.01), ylim=c(0,80),xlim=c(0.8,1.2))
dev.off()

data5=as.matrix(read.table("a_SeqCameleon_5_lettres_Probabilite.txt",header=T))
png("PDB2008NRPROBCameleonFeuillet5lettres.png", width = 800, height = 500)
hist(t(data5[,2]),col="blue",xlab="Probabilité de conformation en feuillet", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 5) à adopter une conformation en feuillet\ndans le jeu de données non redondant (PDB 2008 - Assignation de structure secondaire réalisée par DSSP)",breaks=seq(0.5,1.5,0.01), ylim=c(0,80),xlim=c(0.8,1.2))
dev.off()

data5=as.matrix(read.table("a_SeqCameleon_5_lettres_Probabilite.txt",header=T))
png("PDB2008NRPROBCameleonCoil5lettres.png", width = 800, height = 500)
hist(t(data5[,3]),col="green",xlab="Probabilité de conformation en boucle", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 5) à adopter une conformation en boucle\ndans le jeu de données non redondant (PDB 2008 - Assignation de structure secondaire réalisée par DSSP)",breaks=seq(0.5,1.5,0.01), ylim=c(0,80),xlim=c(0.8,1.2))
dev.off()


data5=as.matrix(read.table("a_SeqCameleon_5_lettres_Probabilite.txt",header=T))
png("PDB2008NRPROBCameleonplot5lettres.png", width = 800, height = 500)
plot((data5[,1]),(data5[,2]), xlab="Probabilité de conformation en hélice", ylab="Probabilité de conformation en feuillet", main="Probabilité de conformation en hélice en fonction de la probabilité de conformation en feuillet\ndes séquences caméléon de longueur 5 (PDB non redondant 2008)",xlim=c(0.90,1.10),ylim=c(0.5,1.7))
dev.off()

















