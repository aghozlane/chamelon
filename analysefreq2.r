#freq
data=as.matrix(read.table("all-jan-07clear_Frequences.txt",header=T))
A=gl(20,1,20,labels=c("L","A","G","V","E","S","D","T","K","I","R","P","N","F","Q","Y","H","M","C","W"))
rownames(data) =levels(A)
barplot(data[,1],names.arg=levels(A))
png("PDB2007Rfreq.png", width = 800, height = 500) 
barplot(t(data),beside=T,ylim=c(0,15),legend.text=colnames(data),col=c("red","blue","green", "yellow"),ylab="Fréquence(%)", xlab="Acides Aminés", main="Fréquence des acides aminés et de leurs structures secondaires dans le jeu de données redondant\n(PDB 2007)")
dev.off()

#propensions
data2=as.matrix(read.table("all-jan-07clear_Propension.txt",header=T))
A=gl(20,1,20,labels=c("L","A","G","V","E","S","D","T","K","I","R","P","N","F","Q","Y","H","M","C","W"))
rownames(data2) =levels(A)
data2[,1]=data2[,1]-1.0
data2[,2]=data2[,2]-1.0
data2[,3]=data2[,3]-1.0
barplot((data2[,1]),names.arg=levels(A))
png("PDB2007Rprop.png", width = 800, height = 500)
barplot(t(data2),beside=T,ylim=c(-1,1),legend.text=colnames(data2),col=c("blue","green","yellow"),ylab="Propension P(aa s2DB)", xlab="Acides Aminés", main="Proprension des acides aminés en hélice, boucle et en feuillet dans le jeu de données redondant\n(PDB 2007)")
dev.off()

#freqglobales
data3=as.matrix(read.table("all-jan-07clear_FrequencesGlobales.txt",header=T))
A=gl(4,1,4,labels=c("1","2","3","4"))
png("PDB2007Rfreqglobal.png", width = 800, height = 500)
barplot(t(data3),ylim=c(0,60),beside=T,names.arg=levels(A),col=c("grey10","grey30","grey50", "grey80"),legend.text=colnames(data3), ylab="Fréquence (%)", main="Fréquence des structures dans le jeu de données redondant\n(PDB 2007)")
dev.off()

#distribution
data4=as.matrix(read.table("all-jan-07clear_Distribution.txt",header=T))
png("PDB2007Rdistribution.png", width = 600, height = 500)
hist4=hist(data4, main="Distribution de la longueur des séquences dans le jeu de données redondant\n(PDB 2007)", xlab="Longueur des séquences (Nombre d'acides aminés)")
dev.off()

#4lettres
data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_4_lettres_Probabilite.txt",header=T))
png("PDB2007RPROBCameleonHelice4lettres.png", width = 800, height = 500)
hist((data5[,1]),col="red",xlab="Probabilité de conformation en hélice",  ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 4) à adopter une conformation en hélice\ndans le jeu de données redondant (PDB 2007 - Assignation de structure secondaire réalisée par DSSP)",ylim=c(0,1000),breaks=seq(0,4,0.01))
dev.off()

data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_4_lettres_Probabilite.txt",header=T))
png("PDB2007RPROBCameleonFeuillet4lettres.png", width = 800, height = 500)
hist((data5[,2]),col="blue",xlab="Probabilité de conformation en feuillet", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 4) à adopter une conformation en feuillet\ndans le jeu de données redondant (PDB 2007 - Assignation de structure secondaire réalisée par DSSP)",ylim=c(0,1000),breaks=seq(0,14,0.01),xlim=c(0,4))
dev.off()

data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_4_lettres_Probabilite.txt",header=T))
png("PDB2007RPROBCameleonCoil4lettres.png", width = 800, height = 500)
hist((data5[,3]),col="green",xlab="Probabilité de conformation en boucle", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 4) à adopter une conformation en boucle\ndans le jeu de données redondant (PDB 2007 - Assignation de structure secondaire réalisée par DSSP)",ylim=c(0,1000),breaks=seq(0,10,0.01),xlim=c(0,4))
dev.off()


data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_4_lettres_Probabilite.txt",header=T))
png("PDB2007RPROBCameleonplot4lettres.png", width = 800, height = 500)
plot((data5[,1]),(data5[,2]), xlab="Probabilité de conformation en hélice", ylab="Probabilité de conformation en feuillet", main="Probabilité de conformation en hélice en fonction de la probabilité de conformation en feuillet\ndes séquences caméléon de longueur 4 (PDB redondant 2007)",ty="p",cex=0.7, pch=18)
dev.off()

#5lettres
data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_5_lettres_Probabilite.txt",header=T))
png("PDB2007RPROBCameleonHelice5lettres.png", width = 800, height = 500)
hist(t(data5[,1]),col="red",xlab="Probabilité de conformation en hélice",  ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 5) à adopter une conformation en hélice\ndans le jeu de données redondant (PDB 2007 - Assignation de structure secondaire réalisée par DSSP)",ylim=c(0,100),breaks=seq(0,5,0.01),xlim=c(0,4))
dev.off()

data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_5_lettres_Probabilite.txt",header=T))
png("PDB2007RPROBCameleonFeuillet5lettres.png", width = 800, height = 500)
hist((data5[,2]),col="blue",xlab="Probabilité de conformation en feuillet", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 5) à adopter une conformation en feuillet\ndans le jeu de données redondant (PDB 2007 - Assignation de structure secondaire réalisée par DSSP)",ylim=c(0,100),breaks=seq(0,15,0.01),xlim=c(0,4))
dev.off()

data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_5_lettres_Probabilite.txt",header=T))
png("PDB2007RPROBCameleonCoil5lettres.png", width = 800, height = 500)
hist((data5[,3]),col="green",xlab="Probabilité de conformation en boucle", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 5) à adopter une conformation en boucle\ndans le jeu de données redondant (PDB 2007 - Assignation de structure secondaire réalisée par DSSP)",ylim=c(0,100),breaks=seq(0,15,0.01),xlim=c(0,4))
dev.off()

data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_5_lettres_Probabilite.txt",header=T))
png("PDB2007RPROBCameleonplot5lettres.png", width = 800, height = 500)
plot((data5[,1]),(data5[,2]), xlab="Probabilité de conformation en hélice", ylab="Probabilité de conformation en feuillet", main="Probabilité de conformation en hélice en fonction de la probabilité de conformation en feuillet\ndes séquences caméléon de longueur 5 (PDB redondant 2007)",ty="p",cex=0.7, pch=18)
dev.off()

#6lettres
data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_6_lettres_Probabilite.txt",header=T))
png("PDB2007RPROBCameleonHelice6lettres.png", width = 800, height = 500)
hist((data5[,1]),col="red",xlab="Probabilité de conformation en hélice",  ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 6) à adopter une conformation en hélice\ndans le jeu de données redondant (PDB 2007 - Assignation de structure secondaire réalisée par DSSP)",breaks=seq(0,4,0.05))
dev.off()

data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_6_lettres_Probabilite.txt",header=T))
png("PDB2007RPROBCameleonFeuillet6lettres.png", width = 800, height = 500)
hist((data5[,2]),col="blue",xlab="Probabilité de conformation en feuillet", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 6) à adopter une conformation en feuillet\ndans le jeu de données redondant (PDB 2007 - Assignation de structure secondaire réalisée par DSSP)")
dev.off()

data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_6_lettres_Probabilite.txt",header=T))
png("PDB2007RPROBCameleonCoil6lettres.png", width = 800, height = 500)
hist((data5[,3]),col="green",xlab="Probabilité de conformation en boucle", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 6) à adopter une conformation en boucle\ndans le jeu de données redondant (PDB 2007 - Assignation de structure secondaire réalisée par DSSP)")
dev.off()

data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_6_lettres_Probabilite.txt",header=T))
png("PDB2007RPROBCameleonplot6lettres.png", width = 800, height = 500)
plot((data5[,1]),(data5[,2]), xlab="Probabilité de conformation en hélice", ylab="Probabilité de conformation en feuillet", main="Probabilité de conformation en hélice en fonction de la probabilité de conformation en feuillet\ndes séquences caméléon de longueur 6 (PDB redondant 2007)")
dev.off()

data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_4_lettres_Probabilite.txt",header=T))
#2 helice 3 feuillet 4boucle
conformation = matrix(0,41,4)
conformation[,1]=seq(0,4,0.1)
#comptage des helices
for (i in 1:dim(data5)[1]){
	for (j in 1:40){
		#print (j)
		#helice	
		if(data5[i,1]<conformation[j+1,1]&&data5[i,1]>=conformation[j,1]){ 
			conformation[j,2]=conformation[j,2]+1			
			break		
		}

	}
	if(data5[i,1]>=conformation[41,1]){
		conformation[41,2]=conformation[41,2]+1	
	}
	#comptage des feuillets
	for (j in 1:40){
		if(data5[i,2]<conformation[j+1,1]&&data5[i,2]>=conformation[j,1]){ 
			conformation[j,3]=conformation[j,3]+1			
			break		
		}
	}
	if(data5[i,2]>=conformation[41,1]){
		conformation[41,3]=conformation[41,3]+1	
	}
	#comptage des coils
	for (j in 1:40){
		if(data5[i,3]<conformation[j+1,1]&&data5[i,3]>=conformation[j,1]){ 
			conformation[j,4]=conformation[j,4]+1			
			break		
		}
	}
	if(data5[i,3]>=conformation[41,1]){
		conformation[41,4]=conformation[41,4]+1	
	}
}
#####################################################inf et sup a 1###############
helice=c(0,0)
feuillet=c(0,0)
coil=c(0,0)
for (i in 1:dim(data5)[1]){
	if(data5[i,1]<1){
		helice[1]=helice[1]+1	
	}
	else{
		helice[2]=helice[2]+1
	}

	if(data5[i,2]<1){
		feuillet[1]=feuillet[1]+1	
	}
	else{
		feuillet[2]=feuillet[2]+1
	}
	if(data5[i,3]<1){
		coil[1]=coil[1]+1	
	}
	else{
		coil[2]=coil[2]+1
	}
}
total=0
for(i in 1:2){
	total=total+helice[i]
	total=total+feuillet[i]
	total=total+coil[i]
}

for(i in 1:2){
	print("helice")
	print(helice[i])
	print("feuillet")	
	print(feuillet[i])
	print("coil")	
	print(coil[i])
}

#####################################################inf et sup a 1###############
conformation = matrix(0,2,2)
coilA=c(0)
coilB=c(0)
coilD=c(0)
coilC=c(0)

for (i in 1:dim(data5)[1]){
	#cas helice <1 et feuillet <1
	if(data5[i,1]<1&&data5[i,2]<1){
		conformation[1,1]=conformation[1,1]+1
		coilC=rbind(coilC,data5[i,3])	
	}
	#cas helice >=1 et feuillet >=1
	else if(data5[i,1]>=1&&data5[i,2]>=1){
		conformation[2,2]=conformation[2,2]+1
		coilB=rbind(coilB,data5[i,3])
	}
	#cas helice <1 et feuillet >=1
	else if(data5[i,1]<1&&data5[i,2]>=1){
		conformation[2,1]=conformation[2,1]+1
		coilA=rbind(coilA,data5[i,3])	
	}
	#cas helice >=1 et feuillet <1
	else if(data5[i,1]>=1&&data5[i,2]<1){
		conformation[1,2]=conformation[1,2]+1
		coilD=rbind(coilD,data5[i,3])
	}
}
coilA=coilA[-1]
coilD=coilD[-1]
coilB=coilB[-1]
coilC=coilC[-1]
png("PDB2007RPROBCameleonCoil4lettres_histoA.png", width = 800, height = 500)
hist(coilA, col="red",ylim=c(0,600),breaks=seq(0,8,0.01),xlim=c(0,3),xlab="Probabilité de conformation en boucle", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 4) à adopter une conformation en boucle\nCAS A (PDB 2007 jeu de données redondant - Assignation de structure secondaire réalisée par DSSP)")
abline(h=0,v=mean(coilA))
dev.off()
png("PDB2007RPROBCameleonCoil4lettres_histoB.png", width = 800, height = 500)
hist(coilB, col="green",ylim=c(0,600),breaks=seq(0,8,0.01),xlim=c(0,3),xlab="Probabilité de conformation en boucle", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 4) à adopter une conformation en boucle\nCAS B (PDB 2007 jeu de données redondant - Assignation de structure secondaire réalisée par DSSP)")
abline(h=0,v=mean(coilB))
dev.off()
png("PDB2007RPROBCameleonCoil4lettres_histoC.png", width = 800, height = 500)
hist(coilC, col="yellow",ylim=c(0,600),breaks=seq(0,8,0.01),xlim=c(0,3),xlab="Probabilité de conformation en boucle", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 4) à adopter une conformation en boucle\nCAS C (PDB 2007 jeu de données redondant - Assignation de structure secondaire réalisée par DSSP)")
abline(h=0,v=mean(coilC))
dev.off()
png("PDB2007RPROBCameleonCoil4lettres_histoD.png", width = 800, height = 500)
hist(coilD, col="blue",ylim=c(0,600),breaks=seq(0,8,0.01),xlim=c(0,3),xlab="Probabilité de conformation en boucle", ylab="Occurence", main="Densité de probabilité des séquences caméléon (longueur 4) à adopter une conformation en boucle\nCAS D (PDB 2007 jeu de données redondant - Assignation de structure secondaire réalisée par DSSP)")
abline(h=0,v=mean(coilD))
dev.off()
####################
data5=as.matrix(read.table("all-jan-07clear_SeqCameleon_4_lettres_Probabilite.txt",header=T))
#2 helice 3 feuillet 4boucle
conformation = matrix(0,41,41)
step=seq(0,4,0.1)
#comptage des helices
for (i in 1:dim(data5)[1]){
	for (j in 1:40){
		#les helices en abscisse
		x=as.integer(data5[i,1]*10)+1
		#les feuillets en ordonnee		
		y=as.integer(data5[i,2]*10)+1
		#print(data5[i,1])
		#print(data5[i,2])		
		#print("x")
		#print(x)
		#print("y")		
		#print(y)
		if(x>40) x=41
		if(y>40) y=41 		
		conformation[y,x]=conformation[y,x]+1	
	}
}

png("PDB2007RScore4lettres.png", width = 800, height = 500)
par(pty="s")
image(x=seq(0,4,0.1),y=seq(0,4,0.1),conformation,col=terrain.colors(100),main="Score des séquences caméléons de longueur 4 (PDB 2007 redondant)",xlab="Probabilité de conformation en hélice",ylab="Probabilité de conformation en feuillet")
dev.off()


image(conformation,axes=F)
lab=c("0-1",">1")
axis(1,c(0,1),lab)
axis(2,c(0,1),lab)

box()






y=c(1,2)
x=conformation[,1]
png("PDB2007RScore4lettres.png", width = 800, height = 500)
image(x,y,conformation[,2:3], col = heat.colors(80),main="Score séquence Caméléon longueur 4 (PDB 2007 redondant)",ylim=c(1,2),xlab="score(probabilité de conformation)(2007 redondant)",ylab="1.0-1.5 helice\t1.5-2.0 feuillet")
y=c(1,2)
contour(x,y,conformation[,2:3],levels = seq(100, 1800, by = 300),add = TRUE, col = "blue")
dev.off()

temp=conformation[,2]/conformation[,3]
x=y=conformation[,1]
image(temp, col = heat.colors(80),main="Score séquence Caméléon longueur 4 (PDB 2007 redondant)",ylim=c(1,2),xlab="score(probabilité de conformation)(2007 redondant)",ylab="1.0-1.5 helice\t1.5-2.0 feuillet")


data5=as.matrix(read.table("all-jan-07clear_bfactor_4_lettres.txt",header=T))
