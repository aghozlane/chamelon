
fichier=open("/home/curry/ghozlane/tests/data/all-2008a","r")
sortie=open("/home/curry/ghozlane/tests/data/all-2008apart1","w")
sortie2=open("/home/curry/ghozlane/tests/data/all-2008apart2","w")
line="rien"
l=0
m=0
t=0
for iline in fichier:
	if(iline[0]!='>'):
		line=iline.split()
		i=len(line)
		#14 pour l'autre
		if(i==17 and t!=1):
			sortie.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\t"+line[9]+"\t"+line[10]+"\t"+line[11]+"\t"+line[12]+"\t"+line[13]+"\n")
			#sortie.write("%s\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n" % line[0] % line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13])
		elif(i==17 and t==1):
			sortie2.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\t"+line[9]+"\t"+line[10]+"\t"+line[11]+"\t"+line[12]+"\t"+line[13]+"\n")		
		else: 
			l=l+1
			print l
			print prot
			sortie.write("#"+iline)
	else:
		prot=iline
		sortie.write(iline)
		m=m+1
	        if(m==126886): 
			t=1
fichier.close()
sortie.close()
sortie2.close()

