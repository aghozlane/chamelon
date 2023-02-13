#!/usr/bin/python
import sys
if len(sys.argv) != 3:
   sys.exit("Usage : ./nettoyage.py fichier nombre de colonnes")

fichier=open(sys.argv[1],"r")
sortie=open(sys.argv[1]+"clear","w")
line="rien"
l=0
m=0
nbcol=int(sys.argv[2])
for iline in fichier:
	if(iline[0]!='>'):
		line=iline.split()
		i=len(line)
		#17 pour 2008
		#14 pour 2007
		if(i==nbcol):
			sortie.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\t"+line[9]+"\t"+line[10]+"\t"+line[11]+"\t"+line[12]+"\t"+line[13]+"\n")
			#sortie.write("%s\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n" % line[0] % line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13])
		else: 
			l=l+1
			print l
			print prot
			sortie.write("#"+iline)
	else:
		prot=iline
		sortie.write(iline)
		m=m+1
	#if(m==49129): 
	#	print iline
fichier.close()
sortie.close()
print m
