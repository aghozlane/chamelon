#include "frequences.h"

/******OBJET CONTENANT LES DONNEES STATISTIQUES******/
//structure contenant les differents domaines d'etude statistique
pstatistique initstruct(pstatistique chamelon,int longueur){
  int i,j;
  char lettres[20]={'L','A','G','V','E','S','D','T','K','I','R','P','N','F','Q','Y','H','M','C','W'};
  
  chamelon=(pstatistique)malloc(sizeof(statistique));
  
  if(longueur!=0){
    chamelon->wordlen=longueur;
    chamelon->seq=(char*)malloc(chamelon->wordlen*sizeof(char));
    chamelon->structseq=(int*)malloc(chamelon->wordlen*sizeof(int));
    chamelon->bfactorlet=(double*)malloc(chamelon->wordlen*sizeof(double));
    chamelon->accesslet=(float*)malloc(chamelon->wordlen*sizeof(float));
    for(i=0;i<chamelon->wordlen;i++){
      chamelon->structseq[i]=0;
      chamelon->bfactorlet[i]=0.0;
      chamelon->accesslet[i]=0.0;
    }
  }  
  chamelon->cas=(int*)malloc(sizeof(int));
  (*chamelon->cas)=0;
  chamelon->nbword=0;
  chamelon->nbworddiff=0;
  chamelon->Nbfragments=0;
  chamelon->longueur=0;
  chamelon->Nblines=0;
  chamelon->Nbaatotaux=0.0;
  chamelon->maillon_suppr=0;
  chamelon->Nbfragdiff2=0;
  chamelon->Nbcameleon=0;
  chamelon->structcameleon[0]=0;chamelon->structcameleon[1]=0;
  chamelon->alphabet=(char*)malloc(20*sizeof(char));
  strcpy(chamelon->alphabet,lettres);
  //initialisation du tableau des acides amines
  for(i=0;i<4;i++) chamelon->FreqStructure[i]=0.0;
  for(i=0;i<21;i++)(chamelon->Nbaa[i])=0.0;
  for(j=0;j<3;j++){
    for(i=0;i<20;i++){
      chamelon->AAstruct[i][j]=0.0;
      chamelon->propension[i][j]=0.0;
    }
  }
  for(i=0;i<7;i++) chamelon->Otherdata[i]=0.0;
  for(i=0;i<3;i++) chamelon->bfactor[i]=0.0;
  chamelon->accessibilite=0.0;
 
  return chamelon;
}
//liberation de la memoire occupee par chatmelon
void freechamelon(pstatistique chamelon){
  free(chamelon->structseq);
  free(chamelon->cas);
  free(chamelon->alphabet);
  free(chamelon->bfactorlet);
  free(chamelon->accesslet);
  free(chamelon);
}


/*************CALCUL DES FREQUENCES**********/
//indentification de la structure de l'aa
int freqStruct(int structure, pstatistique chamelon){
  int sortie=-1;
  //aa en helices
  if(structure<3){ 
    (chamelon->FreqStructure[0])+=1.0;
    sortie=0;
  }
  //aa en boucle
  else if(structure>=3 && structure<=6){ 
    (chamelon->FreqStructure[1])+=1.0;
    sortie=1;
  }
  //aa en feuillet
  else if(structure==7){ 
    (chamelon->FreqStructure[2])+=1.0;
    sortie=2;
  }
  //structure non determinee de aa
  else{ 
    (chamelon->FreqStructure[3])+=1.0;
    sortie=-2;
  }
  return sortie;
}

//calcul de la frequence de l'acide amine
int freqAA(char a, int structure,pstatistique chamelon){
  int conf=freqStruct(structure,chamelon), i,j=0;
  //identification de la lettre
  for(i=0;i<20;i++){
    if(a==(chamelon->alphabet[i])){ 
      (chamelon->Nbaa[i])+=1.0;
      //printf("%d %f\n",i,(chamelon->Nbaa[i]));
      if(conf>=0) chamelon->AAstruct[i][conf]++;
      break;
    }
    else j++;
  }
  if(j==20){
    //acides amines mal identifies
    fprintf(stderr,"%c\n", a);
    chamelon->Nbaa[20]++;
    //fprintf(stderr,"%f",chamelon->Nbaa[20]);
    return NON;
  }
  //la lettre n'a pas structure
  if(conf==-2) return NON;
  return OUI;
}

//accessibilite relative des residus
void calcAccess(char a,pstatistique chamelon){
  fprintf(chamelon->access,"%c\t%f\n",a,chamelon->accessibilite);
}

//Calcul du Bfacteur normalise
void calcBfactor(char a,double bfactoresult,pstatistique chamelon){
  fprintf(chamelon->bfact,"%c\t%lf\n",a,bfactoresult);
}

/***************SORTIE DES RESULTATS DE FREQUENCE************/
//affichage des frequences
void affichFREQ(pstatistique chamelon, char *sortie, char gnufilename[100]){
  FILE *fic1,*fic2,*fic3;
  int i,j;
  
  //determination du nombre d'acide amine total
  for(i=0; i<20;i++){
    chamelon->Nbaatotaux+=chamelon->Nbaa[i];
    //printf("%c %f\n",chamelon->alphabet[i],chamelon->Nbaa[i]);
  }
  
  //determination de la frequence de chacun des acides amines
  sprintf(gnufilename, "%s_Frequences.txt", sortie);
  fic1=fopen(gnufilename,"wt");
  fprintf(fic1,"\"Frequence de l'acide amine\"\"Frequence en helice\"\"Frequence en boucle\"\"Frequence en feuillet\"\n");
  for(i=0; i<20;i++){
    chamelon->Nbaa[i]=chamelon->Nbaa[i]/chamelon->Nbaatotaux*100.0;
    //determination de la frequence des acides amines pour chaque structure
    //for(j=0;j<3;j++) chamelon->AAstruct[i][j]=chamelon->AAstruct[i][j]/chamelon->Nbaatotaux*100.0;
    fprintf(fic1,"%c\t%1.2lf\t%1.2lf\t%1.2lf\t%1.2lf\n",chamelon->alphabet[i], chamelon->Nbaa[i],chamelon->AAstruct[i][0]/chamelon->Nbaatotaux*100.0,chamelon->AAstruct[i][1]/chamelon->Nbaatotaux*100.0,chamelon->AAstruct[i][2]/chamelon->Nbaatotaux*100.0);
    for(j=0;j<3;j++) chamelon->AAstruct[i][j]=chamelon->AAstruct[i][j]/chamelon->FreqStructure[j]*100.0;
  }
  fclose(fic1);
  
  //affichage de la propension
  sprintf(gnufilename, "%s_Propension.txt", sortie);
  fic2=fopen(gnufilename,"wt");
  fprintf(fic2,"\"Propension en helice\"\"Propension en boucle\"\"Propension en feuillet\"\n");
  for(i=0; i<20;i++){
    for(j=0;j<3;j++) chamelon->propension[i][j]=(chamelon->AAstruct[i][j]/chamelon->Nbaa[i]);
    fprintf(fic2,"%c\t%1.2lf\t%1.2lf\t%1.2lf\n", chamelon->alphabet[i],chamelon->propension[i][0],chamelon->propension[i][1],chamelon->propension[i][2]);
  }
  fclose(fic2);
  
  //Frequences globales
  sprintf(gnufilename, "%s_FrequencesGlobales.txt", sortie);
  fic3=fopen(gnufilename,"wt");
  //Affichage de la frequence de chacunes des structures
  fprintf(fic3,"\"Frequence des residus en helices\"\"Frequence des residus en boucles\"\"Frequence des residus en feuillets\"\"Frequence des problemes d'assignation de structure secondaire\"\n");
  fprintf(fic3,"%1.2lf\t%1.2lf\t%1.2lf\t%1.2lf\t\n",(chamelon->FreqStructure[0]/chamelon->Nbaatotaux*100.0),(chamelon->FreqStructure[1]/chamelon->Nbaatotaux*100.0),(chamelon->FreqStructure[2]/chamelon->Nbaatotaux*100.0),(chamelon->FreqStructure[3]/chamelon->Nbaatotaux*100.0));
}

void sortieFREQ(pstatistique chamelon, char *sortie, char gnufilename[100]){
  FILE *fic9;
  int i;
  sprintf(gnufilename, "%s_Chamelon.txt", sortie);
  fic9=fopen(gnufilename,"wt");
  fprintf(fic9,"nbaa =\n");
  for(i=0;i<20;i++){
    fprintf(fic9,"%c\t%.1lf\t%.1lf\t%.1lf\t%.1lf\n",chamelon->alphabet[i],chamelon->Nbaa[i],chamelon->AAstruct[i][0],chamelon->AAstruct[i][1],chamelon->AAstruct[i][2]);
  }
  fprintf(fic9,"Z\t%.1lf\n",chamelon->Nbaa[21]);
  fprintf(fic9,"freqstructure =\n");
  for(i=0;i<4;i++){
    fprintf(fic9,"%.1lf\n",chamelon->FreqStructure[i]);
  }
  fclose(fic9);
}

/**************LECTURE DU FICHIER DE SAUVEGARDE****************/
pstatistique letFREQ(pstatistique chamelon,char *lectfreq){
  FILE *fic8;
  char lettre;
  int nbarguments=0,i;
  if((fic8=fopen(lectfreq,"rt"))!=NULL){
      //la ca doit etre 0 en terme d'arguments lus
    nbarguments+=fscanf(fic8,"nbaa =\n");
    for(i=0;i<20;i++){
      nbarguments+=fscanf(fic8,"%c\t%lf\t%lf\t%lf\t%lf\n",&(lettre),&(chamelon->Nbaa[i]),&(chamelon->AAstruct[i][0]),&(chamelon->AAstruct[i][1]),&(chamelon->AAstruct[i][2]));
      if(lettre!=(chamelon->alphabet[i])){
	fprintf(stderr,"Erreur de lecture des donnees des acides amines du fichier de frequence\n");
	freechamelon(chamelon);
	exit(1);
      }
    }
    nbarguments+=fscanf(fic8,"Z\t%lf\n",&(chamelon->Nbaa[21]));
    //la ca doit etre 0 en terme d'arguments lus
    nbarguments+=fscanf(fic8,"freqstructure =\n");
    for(i=0;i<4;i++){
      nbarguments+=fscanf(fic8,"%lf\n",&(chamelon->FreqStructure[i]));
    }
    fclose(fic8);
    if(nbarguments!=105){
      fprintf(stderr,"Nombre d'arguments lu :%d est diffent du nombre d'arguments total a lire : 110\n",nbarguments);
      freechamelon(chamelon);
      exit(1);
    }
    //allocation des deux tableaux dependant de la longueur des sequences
  }
  else{
    fprintf(stderr,"Erreur d'ouverture du fichier des frequences : %s\n",lectfreq);
    freechamelon(chamelon);
    exit(1);
  }
  return chamelon;
}
