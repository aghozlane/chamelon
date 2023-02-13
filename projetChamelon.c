#include "projetChamelon.h"
#include "options.h"
#include "frequences.h"
#include "arborescence.h"


/**************LECTURE DU FICHIER************/
//construction des mots
parbre wordconstruction(char a, int structure, int lettrecorrect,double bfactoresult,pstatistique chamelon, parbre l){
  int i,j;
  char *seq2;
  //printf("%d\n", structure);
  //char test[4]={'M','L','I','L'};
  //la nouvelle lettre a une structure secondaire
  if(lettrecorrect==OUI){
    //cas ou l'on a pas encore atteint la longueur du mot
    if( (*chamelon->cas)<(chamelon->wordlen-1)){
      chamelon->seq[(*chamelon->cas)]=a;
      chamelon->structseq[(*chamelon->cas)]=structure;
      chamelon->accesslet[(*chamelon->cas)]=chamelon->accessibilite;
      chamelon->bfactorlet[(*chamelon->cas)]=bfactoresult;
      (*chamelon->cas)++;
    }
    //cas ou la nouvelle lettre nous permet d'atteindre la longueur d'un mot
    else if((*chamelon->cas)==(chamelon->wordlen-1)){
      chamelon->seq[(*chamelon->cas)]=a;
      chamelon->structseq[(*chamelon->cas)]=structure;
      chamelon->accesslet[(*chamelon->cas)]=chamelon->accessibilite;
      chamelon->bfactorlet[(*chamelon->cas)]=bfactoresult;
      
      /*for(i=0;i<(chamelon->wordlen);i++){
	printf("accessibilite %f bfactor %lf\n",chamelon->accesslet[i],chamelon->bfactorlet[i]);
	}*/
      
      seq2=(char*)malloc(chamelon->wordlen*sizeof(char));
      strcpy(seq2,chamelon->seq);
      (chamelon->nbword)++;
      
      //cas n'est pas incremente et reste a 3
      //renvoie dans la structure cas 4
      l=buildtree(seq2,chamelon,l);
      //printf("seq buildtree : %s \n", l->seq);
      //decallage du tableau
      for(i=0;i<(chamelon->wordlen);i++){
	chamelon->seq[i]=chamelon->seq[i+1];
	chamelon->structseq[i]=chamelon->structseq[i+1];
	chamelon->accesslet[i]=chamelon->accesslet[i+1];
	chamelon->bfactorlet[i]=chamelon->bfactorlet[i+1];
      }
    }
  }
  //cas ou la lettre n'a pas de structure retour du decompte a 0 
  //le mot est supprime
  else (*chamelon->cas)=0;
  return l;
}

//Traitement de l'acide amine identifie
parbre lancement(char a,char b,FILE *fic,int test,int structure,int lettrecorrect,pstatistique chamelon,parbre l,int option[7]){ 
  int i;
  double bfactoresult=0.0;
  //# cas non lus > nouvelle proteine
  if((a!='#')&&(a!='>')&&(a!=EOF)){
    //DSSP
    if(option[3]==5)test+=fscanf(fic,"\t%c\t%d",&b,&structure);
    //STRIDE
    if(option[3]==6)test+=fscanf(fic,"\t%c\t%d",&b,&structure);
    //KAKSI
    if(option[3]==7)test+=fscanf(fic,"\t%c\t%d",&b,&structure);
    //SEGNO
    if(option[3]==8)test+=fscanf(fic,"\t%c\t%d",&b,&structure);
    
    
    if(test!=2){ 
      fprintf(stderr,"Colonne deficiente, format :  \\tb\\tstructure\\tangles\\tposition\\accessibilite\\bfactor\n");
      exit(1);			
    }
    
    //longueur de la proteine
    chamelon->longueur++;
    
    //calcul des frequences
    lettrecorrect=freqAA(a,structure,chamelon);
    if(lettrecorrect==OUI){
      //lectures des autres donnes de position et d'angle
      for(i=0;i<7;i++){
	test+=fscanf(fic,"\t%f",&chamelon->Otherdata[i]);
      }
      //accessibilite
      test+=fscanf(fic,"\t%f",&chamelon->accessibilite);
      //bfactor
      for(i=0;i<3;i++){
	test+=fscanf(fic,"\t%lf",&chamelon->bfactor[i]);
      }
      
      if(test!=13){ 
	fprintf(stderr,"Colonne deficiente, format :  \\tb\\tstructure\\tangles\\tposition\\accessibilite\\bfactor\n");
	exit(1);			
      }
      //Accessibilite
      //calcAccess(a,chamelon);
      //calcul du Bfactor
      bfactoresult=(chamelon->bfactor[0]-chamelon->bfactor[2])/chamelon->bfactor[1];
      //enregistrement de l'accessibilite
      //calcBfactor(a,bfactoresult,chamelon);
    }
    test=0;
    
    //construction de l'arbre
    if(option[1]==1) l=wordconstruction(a,structure,lettrecorrect,bfactoresult,chamelon,l);
  }
  return l;
}
//Passage de l'entete
void entete(char a,FILE *fic,pstatistique chamelon){
  int i;
  a=fgetc(fic);
  if(a=='>') {
    for(i=0;i<4;i++) chamelon->pdbseq[i]=fgetc(fic);
    chamelon->pdbseq[5]='\0';
  }
  else {
    chamelon->pdbseq[0]=a;
    for(i=1;i<4;i++) chamelon->pdbseq[i]=fgetc(fic);
    chamelon->pdbseq[5]='\0';
  }
  //printf("seq %s",chamelon->pdbseq);
  do    a=fgetc(fic);
  while(a!='\n'&& a!=EOF);
}

//lecture du fragment
parbre lecture(char a,FILE *fic,FILE *fic4,pstatistique chamelon,parbre l,int option[7]){
  char b;
  int test=NON, structure=NON,lettrecorrect=NON;

  //parbre m=NULL;
  entete(a,fic,chamelon);
  //identification de la premiere lettre
  l=lancement(fgetc(fic),b,fic,test,structure,lettrecorrect,chamelon,l,option);
  chamelon->Nbfragments++;
  chamelon->Nblines++;
  //identification des autres lettres
  do{
    a=fgetc(fic);
    if(a=='\n'){
      a=fgetc(fic);
      l=lancement(a,b,fic,test,structure,lettrecorrect,chamelon,l,option);
      chamelon->Nblines++;
    }
  }
  while(a!='>'&& a!=EOF);
  //donnee de distribution au sein des fragments
  fprintf(fic4,"%ld\t%ld\n", chamelon->Nbfragments,chamelon->longueur);
  chamelon->longueur=0;
  if(a=='>') lecture(a,fic,fic4,chamelon,l,option);
  return l;
}

/*****************************FICHIER D'INFORMATION************************************/
void outwork(char gnufilename[100],pstatistique chamelon,poptions prog,FILE *fic6,FILE *fic7){
  //Fichier de sauvegarde sequences cameleons
  sprintf(gnufilename,"%s_SeqCameleon_%d_lettres.txt", prog->sortie,chamelon->wordlen);
  printf("%s\n",gnufilename);
  fic6=fopen(gnufilename,"wt");
  fprintf(fic6,"\"Identifications en helice alpha\"\"Identifications en feuillet beta\"\"PDB\"\n");
  printf("yop\n");
  //fichier de sauvegarde de l'arborescence
  sprintf(gnufilename,"%s_arborescence_%d_lettres.txt", prog->sortie,chamelon->wordlen);
  fic7=fopen(gnufilename,"wt");
  fprintf(fic7,"nbword = %ld\n",chamelon->nbword);
  fprintf(fic7,"nbworddiff = %ld\n",chamelon->nbworddiff);
  fprintf(fic7,"nbfragments = %ld\n",chamelon->Nbfragments);
  fprintf(fic7,"nblines = %ld\n",chamelon->Nblines);
  fprintf(fic7,"wordlen = %d\n",chamelon->wordlen);
  
  //fichier de sauvegarde de l'accessibilite des sequences cameleons
  sprintf(gnufilename,"%s_accessibilite_%d_lettres.txt",prog->sortie,chamelon->wordlen);
  chamelon->access=fopen(gnufilename,"wt");
  fprintf(chamelon->access,"\"Accessibilite cameleon\"\n");
  
  //fichier de sauvegarde du bfactor des sequences cameleons
  sprintf(gnufilename,"%s_bfactor_%d_lettres.txt",prog->sortie,chamelon->wordlen);
  chamelon->bfact=fopen(gnufilename,"wt");
  fprintf(chamelon->bfact,"\"bfactor\"\n");  
}
void informations(char *sortie,char gnufilename[100],pstatistique chamelon){
  FILE *fic5;
  sprintf(gnufilename, "%s_Infolecture_%d_lettres.txt", sortie,chamelon->wordlen);
  fic5=fopen(gnufilename,"wt");
  fprintf(fic5,"Nom du fichier\tNombre de lignes\tNombre de fragments\n");
  fprintf(fic5,"%s\t%ld\t\t\t%ld\n\n",sortie,chamelon->Nblines, chamelon->Nbfragments);
  fprintf(fic5,"Longueur de sequence choisie %d\n",chamelon->wordlen);
  fprintf(fic5,"Nombre d'acide amine compte = %1.0lf\n", chamelon->Nbaatotaux);
  fprintf(fic5,"Nombre d'erreur (Acides amines mal identifies)= %1.0lf\n", chamelon->Nbaa[20]);
  fprintf(fic5,"Nombre de mots identifies = %d\nNombre de mots differents identifies = %d\n",(chamelon->nbword),(chamelon->nbworddiff));
  fprintf(fic5,"Nombre de mots observes plus d'une fois = %d\n",(chamelon->Nbfragdiff2));
  fprintf(fic5,"Nombre de maillon supprime %d\n", (chamelon->maillon_suppr));
  fprintf(fic5,"Nombre de sequences cameleons %d\n",(chamelon->Nbcameleon));
  fprintf(stdout,"Fin de l'analyse\n");
  fclose(fic5);
}

/**************MAIN*********************/
int main(int argc, char**argv){
  FILE *fic,*fic4,*fic6,*fic7;
  char a,gnufilename[100];
  poptions prog=initoption(prog);
  pstatistique chamelon=NULL;
  parbre l=NULL;
  int *nbarguments=(int*)malloc(sizeof(int));
  (*nbarguments)=0;
  
  //Sans choix d'option : affiche les options disponibles
  if(argc==1){
    optiondisp(prog->seq);
    exit(0);
  }
  //un fichier de lecture imposee
  if(argc<3){
    fprintf(stderr,"Usage minimum: ./programme -FicEntre fichier -Longueur chiffre\n");
    optiondisp(prog->seq);
    freeoptions(prog);
    exit(1);
  }
  //determination des options
  defoption(argv,argc,prog);
  //lecture des options et initialisation de chamelon
  
  if(prog->option[0]==2){
    chamelon=initstruct(chamelon,0);
    fprintf(stdout,"Lecture du fichier des frequences...\n");
    chamelon=letFREQ(chamelon,prog->ficfreq);
  }
  else if(prog->option[1]==2){
    if(chamelon==NULL) chamelon=initstruct(chamelon,0);
    fprintf(stdout,"Lecture du fichier d'arborescence...\n");
    l=letarborescence(chamelon,prog->ficarbre);
  }
  else chamelon=initstruct(chamelon,prog->option[2]);
  //s'il y a un fichier d'entree
  if(prog->option[4]==OUI){
    //ouverture du fichier et lecture
    if((fic=fopen(prog->entre,"r"))!=NULL){
      fprintf(stdout,"Analyse du fichier en cours...\n");
      //analyse de la sequence
      
      //distribution des acides amines
      sprintf(gnufilename, "%s_Distribution.txt",prog->sortie);
      fic4=fopen(gnufilename,"wt");
      fprintf(fic4,"\"longueur du fragment\"\n");
      
      //accessibilite
      /*sprintf(gnufilename, "%s_accessibilite.txt",prog->sortie);
      chamelon->access=fopen(gnufilename,"wt");
      fprintf(chamelon->access,"\"Accessibilite\"\n");*/
      
      //B-factor
      /*sprintf(gnufilename, "%s_bfactor.txt",prog->sortie);
      chamelon->bfact=fopen(gnufilename,"wt");
      fprintf(chamelon->bfact,"\"B-facteur normalise\"\n");*/
      
      //lecture du fichier 
      l=lecture(a,fic,fic4,chamelon,l,prog->option);
      
      fclose(fic4);
      //fclose(chamelon->access);
      //fclose(chamelon->bfact);
      fclose(fic);
    }
    else fprintf(stderr,"Il est impossible d'ouvrir ce fichier. Vérifier que ce fichier existe et qu'il n'atteint pas la limite de 2GO du systeme de fichier\n");
  }
  
  //statistiques des fragments
  if(prog->option[0]>NON){
    fprintf(stdout,"Analyse des statistiques du fichier...\n");
    sortieFREQ(chamelon,prog->sortie,gnufilename);
    affichFREQ(chamelon,prog->sortie,gnufilename);
  }
  //analyse de l'arbre
  if(l!=NULL){
    //outwork(gnufilename,chamelon,prog,fic6,fic7);
    
    //Fichier de sauvegarde sequences cameleons
    sprintf(gnufilename,"%s_SeqCameleon_%d_lettres.txt", prog->sortie,chamelon->wordlen);
    fic6=fopen(gnufilename,"wt");
    fprintf(fic6,"\"Identifications en helice alpha\"\"Identifications en feuillet beta\"\"PDB\"\n");
    //fichier de sauvegarde de l'arborescence
    sprintf(gnufilename,"%s_arborescence_%d_lettres.txt", prog->sortie,chamelon->wordlen);
    fic7=fopen(gnufilename,"wt");
    fprintf(fic7,"nbword = %ld\n",chamelon->nbword);
    fprintf(fic7,"nbworddiff = %ld\n",chamelon->nbworddiff);
    fprintf(fic7,"nbfragments = %ld\n",chamelon->Nbfragments);
    fprintf(fic7,"nblines = %ld\n",chamelon->Nblines);
    fprintf(fic7,"wordlen = %d\n",chamelon->wordlen);
    
    //fichier de sauvegarde de l'accessibilite des sequences cameleons
    sprintf(gnufilename,"%s_accessibilite_%d_lettres.txt",prog->sortie,chamelon->wordlen);
    chamelon->access=fopen(gnufilename,"wt");
    fprintf(chamelon->access,"\"Accessibilite cameleon\"\n");
    
    //fichier de sauvegarde du bfactor des sequences cameleons
    sprintf(gnufilename,"%s_bfactor_%d_lettres.txt",prog->sortie,chamelon->wordlen);
    chamelon->bfact=fopen(gnufilename,"wt");
    fprintf(chamelon->bfact,"\"bfactor\"\n");  
    
        
    fprintf(stdout,"Analyse de l'arborescence du fichier...\n");
    analyse_arbre(fic6,fic7,l,chamelon,nbarguments);
    //fprintf(stderr,"out\n");
    fprintf(fic7,"nbarguments = %d\n",((*nbarguments)+1));
    fclose(fic6);
    fclose(fic7);
    free(nbarguments);
    suppr_arbre(l,chamelon);
  }    
  informations(prog->sortie,gnufilename,chamelon);
  freechamelon(chamelon);  
  freeoptions(prog);
  return 0;
}
