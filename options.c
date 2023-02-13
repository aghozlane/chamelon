#include "options.h"

/********************OPTIONS************************/
poptions initoption(poptions prog){
  int i;
  prog=(poptions)malloc(sizeof(options));
  prog->seq=allocseq(prog->seq);
  prog->sortie=NULL;
  prog->entre=NULL;
  prog->ficfreq=NULL;
  prog->ficarbre=NULL;
  for(i=0;i<7;i++){
    prog->option[i]=NON;
  }
  return prog;
}
//allocation des options disponibles
char **allocseq(char **seq){
  seq=(char**)malloc(11*sizeof(char*));
  seq[0]=(char*)malloc(5*sizeof(char));
  seq[1]=(char*)malloc(8*sizeof(char));
  seq[2]=(char*)malloc(5*sizeof(char));
  seq[3]=(char*)malloc(8*sizeof(char));
  seq[4]=(char*)malloc(10*sizeof(char));
  seq[5]=(char*)malloc(5*sizeof(char));
  seq[6]=(char*)malloc(7*sizeof(char));
  seq[7]=(char*)malloc(6*sizeof(char));
  seq[8]=(char*)malloc(6*sizeof(char));
  seq[9]=(char*)malloc(9*sizeof(char));
  seq[10]=(char*)malloc(9*sizeof(char));
  
  strcpy(seq[0],"-Freq");
  strcpy(seq[1],"-LetFreq");
  strcpy(seq[2],"-Arbr");
  strcpy(seq[3],"-LetArbr");
  strcpy(seq[4],"-FicSortie");
  strcpy(seq[5],"-DSSP");
  strcpy(seq[6],"-STRIDE");
  strcpy(seq[7],"-KAKSI");
  strcpy(seq[8],"-SEGNO");
  strcpy(seq[9],"-Longueur");
  strcpy(seq[10],"-FicEntre");
  return seq;
}
//free des differentes donnes d'option
void freeoptions(poptions prog){
  int i;
  for(i=0;i<11;i++){
    free(prog->seq[i]);
  }
  free(prog->seq);
  if(prog->sortie!=NULL) free(prog->sortie);
  if(prog->entre!=NULL) free(prog->entre);
  if(prog->ficfreq!=NULL) free(prog->ficfreq);
  if(prog->ficarbre!=NULL) free(prog->ficarbre);
  free(prog);
}
//affichage des options disponibles
void optiondisp(char **seq){
  fprintf(stdout,"Options disponibles :\n");
  fprintf(stdout,"Options\t\tDépendances\t\t\tEffets\n");
  fprintf(stdout,"1 -rien\t\t%s %s -OPTION\tCalcul des frequences et construction de l'arborescence, sortie des fichiers dans le dossier d'origine des donnees\n",seq[10],seq[9]);
  fprintf(stdout,"2 %s\t\t%s\t\t\tCalcul des frequences\n",seq[0],seq[10]);
  fprintf(stdout,"3 %s\t(%s) et/ou (%s)\tLecture des frequences et poursuite du calcul des frequences\n",seq[1],seq[10],seq[4]);
  fprintf(stdout,"4 %s\t\t%s %s -OPTION\tConstruction de l'arborescence\n",seq[2],seq[10],seq[9]);
  fprintf(stdout,"5 %s\t(%s) et/ou (%s)\tLecture d'un arbre et poursuite de la construction puis analyse de l'arbre\n",seq[3],seq[10],seq[4]);
  fprintf(stdout,"6 %s\t%s %s -OPTION\tChoix du dossier de sortie et du nom du fichier : ex data/all_jan\n",seq[4],seq[10],seq[9]);
  fprintf(stdout,"7 %s\t\t%s %s\t\tPrise en considération des données de DSSP\n",seq[5],seq[10],seq[9],seq[2]);
  fprintf(stdout,"8 %s\t%s %s\t\tPrise en considération des données de STRIDE\n",seq[6],seq[10],seq[9],seq[2]);
  fprintf(stdout,"9 %s\t%s %s\t\tPrise en considération des données de KAKSI\n",seq[7],seq[10],seq[9],seq[2]);
  fprintf(stdout,"10 %s\t%s %s\t\tPrise en considération des données de SEGNO\n",seq[8],seq[10],seq[9],seq[2]);
}
//identifie l'option
int identoption(char **seq,char **argv,int num){
  int i;
  for(i=0;i<11;i++){
    //printf("%s\t%s\n",seq[i],argv[num]);
    if(strcmp(seq[i],argv[num])==0){
      //printf("i %d\n",i);
      return i;
    }
  }
  return (-1);
}
//verification de la frequence proposee
void deflongueur(int longueur){
  //traitement de la longueur de sequence demandee
  if(longueur>0){
    printf("-Longueur de sequence choisie : %d\n",longueur);
    if(longueur>=20) fprintf(stderr,"La longueur de sequence choisie est elevee, risque de saturation de la memoire\n");
  }
  //echec de conversion
  else{
    fprintf(stderr,"Choix de la longueur de la sequence inexact, longueur usuellement comprise entre 4-8\n");
    exit(1);      
  }
}
//verification sommaire du fichier
int compair(char **seq,char **argv,int num){
  int i;
  for(i=0;i<11;i++){
    if(strcmp(seq[i],argv[num])==0) return OUI;
  }
  return NON;
}

//action
void notifoption(int choix,poptions prog,char **argv,int num,int argc){
  //calcul des frequences
  if((choix==0)&&(prog->option[0]==NON)){
    //frequence
    prog->option[0]=1;
    printf("Calcul des frequences\n");
    //pas d'arbre
    prog->option[1]=INTERDIT;
    //longueur
    prog->option[2]=INTERDIT;
    //OPTION
    prog->option[3]=INTERDIT;
  }
  //letfreq
  else if((choix==1)&&(prog->option[0]==NON)){
    if((num+1)==argc){
      fprintf(stderr,"Argument manquant : Fichier de frequence\n");
      freeoptions(prog);
      exit(1);
    }
    else if(compair(prog->seq,argv,(num+1))==OUI){
      fprintf(stderr,"Argument manquant: Fichier de frequence\n");
      freeoptions(prog);
      exit(1);
    }
    prog->ficfreq=(char*)malloc(strlen(argv[num+1])*sizeof(char));
    if(prog->ficfreq==NULL){
      fprintf(stderr,"Erreur de choix du Fichier de frequence\n");
      freeoptions(prog);
      exit(1);
    }
    strcpy(prog->ficfreq,argv[num+1]);
    printf("-Lecture du fichier de frequence : %s\n",prog->ficfreq);
    prog->option[0]=2;
    //longueur
    prog->option[2]=INTERDIT;
  }
  //arbre
  else if((choix==2)&&(prog->option[1]==NON)){
    prog->option[1]=1;
    printf("-Contruction de l'arborescence\n");
    //pas de frequence
    prog->option[0]=INTERDIT;
  }
  //letarbr
  else if((choix==3)&&(prog->option[1]==NON)){
    if((num+1)==argc){
      fprintf(stderr,"Argument manquant : Fichier d'arborescence\n");
      freeoptions(prog);
      exit(1);
    }
    else if(compair(prog->seq,argv,(num+1))==OUI){
      fprintf(stderr,"Argument manquant: Fichier d'arborescence\n");
      freeoptions(prog);
      exit(1);
    }
    //+1 bizzare
    prog->ficarbre=(char*)malloc((strlen(argv[num+1])+1)*sizeof(char));
    if(prog->ficarbre==NULL){
      fprintf(stderr,"Erreur de choix du fichier d'arborescence\n");
      freeoptions(prog);
      exit(1);
    }
    strcpy(prog->ficarbre,argv[num+1]);
    printf("-Lecture du fichier d'arborescence : %s\n",prog->ficarbre);
    //arbre
    prog->option[1]=2;
    //longueur
    prog->option[2]=INTERDIT;
    //OPTION
    prog->option[3]=INTERDIT;
  }
  //longueur
  else if((choix==9)&&(prog->option[2]==NON)){
    if((num+1)==argc){
      fprintf(stderr,"Argument manquant : Longueur\n");
      freeoptions(prog);
      exit(1);
    }
    prog->option[2]=atoi(argv[num+1]);
    deflongueur(prog->option[2]);
  }
  //OPTION : DSSP STRIDE KAKSI SEGNO
  else if((choix>=5&&choix<=8)&&(prog->option[3]==NON)){
    prog->option[3]=choix;
    printf("-Lecture des donnees issues de %s\n",prog->seq[prog->option[3]]);
  }
  //ficEntre
  else if((choix==10)&&(prog->option[4]==NON)){ 
    if((num+1)==argc){
      fprintf(stderr,"Argument manquant : fichier d'entre\n");
      freeoptions(prog);
      exit(1);
    }
    else if(compair(prog->seq,argv,(num+1))==OUI){
      fprintf(stderr,"Argument manquant: fichier d'entre\n");
      freeoptions(prog);
      exit(1);
    }
    prog->entre=(char*)malloc((strlen(argv[num+1])+1)*sizeof(char));
    if(prog->entre==NULL){
      fprintf(stderr,"Erreur de choix du fichier d'entre\n");
      freeoptions(prog);
      exit(1);
    }
    strcpy(prog->entre,argv[num+1]);
    printf("-Fichier d'entree : %s\n",prog->entre);
    prog->option[4]=OUI;
  }
  //ficSortie
  else if((choix==4)&&(prog->option[5]==NON)){
    if((num+1)==argc){
      fprintf(stderr,"Argument manquant : fichier de sortie\n");
      freeoptions(prog);
      exit(1);
    }
    else if(compair(prog->seq,argv,(num+1))==OUI){
      fprintf(stderr,"Argument manquant: fichier de sortie\n");
      freeoptions(prog);
      exit(1);
    }
    prog->sortie=(char*)malloc((strlen(argv[num+1])+1)*sizeof(char));
    if(prog->sortie==NULL){
      fprintf(stderr,"Erreur de choix du dossier : sortie/nomfichier\n");
      freeoptions(prog);
      exit(1);
    }
    strcpy(prog->sortie,argv[num+1]);
    printf("-Fichier de sortie : %s\n",prog->sortie);
    prog->option[5]=OUI;
  }
  else{
    if(choix==9)  fprintf(stderr,"Erreur : choix de la longueur non demandee\n");
    else if(choix==10) fprintf(stderr,"Erreur choix du fichier non demande\n");
    else fprintf(stderr,"Erreur de choix : %d\n",(choix+2));
    optiondisp(prog->seq);
    freeoptions(prog);
    exit(1);
  }
}
//verification des options
void verifoption(poptions prog){
  int i;
  //si c pas letfreq ou letarbr
  if(prog->option[1]!=2&&prog->option[0]!=2){
    for(i=2;i<5;i++){
    //ficentre longueur option doivent avoir ete defini
      if(prog->option[i]==NON){
	fprintf(stderr,"L'option n°%d n'a pas ete selectionne et est necessaire (n°2 -Longueur n°3 -OPTION n°4 -FicEntre)\n",i);
	optiondisp(prog->seq);
	freeoptions(prog);
	exit(1);
      }
    }
    //si l'option de sortie n'a pas ete selectionnee on  donne l'adresse du fichier d'entree pour la sortie
    if(prog->option[5]==NON){
      prog->sortie=(char *)malloc(strlen(prog->entre)*sizeof(char));
      if(prog->sortie==NULL){
	fprintf(stderr,"Erreur d'identification du fichier de sortie");
	exit(1);
      }
      strcpy(prog->sortie,prog->entre);
    }
  }
  //si c le cas 
  else if((prog->option[5]==NON)&&(prog->option[4]==OUI)){
    prog->sortie=(char *)malloc(strlen(prog->entre)*sizeof(char));
    if(prog->sortie==NULL){
      fprintf(stderr,"Erreur d'identification du fichier de sortie");
      exit(1);
    }
    strcpy(prog->sortie,prog->entre);
  }
  //s'il y a pas de fichier d'entre c pas bon
  else if((prog->option[5]==NON)&&(prog->option[4]==NON)){
    fprintf(stderr,"Erreur : Indiquez un fichier de sortie\n");
    optiondisp(prog->seq);
    freeoptions(prog);
    exit(1);
  }
  //cas de l'option rien : activation du calcul des freq et de l'arborescence
  if((prog->option[0]==NON)&&(prog->option[1]==NON)){
    prog->option[0]=OUI;
    prog->option[1]=OUI;
  }
}
//identification des options choisies par l'utilisateur
void defoption(char **argv,int nbarg,poptions prog){
  int i,choix=-2,nbchoix=0; 
  for(i=1;i<nbarg;i++){
    choix=identoption(prog->seq,argv,i);
    if(choix!=-1){
      notifoption(choix,prog,argv,i,nbarg);
      nbchoix++;
    }
  }
  if(nbchoix!=0) verifoption(prog);
  else {
    fprintf(stderr,"Erreur : Aucun choix n'a ete realise\n");
    optiondisp(prog->seq);
    exit(1);
  }
}
