#include "arborescence.h"

/****************** ARBORESCENCE****************/
//alloue la memoire pour l'accessibilite de la sequence
float *newaccess(float *accessibilite,int wordlen){
  int i;
  float *access=(float*)malloc(wordlen*sizeof(float));
  if(!access){
    fprintf(stderr,"Erreur allocation memoire access\n");
    exit(1);
  }
  for(i=0;i<wordlen;i++){
    access[i]=accessibilite[i];
  }
  return access;
}
//alloue la memoire pour le bfacteur
double *newbfacto(double *bfactor, int wordlen){
  int i;
  double *bfact=(double*)malloc(wordlen*sizeof(double));
  if(!bfact){
    fprintf(stderr,"Erreur allocation memoire bfact\n");
    exit(1);
  }
  for(i=0;i<wordlen;i++){
    bfact[i]=bfactor[i];
  }
  return bfact;
}
//alloue la memoire pour les pdb correspondant
char *pdbnewseq(char pdb[5]){
  int i;
  char *pdbseq=(char*)malloc(5*sizeof(char));
  if(!pdbseq){
    fprintf(stderr,"Erreur allocation memoire\n");
    exit(1);
  }
  for(i=0;i<5;i++) pdbseq[i]=pdb[i];
  return pdbseq;
} 

//alloue la memoire pour une nouvelle structure identifiee
int *newstructseq(int *struc, int wordlen){
  int k;
  int *structure=(int*)malloc(wordlen*sizeof(int));
  if(!structure){
    fprintf(stderr,"Erreur allocation memoire\n");
    exit(1);
  }
  for(k=0;k<wordlen;k++){
    structure[k]=struc[k];
  }
  return structure;
}

//CREATION D'UNE FEUILLE DE L'ARBRE
parbre creation(char *seq, int *let, int *struc, int wordlen, char pdb[5], float *accessibilite,double *bfactor){
  parbre ll;
  int i;
  ll=(parbre)malloc(sizeof(arbre));
  ll->gauche=NULL;
  ll->droite=NULL;
  ll->prec=NULL;
  ll->seq=seq;
  //identification du PDB de la sequence
  ll->pdbseq=(char**)malloc(sizeof(char*));
  if(!ll->pdbseq){
    fprintf(stderr,"Erreur allocation memoire pdbseq\n");
    exit(1);
  }
  ll->pdbseq[0]=pdbnewseq(pdb);
  ll->lettre=let;
  
  //indetification de la structure de la sequence
  ll->structseq=(int**)malloc(sizeof(int*));
  if(!ll->structseq){
    fprintf(stderr,"Erreur allocation memoire structseq\n");
    exit(1);
  }
  ll->structseq[0]=newstructseq(struc,wordlen);
  
  //identification de l'accessibilite de la sequence
  ll->access=(float**)malloc(sizeof(float*));
  if(!ll->access){
    fprintf(stderr,"Erreur allocation memoire access\n");
    exit(1);
  }
  ll->access[0]=newaccess(accessibilite,wordlen);
  
  
  //identification du bfacteur de la sequence
  ll->bfacto=(double**)malloc(sizeof(double*));
  if(!ll->bfacto){
    fprintf(stderr,"Erreur allocation memoire bfacto\n");
    exit(1);
  }
  ll->bfacto[0]=newbfacto(bfactor,wordlen);
  /*for(i=0;i<(wordlen);i++){
    printf("accessibilitestruct %f bfactorstruct %lf\n",ll->access[0][i],ll->bfacto[0][i]);
    }*/
  ll->occurence=1;
  return ll;
}
//liberation de la memoire occupe par une structure
/*void freestruct(parbre a){
  if(a->lettre!=NULL) free(a->lettre);
  if(a->seq!=NULL) free(a->seq);
  if(a->structseq[0]!=NULL) free(a->structseq[0]);
  if(a->structseq!=NULL) free(a->structseq);
  if(a->pdbseq[0]!=NULL) free(a->pdbseq[0]);
  if(a->pdbseq!=NULL) free(a->pdbseq);
  free(a);
  }*/

//liberation de la memoire occupe par une structure finaliseee
void freestruct2(parbre a){
  int i;
  if(a->lettre!=NULL) free(a->lettre);
  if(a->seq!=NULL) free(a->seq);
  for(i=0;i<a->occurence;i++){
      if(a->structseq[i]!=NULL) free(a->structseq[i]);
      if(a->pdbseq[i]!=NULL) free(a->pdbseq[i]);
      if(a->bfacto[i]!=NULL) free(a->bfacto[i]);
      if(a->access[i]!=NULL) free(a->access[i]);
  }
  if(a->structseq!=NULL) free(a->structseq);
  if(a->pdbseq!=NULL) free(a->pdbseq);
  if(a->bfacto!=NULL) free(a->bfacto);
  if(a->access!=NULL) free(a->access);
  free(a);
}
//supression d'arbre binaire de recherche selon l'algorithme d'huffman
void suppr_arbre(parbre a, pstatistique chamelon){
  //printf("%s\n", a->seq);
  if(a->gauche!=NULL && a->droite!=NULL){
    suppr_arbre(a->gauche,chamelon);
    suppr_arbre(a->droite,chamelon);
    freestruct2(a);
    (chamelon->maillon_suppr)++;
    a=NULL;
  }
  else if(a->gauche!=NULL){
    suppr_arbre(a->gauche,chamelon);
    freestruct2(a);
    (chamelon->maillon_suppr)++;
    a=NULL;
  }
  else if(a->droite!=NULL){
    suppr_arbre(a->droite,chamelon);
    freestruct2(a);
    a=NULL;
    (chamelon->maillon_suppr)++;
  }
  else{
    freestruct2(a);
    a=NULL;
    (chamelon->maillon_suppr)++;
  }
}

/******************CONSTRUCTION DE L'ARBORESCENCE************/
//Positionnement de la sequence dans l'arbre
void positionnement(parbre l,parbre a,int **newptr,char **newpdbseq,float **newaccessibilite,double **newbfact,int *cas, int i,pstatistique chamelon){
  int k,j;
  //cas 1 ou la lettre est inferieure A face a L par exemple
  if(l!=NULL && a->lettre[i]<l->lettre[i]){
    //s'il n'y a rien a droite
    if(!l->droite){
      //printf("%s a droite %c < %c %d\n", a->seq, a->lettre[i],l->lettre[i], i);      
      l->droite=a;
      a->prec=l;
      assert(l->droite);
      (chamelon->nbworddiff)++;
    }
    //sinon on avance a droite
    else {
      //printf("%s tourne a droite  %c < %c %d\n", a->seq, a->lettre[i],l->lettre[i],i);
      positionnement(l->droite,a,newptr,newpdbseq,newaccessibilite,newbfact,cas,0,chamelon);
    }
  }
  //cas 2 ou la lettre est superieure L face a A
  else if(l!=NULL && a->lettre[i]>l->lettre[i]){
    if(!l->gauche){
      //printf("%s a gauche %c > %c %d\n", a->seq, a->lettre[i],l->lettre[i],i);
      l->gauche=a;
      a->prec=l;
      assert(l->gauche);
      (chamelon->nbworddiff)++;
    }
    else{
      //printf("%s tourne a gauche %c > %c %d\n", a->seq, a->lettre[i],l->lettre[i],i);
      positionnement(l->gauche,a,newptr,newpdbseq,newaccessibilite,newbfact,cas, 0,chamelon);
    }
  }
  //cas de l'egalite des lettre
  else if(l!=NULL && a->lettre[i]==a->lettre[i]){
    //on se deplace a la lettre suivante
    //printf("%s a egalite lettre :%c  == lettre %c, %d\n", a->seq,a->seq[i],l->seq[i],i);
    if(i<(*cas)) positionnement(l,a,newptr,newpdbseq,newaccessibilite,newbfact,cas, i+1,chamelon);
    //cas ou les mots sont identiques
    else {
      l->occurence++;
      //printf("occurence %d\n",l->occurence);
      //printf("l : %s, a : %s,occurence %d\n",l->seq,a->seq,l->occurence);
      //reallocation des donnes de la structure de la sequence
      newptr =(int **)realloc(l->structseq,l->occurence*sizeof(int*));
      if(!newptr){
	fprintf(stderr,"ERREUR REALLOC newptr\n");
	exit(1);
      }
      l->structseq=newptr;
      //-1 comme on est dans un tableau
      //+1 comme le cas s'arrete juste au dernier
      l->structseq[l->occurence-1]=newstructseq(a->structseq[0],chamelon->wordlen);
      
      //reallocation des donnes du PDB de la sequence
      newpdbseq =(char **)realloc(l->pdbseq,l->occurence*sizeof(char*));
      if(!newpdbseq){
	fprintf(stderr,"ERREUR REALLOC newpdbseq\n");
	exit(1);
      }
      l->pdbseq=newpdbseq;
      l->pdbseq[l->occurence-1]=pdbnewseq(a->pdbseq[0]);
      
      //reallocation des donnes pour l'accessibilite
      /*printf("accessibilite avant\n");
      for(k=0;k<l->occurence-1;k++){
	for(j=0;j<chamelon->wordlen;j++){
	  printf("%f ",l->access[k][j]);
	}
	printf("\n");
	}*/
      newaccessibilite=(float **)realloc(l->access,l->occurence*sizeof(float*));
      if(newaccessibilite==NULL){
	fprintf(stderr,"ERREUR REALLOC newaccess\n");
	exit(1);
      }
      l->access=newaccessibilite;
      l->access[l->occurence-1]=newaccess(a->access[0],chamelon->wordlen);
      /*printf("accessibilite apres\n");
      for(k=0;k<l->occurence;k++){
	for(j=0;j<chamelon->wordlen;j++){
	  printf("%f ",l->access[k][j]);
	}
	printf("\n");
	}*/
      //reallocation des donnees pour le bfactor
      /*printf("bfactor avant \n");
      for(k=0;k<l->occurence-1;k++){
	for(j=0;j<chamelon->wordlen;j++){
	  printf("%lf ",l->bfacto[k][j]);
	}
	printf("\n");
	}*/
      newbfact=(double **)realloc(l->bfacto,l->occurence*sizeof(double*));
      if(!newbfact){
	fprintf(stderr,"ERREUR REALLOC newbfact\n");
	exit(1);
      }
      l->bfacto=newbfact;
      l->bfacto[l->occurence-1]=newbfacto(a->bfacto[0],chamelon->wordlen);
      /*printf("bfactor apres\n");
      for(k=0;k<l->occurence;k++){
	for(j=0;j<chamelon->wordlen;j++){
	  printf("%lf ",l->bfacto[k][j]);
	}
	printf("\n");
	}*/
      freestruct2(a);
    }
  }
  else fprintf(stderr,"ca bug dans la structure");
}
//constructeur de la structure
parbre buildtree(char *seq, pstatistique chamelon, parbre l){
  //a modifier
  int i, *word=(int *)malloc(chamelon->wordlen*sizeof(int));
  int **newptr;
  char **newpdbseq;
  float **newaccessibilite;
  double **newbfact;
  parbre p,m=l;
  lettre b;
  
  //printf("seq soumise : %s\n", seq);
  for(i=0;i<chamelon->wordlen;i++){
    b=seq[i];
    word[i]=b;
  }
  if(l==NULL){ 
    l=creation(seq, word, chamelon->structseq, chamelon->wordlen,chamelon->pdbseq,chamelon->accesslet,chamelon->bfactorlet);
    //printf("seq au milieu : %s \n", l->seq);
    (chamelon->nbworddiff)++;
  }  
  else{
    p=creation(seq,word,chamelon->structseq,chamelon->wordlen,chamelon->pdbseq,chamelon->accesslet,chamelon->bfactorlet);
    assert(p);
    positionnement(m,p,newptr,newpdbseq,newaccessibilite,newbfact,chamelon->cas,0,chamelon);
  }
  //printf("MILIEU : %s \n", l->seq);
  return l;
}
/***************PARCOURS EN PROFONDEUR************/
//sequence cameleon
void cameleon(FILE *fic6,parbre a,pstatistique chamelon){
  int i,j,k, alpha=NON, beta=NON,sortie=NON, residus_alpha=0, residus_beta=0;
  //cas simple du alpha et du beta
  for(i=0;i<(a->occurence);i++){
    residus_alpha=0; residus_beta=0;
    for(j=0;j<(chamelon->wordlen);j++){
      if((a->structseq[i][j])<3) residus_alpha++;
      else if((a->structseq[i][j])==7) residus_beta++;
    }
    if(residus_alpha==(chamelon->wordlen)){ 
      alpha=OUI;
      (chamelon->structcameleon[0])++;
    }
    else if(residus_beta==(chamelon->wordlen)){
      beta=OUI;
      (chamelon->structcameleon[1])++;
    }
    if(alpha==OUI &&beta==OUI){
      //pour ne la compter qu'une fois
      if(sortie==NON){
	(chamelon->Nbcameleon)++;
      }      
      sortie=OUI;
    }
  }
  if(sortie==OUI){
    //printf("%d\t%d\n",(chamelon->structcameleon[0]),(chamelon->structcameleon[1]));
    fprintf(fic6,"%s\t%d\t%d\n",(a->seq),(chamelon->structcameleon[0]),(chamelon->structcameleon[1]));
    //ecriture de l'accessibilite et du bfactor
    for(i=0;i<(a->occurence);i++){
      residus_alpha=0; residus_beta=0;
      for(j=0;j<chamelon->wordlen;j++){
	if((a->structseq[i][j])<3) residus_alpha++;
	else if((a->structseq[i][j])==7) residus_beta++;
      }
      if((residus_alpha==chamelon->wordlen)||(residus_beta==chamelon->wordlen)){
	for(k=0;k<chamelon->wordlen;k++){
	  fprintf(chamelon->access,"%c\t%f\n",a->seq[k],a->access[i][k]);
	  fprintf(chamelon->bfact,"%c\t%lf\n",a->seq[k],a->bfacto[i][k]);
	}
      }
    }
  }
  chamelon->structcameleon[0]=0;
  chamelon->structcameleon[1]=0;
}

//parcours des noeuds
void analyse_arbre(FILE *fic6,FILE *fic7,parbre a, pstatistique chamelon,int *nbarguments){
  if(a->gauche!=NULL && a->droite!=NULL){
    analyse_arbre(fic6,fic7,a->gauche,chamelon,nbarguments);
    analyse_arbre(fic6,fic7,a->droite,chamelon,nbarguments);
    if(a->occurence>1){
      (chamelon->Nbfragdiff2)++;
      cameleon(fic6,a,chamelon);
    }
    sortiearborescence(fic7,a,chamelon,nbarguments);
  }
  else if(a->gauche!=NULL){
    analyse_arbre(fic6,fic7,a->gauche,chamelon,nbarguments);
    if(a->occurence>1){
      (chamelon->Nbfragdiff2)++;
      cameleon(fic6,a,chamelon);
    }
    sortiearborescence(fic7,a,chamelon,nbarguments);
  }
  else if(a->droite!=NULL){
    analyse_arbre(fic6,fic7,a->droite,chamelon,nbarguments);
    if(a->occurence>1){
      (chamelon->Nbfragdiff2)++;
      cameleon(fic6,a,chamelon);
    }
    sortiearborescence(fic7,a,chamelon,nbarguments);
  }
  else{
    if(a->occurence>1){
      (chamelon->Nbfragdiff2)++;
      cameleon(fic6,a,chamelon);
    }
    //fprintf(stderr,"Merde 4 av %d\n", (*nbarguments));
    sortiearborescence(fic7,a,chamelon,nbarguments);
    //fprintf(stderr,"Merde 4 ap\n");
  }
}

//fichier de sauvegarde de l'arborescence
void sortiearborescence(FILE *fic7,parbre a,pstatistique chamelon,int *nbarguments){
  int i,j;
  //fprintf(stderr,"Merde 5 ap %s %d\n",a->seq,a->occurence);
  (*nbarguments)+=fprintf(fic7,"%s\t%d\n",a->seq,a->occurence);
  //fprintf(stderr,"Merde 6 ap %d %d\n",a->occurence,(*nbarguments));
  fprintf(fic7,"structseq =\n");
  for(i=0;i<a->occurence;i++){
    for(j=0;j<chamelon->wordlen;j++){
      (*nbarguments)+=fprintf(fic7,"%d\t",a->structseq[i][j]);
      //fprintf(fic7,"%d\t",a->structseq[i][j]);
    }
    //fprintf(stderr,"%s\n",a->pdbseq[i]);
    fprintf(fic7,"%s\n",a->pdbseq[i]);
  } 
  //fprintf(stderr,"Merde 6 ap\n");
  fprintf(fic7,"\n");
  //printf("%s %d\n",a->seq,(*nbarguments));
  //fprintf(stderr,"Merde 7 ap\n");
}


/***********LECTURE DU FICHIER DE SAUVEGARDE************/
//CREATION D'UNE FEUILLE DE L'ARBRE
parbre creation2(char *seq,int occurence, int wordlen){
  parbre ll;
  int i,j;
  lettre b;
  ll=(parbre)malloc(sizeof(arbre));
  ll->gauche=NULL;
  ll->droite=NULL;
  ll->prec=NULL;
  ll->seq=seq;
  ll->lettre=(int*)malloc(wordlen*sizeof(int));
  for(i=0;i<wordlen;i++){
    b=seq[i];
    ll->lettre[i]=b;
  }
  ll->structseq=(int**)malloc(occurence*sizeof(int*));
  for(i=0;i<occurence;i++){
    ll->structseq[i]=(int*)malloc(wordlen*sizeof(int));
    for(j=0;j<wordlen;j++){
      ll->structseq[i][j]=0;
    }
  }
 
  ll->occurence=occurence;
  return ll;
}
//positionnement des maillons
void positionnement2(parbre l,parbre a, int *cas, int i){
    //cas 1 ou la lettre est inferieure A face a L par exemple
  if(l!=NULL && a->lettre[i]<l->lettre[i]){
    //s'il n'y a rien a droite
    if(!l->droite){
      //printf("%s a droite %c < %c %d\n", a->seq, a->lettre[i],l->lettre[i], i);      
      l->droite=a;
      a->prec=l;
      assert(l->droite);
    }
    //sinon on avance a droite
    else {
      //printf("%s tourne a droite  %c < %c %d\n", a->seq, a->lettre[i],l->lettre[i],i);
      positionnement2(l->droite, a, cas,0);
    }
  }
  //cas 2 ou la lettre est superieure L face a A
  else if(l!=NULL && a->lettre[i]>l->lettre[i]){
    if(!l->gauche){
      //printf("%s a gauche %c > %c %d\n", a->seq, a->lettre[i],l->lettre[i],i);
      l->gauche=a;
      a->prec=l;
      assert(l->gauche);
    }
    else{
      //printf("%s tourne a gauche %c > %c %d\n", a->seq, a->lettre[i],l->lettre[i],i);
      positionnement2(l->gauche, a, cas, 0);
    }
  }
  else if(l!=NULL && a->lettre[i]==a->lettre[i]){
    //printf("%s a egalite lettre :%c  == lettre %c, %d\n", a->seq,a->seq[i],l->seq[i],i);
    if(i<(*cas)) positionnement2(l, a,cas, i+1);
    //cas ou les lettres sont identiques normalement impossible
    else {
      fprintf(stderr,"On retrouve deux fois le meme mot, le fichier d'arborescence propose presente une erreur\n");
      exit(1);
    }
  }
}
//lecture de l'arborescence
parbre letarborescence(pstatistique chamelon,char *lectarbre){
  FILE *fic10;
  int nbarguments=0,nbargumentsfic=0,occurence=0,i,j,k;
  char *seq2;
  parbre p=NULL,l=NULL;
  
  if((fic10=fopen(lectarbre,"rt"))!=NULL){
    nbarguments+=fscanf(fic10,"nbword = %ld\n",&(chamelon->nbword));
    nbarguments+=fscanf(fic10,"nbworddiff = %ld\n",&(chamelon->nbworddiff));
    nbarguments+=fscanf(fic10,"nbfragments = %ld\n",&(chamelon->Nbfragments));
    nbarguments+=fscanf(fic10,"nblines = %ld\n",&(chamelon->Nblines));
    nbarguments+=fscanf(fic10,"wordlen = %d\n",&(chamelon->wordlen));
    
    if(nbarguments!=5){
      fprintf(stderr,"Erreur de lecture des donnes globales de frequence de l'arborescence, seulement %d arguments lus\n",nbarguments);
      freechamelon(chamelon);
      exit(1);
    }
    
    else{
      chamelon->seq=(char*)malloc(chamelon->wordlen*sizeof(char));
      chamelon->structseq=(int*)malloc(chamelon->wordlen*sizeof(int));
      (*chamelon->cas)=(chamelon->wordlen)-1;
      nbarguments=0;
    }
    
    for(k=0;k<(chamelon->nbworddiff)&&!feof(fic10);k++){
      //Lecture de la sequence
      for(i=0;i<chamelon->wordlen;i++){
	chamelon->seq[i]=fgetc(fic10);
	nbarguments++;
      }
      nbarguments+=fscanf(fic10,"\t%d\n",(&occurence));
      nbarguments+=2;
      fscanf(fic10,"structseq =\n");
      
      seq2=(char*)malloc(chamelon->wordlen*sizeof(char));
      strcpy(seq2,chamelon->seq);
      p=creation2(seq2,occurence,chamelon->wordlen);
      assert(p);
      for(i=0;i<p->occurence;i++){
	for(j=0;j<chamelon->wordlen;j++){
	  nbarguments+=fscanf(fic10,"%d\t",&(p->structseq[i][j]));
	  nbarguments+=1;
	}
	fscanf(fic10,"\n");
      }
      fscanf(fic10,"\n");
      printf("%s %d\n",p->seq,nbarguments);
      if(l==NULL) l=p;
      else positionnement2(l,p,chamelon->cas,0);
    }
  }
  else{
    fprintf(stderr,"Erreur d'ouverture du fichier d'arborescence : %s\n",lectarbre);
    freechamelon(chamelon);
    exit(1);
  }
  nbarguments+=fscanf(fic10,"nbarguments = %d\n",&(nbargumentsfic));
  fprintf(stderr,"nbargumentsfic %d\n",nbargumentsfic);
  fclose(fic10);
  
  if(nbarguments!=nbargumentsfic){
    fprintf(stderr,"Probleme dans la lecture du fichier d'arborescence\n");
    fprintf(stderr,"Le nombre d'arguments lus : %d est different de ce lui qui a ete ecrit : %d\n",nbarguments,nbargumentsfic);
    //suppr_arbre(l,chamelon);
    //freechamelon(chamelon);
    
    //exit(1);
  }
  return l;
}
 
