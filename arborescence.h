#ifndef _ARBORESCENCE_H
#define _ARBORESCENCE_H

#include "frequences.h"

typedef enum{L,A,G,V,E,S,D,T,K,I,R,P,N,F,Q,Y,H,M,C,W } lettre;

typedef struct arbre{
  int *lettre;
  char *seq,**pdbseq;
  int **structseq;
  float **access;
  double **bfacto;
  //occurence de la sequence
  int occurence;
  //balance du noeud
  //int balance;
  struct arbre *gauche;
  struct arbre *droite;
  struct arbre *prec;
}arbre, *parbre;
extern float *newaccess(float *accessibilite,int wordlen);
extern double *newbfacto(double *bfactor, int wordlen);
extern char *pdbnewseq(char pdb[5]);
extern int* newstructseq(int *struc, int wordlen);
extern parbre creation(char *seq, int *let, int *struc, int wordlen,char pdb[5],float *accessibilite,double *bfactor);
extern void freestruct(parbre a);
extern void freestruct2(parbre a);
extern void suppr_arbre(parbre a,pstatistique chamelon);
extern void positionnement(parbre l,parbre a,int **newptr,char **newpdbseq,float **newaccess,double **newbfact, int *cas, int i,pstatistique chamelon);
extern parbre buildtree(char *seq, pstatistique chamelon, parbre l);
extern void cameleon(FILE *fic6,parbre a,pstatistique chamelon);
extern void analyse_arbre(FILE *fic6,FILE *fic7,parbre a,pstatistique chamelon,int *nbarguments);
extern void sortiearborescence(FILE *fic7,parbre a,pstatistique chamelon,int *nbarguments);
extern parbre creation2(char *seq,int occurence, int wordlen);
extern void positionnement2(parbre l,parbre a, int *cas, int i);
extern parbre letarborescence(pstatistique chamelon,char *lectarbre);

#endif
