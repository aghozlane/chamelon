#ifndef _FREQUENCES_H
#define _FREQUENCES_H

#include "options.h"

typedef struct statistique{
  char *alphabet,*seq, pdbseq[5];
  long int nbword,nbworddiff,Nbfragments,longueur,Nblines;
  int maillon_suppr,Nbfragdiff2,Nbcameleon;
  int wordlen,*structseq,*cas,structcameleon[2];
  double Nbaa[21], FreqStructure[4],AAstruct[20][3],propension[20][3],Nbaatotaux;
  float Otherdata[7],accessibilite,*accesslet;
  double bfactor[3],*bfactorlet;
  FILE *access,*bfact;
}statistique,*pstatistique;
extern pstatistique initstruct(pstatistique chamelon,int longueur);
extern void freechamelon(pstatistique chamelon);
extern int freqStruct(int structure, pstatistique chamelon);
extern int freqAA(char a, int structure,pstatistique chamelon);
extern void calcAccessibilite(char a, pstatistique chamelon);
extern void calcBfactor(char a,double bfactoresult,pstatistique chamelon);
extern void affichFREQ(pstatistique chamelon, char *sortie,char gnufilename[100]);
extern void sortieFREQ(pstatistique chamelon, char *sortie, char gnufilename[100]);
extern pstatistique letFREQ(pstatistique chamelon,char *lectfreq);
#endif
