#ifndef _OPTIONS_H
#define _OPTIONS_H

#include "projetChamelon.h"

typedef struct options{
  char **seq, *sortie,*entre,*ficfreq,*ficarbre;
  int option[7];
}options,*poptions;

extern poptions initoption(poptions prog);
extern char **allocseq(char **seq);
extern void freeoptions(poptions prog);
extern void optiondisp(char **seq);
extern int identoption(char **seq,char **argv,int num);
extern void deflongueur(int longueur);
extern int compair(char **seq,char **argv,int num);
extern void notifoption(int choix,poptions prog,char **argv,int num,int argc);
extern void verifoption(poptions prog);
extern void defoption(char **argv,int nbarg,poptions prog);
#endif
