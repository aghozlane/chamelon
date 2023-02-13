#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char**argv){
  FILE *fic8;
  long int d;
  if((fic8=fopen(argv[1],"rt"))!=NULL){
    fscanf(fic8,"nbword = %ld\n",&(d));
    printf("%ld",d);
    fclose(fic8);
  }
  return 0;
}
