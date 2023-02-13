#Amine GHOZLANE
CC=gcc
CFLAGS=-Wall

#creation de l'executable projetChamelon
all:projetChamelon

projetChamelon:projetChamelon.o options.o frequences.o arborescence.o
	$(CC) -o projetChamelon projetChamelon.o options.o frequences.o arborescence.o $(CFLAGS)

#creation des fichiers binaires utiles pour créer l'executable
projetChamelon.o:projetChamelon.c
	$(CC) -c projetChamelon.c

options.o:options.c options.h
	$(CC) -c options.c

frequences.o:frequences.c frequences.h
	$(CC) -c frequences.c

arborescence.o:arborescence.c arborescence.h
	$(CC) -c arborescence.c

#suppression des fichiers d'extensions .o
clean:
	rm *.o
