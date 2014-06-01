CFLAGS = -O4 -Wall -Wextra

all: fasta2DB DB2fasta DBsplit DBdust Catrack DBread DBstats DBrm simulator

fasta2DB: fasta2DB.c DB.c DB.h
	gcc $(CFLAGS) -o fasta2DB fasta2DB.c DB.c -lm

DB2fasta: DB2fasta.c DB.c DB.h
	gcc $(CFLAGS) -o DB2fasta DB2fasta.c DB.c -lm

DBsplit: DBsplit.c DB.c DB.h
	gcc $(CFLAGS) -o DBsplit DBsplit.c DB.c -lm

DBdust: DBdust.c DB.c DB.h
	gcc $(CFLAGS) -o DBdust DBdust.c DB.c -lm

Catrack: Catrack.c DB.c DB.h
	gcc $(CFLAGS) -o Catrack Catrack.c DB.c -lm

DBread: DBread.c DB.c DB.h
	gcc $(CFLAGS) -o DBread DBread.c DB.c -lm

DBstats: DBstats.c DB.c DB.h
	gcc $(CFLAGS) -o DBstats DBstats.c DB.c -lm

DBrm: DBrm.c DB.c DB.h
	gcc $(CFLAGS) -o DBrm DBrm.c DB.c -lm

simulator: simulator.c DB.c DB.h
	gcc $(CFLAGS) -o simulator simulator.c DB.c -lm

clean:
	rm -f fasta2DB DB2fasta DBsplit DBdust Catrack DBread DBstats DBrm simulator
	rm -f dazz.db.tar.gz

install:
	cp fasta2DB DB2fasta DBsplit DBdust Catrack DBread DBstats DBrm simulator ~/bin

package:
	make clean
	tar -zcf dazz.db.tar.gz README Makefile *.h *.c
