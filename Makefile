DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

ALL = fasta2DB DB2fasta quiva2DB DB2quiva DBsplit DBdust Catrack DBshow DBstats DBrm DBmv \
      simulator fasta2DAM DAM2fasta DBdump rangen arrow2DB DB2arrow DBwipe DBtrim DBa2b DBb2a

all: $(ALL)

fasta2DB: fasta2DB.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o fasta2DB fasta2DB.c DB.c QV.c -lm

DB2fasta: DB2fasta.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DB2fasta DB2fasta.c DB.c QV.c -lm

quiva2DB: quiva2DB.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -DINTERACTIVE -o quiva2DB quiva2DB.c DB.c QV.c -lm

DB2quiva: DB2quiva.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DB2quiva DB2quiva.c DB.c QV.c -lm

DB2arrow: DB2arrow.c DB.c QV.c DB.h QV.h
	gcc $(CFLAGS) -o DB2arrow DB2arrow.c DB.c QV.c -lz

arrow2DB: arrow2DB.c DB.c QV.c DB.h QV.h
	gcc $(CFLAGS) -o arrow2DB arrow2DB.c DB.c QV.c -lz

DBsplit: DBsplit.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBsplit DBsplit.c DB.c QV.c -lm

DBtrim: DBtrim.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBtrim DBtrim.c DB.c QV.c -lm

DBdust: DBdust.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBdust DBdust.c DB.c QV.c -lm

Catrack: Catrack.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o Catrack Catrack.c DB.c QV.c -lm

DBshow: DBshow.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBshow DBshow.c DB.c QV.c -lm

DBdump: DBdump.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBdump DBdump.c DB.c QV.c -lm

DBstats: DBstats.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBstats DBstats.c DB.c QV.c -lm

DBrm: DBrm.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBrm DBrm.c DB.c QV.c -lm

DBmv: DBmv.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBmv DBmv.c DB.c QV.c -lm

simulator: simulator.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o simulator simulator.c DB.c QV.c -lm

rangen: rangen.c
	gcc $(CFLAGS) -o rangen rangen.c

fasta2DAM: fasta2DAM.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o fasta2DAM fasta2DAM.c DB.c QV.c -lm

DAM2fasta: DAM2fasta.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DAM2fasta DAM2fasta.c DB.c QV.c -lm

DBwipe: DBwipe.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBwipe DBwipe.c DB.c QV.c -lm

DBa2b: DBa2b.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBa2b DBa2b.c DB.c QV.c -lm

DBb2a: DBb2a.c DB.c DB.h QV.c QV.h
	gcc $(CFLAGS) -o DBb2a DBb2a.c DB.c QV.c -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f dazz.db.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf dazz.db.tar.gz README.md Makefile *.h *.c
