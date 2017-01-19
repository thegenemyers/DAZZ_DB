DEST_DIR = ~/bin

CFLAGS += -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing
LDLIBS = -lm

ALL = fasta2DB DB2fasta quiva2DB DB2quiva DBsplit DBdust Catrack DBshow DBstats DBrm simulator \
      fasta2DAM DAM2fasta DBdump rangen arrow2DB DB2arrow DBwipe

all: $(ALL)

fasta2DB: fasta2DB.c DB.c DB.h QV.c QV.h

DB2fasta: DB2fasta.c DB.c DB.h QV.c QV.h

quiva2DB: CPPFLAGS += -DINTERACTIVE
quiva2DB: quiva2DB.c DB.c DB.h QV.c QV.h

DB2quiva: DB2quiva.c DB.c DB.h QV.c QV.h

DB2arrow: LDLIBS = -lz
DB2arrow: DB2arrow.c DB.c QV.c DB.h QV.h

arrow2DB: LDLIBS = -lz
arrow2DB: arrow2DB.c DB.c QV.c DB.h QV.h

DBsplit: DBsplit.c DB.c DB.h QV.c QV.h

DBdust: DBdust.c DB.c DB.h QV.c QV.h

Catrack: Catrack.c DB.c DB.h QV.c QV.h

DBshow: DBshow.c DB.c DB.h QV.c QV.h

DBdump: DBdump.c DB.c DB.h QV.c QV.h

DBstats: DBstats.c DB.c DB.h QV.c QV.h

DBrm: DBrm.c DB.c DB.h QV.c QV.h

simulator: simulator.c DB.c DB.h QV.c QV.h

rangen: LDLIBS =
rangen: rangen.c

fasta2DAM: fasta2DAM.c DB.c DB.h QV.c QV.h

DAM2fasta: DAM2fasta.c DB.c DB.h QV.c QV.h

DBwipe: DBwipe.c DB.c DB.h QV.c QV.h


clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f dazz.db.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf dazz.db.tar.gz README.md Makefile *.h *.c
