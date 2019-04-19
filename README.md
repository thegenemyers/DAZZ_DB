# The Dazzler Database Library

## _Author:  Gene Myers_
## _First:   July 17, 2013_
## _Current:   April 19, 2019_

For typeset documentation, examples of use, and design philosophy please go to
my [blog](https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide).

To facilitate the multiple phases of the dazzler assembler, we organize all the read
data into what is effectively a "database" of the reads and their meta-information.
The design goals for this data base are as follows:

1. The database stores the source Pacbio read information in such a way that it can
recreate the original input data, thus permitting a user to remove the
(effectively redundant) source files.  This avoids duplicating the same data,
once in the source file and once in the database.

2. The data base can be built up incrementally, that is new sequence data can be added
to the data base over time.

3. The data base flexibly allows one to store any meta-data desired for reads.  This
is accomplished with the concept of *tracks* that implementors can add as they
need them.

4. The data is held in a compressed form equivalent to the .dexta and .dexqv/.dexar
files of the data extraction module.

5. Quiver or Arrow information can be added separately from the sequence information
and later on if desired, but a database can only hold either Quiver or Arrow information,
but not both.  The Arrow or Quiver information can be removed from the database at any
time leaving a database just containing sequence information.

6. To facilitate job parallel, cluster operation of the phases of our assembler, the
data base has a concept of a *current partitioning* in which all the reads that
are over a given length and optionally unique to a well, are divided up into
*blocks* containing roughly a given number of bases, except possibly the last
block which may have a short count.  Often programs con be run on blocks or
pairs of blocks and each such job is reasonably well balanced as the blocks are
all the same size.  One must be careful about changing the partition during an
assembly as doing so can void the structural validity of any interim
block-based results.

A DB con contain the information needed by Quiver, or by Arrow, or neither, but
not both.  A DB containing neither Quiver or Arrow information is termed a
Sequence-DB (S-DB).  A DB with Quiver information is a Quiver-DB (Q-DB) and
a DB with Arrow information is an Arrow-DB (A-DB). All commands are aware of
the state of a DB and respond to options according to their type. 

A Dazzler DB consists of one named, *visible* file, e.g. FOO.db, and several
*invisible* secondary files encoding various elements of the DB.  The secondary files
are "invisible" to the UNIX OS in the sense that they begin with a "." and hence are
not listed by "ls" unless one specifies the -a flag.  We chose to do this so that when
a user lists the contents of a directory they just see a single name, e.g. FOO.db, that
is used to refer to the DB in commands.  The files associated with a database
named, say FOO,  are as follows:

* "FOO.db": a text file containing
  1. the list of input files added to the database so far, and
  2. how to partition the database into blocks (if the partition
     parameters have been set).

* ".FOO.idx": a binary "index" of all the meta-data about each read allowing, for
  example, one to randomly access a read's sequence (in the store
  ".FOO.bps").  It is 28N + 88 bytes in size where N is the number of
  reads in the database.

* ".FOO.bps": a binary compressed "store" of all the DNA sequences.  It is M/4 bytes
  in size where M is the total number of base pairs in the database.

* ".FOO.qvs": a binary compressed "store" of the 5 Pacbio quality value streams for
  the reads.  Its size is roughly 5/3M bytes depending on the
  compression acheived.  This file only exists if Quiver information has
  been added to the database.

* ".FOO.arw": a binary compressed "store" of the clipped pulse width stream for
  the reads.  Its size is roughly M/4 bytes.  This file only exists if Arrow information has
  been added to the database.

* ".FOO.\<track\>.[anno,data]": a *track* containing customized meta-data for each read.  For
  example, the DBdust command annotates low complexity intervals
  of reads and records the intervals for each read in two files
  .FOO.dust.anno & .FOO.dust.data.  Any kind of information
  about a read can be recorded, such as micro-sats, repeat
  intervals, corrected sequence, etc.  Specific tracks will be
  described as modules that produce them are released.

If one does not like the convention of the secondary files being invisible, then
un-defining the constant HIDE_FILES in DB.h before compiling the library, creates
commands that do not place a prefixing "." before secondary file names, e.g. FOO.idx
instead of .FOO.idx.  One then sees all the files realizing a DB when listing the
contents of a directory with ls.

While a Dazzler DB holds a collection of Pacbio reads, a Dazzler map DB or DAM holds
a collection of contigs from a reference genome assembly.  This special type of DB has
been introduced in order to facilitate the mapping of reads to an assembly and has
been given the suffix .dam to distinguish it from an ordinary DB.  It is structurally
identical to a .db except:

* there is no concept of quality values, and hence no .FOO.qvs or .FOO.arw file.

* every .fasta scaffold (a sequence with runs of N's between contigs estimating the
  length of the gap) is broken into a separate contig sequence in the DB and the
  header for each scaffold is retained in a new .FOO.hdr file.

* the original and first and last pulse fields in the meta-data records held in
  .FOO.idx, hold instead the contig number and the interval of the contig within
  its original scaffold sequence.

A map DB can equally well be the argument of many of the commands below that operate
on normal DBs.  In general, a .dam can be an argument anywhere a .db can, with the
exception of routines or optioned calls to routines that involve quality values, or
the special routines fasta2DAM and DAM2fasta that create a DAM and reverse said,
just like the pair fasta2DB and DB2fasta do for a normal DB.  So in general when we
refer to a database we are referring to either a DB or a DAM.

The command DBsplit sets or resets the current partition for a database which is
determined by 3 parameters: (i) the total number of basepairs to place in each block,
(ii) the minimum read length of reads to include within a block, and (iii) whether or
not to only include the longest (-l) or median (-m) read from a given well or all
reads from a well (-a) (NB: several reads of the same insert in a given well can be
produced by the Pacbio
instrument).  Note that the length and uniqueness parameters effectively select a
subset of the reads that contribute to the size of a block.  We call this subset the
*trimmed* data base.  Some commands operate on the entire database, others on the
trimmed database, and yet others have an option flag that permits them to operate on
either at the users discretion.  Therefore, one should note carefully to which version
of the database a command refers to.  This is especially important for any command that
identifies reads by their index (ordinal position) in the database.

Once the database has been split into blocks, the commands DBshow, DBstats, and DBdust
below and commands yet to come, such as the local alignment finder dalign, can take a
block or blocks as arguments.  On the command line this is indicated by supplying the
name of the DB followed by a period and then a block number, e.g. FOO.3.db or simply
FOO.3, refers to the 3'rd block of DB FOO (assuming of course it has a current
partition and said partition has a 3rd block).  One should note carefully that a block
is a contiguous range of reads such that once it is trimmed has a given size in base
pairs (as set by DBsplit).  Thus like an entire database, a block can be either
untrimmed  or trimmed and one needs to again be careful when giving a read index to
a command such as DBshow.

All programs add suffixes (e.g. .db) as needed.  The commands of the database library
are currently as follows:

```
1. fasta2DB [-v] <path:db> ( -f<file> | -i[<name>] | <input:fasta> ... )
```

Builds an initial data base, or adds to an existing database, either (a) the list of
.fasta files following the database name argument, or (b) the list of .fasta files in
\<file\> if the -f option is used, or (c) entries piped from the standard input if the
-i option is used.  If the DB is being created it is established as a Sequence-DB (S-DB)
otherwise its type is unchanged.  If a faux file name, \<name\>, follows the -i option
then all the
input received is considered to have come from a file by the name of \<name\>.fasta by
DB2fasta, otherwise it will be sent to the standard output by DB2fasta.  The SMRT cells
in a given named input (i.e. all sources other than -i without a name) can only be
added consecutively to the DB (this is checked by the command). The .fasta headers must
be in the "Pacbio" format (i.e. the output of the Pacbio tools or our dextract program)
and the well, pulse interval, and read quality are extracted from the header and kept
with each read record.  The headers may now also be those Pacbio outputs for CCS data
wherein the pulse range is replaced by the string "ccs" and in which case only the
well number is recorded.  If the files are being added to an existing database, and the
partition settings of the DB have already been set (see DBsplit below), then the
partitioning of the database is updated to include the new data.  A file may contain
the data from multiple SMRT cells provided the reads for each SMRT cell are consecutive
in the file.

```
2. DB2fasta [-vU] [-w<int(80)>] <path:db>
```

The set of .fasta files for the given DB are recreated from the DB exactly as they were
input.  That is, this is a perfect inversion, including the reconstitution of the
proper .fasta headers.  Because of this property, one can, if desired, delete the
.fasta source files once they are in the DB as they can always be recreated from it.
Entries imported from the standard input will be place in the faux file name given on
import, or to the standard output if no name was given. 
By default the output sequences are in lower case and 80 chars per line.  The -U option
specifies upper case should be used, and the characters per line, or line width, can be
set to any positive value with the -w option.

```
3. quiva2DB [-vl] <path:db> ( -f<file> | -i | <input:quiva> ... )
```

Adds .quiva streams to an existing DB "path".  The DB must either be an S-DB or a
Q-DB and upon completion the DB is a Q-DB.  The data comes from (a) the given .quiva
files on the command line, or (b) those in the file specified by the -f option, or
(c) the standard input if the -i option is given. The input files can be added incrementally
but must be added in the
same order as the .fasta files were and have the same root names, e.g. FOO.fasta and
FOO.quiva.  This is enforced by the program. With the -l option
set the compression scheme is a bit lossy to get more compression (see the description
of dexqv in the DEXTRACTOR module here).

```
4. DB2quiva [-vU] <path:db>
```

The set of .quiva files within the given Q-DB are recreated from the DB exactly as they
were input.  That is, this is a perfect inversion, including the reconstitution of the
proper .quiva headers.  Because of this property, one can, if desired, delete the
.quiva source files once they are in the DB as they can always be recreated from it.
Entries imported from the standard input will be placed in the faux file name given on
import, or to the standard output if no name was given. 
By .fastq convention each QV vector is output as a line without new-lines, and by
default the Deletion Tag entry is in lower case letters.  The -U option specifies
upper case letters should be used instead.

```
5. arrow2DB [-v] <path:db> ( -f<file> | -i | <input:arrow> ... )
```

Adds .arrow streams to an existing DB "path".  The DB must either be an S-DB or an
A-DB and upon completion the DB is an A-DB.  The data comes from (a) the given .arrow
files on the command line, or (b) those in the file specified by the -f option, or
(c) the standard input if the -i option is given. The input files can be added
incrementally but must be added in the
same order as the .fasta files were and have the same root names, e.g. FOO.fasta and
FOO.quiva.  This is enforced by the program.

```
6. DB2arrow [-v] [-w<int(80)>] <path:db>
```

The set of .arrow files within the given A-DB are recreated from the DB exactly as they
were input.  That is, this is a perfect inversion, including the reconstitution of the
proper .arrow headers.  Because of this property, one can, if desired, delete the
.arrow source files once they are in the DB as they can always be recreated from it.
Entries imported from the standard input will be placed in the faux file name given on
import, or to the standard output if no name was given. 
By default the output sequences are formatted 80 chars per line,
but the characters per line, or line width, can be
set to any positive value with the -w option.

```
7. fasta2DAM [-v] <path:dam> ( -f<file> | -i[<name>] | <input:fasta> ... )
```

Builds an initial map DB or DAM, or adds to an existing DAM, either (a) the list of
.fasta files following the database name argument, or (b) the list of .fasta files in
\<file\> if the -f option is used, or (c) entries piped from the standard input if the -i
option is used.  If a faux file name, \<name\>, follows the -i option then all the input
received is considered to have come from a file by the name of \<name\>.fasta by
DAM2fasta, otherwise it will be sent to the standard output by DAM2fasta.  Any .fasta
entry that has a run of N's in it will be split into separate "contig" entries and the
interval of the contig in the original entry recorded. The header for each .fasta entry
is saved with the contigs created from it.

```
8. DAM2fasta [-vU] [-w<int(80)>] <path:dam>
```

The set of .fasta files for the given map DB or DAM are recreated from the DAM
exactly as they were input. That is, this is a perfect inversion, including the
reconstitution of the proper .fasta headers and the concatenation of contigs with
the proper number of N's between them to recreate scaffolds.
Entries imported from the standard input will be place in the faux file name given on
import, or to the standard output if no name was given.  By default the output
sequences are in lower case and 80 chars per line. The -U option specifies upper case
should be used, and the characters per line, or line width, can be set to any positive
value with the -w option.

```
9. DBsplit [-aflm] [-x<int>] [-s<double(200.)>] <path:db|dam>
```

Divide the database \<path\>.db or \<path\>.dam conceptually into a series of blocks
referable to on the command line as \<path\>.1, \<path\>.2, ...  If the -x option is set
then all reads less than the given length are ignored, and if the -a option is not
set then secondary reads from a given well are also ignored.  The remaining reads,
constituting what we call the trimmed DB, are split amongst the blocks so that each
block is of size -s * 1Mbp except for the last which necessarily contains a smaller
residual.  The default value for -s is 200Mbp because blocks of this size can be
compared by our "overlapper" dalign in roughly 16Gb of memory.  The blocks are very
space efficient in that their sub-index of the master .idx is computed on the fly
when loaded, and the .bps and .qvs files (if a .db) of base pairs and quality values,
respectively, is shared with the master DB.  Any relevant portions of tracks
associated with the DB are also computed on the fly when loading a database block.
If the -f option is set, the split is forced regardless of whether or not the DB in
question has previously bin split, i.e. one is not interactively asked if they wish
to proceed.

By default, the primary read for each well consists of the longest subread of the insert
from the well.  By setting the -m parameter, the subread of median length becomes the
primary read instead.  One can at any later time change this back to the default longest
by splitting again with the -l parameter set.  The setting of the primary reads occurs
regardless of whether the -a parameter is set or not.

```
10. DBtrim [-af] [-x<int>] <path:db|dam>
```

Exactly like DBsplit except that it only resets the trimming parameters (and not the split
partition itself).

```
11. DBdust [-b] [-w<int(64)>] [-t<double(2.)>] [-m<int(10)>] <path:db|dam>
```

Runs the symmetric DUST algorithm over the reads in the untrimmed DB \<path\>.db or
\<path\>.dam producing a track .\<path\>.dust[.anno,.data] that marks all intervals of low
complexity sequence, where the scan window is of size -w, the threshold for being a
low-complexity interval is -t, and only low-complexity intervals of size greater than -m are
recorded.  If the -b option is set then the definition of low complexity takes into
account the frequency of a given base.  The command is incremental if given a DB to
which new data has been added since it was last run on the DB, then it will extend
the track to include the new reads.  It is important to set this flag for genomes with
a strong AT/GC bias, albeit the code is a tad slower.  The dust track, if present,
is understood and used by DBshow, DBstats, and dalign.

DBdust can also be run over an untriimmed DB block in which case it outputs a track
encoding where the trace file names contain the block number, e.g. .FOO.3.dust.anno
and .FOO.3.dust.data, given FOO.3 on the command line.  We call this a *block track*.
This permits job parallelism in block-sized chunks, and the resulting sequence of
block tracks can then be merged into a track for the entire untrimmed DB with Catrack.

```
12. Catrack [-vfd] <path:db|dam> <track:name> ...
```

Find all block tracks of the form .\<path\>.#.\<track\>... and concatenate them into a single
track, .\<path\>.\<track\>..., for the given DB or DAM.  Do so for each track name present on
the command line.    The block track files must all
encode the same kind of track data (this is checked), and the files must exist for
block 1, 2, 3, ... up to the last block number.  If the -f option is set, then the
concatenation takes place regardless of whether or not the single, combined track
already exists or not.  If the -d option is set then every block track is removed after
the successful construction of the combined track.

```
13. DBshow [-unqaUQA] [-w<int(80)>] [-m<mask>]+
                      <path:db|dam> [ <reads:FILE> | <reads:range> ... ]
```

Displays the requested reads in the database \<path\>.db or \<path\>.dam.  By default the
command applies to the trimmed database, but if -u is set then the entire DB is used.
If no read arguments are given then every read in the database or database block is
displayed.  Otherwise the input file or the list of supplied integer ranges give the
ordinal positions in the actively loaded portion of the db.  In the case of a file, it
should simply contain a read index, one per line.  In the other case, a read range is
either a lone integer or the symbol $, in which case the read range consists of just
that read (the last read in the database if $).  One may also give two positive
integers separated by a dash to indicate a range of integers, where again a $
represents the index of the last read in the actively loaded db.  For example,
1 3-5 $ displays reads 1, 3, 4, 5, and the last read in the active db.  As another
example, 1-$ displays every read in the active db (the default).

By default a .fasta file of the read sequences is displayed.  If the -q option is
set and the DB is a Q-DB, then the QV streams are also displayed in a non-standard
modification of the fasta format.
Similarly, if the -a option is set and the DB is an A-DB, then the pulse width stream is
also displayed in a non-standard format.
If the -n option is set then the DNA sequence is *not* displayed.
If the -Q option is set then a .quiva file of the selected reads is displayed and
all other options except -u and -U are ignored.  If the -A option is set then a .arrow
file of the selected reads is displayed and all other options except -u and -w are ignored.

If one or more masks are set with the -m option then the track intervals are also
displayed in an additional header line and the bases within an interval are displayed
in the case opposite that used for all the other bases.  By default the output
sequences are in lower case and 80 chars per line.  The -U option specifies upper
case should be used, and the characters per line, or line width, can be set to any
positive value with the -w option.

The .fasta, .quiva, and .arrow files that are output can be used to build a new DB with
fasta2DB, quiva2D, and arrow2DB, giving one a simple way to make a DB of a subset of
the reads for testing purposes.

```
14. DBdump [-rhsaqif] [-uU] [-m<mask>]+
                      <path:db|dam> [ <reads:FILE> | <reads:range> ... ]
```

Like DBshow, DBdump allows one to display a subset of the reads in the DB and select
which information to show about them including any mask tracks.  The difference is
that the information is written in a very simple "1-code" ASCII format that makes it
easy for one to read and parse the information for further use.  The option flags determine
which items of information are output as follows:

* -r requests that each read number be displayed in an R-line (see below, useful if only a
subset of reads is requested).

* -h requests the header information be output as the source file name on an H-line, the
 If the -d option is set then every block track is removed after the successful construction of the combined track.well # and pulse range on an L-line, and optionally the quality of the read if given on a Q-line.

* -s requests the sequence be output on an S-line.

* -a requests the Arrow information be output as a pulse-width string on an A-line and
the 4 SNR channel values on an N-line,

* -q requests that the 5 Quiver quality streams be output on d-, c-, i-, m-, and s-lines.

* -i requests that the intrinsic quality values be output on an I-line.

* -f requests the source file name is output just before the first read data in the file on a F-line.

* -m\<track\> requests that mask \<track\> be output on a T-line.

Set -u if you want data from the untrimmed database (the default is trimmed) and
set -U if you'd like upper-case letter used in the DNA sequence strings.

The format is very simple.  A requested unit of information occurs on a line.  The
first character of every line is a "1-code" character that tells you what information
to expect on the line.  The rest of the line contains the information where each item is
separated by a single blank space.  Strings are output as first an integer giving the
length of the string, a blank space, and then the string terminated by a new-line.
Intrinsic quality values are between 0 and 50, inclusive, and a vector of said are
displayed as an alphabetic string where 'a' is 0, 'b' is '1', ... 'z' is 25, 'A' is
26, 'B' is 27, ... and 'Y' is 50.
The set of all possible lines is as follows:

```
    R #              - read number
    H # string       - original file name string (header)
    L # # #          - location: well, pulse start, pulse end
    Q #              - quality of read (#/1000)
    N # # # #        - SNR of ACGT channels (#/100)
    Tx #n (#b #e)^#n - x'th track on command line, #n intervals all on same line
    S # string       - sequence string
    A # string       - arrow pulse-width string
    I # string       - intrinsic quality vector (as an ASCII string)
    F # string       - name of source file of following data
    d # string       - Quiva deletion values (as an ASCII string)
    c # string       - Quiva deletion character string
    i # string       - Quiva insertion value string
    m # string       - Quiva merge value string
    s # string       - Quiva substitution value string
    + X #            - Total amount of X (X = H or S or I or F or R or M or T#)
    @ X #            - Maximum amount of X (X = H or S or I or F or T#)
```

1-code lines that begin with + or @ are always the first lines in the output.  They
give size information about what is contained in the output.  That is '+ X #' gives
the number of reads (X=R), the number of masks (X=M), or the total number of
characters in all headers (X=H), sequences (X=S), intrinsic quality vectors (X=I),
file names (X=F), or track (X=T#).  And '@ X #' gives the maximum number of
characters in any header (X=H), sequence (X=S), intrincic quality vector (X=I),
names (X=F), or track (X=T#).  The size numbers for the Quiva strings and
Arrow pulse width strings are identical to that for the sequence as they are all of
the same length for any given entry.

```
15a. DBa2b
15b. DBb2a
```

Pipes (stdin to stdout) that convert an ASCII output produced by DBdump into a compressed
binary representation (DBa2b) and vice verse (DBb2a).  The idea is to save disk space by
keeping the dumps in a more compessed format.


```
16. DBstats [-nu] [-b<int(1000)] [-m<mask>]+ <path:db|dam>
```

Show overview statistics for all the reads in the trimmed data base \<path\>.db or
\<path\>.dam, including a histogram of read lengths where the bucket size is set
with the -b option (default 1000).  If the -u option is given then the untrimmed
database is summarized.  If the -n option is given then the histogran of read lengths
is not displayed.  Any track such as a "dust" track that gives a series of
intervals along the read can be specified with the -m option in which case a summary
and a histogram of the interval lengths is displayed.

```
17. DBrm [-v] <path:db|dam> ...
```

Delete all the files for the given data bases.  Do not use rm to remove a database, as
there are at least two and often several secondary files for each DB including track
files, and all of these are removed by DBrm.
If the -v option is set then every file deleted is listed.

```
18. DBmv [-v] <old:db|dam> <new:db|dam>
```

Rename all the files for the data base old to use the new root.
If the -v option is set then every file move is displayed.

```
19. DBwipe <path:db|dam>
```

Delete any Arrow or Quiver data from the given databases.  This removes the .arw or
.qvs file and resets information in the .idx file containing information for Arrow
or Quiver.  Basically, converts an A-DB or Q-DB back to a simple S-DB.

```
20.  simulator <genome:dam> [-CU] [-m<int(10000)>] [-s<int(2000)>] [-e<double(.15)]
                                  [-c<double(50.)>] [-f<double(.5)>] [-x<int(4000)>]
                                  [-w<int(80)>] [-r<int>] [-M<file>]
```

In addition to the DB commands we include here, somewhat tangentially, a simple
simulator that generates synthetic reads over a given genome reference contained in a
supplied .dam DB.  The simulator first reconstitutes the scaffolds of the reference
genome and fills in their gaps (a run of N's in .fasta format indicating the estimate
gap length) with a random sequence that follows the base distribution of the contigs.
It will then sample reads from these scaffold sequences.

The simulator generates sample reads of mean length -m from a log-normal length
distribution with standard deviation -s, but ignores reads of length less than -x. It
collects enough reads to cover the genome -c times and Introduces -e fraction errors
into each read where the ratio of insertions, deletions, and substitutions are set by
defined constants INS_RATE (default 73%) and DEL_RATE (default 20%) within generate.c.
One can control the rate at which reads are picked from the forward and reverse
strands with the -f option. The -r option seeds the random number generator for the
generation process so that one can reproducibly generate the same dataset. If this
parameter is missing, then the job id of the invocation seeds the random number
generator effectively guaranteeing a different sampling with each invocation.

The output is sent to the standard output (i.e. it is a UNIX pipe). The output is in
Pacbio .fasta format suitable as input to fasta2DB. Uppercase letters are used if the
-U option is given, and the width of each line can be controlled with the -w option.

Finally, the -M option requests that the scaffold and coordinates within said scaffold
from which each read has been sampled are written to the indicated file, one line per
read, ASCII encoded. This "map" file essential tells one where every read belongs in
an assembly and is very useful for debugging and testing purposes. If the map line for
a read is say 's b e' then if b \< e the read is a perturbed copy of s[b,e] in the
forward direction, and a perturbed copy s[e,b] in the reverse direction otherwise.

```
21. rangen <genlen:double> [-U] [-b<double(.5)>] [-w<int(80)>] [-r<int>]
```

Generate a random DNA sequence of length genlen*1Mbp that has an AT-bias of -b.
Output the sequence to the standard output in .fasta format.  Use uppercase letters if
-U is set and -w base pairs per line (default 80).  The result can then be converted
into a .dam DB and given to the simulator to create a read database over a random
synthetic sequence.  The -r option seeds the random number generator for the
generation process so that one can reproducibly generate the same sequence. If this
parameter is missing, then the job id of the invocation seeds the random number
generator effectively guaranteeing a different sequence with each invocation.

Example: A small complete example of most of the commands above. 

```
> rangen 1.0 >R.fasta           //  Generate a randome 1Mbp sequence R.fasta
> fasta2DAM R R.fasta           //  Load it into a .dam DB R.dam
> simulator R -c20. >G.fasta    //  Sample a 20x data sets of the random geneome R
> fasta2DB G G.fasta            //  Create a compressed data base of the reads, G.db
> rm G.fasta                    //  Redundant, recreate any time with "DB2fasta G"
> DBsplit -s11 G                //  Split G into 2 parts of size ~ 11MB each
> DBdust G.1                    //  Produce a "dust" track on each part
> DBdust G.2
> Catrack G dust                //  Create one track for all of the DB
> rm .G.*.dust.*                //  Clean up the sub-tracks
> DBstats -mdust G              //  Take a look at the statistics for the database

Statistics for all reads in the data set

          1,836 reads        out of           1,836  (100.0%)
     20,007,090 base pairs   out of      20,007,090  (100.0%)

         10,897 average read length
          2,192 standard deviation

  Base composition: 0.250(A) 0.250(C) 0.250(G) 0.250(T)

  Distribution of Read Lengths (Bin size = 1,000)

        Bin:      Count  % Reads  % Bases     Average
     22,000:          1      0.1      0.1       22654
     21,000:          0      0.1      0.1       22654
     20,000:          1      0.1      0.2       21355
     19,000:          0      0.1      0.2       21355
     18,000:          4      0.3      0.6       19489
     17,000:          8      0.8      1.3       18374
     16,000:         19      1.8      2.8       17231
     15,000:         43      4.1      6.2       16253
     14,000:         81      8.6     12.0       15341
     13,000:        146     16.5     21.9       14428
     12,000:        200     27.4     34.4       13664
     11,000:        315     44.6     52.4       12824
     10,000:        357     64.0     71.2       12126
      9,000:        306     80.7     85.8       11586
      8,000:        211     92.2     94.8       11208
      7,000:         95     97.3     98.4       11017
      6,000:         43     99.7     99.8       10914
      5,000:          6    100.0    100.0       10897


Statistics for dust-track

  There are 158 intervals totaling 1,820 bases (0.0% of all data)

  Distribution of dust intervals (Bin size = 1,000)

        Bin:      Count  % Intervals  % Bases     Average
          0:        158        100.0    100.0          11

> ls -al
total 66518744
drwxr-xr-x+ 177 myersg  staff        6018 Mar  2 13:28 .
drwxr-xr-x+  20 myersg  staff         680 Feb 26 19:52 ..
-rw-r--r--+   1 myersg  staff     5002464 Mar  2 13:28 .G.bps
-rw-r--r--+   1 myersg  staff       14704 Mar  2 13:28 .G.dust.anno
-rw-r--r--+   1 myersg  staff        1264 Mar  2 13:28 .G.dust.data
-rw-r--r--+   1 myersg  staff       73552 Mar  2 13:28 .G.idx
-rw-r--r--+   1 myersg  staff         162 Mar  2 13:28 G.db
> cat G.db
files =         1
       1836 G Sim
blocks =         2
size =        11 cutoff =         0 all = 0
         0         0
      1011      1011
      1836      1836
```
