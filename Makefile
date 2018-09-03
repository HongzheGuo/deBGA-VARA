CC =	gcc
#CFLAGS = -gstabs+ -Wall -O2 -Wno-unused-variable -Wno-unused-result -Wno-unused-function
#CFLAGS = -g -Wall -O2 -Wno-unused-variable -Wno-unused-result -Wno-unused-function
CFLAGS = -g -Wall -Wno-unused-variable -Wno-unused-result -Wno-unused-function
PROG =	deBGA
#SRCS =	main.c hello.c
SRCS =	$(wildcard *.c)
#OBJS =	$(SRCS:.c=.o)
OBJS =	index_build.o seed_ali_core.o load_input.o bit_operation.o LandauVishkin.o ksw.o extension_tree_alt.o LandauVishkin_tree.o main.o
LIBS =	-lz -lpthread -lgomp -lm
INCLUDES =
SUBDIRS=	.

.SUFFIXES:.c .o .cc

.c.o:
	$(CC) -c $(INCLUDES) $< -o $@
#$(CFLAGS)
all:$(PROG)

#dep:$(DEPS)

deBGA:$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $@

clean:
	@rm -f *.o $(PROG)
	
depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) -- *.c )
	
index_build.o: index_build.h bit_operation.h load_input.h LandauVishkin_tree.h
seed_ali_core.o: binarys_qsort.h LandauVishkin.h load_input.h seed_ali_p.h bit_operation.h extension_tree_alt.h LandauVishkin_tree.h
load_input.o: load_input.h
bit_operation.o: bit_operation.h
LandauVishkin.o: LandauVishkin.h
LandauVishkin_tree.o: LandauVishkin_tree.h
extension_tree_alt.o: extension_tree_alt.h LandauVishkin_tree.h
ksw.o: ksw.h
