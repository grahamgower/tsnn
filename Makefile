CPPFLAGS=-Wall -O2 -g -I kastore/c -I tskit/c
TSK_OBJECTS=kastore.o tskit_tables.o tskit_core.o tskit_trees.o \
	tskit_stats.o tskit_genotypes.o tskit_convert.o

all: kastore tskit train

train: train.o kann.o kautodiff.o libtskit.a
	${CC} ${CPPFLAGS} $^ -o $@ $(LIBS) -lm

libtskit.a: ${TSK_OBJECTS}
	${AR} rcs $@ ${TSK_OBJECTS} 

kastore.o: kastore
	${CC} -c ${CPPFLAGS} kastore/c/kastore.c -o kastore.o

tskit_%.o: tskit
	${CC} -c ${CPPFLAGS} tskit/c/tskit/$*.c -o $@

kastore:
	git clone https://github.com/tskit-dev/kastore.git
	cd kastore && git checkout C_1.1.0

tskit: 
	git clone https://github.com/tskit-dev/tskit.git
	cd tskit && git checkout C_0.99.2

clean:
	rm -f *.a *.o train

mrproper: clean
	rm -fr tskit kastore data
