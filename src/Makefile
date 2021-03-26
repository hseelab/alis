#CC = icc -xSSE4.1 -ansi-alias -openmp
CC = gcc -fopenmp -fopt-info-vec
LINK_FLAGS = -shared -lpng
COMPILE_FLAGS = -O3 -std=c99 -I/usr/include -I/usr/include/hdf5/serial -fpic -c

all: libalis.a libalis.so

libalis.a: world.o object.o shapes.o source.o update.o output.o slices.o fields.o phaser.o farfld.o
	rm -f libalis.a
	ar rcs libalis.a world.o object.o shapes.o source.o update.o output.o slices.o fields.o phaser.o farfld.o

libalis.so: world.o object.o shapes.o source.o update.o output.o slices.o fields.o phaser.o farfld.o
	rm -f libalis.so
	$(CC) $(LINK_FLAGS) -o libalis.so world.o object.o shapes.o source.o update.o output.o slices.o fields.o phaser.o farfld.o

world.o: alis_c.h world.c
	$(CC) $(COMPILE_FLAGS) -o world.o world.c

object.o: alis_c.h object.c
	$(CC) $(COMPILE_FLAGS) -o object.o object.c

shapes.o: alis_c.h shapes.c
	$(CC) $(COMPILE_FLAGS) -o shapes.o shapes.c

source.o: alis_c.h source.c
	$(CC) $(COMPILE_FLAGS) -o source.o source.c

update.o: alis_c.h update.c
	$(CC) $(COMPILE_FLAGS) -o update.o update.c

output.o: alis_c.h output.c
	$(CC) $(COMPILE_FLAGS) -o output.o output.c

slices.o: alis_c.h slices.c
	$(CC) $(COMPILE_FLAGS) -o slices.o slices.c

fields.o: alis_c.h fields.c
	$(CC) $(COMPILE_FLAGS) -o fields.o fields.c

phaser.o: alis_c.h phaser.c
	$(CC) $(COMPILE_FLAGS) -o phaser.o phaser.c

farfld.o: alis_c.h farfld.c
	$(CC) $(COMPILE_FLAGS) -o farfld.o farfld.c

clean:
	rm -f *.a *.o *.so
