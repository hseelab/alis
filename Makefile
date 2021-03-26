prefix=/usr

all:
	cd src && make && cd ..
	gcc -std=c99 -fopenmp -O3 -o bin/alispp bin/alispp.c
	gcc -std=c99 -fopenmp -O3 -o bin/abench bin/abench.c

install:
	install -m755 -d $(prefix)/share/doc/alis/examples
	install -m644 doc/changelog $(prefix)/share/doc/alis/
	install -m644 doc/copyright $(prefix)/share/doc/alis/
	install -m644 doc/DFFit.nb $(prefix)/share/doc/alis/
	install -m644 examples/*.c $(prefix)/share/doc/alis/examples/
	install -m644 src/libalis.* $(prefix)/lib/
	install -m644 src/alis*.h $(prefix)/include/
	install -m755 bin/alis bin/alisp bin/alispp $(prefix)/bin/

uninstall:
	rm -f $(prefix)/bin/alis*
	rm -f $(prefix)/include/alis*.h
	rm -f $(prefix)/lib/libalis.*
	rm -fr $(prefix)/share/doc/alis
	
clean:
	cd src && make clean && cd ..
	rm -f bin/alispp
	rm -fr examples/*[\ -bd-~] examples/*[\ -\-0-~]c

cleanexamples:
	rm -fr examples/*[\ -bd-~] examples/*[\ -\-0-~]c
