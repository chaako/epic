all: Cuba-3.0/libcuba.a

Cuba-3.0/libcuba.a:
	tar -xf Cuba-3.0.tar.gz
	cd Cuba-3.0 && patch < ../configure.patch && autoreconf && \
	./configure CFLAGS="-DNDIM=3 -DNCOMP=1" && make lib
# The patch disables fork detection and thus threading in cuba...
# setting CUBACORES didn't do the trick
#	cd Cuba-3.0 && ./configure CUBACORES=1 CFLAGS="-DNDIM=3 -DNCOMP=1" && make lib

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v Cuba-3.0.tar.gz | \
	grep -v configure.patch)' > clean.sh
	chmod +x clean.sh
	./clean.sh

