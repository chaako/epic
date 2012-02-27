all: Cuba-3.0/libcuba.a

Cuba-3.0/libcuba.a:
	tar -xf Cuba-3.0.tar.gz
	cd Cuba-3.0 && ./configure && make lib

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v Cuba-3.0.tar.gz)' > clean.sh
	chmod +x clean.sh
	./clean.sh

