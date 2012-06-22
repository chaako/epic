all: cmake-2.8.8/bin/cmake

cmake-2.8.8/bin/cmake:
	tar -xf cmake-2.8.8.tar.gz
	cd cmake-2.8.8 && ./bootstrap && make

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v cmake-2.8.8.tar.gz)' > clean.sh
	chmod +x clean.sh
	./clean.sh

 