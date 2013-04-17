all: Python-2.7.3/python

Python-2.7.3/python:
	tar -xf Python-2.7.3.tgz
	cd Python-2.7.3 && ./configure --prefix=@abs_top_builddir@ && make && make install

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v Python-2.7.3.tgz)' > clean.sh
	chmod +x clean.sh
	./clean.sh

