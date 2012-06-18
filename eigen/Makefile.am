all: eigen-eigen-13a11181fc5a/Eigen

eigen-eigen-13a11181fc5a/Eigen:
	tar -xf eigen-eigen-13a11181fc5a.tar.gz

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v eigen-eigen-13a11181fc5a.tar.gz | grep -v .py)' > clean.sh
	chmod +x clean.sh
	./clean.sh

