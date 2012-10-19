all: eigen-eigen-43d9075b23ef/Eigen

eigen-eigen-43d9075b23ef/Eigen:
	tar -xf eigen-eigen-43d9075b23ef.tar.gz

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v eigen-eigen-43d9075b23ef.tar.gz | grep -v .py)' > clean.sh
	chmod +x clean.sh
	./clean.sh

