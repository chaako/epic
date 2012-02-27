all: cubature-20101018/cubature.c

cubature-20101018/cubature.c:
	tar -xf cubature-20101018.tgz

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v cubature-20101018.tgz)' > clean.sh
	chmod +x clean.sh
	./clean.sh

