noinst_LIBRARIES = libcubature.a
libcubature_a_SOURCES = cubature-20101018/cubature.c

cubature-20101018/cubature.c:
	tar -xf cubature-20101018.tgz
	cp cubature-20101018/cubature.h @abs_top_builddir@/include/

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v cubature-20101018.tgz)' > clean.sh
	chmod +x clean.sh
	./clean.sh

