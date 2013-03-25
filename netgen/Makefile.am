all: @abs_top_builddir@/lib/libnglib.so

@abs_top_builddir@/lib/libnglib.so:
	tar -xf netgen-5.0.0.tar.gz
	cd netgen-5.0.0 && ./configure --disable-gui --prefix=@abs_top_builddir@ \
	&& make -j8 && make install

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v netgen-5.0.0.tar.gz)' > clean.sh
	chmod +x clean.sh
	./clean.sh
