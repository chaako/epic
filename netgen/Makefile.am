all: opt/netgen/lib/libnglib.so

opt/netgen/lib/libnglib.so:
	tar -xf netgen-5.0.0.tar.gz
	cd netgen-5.0.0 && ./configure --disable-gui --prefix=$$PWD/../opt/netgen/ \
	&& make -j32 && make install

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v netgen-5.0.0.tar.gz)' > clean.sh
	chmod +x clean.sh
	./clean.sh
