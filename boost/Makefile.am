all: boost_1_49_0/boost/type_traits.hpp

boost_1_49_0/boost/type_traits.hpp:
	tar -xf boost_1_49_0.tar.gz
	cd boost_1_49_0 && ./bootstrap.sh --prefix=@abs_top_builddir@ --with-libraries=python,mpi,program_options \
	--with-python=@abs_top_builddir@/bin/python && echo "using mpi ;" >> project-config.jam && ./b2 install

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v boost_1_49_0.tar.gz)' > clean.sh
	chmod +x clean.sh
	./clean.sh

