all: boost_1_49_0/boost/type_traits.hpp

boost_1_49_0/boost/type_traits.hpp:
	tar -xf boost_1_49_0.tar.gz

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v boost_1_49_0.tar.gz)' > clean.sh
	chmod +x clean.sh
	./clean.sh

