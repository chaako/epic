all: boostable/boost/numeric/quadrature/adaptive.hpp

boostable/boost/numeric/quadrature/adaptive.hpp:
	tar -xf quadrature.tgz
	cd boostable/boost/numeric/quadrature && \
	find -name "*.hpp" -print | xargs sed -i 's/F integrand/F\& integrand/g'
	cp -r boostable/boost/numeric/quadrature @abs_top_builddir@/include/boost/numeric/

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v quadrature.tgz)' > clean.sh
	chmod +x clean.sh
	./clean.sh
