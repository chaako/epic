all: headmyshoulder-odeint-v2-47dbdac/boost/numeric/odeint.hpp

headmyshoulder-odeint-v2-47dbdac/boost/numeric/odeint.hpp:
	tar -xf headmyshoulder-odeint-v2-v2.1-128-g47dbdac.tar.gz

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v \
	headmyshoulder-odeint-v2-v2.1-128-g47dbdac.tar.gz)' > clean.sh
	chmod +x clean.sh
	./clean.sh
