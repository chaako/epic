all: VTK/install/lib/vtk-5.10/libvtkFiltering.so

VTK/install/lib/vtk-5.10/libvtkFiltering.so:
	tar -xf vtk-5.10.0.tar.gz
	cd VTK && mkdir install && \
	../../cmake/cmake-2.8.8/bin/cmake -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=./install/ \
	-DBUILD_SHARED_LIBS=ON . \
	&& make && make install
#	-DBUILD_EXAMPLES=OFF -DBUILD_SHARED_LIBS=ON \
#	-DCMAKE_BUILD_TYPE=RelWithDebInfo \
#	-DVTK_USE_CHARTS=OFF -DVTK_USE_CHEMISTRY=OFF \
#	-DVTK_USE_GEOVIS=OFF -DVTK_USE_INFOVIS=OFF \
#	-DVTK_USE_N_WAY_ARRAYS=OFF -DVTK_USE_PARALLEL=OFF \
#	-DVTK_USE_QT=OFF -DVTK_USE_RENDERING=OFF \
#	-DVTK_USE_TEXT_ANALYSIS=OFF -DVTK_USE_VIEWS=OFF \
#	-DVTK_WRAP_JAVA=OFF -DVTK_WRAP_PYTHON=OFF \
#	-DVTK_WRAP_PYTHON_SIP=OFF -DVTK_WRAP_TCL=OFF . \
#

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v vtk-5.10.0.tar.gz)' > clean.sh
	chmod +x clean.sh
	./clean.sh
