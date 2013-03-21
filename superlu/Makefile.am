if USE_MPI
all: SuperLU_4.3/lib/libsuperlu_4.3.a SuperLU_DIST_3.2/lib/libsuperlu_dist_3.2.a
else
all: SuperLU_4.3/lib/libsuperlu_4.3.a
endif

SuperLUrootVAL = $(PWD)/SuperLU_4.3

SuperLU_4.3/lib/libsuperlu_4.3.a:
	tar -xf superlu_4.3.tar.gz
	echo "PLAT = _linux" > make.inc
	echo "SuperLUroot   	= $(SuperLUrootVAL)" >> make.inc
	echo 'SUPERLULIB   	= $$(SuperLUroot)/lib/libsuperlu_4.3.a' >> make.inc
	echo 'BLASLIB   	= $$(SuperLUroot)/lib/libblas.a' >> make.inc
	echo "TMGLIB       	= libtmglib.a" >> make.inc
	echo 'LIBS		= $$(SUPERLULIB) $$(BLASLIB)' >> make.inc
	echo "ARCH         = ar" >> make.inc
	echo "ARCHFLAGS    = cr" >> make.inc
	echo "RANLIB       = ranlib" >> make.inc
	echo "CC           = $(CC)" >> make.inc
	echo "CFLAGS       = -O3" >> make.inc
	echo "NOOPTS       = " >> make.inc
	echo "FORTRAN	     = $(F77)" >> make.inc
	echo "FFLAGS       = -O2" >> make.inc
	echo "LOADER       = $(CC)" >> make.inc
	echo "LOADOPTS     =" >> make.inc
	echo "CDEFS        = -DAdd_" >> make.inc
	cp make.inc SuperLU_4.3/make.inc
	cd SuperLU_4.3 && make -j 8 blaslib && make install && make -j 8 lib && make testing

DSuperLUrootVAL = $(PWD)/SuperLU_DIST_3.2

SuperLU_DIST_3.2/lib/libsuperlu_dist_3.2.a:
	tar -xf superlu_dist_3.2.tar.gz
	echo "PLAT = _linux" > make.inc
	echo "DSuperLUroot   	= $(DSuperLUrootVAL)" >> make.inc
	echo 'DSUPERLULIB   	= $$(DSuperLUroot)/lib/libsuperlu_dist_3.2.a' >> make.inc
	echo 'BLASLIB   	= $$(DSuperLUroot)/lib/libblas.a' >> make.inc
#	echo 'BLASLIB   	= ../lib/libblas$(PLAT).a' >> make.inc
	echo "TMGLIB       	= libtmglib.a" >> make.inc
	echo 'LIBS		= $$(DSUPERLULIB) $$(BLASLIB)' >> make.inc
	echo "ARCH         = ar" >> make.inc
	echo "ARCHFLAGS    = cr" >> make.inc
	echo "RANLIB       = ranlib" >> make.inc
	echo "CC           = mpicc" >> make.inc
	echo "CFLAGS       = -O3" >> make.inc
#	echo "CFLAGS       = -O3 -DDEBUGlevel=0 -DPRNTlevel=1 -DPROFlevel=0 -D_LONGINT" >> make.inc
	echo "NOOPTS       = -O0" >> make.inc
	echo "FORTRAN	     = mpif90" >> make.inc
	echo "FFLAGS       = -O2" >> make.inc
	echo "LOADER       = mpicc" >> make.inc
	echo "LOADOPTS     =" >> make.inc
	echo "CDEFS        = -DAdd_" >> make.inc
	cp make.inc SuperLU_DIST_3.2/make.inc
	patch SuperLU_DIST_3.2/SRC/get_perm_c.c get_perm_c.patch
	patch SuperLU_DIST_3.2/SRC/get_perm_c_parmetis.c get_perm_c_parmetis.patch
	cd SuperLU_DIST_3.2 && make -j 8 blaslib && make -j 8 lib && make -j 8 install

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v superlu_4.3.tar.gz | grep -v superlu_dist_3.2.tar.gz | grep -v .patch)' > clean.sh
	chmod +x clean.sh
	./clean.sh

