all: SuperLU_4.3/lib/libsuperlu_4.3.a

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
	cd SuperLU_4.3 && make blaslib && make

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v superlu_4.3.tar.gz)' > clean.sh
	chmod +x clean.sh
	./clean.sh

