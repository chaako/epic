all: fmdbSerial/lib/iMesh-Defs.inc

fmdbSerial/lib/iMesh-Defs.inc:
	./build_Serial_FMDB_GMI.sh

clean:
	echo 'rm -rf $$(ls | grep -v Makefile | grep -v build_Serial_FMDB_GMI.sh | grep -v tar | grep -v buildFMDB | grep -v SCOREC_Software); cd SCOREC_Software && rm -rf $$(ls | grep -v tar)' > clean.sh
	chmod +x clean.sh
	./clean.sh

