#!/bin/bash -ex


export SCRIPT_HOME=$PWD

source $SCRIPT_HOME/auxillaryBuildScripts/setInstallDir.sh
set_install_dir $INSTALL_TARGET

if [ "x$CREATE_TARBALLS" == "x1" ]
then
  echo "CREATE_TARBALLS specified as: $CREATE_TARBALLS"
  echo "creating tarballs for distribution"
fi

source $SCRIPT_HOME/auxillaryBuildScripts/disableSharedLibs.sh
disable_shared_libs $DISABLE_SHARED_LIBS

if [ "$SCOREC_SOFTWARE" ]
then
  echo "SCOREC_SOFTWARE specified as: $SCOREC_SOFTWARE"
elif [ ! -e SCOREC_Software ]
then
  echo "Creating a directory for SCOREC software"
#  mkdir SCOREC_Software
  cd SCOREC_Software
  export SCOREC_SOFTWARE=$PWD
else 
  echo "directory for SCOREC software exists"
  cd SCOREC_Software
  export SCOREC_SOFTWARE=$PWD
fi

if [ "x$SCOREC_SVN" == "x1" ];
then
  echo "downloading the source trees from the SCOREC SVN to the SCOREC_Software directory"
  svn co http://redmine.scorec.rpi.edu/anonsvn/fmdb/software/trunk/FMDB/FMDB fmdb
  svn co http://redmine.scorec.rpi.edu/anonsvn/gmi/trunk gmi
  svn co http://redmine.scorec.rpi.edu/anonsvn/fmdb/software/trunk/SCORECUtil/SCORECUtil scorecutil
  svn co http://redmine.scorec.rpi.edu/anonsvn/siter/trunk siter
  svn co http://redmine.scorec.rpi.edu/anonsvn/buildutil/trunk/GNUautoTools/m4Macros m4

  echo "setting variables for FMDB, GMI, and Utilities"
  M4=$SCOREC_SOFTWARE/m4
  SITER=$SCOREC_SOFTWARE/siter
  SCORECUTIL_SRC=$SCOREC_SOFTWARE/scorecutil
  SCORECUTIL=$INSTALL_TARGET
  GMI_SRC=$SCOREC_SOFTWARE/gmi
  GMI=$INSTALL_TARGET
  FMDB_SRC=$SCOREC_SOFTWARE/fmdb
  FMDB=$INSTALL_TARGET
elif [ "x$SCOREC_SVN" == "xJENKINS" ]; then
  # Jenkins downloads source from SVN 
  echo "setting variables for FMDB, GMI, and Utilities"
  M4=$SCOREC_SOFTWARE/m4
  SITER=$SCOREC_SOFTWARE/siter
  SCORECUTIL_SRC=$SCOREC_SOFTWARE/scorecutil
  SCORECUTIL=$INSTALL_TARGET
  GMI_SRC=$SCOREC_SOFTWARE/gmi
  GMI=$INSTALL_TARGET
  FMDB_SRC=$SCOREC_SOFTWARE/fmdb
  FMDB=$INSTALL_TARGET
else
  if [ "x$HTTP_DOWNLOAD" == "x" ]; then
    HTTP_DOWNLOAD=http://www.scorec.rpi.edu/FMDB/source
  fi
  echo "downloading the tarballs to the SCOREC_Software directory from $HTTP_DOWNLOAD"
#  wget $HTTP_DOWNLOAD/FMDB-1.3.8.tar.gz
#  wget $HTTP_DOWNLOAD/GMI-1.0.1.tar.gz
#  wget $HTTP_DOWNLOAD/SCUtil.tar.gz

  echo "extracting the source from the tarballs and setting install variables"
  tar -xzf SCUtil.tar.gz
  cd SCUtil
  UTILITIES=$PWD

  cd $UTILITIES
  tar -xzf siter.tar.gz
  cd siter
  SITER=$PWD

  cd $UTILITIES
  tar -xzf SCORECUtil-0.1.tar.gz
  cd SCORECUtil-0.1/
  SCORECUTIL_SRC=$PWD
  SCORECUTIL=$INSTALL_TARGET

  cd $SCOREC_SOFTWARE
  tar -xzf GMI-1.0.1.tar.gz
  cd GMI-1.0.1
  GMI_SRC=$PWD
  GMI=$INSTALL_TARGET

  cd $SCOREC_SOFTWARE
  tar -xzf FMDB-1.3.8.tar.gz
  cd FMDB-1.3.8
  FMDB_SRC=$PWD
  FMDB=$INSTALL_TARGET

fi

source $SCRIPT_HOME/auxillaryBuildScripts/runAutoreconf.sh
source $SCRIPT_HOME/auxillaryBuildScripts/createTarballs.sh

if [ "x$ENABLE_ITAPS" == "x1" ]
then
   IMESH="--enable-imesh"
   IGEOM="--enable-igeom --enable-meshmodel"
   HAVE_FMDB="--with-fmdb=$FMDB_SRC"
   echo "IMESH specified as $IMESH and IGEOM specified as $IGEOM"
fi

echo "Building SCORECUtil"
cd $SCORECUTIL_SRC
run_autoreconf $SCOREC_SVN
./configure $SHARED_LIBS --with-gmi=$GMI_SRC --with-fmdb=$FMDB_SRC --with-iterators=$SITER --prefix=$SCORECUTIL
make -j 4
make install
create_tarballs $CREATE_TARBALLS "SCORECUtil"
make check


echo "Building GMI"
cd $GMI_SRC
run_autoreconf $SCOREC_SVN 
./configure  $SHARED_LIBS $IGEOM $HAVE_FMDB --with-scorecutil=$SCORECUTIL --with-iterators=$SITER --prefix=$GMI
make -j 4
make install
create_tarballs $CREATE_TARBALLS "GMI"
make check

echo "Building FMDB"
cd $FMDB_SRC
run_autoreconf $SCOREC_SVN
./configure $SHARED_LIBS $IMESH --with-gmi=$GMI --with-scorecutil=$SCORECUTIL --with-iterators=$SITER --prefix=$FMDB
make -j 4
make install
create_tarballs $CREATE_TARBALLS "FMDB"
make check

if [ "x$CREATE_TARBALLS" == "x1" ]; then
  echo "Creating SCUtil Tarball"
  cd $SCOREC_SOFTWARE
  if [ "$SCOREC_SVN" ]; then 
    tar --exclude=".svn" -czf siter.tar.gz siter
  else 
    tar -czf siter.tar.gz SCUtil/siter
  fi

  mkdir SCUtil_dist
  mv $SCORECUTIL_SRC/*.tar.gz SCUtil_dist
  mv $SCOREC_SOFTWARE/siter.tar.gz SCUtil_dist
  if [ "x$SCOREC_SVN" != "x1" ]; then 
    mv SCUtil SCUtil_build
  fi

  mv SCUtil_dist SCUtil
  tar -czf SCUtil.tar.gz SCUtil

  if [ "x$SCOREC_SVN" != "x1" ]; then 
    mv SCUtil_build SCUtil
  fi

fi
