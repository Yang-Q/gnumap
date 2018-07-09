#!/bin/sh -l

#
# Arguments check and suage information
#

INSTALL_DIR=$(pwd)
DEP_DIR="$INSTALL_DIR/dep"

QCONTROL="$1"
[ "$QCONTROL" ] || QCONTROL="null"

BLAT="$2"
[ "$BLAT" ] || BLAT="null"

MPI="$3"
[ "$MPI" ] || MPI="null"

MACHTYPE="$4"
[ "$MACHTYPE" ] || MACHTYPE="null"

echo "=========================================="
echo "GNUMAPs PIPELIE (v0.1) installation script"
echo "=========================================="

echo "$INSTALL_DIR"
#
# Sanity checks for the specified galaxy directory
#
if [ "$QCONTROL" = "qc" ] ; then

	PRINSEQD="$DEP_DIR/prinseq-lite-0.20.4"
	echo -e "Checking prinseq directory...\n"
	[ -d "$PRINSEQD" ] || { echo "Error: directory '$PRINSEQD' does not exist!" ; exit 1 ; }
	
	CUTADPT="$DEP_DIR/cutadapt-1.2.1"
	echo -e "Checking cutadpt directory...\n"
	[ -d "$CUTADPT" ] || { echo "Error: directory '$CUTADPT' does not exist!" ; exit 1 ; }
	
# 	FASTQC="$INSTALL_DIR/FastQC"
# 	echo -e "Checking fastqc directory...\n"
# 	[ -d "$FASTQC" ] || { echo "Error: directory '$FASTQC' does not exist!" ; exit 1 ; }	
fi

SAMTOOLS="$DEP_DIR/samtools"
echo -e "Checking samtools directory...\n"
[ -d "$SAMTOOLS" ] || { echo "Error: directory '$SAMTOOLS' does not exist!" ; exit 1 ; }	

GNUMAP="$DEP_DIR/gnumap_MPI_opt"
echo -e "Checking gnumap_MPI_opt directory...\n"
[ -d "$GNUMAP" ] || { echo "Error: directory '$GNUMAP' does not exist!" ; exit 1 ; }
echo "ok"

if [ "$BLAT" = "blat" ] ; then
	PBLAT="$DEP_DIR/pblat"
	echo -e "Checking pblat directory...\n"
	[ -d "$PBLAT" ] || { echo "Error: directory '$PBLAT' does not exist!" ; exit 1 ; }
	echo "ok"

	BLAT1="$DEP_DIR/blatSrc"
	echo -e "Checking blat directory...\n"
	[ -d "$BLAT1" ] || { echo "Error: directory '$BLAT1' does not exist!" ; exit 1 ; }
	echo "ok"
fi

USEQ="$DEP_DIR/USeq_8.6.4"
echo -e "Checking USeq directory...\n"
[ -d "$USEQ" ] || { echo "Error: directory '$USEQ' does not exist!" ; exit 1 ; }
echo "ok"

SCRIPTS="$INSTALL_DIR/scripts"
echo -e "Checking gnumap_wrapper script directory...\n"
[ -d "$SCRIPTS" ] || { echo "Error: directory '$SCRIPTS' does not exist!" ; exit 1 ; }
echo "ok"

#
# Compile each package
#
# 
if [ "$QCONTROL" = "qc" ] ; then
	echo -e "Installing FASTQ quality control programs...\n"
	
	cd "$CUTADPT"
	export PYTHONPATH=$CUTADPT/lib/python2.7/site-packages
	python setup.py install --prefix=$CUTADPT

	echo -e "Installing cutadapt...\n"
	if [ -f "$DEP_DIR/cutadapt" ] ; then
		unlink "$DEP_DIR/cutadapt"
	fi
	ln -s "$CUTADPT" "$DEP_DIR/cutadapt"
	cd "$INSTALL_DIR"

	#install prinseq
	echo -e "Installing prinseq...\n"

	if [ -d "$DEP_DIR/prinseq-lite" ] ; then
		unlink "$DEP_DIR/prinseq-lite"
	fi
	ln -s "$PRINSEQD" "$DEP_DIR/prinseq-lite"
	cd "$INSTALL_DIR"
	echo -e "ok"
fi

#TODO: compile faToBit from Kent's original source

echo -e "Installing samtools...\n"
cd "$DEP_DIR/samtools"
make
cd "$INSTALL_DIR"
echo -e "ok"

if [ -d "$DEP_DIR/USeq" ] ; then
	unlink "$DEP_DIR/USeq"
fi
ln -s "$USEQ" "$DEP_DIR/USeq"
echo -e "ok"

echo -e "Installing gnumap_MPI_opt...\n"
cd "$GNUMAP"
if [ -f "$GNUMAP/Makefile" ] ; then
		unlink "$GNUMAP/Makefile"
fi

if [ "$MPI" = "mpi" ] ; then
	ln -s ${GNUMAP}/Makefile.mpi ${GNUMAP}/Makefile
else
	ln -s ${GNUMAP}/Makefile.pthread ${GNUMAP}/Makefile
fi

make all
echo -e "ok"


cd "$INSTALL_DIR"

if [ "$BLAT" = "blat" ] ; then
	echo -e "Installing BLAT (MACHTYPE=${MACHTYPE})...\n"
	if [ ! -d ${HOME}/bin/${MACHTYPE} ] ; then
		echo -e "make directory [${HOME}/bin/${MACHTYPE}] ...\n"
		mkdir -p ${HOME}/bin/${MACHTYPE}
	fi
	
	if [ ! -d ${BLAT1}/lib/${MACHTYPE} ] ; then
		echo -e "make directory [${BLAT1}/lib/${MACHTYPE}] ...\n"
		mkdir -p ${BLAT1}/lib/${MACHTYPE}
	fi
	
	cd $BLAT1
	echo -e "make MACHTYPE=${MACHTYPE}"
	make "MACHTYPE=${MACHTYPE}"
	echo -e "ok"

	echo -e "Installing PBLAT...\n"
	cd $PBLAT
	make
	
	if [ -f ${PBLAT}/faToTwoBit ] ; then
		unlink ${PBLAT}/faToTwoBit
	fi
	ln -s ${HOME}/bin/${MACHTYPE}/faToTwoBit ${PBLAT}/faToTwoBit
	cd $INSTALL_DIR
	echo -e "ok"
fi

clear

echo -e "GNUMAPS succesfully is installed\n"
echo -e "The installation directory is $INSTALL_DIR\n"
echo -e "Add this path to your shell script environment file. (e.g, export GNUMAPS=$INSTALL_DIR for bash) \n"

ls $INSTALL_DIR/dep/gnumap_MPI_opt/bin/gnumap-stl
ls $INSTALL_DIR/dep/pblat/pblat
ls $INSTALL_DIR/dep/pblat/faToTwoBit
