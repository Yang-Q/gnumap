#!/bin/sh -l

#
# Arguments check and suage information
#
echo "=========================================="
echo "GNUMAPs PIPELIE (v0.1) uninstallation script"
echo "=========================================="

INSTALL_DIR=$(pwd)

DEP_DIR="$INSTALL_DIR/dep"
CUTADPTD="$DEP_DIR/cutadapt"
if [ -d $CUTADPTD ] ; then
		unlink $CUTADPTD
fi
	

PRINSEQD="$DEP_DIR/prinseq-lite"
if [ -d $PRINSEQD ] ; then
		unlink $PRINSEQD
fi


USEQD="$DEP_DIR/USeq"
if [ -d $USEQD ] ; then
		unlink $USEQD
fi


SAMTOOLSD="$DEP_DIR/samtools"
cd $SAMTOOLSD
make clean

GNUMAP="$DEP_DIR/gnumap_MPI_opt"
cd $GNUMAP
make deep-clean


BLAT="$DEP_DIR/blatSrc"
cd $BLAT
make clean


PBLAT="$DEP_DIR/pblat"
cd $PBLAT
make clean
if [ -f "$PBLAT/faToTwoBit" ] ; then
	unlink "$PBLAT/faToTwoBit"
fi

cd "$INSTALL_DIR"

clear
echo -e "GNUMAPs succesfully is uninstalled\n"
echo -e "The installation directory is $INSTALL_DIR\n"
echo -e "Remove this path from your shell script environment file. (e.g, export GNUMAPS=$INSTALL_DIR for bash) \n"