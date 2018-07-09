#! /bin/bash -l
blat "$1" -minIdentity=93 -out=pslx -noHead -q=rna -repMatch=256 "$2" "$3"
