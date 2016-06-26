# MIDAS_treemaker

Auxiliary script to help with some tree making tasks in MIDAS.
For details on MIDAS, see https://github.com/snayfach/MIDAS/
Tested in Linux, specifically Ubuntu Xenial.

Needs:
* FastTree http://www.microbesonline.org/fasttree/
** (assumes a binary called FastTree is in the execution path)
* Muscle http://www.drive5.com/muscle/
** (assumes a binary called muscle is in the execution path)
* BioPython http://biopython.org/wiki/Download
** imports SeqIO, AlignIO, Phylo

## Installation

Assuming you have a folder called "~/bin" in your path, just type "make"
and it will copy the script to ~/bin/midas_treemaker.py and change
its permissions to executable.
