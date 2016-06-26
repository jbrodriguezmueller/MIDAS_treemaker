# MIDAS_treemaker

Auxiliary script to help with some tree making tasks in MIDAS.
For details on MIDAS, see https://github.com/snayfach/MIDAS/

Tested in Linux, specifically Ubuntu Xenial. Should work in most
Linuxes. *Might* work in macOS, if it doesn't, should only need
trivial changes. I don't have an easy, reliable way to test it on
macOS, so YMMV there.

Needs:
* FastTree http://www.microbesonline.org/fasttree/
** (assumes a binary called FastTree is in the execution path)
* Muscle http://www.drive5.com/muscle/
** (assumes a binary called muscle is in the execution path)
* BioPython http://biopython.org/wiki/Download
** imports SeqIO, AlignIO, Phylo

## Installation

Assuming you have a folder called "~/bin" in your path, just type "make", and it will copy the script to ~/bin/midas_treemaker.py and change its permissions to executable.
