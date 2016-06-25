all:
	cp FastTree ~/bin
	chmod 755 ~/bin/FastTree
	cp FastTreeMP ~/bin
	chmod 755 ~/bin/FastTreeMP
	cp midas_treemaker.py ~/bin
	chmod 755 ~/bin/midas_treemaker.py

delete_install:
	rm ~/bin/FastTree
	rm ~/bin/midas_treemaker.py
