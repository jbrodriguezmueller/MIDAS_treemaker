all:
	cp midas_treemaker.py ~/bin
	chmod 755 ~/bin/midas_treemaker.py

localdepshere:
	cp FastTree ~/bin
	chmod 755 ~/bin/FastTree
	cp FastTreeMP ~/bin
	chmod 755 ~/bin/FastTreeMP

delete_install:
	rm ~/bin/FastTree
	rm ~/bin/midas_treemaker.py
