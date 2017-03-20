# Use this file to build debian packages
# Create a file called stdb.cfg in this directory with the following
# contents, where # "precise" is your linux version:
# [DEFAULT]
# Suite: precise
#
default: nothing

nothing:
	@ echo "This makefile does not build wxgen, use setup.py"

VERSION=$(shell grep __version__ wxgen/version.py | cut -d"=" -f2 | sed s"/ //g" | sed s"/'//g")
coverage:
	#nosetests --with-coverage --cover-erase --cover-package=wxgen --cover-html --cover-branches
	nosetests --with-coverage --cover-erase --cover-package=wxgen --cover-html

test:
	nosetests

deb_dist: makefile
	echo $(VERSION)
	rm -rf deb_dist
	python setup.py --command-packages=stdeb.command bdist_deb
	cd deb_dist/wxgen-$(VERSION)/ || exit; debuild -S -sa

clean:
	python setup.py clean
	rm -rf build/
	find . -name '*.pyc' -delete
	rm -rf deb_dist
	rm -rf wxgen.egg-info

lint:
	python wxgen/tests/pep8_test.py

count:
	@wc -l wxgen/*.py | tail -1
