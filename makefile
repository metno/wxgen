# Use this file to build debian packages
# Create a file called stdb.cfg in this directory with the following
# contents, where # "precise" is your linux version:
# [DEFAULT]
# Suite: precise
#
.PHONY: dist

default: nothing

nothing:
	@ echo "This makefile does not build wxgen, use setup.py"

VERSION=$(shell grep __version__ wxgen/version.py | cut -d"=" -f2 | sed s"/ //g" | sed s"/'//g")
coverage:
	#nosetests --with-coverage --cover-erase --cover-package=wxgen --cover-html --cover-branches
	nosetests --with-coverage --cover-erase --cover-package=wxgen --cover-html

test:
	nosetests

dist: makefile
	echo $(VERSION)
	rm -rf dist
	python setup.py sdist
	python setup.py bdist_wheel
	@ echo "Next, run 'twine upload dist/*'"

clean:
	python setup.py clean
	rm -rf build/
	find . -name '*.pyc' -delete
	rm -rf dist
	rm -rf wxgen.egg-info

lint:
	python wxgen/tests/pep8_test.py

count:
	@wc -l wxgen/*.py | tail -1
