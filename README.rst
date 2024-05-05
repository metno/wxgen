Weather generator software
==========================

.. image:: https://travis-ci.org/metno/wxgen.svg?branch=master
  :target: https://travis-ci.org/metno/wxgen
.. image:: https://coveralls.io/repos/metno/wxgen/badge.svg?branch=master&service=github
  :target: https://coveralls.io/github/metno/wxgen?branch=master

``wxgen`` is a command-line tool for generating **arbitrarily long** weather time-series. The generator
produces **gridded** output for **multiple variables** (e.g. temperature, precipitation) and aims to
have realistic covariances in space, time, and across variables.

The generator uses a database of past weather model simulations (e.g 15 day forecasts) and combines
the segments stochastically. Longer trajectories are created by concatenating the shorter
trajectories from the database.  This is done by matching the end state of one trajectory with the
beginning state of another. The matching is done using a specified metric, such as the sum of the
square differences between states (with some kind of normalization strategy as each atmospheric
variable has different variances).

Documentation
-------------

For information about how to use wxgen, check out the wiki page at https://github.com/metno/wxgen/wiki

Features
--------

.. image:: examples/example.gif
    :alt: Example plot
    :width: 400
    :align: center

Installing on Ubuntu
----------------------

**Prerequisites**

Install the required pacakges:

.. code-block:: bash

  sudo apt-get update
  sudo apt-get install netcdf-bin libnetcdf-dev libhdf5-serial-dev
  sudo apt-get install python-setuptools python-pip
  sudo apt-get install python-numpy python-scipy python-matplotlib python-mpltoolkits.basemap

**Installing using pip**

After this, the easiest is to install the latest version of ``wxgen`` using pip:

.. code-block:: bash

  sudo pip install wxgen


This will create the executable ``/usr/local/bin/wxgen``. If ``/usr/local/bin`` is not in your PATH
environment variable, then add it (i.e add ``export PATH=/usr/local/bin/:$PATH`` to ``~/.bashrc``).
If you do not have sudo privileges do:

.. code-block:: bash

  pip install wxgen --user

This will create the executable ``~/.local/bin/wxgen``. Add this to your PATH environment
variable if necessary (i.e add ``export PATH=$PATH:~/.local/bin`` to ``~/.bashrc``).

If at a later time you wish to update to a newer release, do:

.. code-block:: bash

   pip install wxgen --upgrade

(adding ``--user`` if appropriate)

**Installing from source**

Alternatively, to install from source, download the source code of the latest version:
https://github.com/metno/wxgen/releases/. Unzip the file and navigate into the extracted folder.

Then install ``wxgen`` by executing the following inside the extracted folder:

.. code-block:: bash

  sudo pip install -r requirements.txt
  sudo python setup.py install

If you are working on the code, the -e flag ensures that you do not need to rerun pip install every
time you make changes to the code:

.. code-block:: bash

  sudo pip install -e .

Again, if you do not have sudo privileges, remove ``sudo`` and append ``--user`` to the previous pip
and python commands.

Copyright and license
---------------------
Copyright (C) 2017 MET Norway. Wxgen is licensed under `LGPL version 3
<https://github.com/metno/wxgen/blob/master/LICENSE>`_ or (at your option) any later version.
