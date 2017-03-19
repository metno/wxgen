Weather generator software
==========================

.. image:: https://travis-ci.org/metno/wxgen.svg?branch=master
  :target: https://travis-ci.org/metno/wxgen
.. image:: https://coveralls.io/repos/metno/wxgen/badge.svg?branch=master&service=github
  :target: https://coveralls.io/github/metno/wxgen?branch=master

``wxgen`` is a command-line tool for generating arbitrarily long weather time-series. The generator
produces **gridded** output for **multiple variables** (e.g. temperature, precipitation) and aims to
have realistic covariances in space, time, and across variables.

The generator uses a database of past weather model simulations (e.g 15 day forecasts) and combines the segments 
stochastically. Longer trajectories are created by concatenating the shorter trajectories from the database.
This is done by matching the end state of one trajectory with the beginning state of another. The
matching is done using a specified metric, such as the sum of the square differences between states
(with some kind of normalization strategy as each atmospheric variable has different variances).

Installation
------------

**Ubuntu**

Install the required pacakges:

.. code-block:: bash

  sudo apt-get update
  sudo apt-get install netcdf-bin libnetcdf-dev libhdf5-serial-dev
  sudo apt-get install python-setuptools python-pip
  sudo apt-get install python-numpy python-scipy python-matplotlib python-mpltoolkits.basemap

Then install ``wxgen`` by executing the following inside the extracted folder:

.. code-block:: bash

  sudo pip install .

This will create the executable ``/usr/local/bin/wxgen``. If ``/user/local/bin`` is not in your PATH
environment variable, then add it (i.e add ``export PATH=/usr/local/bin/:$PATH`` to ``~/.bashrc``).
If you do not have sudo privileges do:


.. code-block:: bash

  pip install .

This will create the executable ``~/.local/bin/wxgen``. Add the folder to your PATH environment
variable (if necessary).

If you are working on the code, the -e flag ensures that you do not need to rerun pip install every
time you make changes to the code:

.. code-block:: bash

  sudo pip install -e .

Example use
-----------

Wxgen has three commands. The first simulates sequences of weather as follows:

.. code-block:: bash

   wxgen sim -db examples/database.nc -n 10 -t 1000 -o output.nc

This command uses a NetCDF database file and creates 10 trajectories that are 100 days long. Results
are stored in output.nc. In order to evaluate the performance of the generated data, a "truth" file
can be created using the following command:

.. code-block:: bash

   wxgen truth -db examples/database.nc -o truth.nc

This command uses the day 1 forecasts from the database as the truth, and joins these together into
a sequence that spans the dates in the dataset.

Finally, the generated sequences can be evaluated using th verif command:

.. code-block:: bash

   wxgen verif -truth truth.nc output.nc -m timeseries
   wxgen verif -truth truth.nc output.nc -m variance

The -m switch selects the verification metric, in this case a timeseries.

.. image:: examples/example.gif
    :alt: Example plot
    :width: 400
    :align: center

Copyright and license
---------------------
Copyright (C) 2017 MET Norway. Wxgen is licensed under `LGPL version 3
<https://github.com/metno/wxgen/blob/master/LICENSE>`_ or (at your option) any later version.
