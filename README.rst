Hybrid weather generator
========================

.. image:: https://travis-ci.org/tnipen/wxgen.svg?branch=master
  :target: https://travis-ci.org/tnipen/wxgen
.. image:: https://coveralls.io/repos/tnipen/wxgen/badge.svg?branch=master&service=github
  :target: https://coveralls.io/github/tnipen/wxgen?branch=master

``wxgen`` is a command-line tool that generates arbitrarily long trajectories (time-series) of
weather sequences. It is a hybrid approach as it uses data from weather models and combines
sequences in a stochastic way. Weather model data is stored in a database of short weather sequences
(e.g. 15 days). Longer trajectories are created by concatenating the shorter trajectories from the
database. This is done by matching the end state of one trajectory with the beginning state of
another. The matching is done using a specified metric, such as the sum of the square differences
between states (with some kind of normalization strategy as each atmospheric variable has different
variances).

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

Example use
-----------

Wxgen has three commands. The first simulates sequences of weather as follows:

.. code-block:: bash

   wxgen sim -db database.nc -n 10 -t 100 -o output.nc

This command uses a NetCDF database file and creates 10 trajectories that are 100 days long. Results
are stored in output.nc. In order to evaluate the performance of the generated data, a "truth" file
can be created using the following command:

.. code-block:: bash

   wxgen truth -db database.nc -o truth.nc

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
