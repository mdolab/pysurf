.. _installation:

Installation
============

Before installing pySurf, make sure you have CGNS Library and `Complexify <https://github.com/mdolab/complexify>`__ installed.

pySurf follows the standard MDO Lab build procedure.
First, find a configuration file close to your current setup in ``config/defaults`` and copy it to ``config/config.mk``.
For example:

.. prompt:: bash

    cp config/defaults/config_LINUX_GFORTRAN.mk config/config.mk

Once you have copied the config file, compile pySurf by running:

.. prompt:: bash

    make

If the real build was successful, the following lines will be printed to
the screen near the middle of the output::

   Testing if module cgnsAPI can be imported...
   Module cgnsAPI was successfully imported.

   Testing if module adtAPI can be imported...
   Module adtAPI was successfully imported.

   Testing if module curveSearchAPI can be imported...
   Module curveSearchAPI was successfully imported.

   Testing if module intersectionAPI can be imported...
   Module intersectionAPI was successfully imported.

   Testing if module utilitiesAPI can be imported...
   Module utilitiesAPI was successfully imported.

Similarly, if the complex build was successful, the following lines will be
printed near the end of the output::

    Testing if module cgnsAPI_cs can be imported...
    Module cgnsAPI_cs was successfully imported.

    Testing if module adtAPI_cs can be imported...
    Module adtAPI_cs was successfully imported.

    Testing if module curveSearchAPI_cs can be imported...
    Module curveSearchAPI_cs was successfully imported.

    Testing if module intersectionAPI_cs can be imported...
    Module intersectionAPI_cs was successfully imported.

    Testing if module utilitiesAPI_cs can be imported...
    Module utilitiesAPI_cs was successfully imported.

Finally, install the Python interface with:

.. prompt:: bash

    pip install .

Testing Your Installation
-------------------------

You can run the tests to verify your installation.

First, download the input files:

.. prompt:: bash

    ./input_files/get-input-files.sh

Then, in the root directory, run:

.. prompt:: bash

    testflo -v
