uvrct
======

``uvrct`` is a python package for UV photochemical reactor modelling.

Installation
------------

These instructions will get you a copy of the project up and running on your local machine for usage, development and testing purposes. **Please note** that only Linux environments have been tested in the current implementation but the code should work independently of the OS.

Open a terminal window and clone this repository by writing:

.. code-block:: bash

    git clone https://github.com/marrov/uvrct

In order to use ``uvrct`` several `Python 3 <https://www.python.org/>`__ packages are required. Creating a brand new `Conda <https://docs.conda.io/en/latest/>`__ environment for this is recommended. This can be done easily with the provided ``environment.yml`` file as follows:

.. code-block:: bash

    conda env create -f environment.yml
    conda activate uvrct

After executing these commands a new Conda environment named ``uvrct`` will be created with all necessary packages. The environment is self-contained so as to not influence other local python installations and avoid conflicts with previously installed packages. 

To deactivate this environment simply type:

.. code-block:: bash

    conda deactivate

Authors
-------

-  `Marc Rovira <https://github.com/marrov>`__

