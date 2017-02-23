Getting started
===============

This is where you describe how to get set up on a clean install, including the
commands necessary to get the raw data (using the `sync_data_from_s3` command,
for example), and then how to make the cleaned, final data sets.

Install
=======

Yet to come


Simulations type
================


3D simulations
--------------

To launch a simulation: src/data/make_3D_simulation.py
example: 
  smt run -m src/data/make_3D_simulation.py src/data/replication/coarsed.json p_inte=0.01 --label="yeast-5kb/traj1"

the main file coupling molecular dynamic and motion of the diffusing elements
on the polymer is src/data/replication/simulate.py

The file that take cake of the 1D sequence is src/data/replication/PMotion.py
