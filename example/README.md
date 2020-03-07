Using cetoxa
========

To try out for yourself, once qvina02 is installed, run:

> ../cetoxa.py -f ../data/Compounds/Toxcast/pdbqt/20006.pdbqt -ncpu 2
> A Computational Ecotoxicity Assay
> https://doi.org/10.26434/chemrxiv.11944371.v1
> A summary of results is in results.csv.

This will produce an output file **results.csv** similar to the one in the github repository. Please keep in mind that the docking algorithm is not reproducible since it uses random numbers. That means that you will not obtain exactly the same results as the example.csv in the repository.

The output provided in **example.csv** is sorted according to protein family first and score then. Further subtleties like class average and standard deviation may be added in the future.
