corrMCM
=======

Julia script to correct the MCM for duplicate species definitions (i.e. species
with the same structure, but different MCM names).

RUN with `julia corrMCM.jl <kpp file>`, where `kpp file` is the directory and file
name of the kpp file to be corrected.  
A variable default directory is assumed, depending on your current location:

- `./`, if current folder is "mechanisms"
- `../..`, else

If you want to look for files in the current directory, when you are not in
`.\mechanisms`, use `./<kpp file>` as script argument.


Version history
===============

Version 1.0
-----------
- Correcting C4CONO3CO & CO2N3CHO
- Correcting COHM2/CHOMOH system
- Correcting NC3OO
