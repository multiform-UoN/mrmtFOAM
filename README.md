# mrmtFOAM
Multi-Rate Mass Transfer (MRMT) model for mobile-immobile transport in porous media using the opensource libary OpenFOAM.

Installation
------------

With OpenFOAM-7 installed simply run the _Allwmake_ script. Use the _Allwclean_ script to clean (uninstall) the library.

# Library structure

Applications
------------

This folder contains applications created using the _multiContinuumModels_ library and OpenFOAM.

Currently, two applications are available:

* __simpleDarcyFoam__: Creates a velocity field solving the Darcy equation.
* __multiRateScalarTransportFoam__: Solves the scalar transport equation using the multi-rate model.

Etc
------------

This folder contains some useful functions.


Src
------------

This folder contains the _multiContinuumModels_ library.




Tutorials
--------

This folder contains the cases analysed in the paper (from now on referred as [1]):

"Heterogeneous Multi-Rate mass transfer models in OPENFOAM"
by Federico Municchi,Nicodemo di Pasquale,Marco Dentz,Matteo Icardi


as tutorials ready to be run.

Each directory contains all the files needed to generate the mesh, run the simulations and the postprocessing (such as calculations of breakthrough curves and averages). All the files should allow to rerun the same calculations shown in [1].

There are four main directories:

0D. Simpleremediation: cases analogue to the one presented in [2]. It employes the same beta and omega reported in Tab. 2 in [2] 

1D. spherical_immobile_region: single immobile region (Sphere)

2D. Two dimensional cases as discussed in [1]. 

    Generation of the velocity field:
        First, the flow field must be calculated by using < simpleDarcyFoam > in the simpleDarcyFoam directory. The permeability can be chosen constant or  defined as a random field. An example of a random permeability field (0/K) the same used to generate the results in [1] is included.
 
    MultiRate models:
        7Sp_Krand: 7 Spheres model taken from [2] and indicated as 7Sp in [1]
        Comp_Krand: Composite model taken from [2] and indicated as Comp in [1]
        R1Sp_KrandOmrandBetaRand: single sphere model, where omega and beta are defined as random fields derived from K. The relevant functionobjects needed to generate beta and gamma are included in the /constant/ directory

3D. Three dimensional case as discussed in [1] (production of the results is similar to 2D case).

    The flow field can be generated using the command < simpleDarcyFoam > in the simpleDarcyFoam directory. 
    7Sp_Krand: 7 Spheres models as discussed in [1]. Coeffients for omega and beta are taken from Tab. 2 in [1]






REFERENCES:

[1] Federico Municchi, Nicodemo di Pasqualea, Marco Dentzband Matteo Icardi, 2020. Heterogeneous Multi-Rate mass transfer models in OPENFOAM. https://arxiv.org/abs/2006.02704

[2] Municchi F, Icardi M. Generalized multirate models for conjugate transfer in heterogeneous materials. Phys Rev Res. 2020;2(1):https://link.aps.org/doi/10.1103/PhysRevResearch.2.013041

Aknowledgements
---------------

This work has been funded by the European Union's Horizon 2020
research and innovation programme, grant agreement number 764531, "SECURe -- Subsurface Evaluation of Carbon capture and storage and Unconventional risks".
