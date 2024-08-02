# MINTpy
2D Magnetotelluric inversion using linear finite element methods and a discretize-last strategy with first and second-order anisotropic regularization.  Uses esys-escript 5.10 https://doi.org/10.48610/6a7f0c5

A newer version of esys-escript has been released and the code could be adapted for this, see https://github.com/LutzGross/esys-escript.github.io/.

Commeni 4 example: uses a simpe mesh created with gmsh.  Original geo code is comm4Mesh.geo, converted to a .msh file then to a .fly file: comm4Mesh.fly.

![image](https://github.com/user-attachments/assets/1380311d-6ad9-4ef9-82c5-848240a522af)

Simple first order and second order codes simply run the forward problem first to make the data, then the inversion is run.



