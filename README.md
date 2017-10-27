AP-Cloud method code for cosmology simulation
Author:Xingyu Wang <xingyuwangcs@gmail.com>
For the detailed description of the method, see <https://arxiv.org/pdf/1603.04039.pdf>

cosmology_main.cpp: Main function

cosmology_controller.cpp/h: A wrap to call cosmology_solver, cosmology_mover and cosmology_viewer

cosmology_solver.cpp/h: Density estimation, Poisson solver for potential, differentiation for gravity

cosmology_mover.cpp/h: Selection of computational particles, interpolation from computational to physical particles, particle mover

physical_particle.cpp/h: Physical particles should be storaged in an array since their size is fixed. Data structure: Position, velocity, acceleration, mass, type

computational_particle.cpp/h: Computational particles should be storaged in vector since their size is unknown in advance and there are frequent random accessing. Data structure: Position, size, number of physical particles, density, potential, gravity field

initialization.cpp/h: Read initial condition

read_parameter.cpp/h: Read parameter file

cosmology_viewer.cpp/h: Output physical/computational particles into .vtk files

octree.cpp/h: Build octree for physical particles, select computational particles


Governing Equations

dp/da=-F(a)\Nabla \phi, dx/da=F(a)p/a^2

\Laplacian \phi=3\Omega_0/2/a(\rho-1)

F(a)=((\Omega_0+\Omega_{curv}a+\Omega_{\Lambda,0}a^3)/a)^{-1/2}
