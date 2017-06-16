# Stochastic Reconstruction using X-Ray Tomography Projections
> In this repo, I will present the stochastic optimization algorithms for reconstructing heterogeneous materials using X-ray tomographic projections

## Table of Content
- X-ray tomography
    - Parallel-beam geometry
    - Cone-beam geometry
    - Generating X-ray projections
- 3D reconstruction
    - Stochastic optimization
    - Simulated annealing
    - Multi-modal reconstruction
    
## X-ray tomography
X-ray tomography microscopy, when properly combined with in situ experiments, is an extremely attractive, non-destructive technique for characterizing microstructure in 3D and 4D. The use of high brilliance and partially coherent synchrotron light allows one to image multi-component materials from the sub-micrometer to nanometer range. X-ray tomography can be conducted in imaging modes based on absorption or phase contrast. The technique can also be used using lab-scale systems. In x-ray tomography, 2D projections are usually obtained at small angular increments. Given a sufficiently large number of such 2D projections, different reconstruction algorithms can be allpied to reconstruction the 3D representation of the target microstructure.

<img height="350px" src="/images/x-ray-tomography.png?raw=true">

In x-ray tomographic microscopy, an x-ray beam is generated from an x-ray source and then the rays pass through the material sample that is being examined and eventually being captured on a screen, or a detector to generate the tomographic radiographs. As an x-ray passes through the material, some of photons will be scattered and/or absorbed by the materials via different light-matter interactions. Without going into details of such interactions, one can consider that the intensity of the x-ray is attenuated as it passes through the material. 
Based on the x-ray source, the x-ray tomgoraphy can be categorized as parallel beam geometry or cone beam geometry. For the parallel-beam projection geometry, a set of parallel rays are sent through the material sample from different projection angles and the attenuated intensities of the x-rays passing through the sample are recorded via a detector behind the sample.

<img height="350px" align="right" src="/images/parallel-geometry.png?raw=true">

Different from the parallel-beam projection geometry, the cone-beam geometry has a point x-ray source, and the x-rays are radically emitted from the source to form a cone pattern. A consequence of such projection geometry is that the spatial density of incident rays is higher at the front side of the material sample which is close to the source and is lower at rear side of the sample. Therefore, each voxel may affect the total attenuation values in multiple detector bins and thus, a detector-bin-based projection method is used to compute the contribution of each voxel to the total attenuation data.

<img height="350px" align="right" src="/images/cone-geometry.png?raw=true">

## Stochastic reconstruction
Once we obtain the projections, the next step is to develop the algorithm for 3D reconstruction. In the current project, the reconstruction algorithm is followed the simulated annealing procedure, which considers the reconstruction problem as an inverse optimization problem and allows one to generate virtual microstructures compatible with prescribed structural information and statistics. In particular, we formulate the reconstruction problem as an “energy” minimization problem, where the energy functional 'E' is defined as the square difference between the target data set 'D' and the corresponding data set 'D''
