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

In x-ray tomographic microscopy, an x-ray beam is generated from an x-ray source and then the rays pass through the material sample that is being examined and eventually being captured on a screen, or a detector to generate the tomographic radiographs. As an x-ray passes through the material, some of photons will be scattered and/or absorbed by the materials via different light-matter interactions. Without going into details of such interactions, one can consider that the intensity of the x-ray is attenuated as it passes through the material. The attenuated intensity is then given by, 
'I=I_0 e^(-μL)'
