# Dendritic-spine
This Matlab code has been developed for solving the membrane shape equations for dendritic spines under the action of actin-mediated forces and induced spontaneous deviatoric curvatures (Alimohamadi, Haleh, et al. "Mechanical principles governing the shapes of dendritic spines. Frontiers in Physiology 12 (2021).). 

We used the bvp4c solver in Matlab to solve the system of ordinary differential equations with prescribed boundary conditions.

The Matlab code includes 4 main functions: 
1- Init.m. This function initializes the shape of the membrane. You can start from a flat disk where the radius (r) is proportional to the arc-length, the height (z) is equal to zero, the membrane tension (𝛌) is equal to prescribed boundary tension, the tangent angle (ᴪ), the mean curvature (H) and the defined parameter (M) are equal to zero. Init.m has four inputs; alpha is the dimensionless patch area. mesh is the discretized domain runs from 0 to 1, i.e. 0:0.01:1. lambda is membrane tension at the boundary in units of pN/nm. k0 is the bending rigidity of bare membrane in units of pN.nm, and R0 is non-dimensionalization length that can change depending on the problem. Init.m returns the initial solution for the membrane profile.
2- Loop.m. We used this function to loop over a range of pole heights to extract the force vs. membrane displacement curve. Loop.m has twelve inputs; alpha is the dimensionless patch area. mesh is the discretized domain runs from 0 to 1, i.e. 0:0.01:1. lambda is membrane tension at the boundary in units of pN/nm. acoat is the domain of the applied deviatoric curvature. rF is the domain of the applied force. zpRng is the height of the membrane at the pole. f0 represents the magnitude of applied force. k0 is the bending rigidity of bare membrane in units of pN.nm. gamma represents the slope of applied force using a hyperbolic tangent function. C0 is the spontaneous deviatoric curvature. R0 is non-dimensionalization length and initSol is the initial solution which is the output of Init.m function. Loop.m returns the height vs force, the solution of ODE, and the height of the membrane.
3- PullMembrane.m. In this function, we used bvp4c to solve the coupled ordinary differential equations with a pulling force. Here, we have an extra boundary condition (the height of the membrane at the pole). This will allow us to extract one extra unknown parameter (applied force) along with solving the ODEs equations. PullMembrane.m has twelve inputs; alpha is the dimensionless patch area. mesh is the discretized domain runs from 0 to 1, i.e. 0:0.01:1. lambda is membrane tension at the boundary in units of pN/nm. Alpha0 is the domain of the applied deviatoric curvature. rF is the domain of the applied force. Zp is the height of the membrane at the pole. f0 represents the magnitude of applied force. k0 is the bending rigidity of bare membrane in units of pN.nm. gamma represents the slope of applied force using a hyperbolic tangent function. C0 is the spontaneous deviatoric curvature. R0 is non-dimensionalization length and initSol is the initial solution which is the output of Init.m function. PullMembrane.m returns the domain discretization, the solution of ODE, and the magnitude of applied force.
4- plotMembraneProfile.m. This function plots the membrane profile according to the applied force for each prescribed membrane height. PullMembrane.m has ten inputs; Sol is the solution of ODE. t is membrane mesh points. R0 is non-dimensionalization length.coatArea is the domain of applied spontaneous deviatoric curvature.actArea is the domain of applied force, limit of x and y axes, and the title of the plot.

Instruction for running the Matlab code for each figure:
Figure2a: Fix the height of membrane and run the Matlab code to calculate the magnitude of applied axial force for different boundary tension.
Figure2b: Fix the tension of membrane at 10 pN/μm and run the Matlab code to calculate 
the magnitude of applied axial force for different membrane heights.
Figure2d: Fix the height of membrane and run the Matlab code to calculate the radius for different boundary tension.
Figure2e: Fix the height of membrane and run the Matlab code to calculate the magnitude of applied axial force for different boundary tension.
Figure3a: Fix the height of membrane and run the Matlab code to calculate the radius for different boundary tension.
Figure3b: Fix the height of membrane and run the Matlab code to calculate the magnitude of applied normal force density for different boundary tension.
Figure3c: Fix the height of membrane and membrane tension at 10 pN/μm and run the Matlab code to calculate the magnitude of applied normal force density for the different domains of applied force.
Figure4b: Fix the height of membrane and run the Matlab code to calculate the volume of spine head for different boundary tension.
Figure4C: Fix the height of membrane and run the Matlab code to calculate the magnitude of normal force density in the head and PSD domains for different boundary tension.
Figure5b: Fix the height of membrane and membrane tension and run the Matlab code to calculate the radius for different effective tension.
Figure5d: Fix the height of membrane and membrane tension and run the Matlab code to calculate the effective axial force for different effective tension.
Figure6a: Fix the height of membrane and run the Matlab code to calculate the normal force density for different effective tension.
Figure6b: Fix the height of membrane and run the Matlab code to calculate the axial force for different effective tension.
Figure6c: Fix the height of membrane and run the Matlab code to calculate the normal force density in the spine head and PSD for different effective tension.



