# Tiltrotor Morphology Optimization Tool
Morphology optimization of a tiltrotor MAV in MATLAB

This package is based on the Masters Thesis work of Luca Rinsoz at the Autonomous Systems Lab, ETH Zurich. 

Abstract: Selecting an optimal morphology of an omni-directional flying platform with changing propeller axes is not a straightforward problem. Various factors influence the desired morphology, such as flight efficiency, omni-directionality, and control authority. In the present work the development of a tool that solves the design optimization problem for a tilt-rotor MAV is presented. Different optimal morphologies are then acquired.

## Citation
Please check back for citation reference, coming soon.

## Assumptions
* The MAV is composed of n+1 rigidly attached bodies: body core + n propellers.
* Thrust units can instantaneously achieve a desired force.
* Second order aerodynamic effects and disturbances are considered negligible.
* Airflow interactions between the rotors are not considered.

## System Requirements
* MATLAB installation, 2018a or newer
* Powerful computer will significantly improve computation time.

## Usage
Open MATLAB and set path to the main directory of tiltrotor_morphology_optimization. In the MATLAB command line, open the optimization GUI:
```
Mav_GUI
```
Select the parameters over which to optimize, and the desired optimization function.

## System Modeling
* <img src="https://latex.codecogs.com/gif.latex?n=\text{%20number%20of%20propeller%20groups}" />
* <img src="https://latex.codecogs.com/gif.latex?L=\text{%20arm%20length%20from%20body%20origin%20to%20propelle%20group%20origin}" />
* <img src="https://latex.codecogs.com/gif.latex?{\beta}_i=\text{%20angle%20of%20declination%20from%20the%20horizontal%20plane}" />
* <img src="https://latex.codecogs.com/gif.latex?{\alpha}_i=\text{%20angle%20of%20deviation%20from%20equally%20spaced%20arm%20position%20within%20the%20horizontal%20plane}" />

## Optimization
The tool uses MATLAB's fmincon optimization function to find the following optimal parameters:
* <img src="https://latex.codecogs.com/gif.latex?{\beta}_i\text{%20for%20}i\in\{1{\hdots}n\}" /> (default)
* <img src="https://latex.codecogs.com/gif.latex?{\alpha}_i\text{%20for%20}i\in\{1{\hdots}n\}" /> (optional)
* <img src="https://latex.codecogs.com/gif.latex?L" /> (optional)
* <img src="https://latex.codecogs.com/gif.latex?n" /> (optional)

## Cost Functions
Various optimization functions can be selected in the GUI dropdown. For example:

### Maximize the minimal attainable force and torque that the MAV can produce in any direction.
Omnidirectionality is defined as the capacity to accelerate in every direction. This requires a high minimal attainable force and torque. The minimal attainable forces and torques for the drone are assumed to be in the directions where one of the propellers can’t apply thrust, i.e. along the arm axis. Therefore, this optimization is computationally efficient, since it is enough to optimize the force and torque in n directions.

### Maximize minimal force/torque and minimize the system inertia.
The last term is added in order to penalize increasing arm lengths when the arm length is an optimization argument.

### Maximize the volume of the reachable force and torque space.
The force and torque spaces are two polyhedra formed by the drone’s attainable forces and torques in every direction. The idea behind this cost function is to have the largest task space for the drone and hence increase the MAV’s ability to navigate in 6 DOF. This cost function is computationally heavy for the solver, because in order to generate precise polyhedra, the forces and torques should be computed in at least 578 directions (number found empirically).

### Maximize the force, torque and the hover efficiency in all directions.
The aim of this cost function is to maximize the agility of the MAV for good disturbance rejection. Moreover, the term that maximizes the hover efficiency is designed to give the drone the ability to efficiently perform aerial interaction in any orientation. Solving the optimization problem for this cost function can also be computationally heavy depending on the number of directions computed.

### Maximize the force and the torque in one defined direction. 
This objective is computationally light and, given specific directions, the optimal design is evident. For instance, if you maximize the force in the Z direction for a 4-rotor MAV, the expected optimal solution would be a traditional quad-copter.

## Metrics
The following metrics are computed to drive the optimization function:
* Volumetric envelope of maximum forces in the body fame, under 0 torque.
* Volumetric envelope of maximum torques in the body fame, under 0 force.
* Volumetric envelope of hover efficiency in the body fame.
* System mass and inertia.

## Results
Optimization results fall into two categories: preferred efficiency, and preferred omnidirectionality. A folder named Results contains MATLAB figures of various optimized solutions.

### Preferred efficiency
Results for preferred efficiency reduce to a tilt-rotor version of a standard underactuated MAV with n arms equally distributed about the body Z axis.

### Preferred omnidirectionality 
Results for preferred omnidirectionality generate certain offsets in beta. In the case of even numbered rotor groups these find equally distributed points on a sphere, at the vertices of an n-vertex platonic solid. In the case of odd numbered rotor groups, arms deviate from the body v plane to lie evenly distributed along a conic surface, or find a relatively even distribution about a sphere.

## Feasibility Verification
As an extension to this package, the morphologies have been tested in the RotorS Gazebo simulation environment. A representation of the 6-rotor hexacopter model is shown here:

