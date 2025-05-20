# Collinear-Lagrange-Points

A Fortran collinear Lagrange point calculator.

- forclap.f90 &nbsp;&nbsp; Modern Fortran source code.<br />

ForCLaP (Fortran Collinear Lagrange Points Calculator) is a Fortran command line application demonstrating Newton-Raphson convergence on collinear solutions to the circular restricted three-body problem. Optional displacement of the collinear Lagrange points due to radiation pressure, relevant to small bodies such as asteroidal dust particles, is incorporated.

ForCLaP assumes a traditional three-body system with a large, central primary (m<sub>1</sub>), a smaller secondary (m<sub>2</sub>) and an infinitesimal tertiary mass (m<sub>3</sub>). ForCLaP accepts a user-input value for the m<sub>1</sub>/m<sub>2</sub> mass ratio of the system. Use of a mass ratio within the Routh value limit, is enforced however.

ForCLaP also accepts a user-input value for a ratio of solar radiation pressure to gravitational force experienced by m<sub>3</sub>, &beta;. This allows simulation of a luminous primary such as the Sun, but ForCLaP can just as easily be used for simulation of a non-luminous primary, such as planet, with an input ratio, &beta;, of zero.

Collinear Lagrange point solution calculation is done with a Newton-Raphson algorithm. The program therefore inherits limitations of that algorithm. The algorithm can suffer runaway solution divergence, for example, but ForCLaP explicitly flags this and provides graceful algorithm termination. For algorithm convergence, iteration continues until the solution and the Newton derivative are both stable to four significant figures.

Program source code requires compilation on the target system with a Fortran compiler. E.g., via:

    gfortran forclap.f90 -o forclap

ForCLaP has been tested with the GNU Fortran (gfortran) compiler, version 15.1.0.

Reference: Stenborg, TN 2008, "[Collinear Lagrange Point Solutions in the Circular Restricted Three-Body Problem with Radiation Pressure using Fortran](https://aspbooks.org/custom/publications/paper/394-0734.html)", in RW Argyle, PS Bunclark & JR Lewis (eds), Astronomical Data Analysis Software and Systems XVII, Astronomical Society of the Pacific Conference Series, vol 394, pp. 734-737.
