# MOredox23

MATLAB code and data files for manuscript "Impact-driven redox stratification of Earth's mantle"

Estefania Larsen Henningsen, Yale University -- estefania.larsen@yale.edu -- Summer 2022-Fall 2024

This code produces all data and figures from the manuscript and is a standalone repository. This code was originally written with MATLAB 2021a, but also works in v2023b. MATLAB's Optimization Toolbox may be necessary for some scripts to run.

This code is organized such that any .m script files in the main will work alone, but may access either the /functions folder or the /db folder.
Scripts without ...fig ending will generate the data necessary, with a line to write to file. There is sometimes a figure output in these files for testing purposes.
To write data, I used .xlsx by default. Scripts with a ...fig ending read the necessary data file and generate the figure corresponding with it. 
For instance, EffvsD.m calculates the magma ocean Fe3+/sumFe after a giant impact with given parameters (initial ratio, mass of GI) as a function of degree of
equilibration ('efficiency') and depth of melting. This script writes to a .xlsx file the matrix of values for both pre-mix and post-mixing with underlying mantle.
The script EffvsD_fig.m reads this data and produces the corresponding figure in the manuscript.

Original geotherms were calculated using scripts provided by Jun Korenaga. Please ask if these scripts are desired.

The dVvsP.xlsx data is the difference in molar volumes as a function of pressure, and I used Jie Deng's original Python code to output this. His github repository can
be found here: https://github.com/neojie/oxidation_lite
From this, I calculated the PV integration for oxygen fugacity in /testing scripts/PVintegration.m. This produces the output file PVcalc.xlsx.

This code was created by me for my needs on the project. I have attempted, to the best of my ability, to clean up the code and comment for my own sake. I hope that 
someone with sufficient MATLAB experience (or general coding experience) will be able to navigate the code. However, I cannot obviously guarantee that. Any questions
regarding the code should be directed to me (E. Larsen).

"We stand on the shoulders of giants."
Special thanks to Jie Deng and Marc Hirschmann, whose work before me was instrumental to this project.
Also, thank you to a vast group of other scientists involved in magma ocean geodynamics and geochemistry... including, but not limited to, D. C. Rubie, D. J. Frost, 
S. Aulbach, V. Stagno, R. Degueun, M. Landeau, P. Olson, C. McCammon, H. StC. O'Neill, V. Solomatov, and D. J. Stevenson.
