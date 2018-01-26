# InverseBarbourModel
This repository contains two codes (one in Matlab and one in Excel) that enable the back-calculation of isotopic ratios (d18O) of water sources used by plants based on measured d18O values from plant cellulose, climate variables, and tree physiology variables. Both codes are based on the Barbour model for isotopic fractionation in plants (Barbour et al, 2004), but they have been modified to invert the problem. The codes use optimization to arrive at a predicted value source water d18O. The Excel version uses the Solver AddIn to minimize the difference between the observed and predicted values of d18O in tree ring cellulose. To calculate water sources, run the Solver on the values in Column V with the specification that Column V values should be minimized based on adjusting values in Column C. The Matlab version does the optimization explicitly through looping, but also includes a Monte Carlo error estimate for each water source value based on repeated sampling from ranges of input parameters. 

Note 1: There are numerous assumptions that go into the Barbour theory, but also assumptions that are used to determine the input variables and these are documented in the codes. Model sensitivity to input parameters is discussed in the companion paper: Sargeant, C.I., C. Vallet-Coulomb, and M.B. Singer. In Review. A toolkit for detecting historical water use by forest trees, Water Resources Research.

Note 2: Back-calculated water sources should be carefully evaluated against best available information about water sources that are likely contributing to the root zone. This is also discussed in the companion paper.

This code can be cited as: 

Singer, MB, Sargeant, CI, Evans, C, A tool for back-calculating isotopic signatures of water sources used by vegetation, doi:10.5281/zenodo/1161221 <a href="https://zenodo.org/badge/latestdoi/119047089"><img src="https://zenodo.org/badge/119047089.svg" alt="DOI"></a>
