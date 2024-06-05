# ASRE_Timoshenko
## Introduction
This repository contains the Timoshenko beam soil-structural interaction model used in the open-source GIS-based tool [REMEDY_GIS_RiskTool](https://github.com/norwegian-geotechnical-institute/REMEDY_GIS_RiskTool). The implementation is created to realize the Probabilistic Performance Based Assessment approach for tunneling- and excavation-induced surface building damage proposed by Zhao et al. (2021) and Zhao et al. (2022). The implementation is extended from the original ASRE model proposed by Franza et al. (2020). 

## Compile
The source codes are in the folder ASRE_Timo. The main file is `dllmain.cpp` and the dependencies are the `ASRE_Timo_functions.h`, `Eigen`, and `nlohmann`. The codes will compile into a dynamically linked library `ASRE_Timo.dll`. Which can be linked to another C program or be used in a Python script through the python package `ctype`. To compile and use the code on Windows computers,  open the solution `ASRE_Timo.sln` with Visual Studio and click build in the top menu bar. Compiling and using the code on other operation systems is also possible and see [this link](https://stackoverflow.com/questions/496664/c-dynamic-shared-library-on-linux) for some hints. 

## Feedback and Contribute
Have you identified a problem with the code? Run into any difficulty using it? We want to hear about it!
Submit an issue or start a discussion!

## References
1. Franza, A., Acikgoz, S. and DeJong, M.J., 2020. Timoshenko beam models for the coupled analysis of building response to tunneling. _Tunnelling and Underground Space Technology_, _96_, p.103160.
2. Zhao, J., Franza, A. and DeJong, M.J., 2021. Method for probabilistic assessment of tunneling-induced damage to surface structures considering soil-structure interaction effects. _ASCE-ASME Journal of Risk and Uncertainty in Engineering Systems, Part A: Civil Engineering_, _7_(4), p.04021055.
3. Zhao, J., Ritter, S. and DeJong, M.J., 2022. Early-stage assessment of structural damage caused by braced excavations: Uncertainty quantification and a probabilistic analysis approach. _Tunnelling and Underground Space Technology_, _125_, p.104499.