## Introduction
Subchannel CFD (SubChCFD) is a CFD-based Subchannel analysis tool. This tool is designed to perform thermal hydraulics analyses of full-scale engineering geometry of reactor systems (e.g. a full nuclear reactor core), while maintaining reasonable computational cost. Remarkably, it has been demonstrated that SubChCFD is 2 to 3 orders of magnitude less computationally expensive than traditional CFD tools. This significant reduction in computational cost enables routine engineering design calculations and safety analyses of nuclear reactors. With years of development, SubChCFD can be applied to various types of fission reactors, including Pressurised Water Reactor (PWR), Advanced Gas-cooled Reactor (AGR), Super-critical Water Reactor (SCWR) and prismatic High Temperature Gas-cooled Reactor (HTGR).

SubChCFD is built on EDF’s open-source general purpose CFD software [code_saturne](https://github.com/code-saturne/code_saturne). To use SubChCFD, the same version of [code_saturne](https://github.com/code-saturne/code_saturne) should be installed as prerequisite. 

## Usage
1. Install the correct version of [code_saturne](https://github.com/code-saturne/code_saturne)  
2. Clone the repository:  
```
git clone git@github.com:lbyuan0831/SubchannelCFD.git
```
3. Create your model for a fuel assembly of a specific type of reactor. For starters, this can be done by adapting the examples provided  
4. Change the directory of the binary excutable of [code_saturne](https://github.com/code-saturne/code_saturne) based on your own installation in the `runcase` script and specify the number of CPUs to be used within the script
5. Run the simulation simply by
```  
cd <directory_to_your_case>/DATA
./runcase
```

## Contact
* Bo Liu, SC STFC Daresbury Laboratory, UK (bo.liu@stfc.ac.uk; bo.liu@sheffield.ac.uk)
* Shuisheng He, The University of Sheffield, UK (s.he@sheffield.ac.uk)

## References

For more details of SubChCFD, please refer to our publications:


[1] B. Liu, et al. Sub-channel CFD for nuclear fuel bundles, Nuclear Engineering and Design, 355(2019), 110318  
[2] B. Liu, et al. A coupling approach between resolved and coarse-grid sub-channel CFD, Nuclear Engineering and Design, 377(2021), 111124  
[3] B. Liu, et al. Coupled porous media appraoches in sub-channel CFD, Nuclear Engineering and Design, 377(2021), 111159  
[4] B. Liu, et al. A multiscale model of a rod bundle using subchannel CFD, Nuclear Engineering and Design, 393(2022), 111793  
