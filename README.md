## Sample analysis scripts for molecular dynamics simulations of A6 T-cell receptor

Authors: Ana C. Chang-Gonzalez & Wonmuk Hwang ( hwm@tamu.edu ). see LICENSE.

**Primary reference:**

Ana C. Chang-Gonzalez, Robert J. Mallis, Matthew J. Lang, Ellis L. Reinherz, and Wonmuk Hwang

**Asymmetric framework motion of TCR&#x03B1;&#x03B2; controls load-dependent peptide discrimination**

*eLife* **13**:e91881 (2024). https://doi.org/10.7554/eLife.91881

* This README.md contains general notes. See the README.md file in each folder for further details.
* See our paper (the primary reference above) for figure numbers mentioned in README.md in each folder.
* MD trajectory is available upon request. Not included here due to size. 
* Each folder contains:
  ./source/: Scripts that read & process the MD trajectory. 
  ./data/: Output from ./source/ scripts.
  ./anal/: Analysis files. Input data files necessary to execute these scripts are in ./data/ folder.

- Using the scripts on a different system
  - In a script file, parts to be modified by the user are marked with "MODIFIABLE SECTION" 
  - For reading trajectories into CHARMM scripts, dcd_header*.str and input0.str files are included. They contain the file path for the trajectory data that needs to be adjusted if you want to read your own trajectories.
- The label "1ao712" in file names indicate WThigh. The two labels are used interchangeably.

