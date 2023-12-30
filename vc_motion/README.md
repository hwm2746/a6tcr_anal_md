**\./vc_motion/ : Files for analyzing V-C motion.**

- Several files are the same as those in ../vv_motion/, which we include here as symbolic links.

************

* ./source/ : Source scripts.

[0] Select residues for triad assignment.
    - See 'V-C BOC and PCA' in the Trajectory Analysis section of our paper.
      List selected residues under "BRESI SELECTION" of orie.inp
    - Decide on reference structure. Use the same one as for ../vv_motion/.
      Set path to this structure in 'set rcor' line of orie.inp
    
[1] Align trajectories to the reference frame using selected residues from step [0]
    - For space, save only backbone residues. 
    - orie.inp streams dcd_header.str, which defines the location of the input trajectory.
      charmm -chsize 100000 < orie.inp > out/out_orie.dat
    - Generate backbone-only psf. This only needs to be done once, so if done during ../vv_motion/, skip this step.
      charmm -chsize 100000 < gen_bb.inp > out/out_bb.dat    
    - Can be used to view aligned trajectory in VMD.

[2] Write bresi\_\*.str files listing residues from step [0]
    - bresi\_\*.str files also include V-module residues for triads. 

[3] Compile pca analysis scripts
    f95 pca_boc.f95 matvec.f95 -llapack -lblas -o te

--> Steps [4--7] are similar to those in ../vv_motion/. Difference is in using "cab" for domain name instead of "vab". Note also the trajetory aligned to the C-module from step [0] is used to generate the BOC, instead of the V-module aligned trajectory in ../vv_motion/. 

[4] Generate triads & BOC
      mkdir -p data_cab
      bash action1.sh cab

    - Value "cab" is the domain name to use.
    - action1.sh executes two files:
      action_triad.sh = for triads, streams get_triad.str
      action_cm.sh    = for BOC, streams get_cm.str
    - Both stream dcd_header_orie.str, which defines the location of the aligned trajectory from step [1].
    
    - Output: data_cab/out_temp_*, raw_*.dat, cm_*.dat.
    - Check any errors or warnings in the out*dat files. 
      grep -a 'MOST SEVERE' data_cab/out*dat      
    - If all warnings are at level 0: rm -f data_cab/out*.dat

    - Here, ./source/data_cab/ contains output files for WThigh. No files were deleted.
    - The label "1ao712" is equivalent to "WThigh"

[5] Process triad, build BOC.
      python3 proc_triad.py -d cab

    - Value "cab" is the domain name to use.      
    - Output: ./source/data_cab/boc_*.{psf,pdb,dat}, data_cab/v*.{psf,pdb}

[6] Collect measurements from triads & BOC.
      python3 anal_triad.py -d cab

    - Value "cab" is the domain name to use.
    - Output: data_*/meas_*.dat 

[7] Perform PCA on triads.
    mkdir -p data_cab/pca3
    ./te -d data_cab -i 1ao712 -o 3 -f0 25000

    - Flags: -d=directory, -i=input, -o=output, -f0=start frame for PCA. 
    - Can also set the last frame using -f1 flag. 
    - Output: Files in ./source/data_cab/pca3/.
    - The above will perform PCA of system 1ao712 (WThigh) using frames [25000,end], which is 500 ns to the end of simulation.
    - It will use the data in ./source/data_cab, and prints output to ./source/data_cab/pca3.
    - See note in ../vv_motion/README.md about -o flag. 

**********
* ./data/ : 

   - To view Fig 5A:
  cd ./fig5a
  vmd -e view_fig5a.vmd

   - Each folder contains the following data for the respective trajectory:

     - Hinge angles (Output from step [6] above): meas_*.dat 

     - PC amplitude (Output from step [7] above):
    dir_cm.dat = PC 18-dim direction vector
    std_cm.dat = PC amplitude (these valus are plot in Fig 5 Supp 1A)
    
    
**********
* ./anal/ : 

  - Fig 5C: Difference in PC amplitude between alpha and beta chains.
  python3 diff_pc_amp.py
  Also plots the PC amplitude values used to get the difference in Fig 5C. 

  - Fig 5D: Distribution of hinge angles.
  python3 hinge_angle.py

  - Fig 5E: CDR3 distance versus triad arm angle.
  python3 cdr3_v_hinge_angle.py
  - CDR3 distance files are in ../vv_motion/data/. 

**********
* ./dotp/ :  Files to plot dot products in Fig 5B. 

  - Data file (18-dim PCA vectors): \*\_pca.txt
   These are part of the printed output from BOC PCA. 
   To create these files, execute the command from step [7].
    ./te -d data_cab -i 1ao712 -o 3 -f0 25000 > temp_pca.txt
  - Manually modify temp_pca.txt so it looks like *_pca.txt. 
    This involves deleting lines that start with "Reading:", "Writing:", and first set of " pc            1" to " pc           18". The first set pertains to PCA of triads. 

  - Calculate dot product for various systems.
   List systems in setSys() function in pca_dotp.cpp. 
    g++ -std=c++11 pca_dotp.cpp -o dotp
    ./dotp > out_dotp.txt
    (pca_dotp.cpp is  thesame as in ../vv_motion/dotp, except for systems list in setSys()). 

   - Plot dot product output.
      python3 plot_dotp.py
      List of names in MODIFIABLE SECTION of plot_dotp.py need to be in the same order as setSys mols list in pca_dotp.cpp. 

