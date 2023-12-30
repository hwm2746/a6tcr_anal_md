**\./vv\_motion/ : Files for analyzing V&#x03B1;-V&#x03B2; motion.**

**********

* ./source/ : Source scripts.

[0] Select residues for triad assignment.
    - See 'Variable domain triads and PCA' in the Trajectory Analysis section of our paper.
      List selected residues under "BRESI SELECTION" of orie.inp
    - Decide on the reference structure.
      Set path to this structure in 'set rcor' line of orie.inp
    
[1] Align trajectories to the reference frame using selected residues from step [0]
    - For space, save only backbone residues. 
    - orie.inp streams dcd_header.str, which defines the location of the input trajectory.
      charmm -chsize 100000 < orie.inp > out/out_orie.dat
    - Generate backbone-only psf:
      charmm -chsize 100000 < gen_bb.inp > out/out_bb.dat    
    - Can be used to view aligned trajectory in VMD.

[2] Write bresi\_\*.str files listing residues from step [0]
    - bresi\_\*.str files also include hinge and C-module residues for BOC. 

[3] Compile pca analysis scripts
    f95 pca_boc.f95 matvec.f95 -llapack -lblas -o te

--> Steps [4--7] are similar to those in ../vc_motion/.

[4] Generate triads & BOC
      mkdir -p data_vab
      bash action1.sh vab

    - Value "vab" is the domain name to use.
    - action1.sh executes two files:
      action_triad.sh = for triads, streams get_triad.str
      action_cm.sh    = for BOC, streams get_cm.str
    - Both stream dcd_header_orie.str, which defines the location of the aligned trajectory from step [1].
    
    - Output: data_vab/out_temp_*, raw_*.dat, cm_*.dat.
    - Check any errors or warnings in the out*dat files. 
      grep -a 'MOST SEVERE' data_vab/out*dat      
    - If all warnings are at level 0: rm -f data_vab/out*.dat

    - Here, ./data_vab/ contains output files for WThigh. No files were deleted.
    - The label "1ao712" is equivalent to "WThigh". This is just how "WThigh" was named in my personal files. 

[5] Process triad, build BOC.
      python3 proc_triad.py -d vab
      
    - Value "vab" is the domain name to use.      
    - Output: data_vab/boc_*.{psf,pdb,dat}, data_vab/v*.{psf,pdb}

[6] Collect measurements from triads & BOC.
      python3 anal_triad.py -d vab

    - Value "vab" is the domain name to use.
    - Output: data_*/meas_*.dat 

[7] Perform PCA on triads.
    mkdir -p data_vab/pca3
    ./te -d data_vab -i 1ao712 -o 3 -f0 25000

    - Flags: -d=directory, -i=input, -o=output, -f0=start frame for PCA. 
    - Can also set last frame using -f1 flag. 
    - Output: Files in ./source/data_vab/pca3/.
    - The above will perform PCA of system 1ao712 (WThigh) using frames [25000,end], which is 500 ns to end of simulation.
    - It will use the data in data_vab trajectory, and print output to ./source/data_vab/pca3/.
    - The "-o 3" is arbitrary. It can be used to save PCA of other frame intervals. For example, 
      ./te -d data_vab -i 1ao712 -o 1 -f0 10000 -f1 25000
      This will perform PCA using frame range [10000,25000], and save the output to ./source/data_vab/pca1/. 

[8] Calculate angles between triad arms.
    - Compile: g++ -std=c++11 calc_triad_angle.cpp
    - Execute: ./a.out

    - Output: phi*.txt. We placed these in ./data/*/ folders.
    - No user modifications necessary, as long as output files from step [5] are in ./source/data_vab/. 
    - Search for "if0" in calc_triad_angle.cpp to see input file definition.
    
[9] Calculate CDR3a-CDR3b distance.
    charmm -chsize 200000 < cdr3.inp > out_cdr3.dat
    grep -a 'DR <- ' out_cdr3.dat | sed -e 's/"//g' | gawk '{print $4}' > cdr3.dat

    - Output: out_cdr3.dat and cdr3.dat. We placed cdr3.dat in ./data/*/.
    - CDR3 distanced is calculated from trajectory file prior to orientation in step [1].
    - Script streams dcd_header.str.  
    - User will need to adjust cdr3 residues to use for distance measure. 

**********
* ./data/: 

  - To view Fig 4B:
    cd fig4b/
    vmd -e view_fig4b.vmd

  - Each folder contains the following data for the respective trajectory:

  - Triad arm angles (from step [8] above):
    phi1.txt = e1-e1 angle
    phi2.txt = e2-e2 angle
    phi3.txt = e3-e3 angle
    
  - CDR3a-CDR3b distance (from step [9] above): cdr3.dat

  - Projection of Va-Vb triads at each frame onto first 3 PC directions (from step [7] above):
    traj01.dat = PC mode 1
    traj02.dat = PC mode 2
    traj03.dat = PC mode 3

**********
* ./anal/ : 

  - Fig 4C: Distribution of e1-e1, e2-e2, and e3-e3 angles.
    python3 va-vb_angle.py
  - Also plots Fig 4 Supp 2D--F in one plot. 

  - Fig 4D: CDR3 distance versus triad arm angle.
    python3 cdr3_v_angle.py

  - Fig 4 Supp 1D: PC versus CDR3 distance.
    python3 proj_v_cdr3.py

**********
* ./dotp/ : Files to plot dot products in Fig 4 Supp 1C.

  - Data files:  *_pca.txt
    These are part of the printed output from the triad PCA. 
    To create them, follow step [7]: ./te -d data_vab -i 1ao712 -o 3 -f0 25000 > temp_pca.txt
  - Output: temp_pca.txt, which prints the 18-dim PCA vectors. 
  - Manually modify temp_pca.txt so it looks like *_pca.txt. 
    This involves deleting lines that start with "Reading:" and everything after "Writing:."

  - Calculate dot product for various systems.
    List systems in setSys() function in pca_dotp.cpp. 
    g++ -std=c++11 pca_dotp.cpp -o dotp
    ./dotp > out_dotp.txt

  - Plot dot product output.
    python3 plot_dotp.py

   - List of names in MODIFIABLE SECTION of plot_dotp.py need to be in same order as setSys mols list in pca_dotp.cpp. 

