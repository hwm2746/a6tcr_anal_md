**\./contact/ : Files for analyzing contacts from trajectory.**

- Note that we have not included CHARMM trajectory (*.dcd) files due to size. Use below as a guide to apply to your own trajectory. For this, you need to modify file paths, e.g., in ./source/input0.str

************
* ./source/ : Source scripts.
  	      Executing these will generate the files in ./data.

[0] Compile C++ codes to write executable file (te) to get various contacts.
    g++ get_bond.cpp ftn_bond.cpp -O3 -o te
  
[1] Create directories. Since we will get contacts for MHC-Valpha, MHC-Vbeta, peptide-Valpha, and peptide-Vbeta, we create four folders to save data. 
    mkdir data_ma data_mb data_pa data_pb
    (These folders are already present in subfolders of ./data)

[2] Copy the executable "te" created in step [0] in each directory, or make symbolic links to it.

[3] action_*.sh = lists the interfaces for which we want to get contacts.

- Execute the script, e.g.: nohup bash action_pmhc.sh

[4] Check output
    grep -a 'SEVERE' data*/out*.dat                   
    - If no errors, remove temporary files:
      rm -f data\*/{\*temp\*,out\*dat,\*inp}  
      rm -f ./{\*temp\*,out\*dat}  

**********
* ./contact/data/ : Files generated from source scripts. 

  - Each directory contains four subdirectories: data_ma, data_mb, data_pa, data_pb. 
  - subdirectories contain data files used in analysis scripts.
    - hb = hydrogen bond; np = nonpolar contacts. 

    - {hb,np}\_*_occ.dat	   = average contact occupancy.
    - {hb,np}\_*_occ_traj.dat = occupancy trajectory.
    - {hb,np}\_*.dat 	   =  Raw contact trajectory (0=no, 1=yes).

**********
* ./contact/anal/ : Analysis scripts and outputs.
  Calls files in ./contact/data/. 

  ./contact/anal/out/ : Output files from included scripts.
  These files are used to generate figures in our paper.

  - Fig 3A: Number of contacts with pMHC.

  - Print avg and std: 
    python3 count_all.py -ocut .8 -ocut1 .5 -ofn out/ncont_
    
    -ocut : instantaneous occupancy
    -ocut1 : overall average occupancy
    -ofn : output filename prefix

   - Plot: python3 plot_count.py
  
  - Fig 3B: Hamming distance.

    - Calculate the Hamming distance:
    g++ hamming_dist.cpp ftn_bond.cpp -O3 -o te1   
    ./te1 input_hamming_WThigh.dat
    ./te1 input_hamming_WTlow.dat    

    - te1 is the executable file. 
    - Output: ./out/hamming_WThigh.dat, ./out/hamming_WTlow.dat. 

    - Plot: python3 hamming_dist.py 
  
- Fig 3C: WT high contact occupancy heat maps.

  - Plot heat maps: 
    python3 occ_traj.py -o ./out/WThigh_map_ma -n namedef_ma ;
    python3 occ_traj.py -o ./out/WThigh_map_mb -n namedef_mb ;
    python3 occ_traj.py -o ./out/WThigh_map_pa -n namedef_pa ;
    python3 occ_traj.py -o ./out/WThigh_map_pb -n namedef_pb ;

  - Flags: -o=output file prefix, -n=input file 
  - Output: ./out/WThigh_map\_{ma,mb,pa,pb}\_{hb,np}_occ.pdf
