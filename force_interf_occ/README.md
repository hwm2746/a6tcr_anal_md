**\./force_interf_occ/ : Files for analyzing trajectory force, interface, and contact occupancy.**

- Note that we have not included CHARMM trajectory (*.dcd) files due to size. Use below as a guide to apply to your own trajectory. For this, you need to modify file paths, e.g., in dcd_header.str and include.str

************
* ./source/ : Source scripts.
  	      Executing these will generate the files in ../contact/data.

[0] Some input files necessary for these scripts were generated in previous steps.
    - These files are in ./source/in_data/.

[1] Calculate instantaneous force.
    charmm -chsize 200000 < get_pos.inp > out_get_pos.dat
    grep -ia 'Parameter: [XYZQ]' out/out_get_pos.dat | sed -e 's/ Parameter: //g' | sed -e 's/"//g' | sed -e 's/<-//g' > pos.dat
    python3 tows.py -o f_WThigh_40_0 -w 2000 -b 10000
    python3 tows.py -o f_WThigh_40_1 -w 2000 -b 11000

    - Files get_pos.inp and tows.py need to be modified for each structure. 
    - Here, we include the scripts for WThigh. 
    - Streamed file domain_1ao7.str will need to be modified for a different structure. The example included here pertains to A6 TCR. 
    - In get_pos.inp, modify MODIFIABLE SECTION.
    - In tows.py, modify MODIFIABLE SECTION.
    - Output of the above is placed in ./source/in_data/. 

[2] Calculate contact occupancy.
   - Compile.
     g++ -Wall -std=c++11 calc_fvocc_noload.cpp -o fvocc_nl   
     g++ -Wall -std=c++11 calc_fvocc.cpp -o fvocc
   - Execute.      
     ./fvocc_nl WT0 40 pmhc   " "  .5 .8 
     ./fvocc_nl WT0 40 intranoc   " "  .5 .8 
     ./fvocc WTlow 40 pmhc   " "  .5 .8 
     ./fvocc WTlow 40 intranoc   " "  .5 .8 
     ./fvocc WThigh 40 pmhc   " "  .5 .8 
     ./fvocc WThigh 40 intranoc   " "  .5 .8 

   - Search for "ifname=" in calc_fvocc\*.cpp. calc_fvocc*.cpp needs:
     - ./source/in_data/\*{hb,np}_*.dat     =  contact prescence
       - Created in ../contact/ during step [3].
       
     - ./source/in_data/\*/f\_\*\_\*\_\*\_interval.dat  =  instantaneous force (WTlow and WThigh only)
       - Created during step [0] above.

   - In "Define domains" section in \*.cpp files, there are several domain combinations built-in for which to calculate total occupancy.
   - We use "pmhc" and "intranoc" in this example.
   - Other domain combinations can be selected during the execute step. 
   - Output messages from executing the code help track how total occupancy is determined.
   - May need to update "dt" variable = occupancy data saving rate (search for "dt=" in calc_fvocc.cpp)
   - Frame rates (frate, fstart) are set in *.cpp function. These should be adjusted based on frame rate from step [1].

[3] Get the least-square fit plane of peptide.
    - Need the V-module oriented trajectory and backbone-only psf from ../vv_motion/source step [1].
    charmm < calc_pep_lsqp.inp > out_calc.dat
    bash get_lsqp.sh
    rm out_calc.dat

**********
* ./force_interf_occ/data/ : Output data generated from source scripts. 

  - Each folder contains the following data for the respective trajectory:

  - For force vs. time with contact occupancy:
    fvocc_pmhc.txt       = TCRab-pMHC occupancy 
    fvocc_intranoc.txt	 = Intra-TCRab occupancy, exclude Ca-Cb contacts
    - Outputs from step [2] above. 

  - For peptide angle:
    pep_lsqp.txt
    - Output from step [3] above. 

**********
* ./force_interf_occ/anal/

- Fig 8A,B: Instantaneous force versus time with contact occupancy.
  python3 force_occ.py 

- Fig 8C: Peptide angle.
  python3 angle_pep.py
  - Calls va.pdb and vb.pdb from ../vv_motion/data/.
  - These files are generated during ../vv_motion step [5]. 
