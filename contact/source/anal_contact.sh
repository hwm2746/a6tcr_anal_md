# Get contacts using the whole frames. 

d=$1 ; f0=$2 ; t0=$3 ; t1=$4 ; iseg=$5

# d: output dir
# f0: filename prefix (vab, cab, tcra, tcrb, fg)
# t0: segid for the domain whose contact with TCR will be traced (domain1) 
# t1: TCR subunit name (domain2)
# iseg: segid for TCR subunit (E000,E001)
# frm_{ini/fin}: Initial/last frame. Set to -1 for default values

j=0
  f=${f0}
  echo $f
  
  ### hbond ###
  cd $d
  fn0=temp_hb$f.inp   # charmm input file
  echo "* created by anal_contact.sh" > $fn0
  echo "*"  >> $fn0
  echo "set s $s" >> $fn0
  echo "stream ../input0.str" >> $fn0
  echo "define t0 sele ${t0} end" >> $fn0
  echo "define t1 sele ${t1} end" >> $fn0
  echo "" >> $fn0
  echo "stream ../contact_hb.str" >> $fn0

  charmm < $fn0 > out_temp_hb$f.dat
  
  fn1=input_hb$f.dat  # contact analysis file
  echo "ifname out_temp_hb$f.dat" > $fn1
  cat ../input_hb.str >> $fn1
  echo "segname $iseg" >> $fn1
  echo "ofname hb_$f " >> $fn1
  ../te $fn1

  ### nonpolar ###
  fn0=temp_np$f.inp   # charmm input file
  echo "* created by anal_contact.sh" > $fn0
  echo "*"  >> $fn0
  echo "set s $s" >> $fn0
  echo "stream ../input0.str" >> $fn0
  echo "define t0 sele ${t0} .and. (property abs CHARGE .lt. 0.30) end" \
       >> $fn0
  echo "define t1 sele ${t1} .and. (property abs CHARGE .lt. 0.30) end" \
       >> $fn0
  echo "" >> $fn0
  echo "stream ../contact_np.str" >> $fn0
  charmm < $fn0 > out_temp_np$f.dat

  fn1=input_np$f.dat  # contact analysis file
  echo "ifname out_temp_np$f.dat" > $fn1
  cat ../input_np.str >> $fn1
  echo "segname $iseg" >> $fn1
  echo "ofname np_$f " >> $fn1
  ../te $fn1

  cd - 

#done
