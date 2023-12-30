# Get mhc - V{alpha,beta} contacts and pep-Vab contacts.
# MODIFIABLE FILE. 

bash anal_contact.sh data_ma ma MHC0 Valpha MHC0 > out_contact_ma.dat
bash anal_contact.sh data_mb mb MHC0 Vbeta MHC0 > out_contact_mb.dat
bash anal_contact.sh data_pa pa peptide Valpha C000 > out_contact_pa.dat
bash anal_contact.sh data_pb pb peptide Vbeta  C000 > out_contact_pb.dat
