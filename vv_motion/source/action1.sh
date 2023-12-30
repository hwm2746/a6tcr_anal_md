dname=$1
bash action_triad.sh va ${dname}
bash action_triad.sh vb ${dname}

bash action_cm.sh ca ${dname} &
bash action_cm.sh cb ${dname} &
bash action_cm.sh ha ${dname} &
bash action_cm.sh hb ${dname}
