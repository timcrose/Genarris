rm Genarris*
rm ZELDO_4_volume_estimate/*lock
rm ZELDO_4_volume_estimate/*failed
rm -r relax_tmp_ZELDO_4/*
python ../../../../genarris/src/genarris_master.py relaxation.conf
