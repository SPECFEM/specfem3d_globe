------------------

Reference profiles

------------------

These depth profiles have been created by

- based on Par_file setting from s362ani:
  > cp EXAMPLES/global_s362ani_shakemovie/DATA/Par_file ./DATA/Par_file

- modified model name:
  > sed -i "s:^MODEL .*:MODEL   = 1d_1066a:" DATA/Par_file

- at single location lat/lon = 47.5 / -22.0
  > ./bin/xwrite_profile 1 1 47.5 -22.0

The created output file OUTPUT_FILES/CARDS_th0042_ph0338 was copied into this folder
as reference profile to test the consistency of the model creation routines.


