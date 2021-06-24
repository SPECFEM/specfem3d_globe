----------------------------
Mars - Topography
----------------------------

topography:

- to download and store in SPECFEM binary format, type:
  ./run_create_topo_bathy_file.py mars4

  uses smoothed version of original MOLA data set with 32-pixel per degree resolution
  original topography min/max = -8206 / 21181 m


Mars Orbital Laser Altimeter (MOLA) data:
  https://pds-geosciences.wustl.edu/mgs/mgs-m-mola-5-megdr-l3-v1/mgsl_300x/aareadme.txt

  Mission Experiment Gridded Data Record (MEGDR)
  32-pixel per degree data set:
  - megt90n000fb.img      - data
    https://pds-geosciences.wustl.edu/mgs/mgs-m-mola-5-megdr-l3-v1/mgsl_300x/meg032/megt90n000fb.img

  - megt90n000fb.lbl.txt  - detailed table info
    https://pds-geosciences.wustl.edu/mgs/mgs-m-mola-5-megdr-l3-v1/mgsl_300x/meg032/megt90n000fb.lbl


citation:
  The MOLA MEGDR data set may be cited in published literature using the
  following reference:

     Smith, D., G. Neumann, R. E. Arvidson, E. A. Guinness,
     and S. Slavney, "Mars Global Surveyor Laser Altimeter Mission
     Experiment Gridded Data Record", NASA Planetary Data System,
     MGS-M-MOLA-5-MEGDR-L3-V1.0, 2003.



