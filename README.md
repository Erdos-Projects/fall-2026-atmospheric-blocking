# fall-2026-atmospheric-blocking
Team project: fall-2026-atmospheric-blocking

The data_excerpt directory contains 4 years worth of data. This will be kept within the 
Github repo for easy testing of files. The full 40 years of data (including these 4)
will be stored on Google Drive. 

To visualize the potential temperature theta, run this from the root (fall-2026-atmospheric-blocking)
directory: 

python plot_era40.py --theta data_excerpt/theta/theta2pvu_199812.grib --time 1998-12-15T12

The above command displays a heatmap showing the potential temperature at 2PVU at 12:00 on 15 December 1998. 
The below command adds arrows to the heatmap showing the direction of the wind.

python plot_era40.py --theta data_excerpt/theta/theta2pvu_199812.grib --uv data_excerpt/uv/uv300_199812.grib --time 1998-12-15T12
