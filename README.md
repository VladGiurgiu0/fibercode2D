# fibercode2D
Matlab codes to process 2D images containing tracers and fibers to measure fiber and flow velocities.

For an example of use of this code see: https://doi.org/10.1016/j.ijmultiphaseflow.2024.105021

Procedure:
1. Acquire images of fibers and spherical tracers in a 2D-PIV configuration and save them with the naming format "A0001.tiff".
2. Use main.m to discriminate the fibers from the tracers. The "compute" section saves images containing only fibers separetely from ones containing only tracers. Note: "loop" means a subfolder in the root_folder and main.m can be used to run the same operations over multiple subfolders in sequence and this is useful when multiple statistically independent sets ("loops") of time-resolved images have been acquired. 
3. Use PaIRS (https://pairs.unina.it) to process the images containing only tracers and compute the PIV velocity fields.
4. Use main.m to process the velocity fields. The "AFTER PIV: compute flow and quantities" section filters and interpolates the velocity fields and computes vorticity, shear rate, and swirling strength.
5. Use fiber_statistics_v2.m to extract fiber computed quantities, e.g. fiber velocity, and flow quantities at the fiber location, e.g. flow velocity, from all the "loops" and save them in single matrices, one for every quantity.

Use the data in "example_data" folder for a processing example.


Things to do:
1. Check in "fct_track_fibers.m" section "%% correct angles" if the logic of fixing the orientation angles works well for computing rotation rates. 
