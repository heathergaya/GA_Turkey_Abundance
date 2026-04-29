Code for manuscript: Spatial mark-resight models reveal declining wild turkey population in Georgia, USA

There are 3 folders in this repository. 

Simulation
-

- Abundance_Simulation.R
  + Simulation of data to evaluate precision and accuracy of abundance estimates
  + Includes code to evaluate simulation results
- NimbleDists.R
  + Nimble functions used in analysis 
- Sim_Models.R
  + Code to run simulation analysis
  + Requires runMCMCbites2.R in CaseStudy folder
- sim.res.csv
  + Output from simulation

CaseStudy
- 

- NimbleDists.R
  + Repeat of file in simulation folder
- TwoStageSMR.R
  + Code to run the 2 stage transmitter SMR model
  + uses Stage1_inputs.rds and Stage2_Winter.rds
  + To just run stage 2, use Feb20_Stage1_res.rds object
- SMR_BandsOnly.R
  + Code to run the band SMR model
  + uses SMR2_inputs.rds 
- runMCMCbites2.R
  + code to have nimble output results every 500 iterations 

RawData
-

  - CamTaggedTurkeys_20_25.csv
    + All turkeys with transmitters/bands detected on game cameras
  - Camera_effort_Feb14_Marc14.csv
    + Days cameras were on/off
  - WinterCams_20_25.csv
    + All turkey detections (marked or unmarked)
  - WinterCamera_Locations.csv
    + Locations of camera traps
    + Used in both case study and simulation
  - bf grant boundary.shp
    + Boundary of BF Grant property
    + Used to determine spatial area
  - cedarcreek_boundary2014.shp
    + Boundary of Cedar Creek WMA
    + Used to determine spatial area
    
