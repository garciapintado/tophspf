# tophspf
A F90 version of the PWATER module in HSPF with distributed overland flow and Data Assimilation capability. A grid-to-grid kinematic routing and a 2D shallow equations (and adaptation of Geoclaw) are available as overland flow routing options, which are strongly coupled to the hydrologic model through the state variable considering the overland stored water. A generic interface serves to switch between the routing options. This version has been successfully tested for assimilation of satellite Synthetic Aperture Radar. 

<<<<<<< HEAD
The original assimilation software is written in R, which is also the scripting language to organize the input/output and model pre-/postprocessing. However, the original assimilation code in the above mentioned tests has been restructured as a core R package (rDAF) [available in Github] plus additional functions for the MPI parallelization and localization functions in the package rPDAF [available in Bitbucket under personal request to me]. This version need to be slightly udpdated to be made compatible with rDAF & rPDAF.
=======
The original assimilation software is written in R, which is also the scripting language to organize the input/outpu and pre-/postprocessing of the model. However, the original assimilation code in the above mentioned tests has been restructured as an R package (rDAF), and this version need to be udpdated to be made compatible with the rDAF, for which the MPI parallelization and the localization capabilites (available in this original version) are currently under preparation. The modules written as replacement for the original geoclaw still need to be submitted to Github. Please contact the author for details on progress or requests.
>>>>>>> 460dca915c12160321642b4ec6db0ca513ef6201
 
