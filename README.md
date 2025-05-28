# iSTARMOD
python code implementing Spectral Subtraction Technique to study chromospheric activity in Cool Stars

SPECTRAL SUBTRACTION TECHNIQUE

Synthetic Spectrum:

  •	Reference star spectra:
  
	    Inactive (no emission in the activity indicator considered)
     
    	    Low projected rotational velocity (vsini)
	 
	    Same effective temperature (Teff) and luminosity class
     
  •	Doppler shift is applied
  
  •	Rotational broadening (vsini)
  
  •	Weighting based on the relative intensity in the Halpha and continuum regions for the primary and the secondary component
  

Subtracted Spectrum:

  •	Computed as: *Observed Spectrum – Synthetic Spectrum*, thus isolating the chromospheric contribution
  


iSTARMOD Program

This program allows for the creation of a synthetic spectrum composed of one or two stars, each with specified relative contributions or weights. Once constructed, the synthetic spectrum is subtracted from the original observed spectrum to yield the subtracted spectrum.
  •	The command is executed as follows:

starmod(file=filename.sm)

  •	The file must have the .sm extension. Upon editing or running this file, the user provides stellar parameters including:
    o Filename of the .fits file of the object
    o Number of iterations
    o	Wavelenght range to consider
    o	Optionally, wavelength intervals to exclude (as shown in the example)
    o	Output filenames for the synthetic and subtracted spectra
    
    o	Parameters for the primary and secondary stars separately, including:
    
        - Reference star (of the same spectral type and luminosity class): filename of the fits file
	
        - Radial velocity
	
        - Projected rotational velocity (vsini)
	
        - Spectral weight  
	
        - Whether each parameter is fixed or variable during the iteration (denoted by *fix* or *var*)
	
        - The aperture number (spectral order) is also included
        
Once all parameters are set, execution begins. This generates the synthetic and subtracted spectra and prints the results of the iterations along with the final fit for both components.

As noted, the command starmod("filename.sm") must be executed from the same directory where the executable is located. From the files included here:



from istarmod_automat_reloaded import *

starmod("pwand_n3_cah_34_wvl.sm")

