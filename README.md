![](https://img.shields.io/github/v/release/CNIC-Proteomics/RefMod)

![RefMod logo](/assets/RefMod_logo_text_white.svg)

# Description
 RefMod is an implementation of the ReCom concept ([Laguillo-Gómez et al., 2023](https://doi.org/10.1016/j.jprot.2023.104968)) designed to be run as a post-processing step after an open [MSFragger](https://msfragger.nesvilab.org/) search. It is compatible with both DDA and DIA data. When using RefMod, DIA data can be searched in a “pseudo-DDA” workflow, using a curated list of theoretical Δmass values to correct errors caused by the uncertainty in precursor masses contained within the same fragmentation window.
 
 RefMod has been developed at the [Cardiovascular Proteomics Lab](https://www.cnic.es/en/investigacion/cardiovascular-proteomics) / [Proteomics Unit](https://www.cnic.es/en/investigacion/proteomics) at CNIC (Spanish National Centre for Cardiovascular Research).

# Quick Start Guide
This section provides a brief overview of the basic requirements to run RefMod from the command line. A GUI is also available and can be launched using the `RefMod.bat` file located in the root directory. If this is your first time using RefMod, or if you’re looking for an in-depth explanation of any step, please refer to the [full documentation](/docs/docs.pdf).
1. Download the latest release from the GitHub page.
2. Install Python 3.11 and the packages from `requirements.txt`.
3. Prepare your configuration file and run RefMod from the GUI or from the command line:
`$ python RefMod.py -i [MSFragger files] -r [MGF or mzML files] -d [Theoretical Δmasses file] -c [Configuration file]`
5. Find the output files in a newly created `refmod` directory in the same path specified by the -i option.
          
          
