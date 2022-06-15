# METIS
* Publication: Pandi A*‡, Diehl C*, Kharrazi AY, Faure L, Scholz SA, Nattermann M, Adam D, Chapin N, Foroughijabbari Y, Moritz C, Paczia N, NS Cortina, Faulon J-L, Erb TJ‡. A versatile active learning workflow for optimization of genetic and metabolic networks. bioRxiv (2021). 10.1101/2021.12.28.474323

## Full Performance Can Be Achieved On Google Colab Online
* [Optimization NoteBook](https://colab.research.google.com/github/amirpandi/METIS/blob/main/METIS_Optimization_Notebook.ipynb)
* [Prediction NoteBook](https://colab.research.google.com/github/amirpandi/METIS/blob/main/METIS_Prediction_Notebook.ipynb)

## Guide
* Apart from the two Google Colab notebooks above that can be tailored for custom applications. All example notebooks used in this work along with their data (combinations and yields) are provided in the Applications folder of the METIS repository. 

* To tailor METIS for a custom application, check Fig. 2 of the above paper (see below), and a detailed step-by-step guide in Supplementary Note 2 and Supplementary Fig. 4-6.


![Fig2  METIS](https://user-images.githubusercontent.com/55136474/173783016-43c756bd-0f14-4a66-b00e-743660e4bdba.png)



## Examples
This folder contains METIS applications as well as other scripts used in the paper:

## Folders Description:
### Cell_Free_System:
* Data: Data for 10 Days of active learning optimization (Fig. 1e).
* Code Description: This folder uses old version of METIS algorithms and codes. The adapted METIS notebook for lysate cell-free systems is also available under Code_METIS.

### Cell_Free_System_Simulation:
* Data: Data for the simulation on Borkowski et al. (2020) data (Fig. 1b).
* Code Description: This folder uses an old version of METIS algorithms and codes.

### Enzyme_Engineering:
* Data: Data for 800 enzyme mutants from Nattermann et al (2021) (Supplementary Fig. 19).
* Codes: METIS optimization and simulation notebooks. 

### LacI (Gene Circuit):
* Data: Data for 10 Days of active learning optimization (Fig. 3c) and 2 rounds of "20 Most Informative" experiments (Fig. 3k).
* Codes: METIS optimization notebook. 

### PURE (Cell-Free System:
* Data: Data for 10 Days of simulation (Supplementary Fig. 19).
* Codes: METIS prediction notebook.
* 
### TTU (Transcription & Translation Unit):
* Data: Data for 4 Days of active learning optimization (Fig. 4c) and 2 rounds of "20 Most Informative" experiments (Fig. 4e).
* Codes: METIS optimization notebook. 

### CETCH
* Data: Data for 5 Days of "yield" active learning optimization (Fig. 5c) and 5 Days of "efficiency" active learning optimization (Fig. 5e).
* * Codes: METIS optimization notebook.


## Documentation
* All functions have Help as DocString (i.e. Help(Name of Function))

### Support
* Feel free to open an issue.
