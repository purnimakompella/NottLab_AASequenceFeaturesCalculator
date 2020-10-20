# NottLab_AASequenceFeaturesCalculator

Author: Purnima Kompella
Email: purnima.kompella@bioch.ox.ac.uk
Last updated: September 28, 2020



#What you need to run it:
*This script accepts a directory of csv files with amino acid sequences. Each csv file should have the sequences in a column labeled "Residue". The csv file name should contain the name of the gene.

*The script requires the following packages: stringr and Peptides.

#What it calculates: 
* **Net Charge** at user-inputted pH which takes into account contributions from the following amino acids: D, E, C, Y, K, R, and H, and also the N- and C-termini. It also calculates "Simple" Net Charge which is not pH dependent, and only takes into account amino acids D, E, K and R.
* **Isoelectric Point**
* **Number of times is residue occurs in a sequence**
* **Fraction of each residue in the sequence**
* **Fraction of Charged Residues** ("Charge Density" for both Net Charge and "Simple" Net Charge)*8Fraction of Aromatic Residues


#Output
*The output is stored in a subdirectory in the input directory called "Analysis" (this is created if not already present).
*The output is one csv file containing all the calculations listed above for each csv file in the input directory
