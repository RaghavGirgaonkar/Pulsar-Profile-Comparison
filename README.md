# Pulsar-Profile-Comparison

A script to compare two pulsar profiles using the Fourier Domain Alignment Method proposed by Taylor (1992)

# Usage

#### python3 profilecomparison.py profile1.txt profile2.txt

The above command compares profiles 1 and 2 by taking profile 1 as the reference profile and scaling and aligning profile 2 to profile 1.

If you want to compare a single reference profile to multiple profiles,

#### python3 profilecomparison.py profile1.txt profile2.txt profile3.txt profile4.txt

The above command compares profiles 2 to 4 to profile 1 by taking profile 1 as the reference profile.

# Output

The code outputs "Difference" files for all of the profiles comapred. It outputs the chisquare value for each file along with the peak deviation value. 

Along with this, a comparison plot with the two profiles overlaid and their differences is plotted for all pairs of files compared.

# Note

1. The Chi-square value is calculated using the off pulse region of the reference profile (i.e template), the code is still the same but your off pulse region may be different depending on which profile you use for profile1.txt. This will have to be manually changed in the get_chisqrs() function in the script if needed.
2. Please ensure that all txt files have the same number of bins, other wise the script will skip those files.
3. Please also ensure that the txt files have only 1 column, i.e the values/counts of the profile. This can be ensured by selecting the correct column while doing pdv -t
4. Also ensure that the fits files are frequency and time scrunched and also de dispersed. (And the nbins of all are the same) before their text files are made.
