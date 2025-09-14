# User Input Options

## User Input I: Steps 1-4 below must be changed to according to a specific dataset on a local computer.

Step 1:  Select the working directory. Indicate where the input file resides. The same folder will be the location for the output files.

Step 2: Select the input file. The input file should be a CSV file correctly formatted for the R script (Figure 2B). Specifically, the top row has the time series, and the first column contains individual sample positions on a 96-well plate. An example of such an input CSV file can be found in the supplemental material, NO7.csv.

Step 3: Name the samples and/or treatments. Please note this analysis is used for data obtained from 96-well settings with time course. The experiments can be designed as 8 replicates per treatment for a total of up to 12 treatments per plate, or as 12 replicates per treatment for a total of up to 8 treatments per plate, based on rows or columns of the 96-well plate. It is important to enter the correct number of samples, either 8 or 12, because the sample number counts. If the sample count is not at 8 or 12, the code will not run. However, if there are rows of empty wells, just list them as such in the space for treatment names, e.g. naming as empty 1, empty 2... It is worth noting that for sample names, one should avoid using backslashes or similar symbols as they may cause file name issues. Additionally, one should also avoid using "NA" as a label when selecting to use ANOVA in User Input II. Any treatments with the exact same name will be combined into one large treatment group.

Step 4: Indicate the relative start time for the assay. For a given day of 24 h, subjective dawn is when the chamber light is turned on during a normal light/dark cycle. For example, if light is on at 7 am, then this time is regarded as the circadian time 0. If a luciferase assay begins at 9 am, the assay time is at the circadian time 2. The user would input 2 as the relative start time. This number must be a whole number and cannot be negative. This will affect the phase reading of the samples.

## User Input II: Options 1-67 can be changed per the user's specific needs.  

Option 1: Include graphs. These are luminescence curves, plots for the period, phase, and amplitude per genotype and treatment. 

Option 2: Using the ANOVA test with the Tukey HSD. This will compare treatments on their period, phase, and amplitude.

Option 2 addition 1: Choose a control for the ANOVA test output; the user can specify the control for the ANOVA test or leave the option empty to use each treatment as a control in the pairwise comparison. 

Option 2 addition 2: Choose whether to number the ANOVA output files for quicker reference and organization. 

Option 3: Using a t-test to compare the data. This will compare treatments on their period, phase, and amplitude. This should only be used to compare two treatments at once, such as two lines of the same genotype with either SA or Mock added. The p-values shown are not adjusted to account for multiple comparisons with the same line.

Option 3 addition 1: Choose whether the t-test should be a pairwise comparison.

Option 3 addition 2: If you selected a pairwise t-test, select whether the data is paired.

Option 3 addition 3: Choose a control for the t-test; the user can specify the control for the t-test or leave the option empty to use each treatment as a control.

Option 3 addition 4: Choose whether to number the t-test output files for quicker reference and organization.

Option 4: Shifting the phase values. This allows the user to avoid values outside of the 24-hour cycle. To get low values or negative values instead of high values, the original phase value n can be subtracted by XX (i.e. 24) if the value of the original phase value n is greater than YY (i.e. 16). The user would need to set the values of XX and YY to use this option. 

Option 5: Choose whether to round the time points. This code calculates the period, phase, and amplitude using the ARS method 23. This method requires a consistent time series in hours to run; if the time points are off by only a minute or two, the user may round to the nearest hour and use this analysis. 

Option 6: Depending on how the plate reader records wells, the user may change the way the input is read. There are two options, one for the standard data listing of wells by A1, A2, A3... and one for wells listed by A1, B1, C1…
