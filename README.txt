TFBS_Visualization
==================

---------------------------------------------------------------------------
R_function--CB-plot_TFBS-landscape_TFBS-bi-motif_Dinucleotide-environment.R
---------------------------------------------------------------------------

If R is not yet installed, visit http://cran.r-project.org

To load the functions within the file R_function--CB-plot_TFBS-landscape_TFBS-bi-motif_Dinucleotide-environment.R, use the command source() e.g. source("/dir1/dir2/R_function--CB-plot_TFBS-landscape_TFBS-bi-motif_Dinucleotide-environment.R"), where dir1 and dir2 represent the location (absolute path) of the R code.

To view the usage of the four plotting functions, use the command view.help(), which is the last function provided in the R file.


Composition-bias plot 

function call: CB.plot(data.file, col.X=F, col.Y=F, headerline=T, ...)
Generates a plot with PFM GC nucleotide composition on the x-axis and motif enrichment score on the y-axis. The PFM GC nucleotide composition must be provided by the user if the motif enrichment software does not provide it. We only know of one motif enrichment software that provides PFM GC composition: oPOSSUM 3.0.

TFBS-landscape view

function call: landscape.view( data.file, col.X=F, col.Y=F, score.type="relative", headerline=T, ...)
Requires data for TF motif distance to the peakMax or peak centre of ChIP-Seq sequences (x-axis) and the associated PWM motif scores (y-axis).  PWM motif scores can be either the raw score, or the relative score, but the code will not handle p-values. If PWM raw scores are used, they may produce different results than the relative score. Relative scores allow single value hard-coded parameters to be used in the code (such as cutoff values) across all PWMs, but because the code could not assume access to the full range of each PWM's raw scores, it could not convert the hard coded parameter values to the equivalent raw scores specific to each PWM. The hard-coded parameters were approximated but will not be a perfect match for each PWM. Use the optional heuristic decision plots to confirm results.

When the optional arguments 'threshold' or 'heuristic.plot' are set to T (true), the enrichment zone boundary and the motif score threshold are returned as a vector. They can be stored in a variable, by setting a variable equal to the function call e.g. my_result = landscape.view( data.file, col.X=2, col.Y=5, score.type="relative", headerline=T, threshold=T)
To view the results type the variable name and hit enter/return. If a variable is not set equal to the function, the results will be printed to the screen.

TFBS-bi-motif view

function call: bimotif.view( data.file, col.X=F, col.Y=F, headerline=T, ...)
While col.X and col.Y provide the column identification for the motif_1 distance to peakMax and motif_2 distance to peakMax respectively, the R code will convert the distances to plot the spacing between motif_2 and motif_1 on the x-axis and the distance from motif_1 to peakMax on the y-axis.

Dinucloetide-environment view

function call: dinucleotide.view( data.file, ...)
Plotting this view first requires the use of the provided perl code: dinucleotide_align-on-position_and_count.pl
The perl code generates an output file with extension *.dinuc.align which is submitted to dinucleotide.view()
If the plot provides a motif alignment that is confusing, check the input into the perl code.


---------------------------------------------------------------------------
dinucleotide_align-on-position_and_count.pl
---------------------------------------------------------------------------

The dinucleotide alignment code, requires that the user have a fasta file for the sequences of interest, and know the location of the feature of interest (e.g. TFBS motif) relative to the start of each sequence. If no strand is provided, all feature and sequences are assumed to be in a 5'-3' orientation (also indicated with +ve sign). If strand is provided, the assumption is that a -ve sign indicates the feature is in a 3'-5' orientation in the sequence and thus the code will reverse compliment the sequence to orient the feature in a 5'-3' direction to align with other 5'-3' features.  

The perl code requires four arguments and has two optional arguments. Required is a fasta file of sequences, a seperate meta file with at least one column providing the offset of the feature relative to the sequence start (start = 1), or at least two columns if the strand/orientation of the feature is provided. There is no expected format for the meta file, therefore the columns containing the one (or two) types of information of interest must be specified. The last required argument is an output name for the results file, to which an extension of *.dinuc.align will be appended.

To view the options, type on the command line:  perl dinucleotide_align-on-position_and_count.pl
Options:
-f  = fasta seq file
-m  = meta file: motif offset RELATIVE TO START OF SEQ and strand if using it
-p  = column number containing motif offset (a column in -m metafile). First column is 1. 
-o  = output name  (output extension will be: *.dinuc.align)
-s  = (optional) column number containing strand +/- (a column in -m meta file). First columnn is 1.
-w  = (optional) motif 'width' i.e. length of PFM. Necessary for alignment adjustment between +/- strands.



