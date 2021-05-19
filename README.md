# Scripts

## BedGraphToBigWig.py
Script that takes as input a BedGraph file, sorts it, and obtains a BigWig file. It uses bedGraphToBigWig from UCSC.

## BinGenome.py
Script that divides mouse genome (mm10) in same size bins.

## H3K27me3_LoessNormalization_example.R
Example of loess normalization of ChIP-seq levels across differentiation time points. Loess is fitted in the whole genome (in 2KB bins), and the corrections are aplied to the regions of interest (BP and PE). Input files are outputs of function recoverChIPlevels from SeqCode.

## LM_Enhancer_ggplot.R
Script that learns gene expression predictive models from ChIP-seq levels. The input is a file that contains gene expression (FPKMs) and ChIP-seq levels for several histone modifications. The model is trained in a subset of genes (training set), and tested in another one (test set). It uses a 10X cross-validation repeated three times to learn the model. It also trains a random predictive model which randomizes gene expression of the training set. The predictive model and the random predictive model are both evaluated in the same training set. It outputs text files with summaries of the models and scatter plots which represent predicted vs. measured expression.

## LM_Enhancer_ggplot_rd.R
Similar to LM_Enhancer_ggplot.R, but in this case models are trained in one cell type, and evaluated in a different cell type.

## Merge2bams.py
Script that merges two bam files into one, and produces a BedGraph with the profile. It calls samtools and function BuildChIPprofile from SeqCode.

## MergePeaks.py
Script that merges two lists of peaks and homogenizes coordinates of overlapping peaks. It outputs the merged list of peaks, and the specific peaks of each list after homogenizing coordinates.
