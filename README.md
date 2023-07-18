# operon_analyzer
operon_analyzer is a Python-based tool for analyzing operon structures in bacterial genomes. This repository provides a collection of scripts and functions to process and analyze gene annotations from ptt and gff files, with a focus on identifying operons based on gene proximity and strand orientation.

Programming Language and Version: Python 3.9.13

Required packages/Libraries: 
1, Pandas
- To install the Pandas package, type the following command in the command prompt
	py -m pip install "pandas"
- To import the package to start using it, use the following command
	import pandas as pd

Required input files: 
1,B_subtilis_168.ptt
2,E_coli_K12_MG1655.ptt
3,Halobacterium_NRC1.ptt
4,Synechocystis_PCC6803_uid159873.ptt
5,2088090036.gff

Description: 
-The read.csv function in the code reads the ptt and gff files as a pandas data frame.
-The defined function "process_dataframe" takes one argument: data frame. It makes necessary changes to the data frames derived from the ptt files.
-The defined function "get_operons" takes two arguments: the data frame and the column whose values are used to construct operons. The function calculates the total number of operons present in each data frame by checking if the distance between 2 genes is less than 50 nucleotides and if the genes are on the same strand. 

Execution: 
1, Open Command Prompt.
2, Change the working directory to the location of the file.
	cd path\to\the\file\location
3, Run the script using the following command.
	operon_analyzer.py
4, Once execution is complete check output.
	
Output Files:
None

Author: Amulya Saini
Date: 03/24/2023
