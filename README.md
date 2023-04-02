# AptaZ

Repository of code associated with article "A high-dimensional microfluidic approach for
selection of aptamers with programmable binding affinities" published in Nature Chemistry

AptaZ is an algorithm that generates aptamers with desired binding affinity on the basis of high-content sequencing dataset.

(*We are working to publish a protocol/technical note paper regarding Apta Z, more complied code and multiple datasets for testing will be released with the technical note upon its acceptance) 

# Installation
•	Make sure your computer has at least 8GB RAM
•	Install MATLAB 2021a or later version

# Running
•	Convert the raw sequencing data (fastq) into a txt file with the following format using USEARCH or other search algorithms. Make sure to use semicolon ; as the spacer

#sequence_name;count;sequence

For example, #seq1;size=7;TGTAGCAGCACAGAGGTCAGATGTGTAGCAGCACAGAGGTCAGATG

•	Deposit all sequence txt files from the SORTED samples into a folder, index the files following the following format. Make sure the folder only contains the txt of SORTED samples.

Flow_rate –target_concentration–zone_number.txt

For example, 16-10p-z1.txt.
This name indicates the sequences were sorted under the condition of 16 mL/hr, 10p target concentration, and recovered from the first zone (zone 1).

•	Deposit the reference (unsort) txt file in any folders other than the SORTED samples.

•	Run ‘Z_score_calculation.m’, follow the instructions and select the reference file in txt format and the folder containing the txt files of SORTED samples.

The code will create a new folder named ‘Z-results’ under the current location ‘Z_score_caultion.m’ and store the calculated Z scores per condition there in the format of mat

The calculation typically takes 1 – 7 days depending on the size of dataset and consumes ~ 4GB RAMs during calculation.

•	Run ’Sum_Z_calculation.m’ and select the ‘Z-results’ or any renamed folder

The code will create a new folder named ‘Sum-Z-results’ under the current location ‘Sum_Z_calculation.m’ and store the SumZ scores per sequences there in the format of csv

•	Interpret the data using excel, ultraedit or other softwares.

# Contact
Please send your inquiry to zongjie.wang@northwestern.edu and cc shana.kelley@northwestern.edu
