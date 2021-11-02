# Vi-Explore
Virus Identification Pipeline for NParks/CAVS

Vi Explore is a tool for quickly and easily extracting accurate, meaningful information from Next Generation Sequencing data.

Please also have generate_html.py and style.css, when running generate_illumina_virus_detection_pipeline.py for best results.

*********************************************
Arguments:

-h, --help            show this help message and exit

-o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                      Takes in path of directory that you want to create to
                      store the output
                      
-v, --version         Tells you the version of script

-s SAMPLE_FILE, --sample_file SAMPLE_FILE
                      Takes in the data that you wish to use to run FastQC
                      on

*********************************************

Example format when running generate_illumina_virus_detection_pipeline.py:

./generate_illumina_virus_detection_pipeline.py -v -o /home/melody/testdir/vnnv_samples -s /home/melody/testdir/vnnv_samples.sa

*********************************************

Running generate_illumina_virus_detection_pipeline.py creates a directory to store the output as well as a makefile to generate the output.


  To run the makefile:
  
  make -f illumina_virus_detection_pipeline.mk
  
*********************************************

Multiple paired end reads can be run when using this tool.

Please store your sample file as a .sa file in the following format: 
{Sample name} {first read} {second read}
  
Ensure that the contents are separated by tabs and not spaces.
