# RPPA
Script for normalizing RPPA data

Installation:

1) install anaconda with python 3.10

2) clone this github repository

3) run the following in an anaconda prompt:
       pip install -r requirements.txt

4) copy input files into the folder. you will need at leas one of each below:
     proten stain files for normalization.  Theese need to end with "PROTEIN STAIN.xlsx"
     antibody files need to end with " <antibody_name>.xlsx", e.g. "... VEGFR2.xlsx"
     files grouping samples for plotting need to end with "GROUP.xlsx"
   
       examples for these files are in this repository
   
5) run the script
   A folder for protein fits will be created as well as one for the graphs along with .csv files with the data of each graph
