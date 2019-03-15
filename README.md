# genindToolbox

Some tools to export genind object (adegenet R package) into other formats (e.g. STRUCTURE input files, GENEPOP format).

 ## DESCRIPTION

 R function to export to STRUCTURE format from genind object. 'genind' objects are created in the R package adegenet.

 STRUCTURE input files are described as following by the STRUCTURE user manual:
 
 Each row of individual data contains the following elements. These form columns in the data file.
 
 1. Label (Optional; string) A string of integers or characters used to designate each individual in the sample.
 
 2. PopData (Optional; integer) An integer designating a user-defined population.
 
 3. PopFlag (Optional; 0 or 1) A Boolean flag which indicates whether to use the PopData when using learning samples (Note: A Boolean variable (flag) is a variable which takes the values TRUE or FALSE, which are designated here by the integers 1 (use PopData) and 0 (don’t use PopData), respectively).
 
 4. LocData (Optional; integer) An integer designating a user-defined sampling location.
 
 5. Phenotype (Optional; integer) An integer designating the value of a phenotype of interest.
 
 6. Extra Columns (Optional; string) It may be convenient for the user to include additional data in the input file which are ignored by the program.
 
 7. Genotype Data (Required; integer) Each allele at a given locus should be coded by a unique integer (eg microsatellite repeat score).

 ## USAGE
 genind2structure(obj, file="", pops=FALSE, missingdata.label=-9, popflag=NA, locdata=NA, phenotype=NA, extra=NA)

## ARGUMENTS
 obj: a genind object to convert
 
 file: file name to write
 
 pops: whether to include population info in the file, population names are taken in the genind object
 
 missingdata.label: value is the expected value for missing data
 
 popflag: if popflag is not NA, then it is a vector of popflags (numeric, 0 or 1) of the same size as the number of individuals
 
 locdata: if locdata is not NA, then it is a a numeric vector of user-defined sampling locations
 
 phenotype: if phenotype is not NA, then it is a numeric vector of phenotypes
 
 extra: if extra is not NA, then it is a character vector with extra informations

 example use:
 
 data(nancycats)
 
 genind2structure(nancycats, file="nancy_structure.txt", pops=TRUE)

 ## VALUE

 Return either a data.frame if file="", or write the data.frame into the file given in argument.
