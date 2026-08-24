# Creates a data sheet calculating total protein intensity normalized to the number of nuclei per cyst
# In order to run this code, you MUST have csv files acquired from running the "Protein_Intensity.ijm" macro in Fiji/ImageJ
# You ALSO MUST HAVE nuclei quantification csv files acquired from Napari-zelda to run this code
# Author: Dr. Madeline Lovejoy
# Date: 2026-08-21
# R version 4.1.2 (2021-11-01)

# dplyr version 1.1.2
# tidyr version 1.3.0
# stringr version 1.5.0

# load packages
# if you have never used these packages in your R console, you will need to install them first with the code: install packages ("package_name")
library(dplyr)
library(tidyr)
library(stringr)

# ***BEFORE RUNNING THIS LINE, MANUALLY SET THE WORKING DIRECTORY TO YOUR QUANTIFICATION FOLDER BYU GOING TO SESSION -> SET WORKING DIRECTORY -> CHOOSE DIRECTORY***
# defines working directory as "dir" data frame
dir<-getwd()

# Places Excel sheets ending in "E-cadInt.csv" in the "FilesE" data frame
FilesE <- list.files(pattern = "E-cadInt.csv", recursive = TRUE)
CSVlistE <- lapply(FilesE, read.csv)
names(CSVlistE) <- FilesE
# connects the data between columns to each other
data <-  bind_rows(CSVlistE, .id = "filename")

# separates the file name into different variables based on the placement of underscores ( _ )
data <- data %>% separate(filename, c( NA, NA, "Condition", NA, NA, NA, NA, NA, NA, "IMS", "ROI"), sep="_")

# Discards values within the protein intensity csv files that are "NA"
data_filtered <- na.omit(data)

# Adds the average integrated density values calculated from each z slice
sum_data <- data_filtered %>% group_by(Condition, IMS, ROI) %>% summarize(sum_area = sum(Area), sum_IntDen = sum(IntDen))

# Nuclei files from Napari
FilesN <- list.files(pattern = "Nuclei.csv", recursive = TRUE)
CSVlistN <- lapply(FilesN, read.csv)
names(CSVlistN) <- FilesN
CSV <-  bind_rows(CSVlistN, .id = "filename")

Nuc_Tally<-CSV %>% separate(filename, c( NA, NA, "Condition", NA, NA, NA, NA, NA, NA,  "IMS", "ROI"), sep="_")%>%
  # Discards any artifacts that are too small to be nuclei (less than 40 microns^3)
  filter(Volume > 40) %>%
  group_by(Condition, IMS, ROI) %>%
  # Counts the number of rows in each nuclei csv - each row should represent one nucleus
  tally() %>%
  # Eliminates the "Nuclei.csv" suffix from the end of the ROI, so different values can be grouped based on their ROI
  mutate(ROI =  str_replace(ROI, "Nuclei.csv", ""))

sum_data <- sum_data %>% mutate(ROI =  str_replace(ROI, "E-cadInt.csv", ""))

# Combines the sum integrated density and nuclei number (n) in the same data frame
combined_data <- left_join(sum_data, Nuc_Tally)%>%
# Divides sume integrated density by nuclei number
mutate(Int_per_Cell = sum_IntDen / n)

# Saves the combined data frame as a csv files
write.csv(combined_data,file = "Ecad_Nuc_Quant.csv" )
