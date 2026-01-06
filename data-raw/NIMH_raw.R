library(dplyr)
library(ggpubr)
library(huge)

###############
#Data selection
# ###############
#
## 'Biology' data block
#Load data
chem_panel <- read.table('acute_care.txt', sep = "\t", header = TRUE)
hep_panel <- read.table('hepatic.txt', sep = "\t", header = TRUE)
blood_count <- read.table('cbc.txt', sep = "\t", header = TRUE)
blood_chem <- read.table('blood_chem.txt', sep = "\t", header = TRUE)

#Remove variables that are full of missing data or not continious
hep_panel <- hep_panel %>% select(-3, -6)
blood_count <- blood_count %>% select(-14, -15, -16, -17, -19, -20)
blood_chem <- blood_chem %>% select(1, 6, 7, 26)

#merge the datasets by participant id
stepa <- merge(chem_panel, hep_panel, by = "participant_id")
stepb <- merge(stepa, blood_count, by = "participant_id")
stepc <- merge(stepb, blood_chem, by = "participant_id")

#replace -999 and non-values with NA
stepd <- replace(stepc, stepc<0, NA)
stepd[105, 11] <- NA
stepd[21, 11] <- NA

#remove all cases with one or more NA
bio <- stepd[complete.cases(stepd),]

## 'Psychology' data block
#Load data
alcohol <- read.table('audit.txt', sep = '\t', header = TRUE)
cog <- read.table('intel_vas.txt', sep = '\t', header = TRUE)
#Retain possibly useful variables
alcohol <- alcohol %>% select(1, 2)
cog <- cog %>% select(1, 5, 10, 15, 16, 17)
#merge datasets by matching participant id
psy <- merge (alcohol, cog, by = "participant_id")

## Combine bio and psy block
data <- merge (bio, psy, by = "participant_id")
#n = 183
#Remove pp id and make variables numeric
data <- data[,-1]
data <- apply(data, 2, as.numeric)
#remove variables that are cannot be transformed into normal distribution
data <- data[,-31]
data <- data[,-16]
NIMH_raw <- data

write(NIMH_raw)
usethis::use_data(NIMH_raw, overwrite = TRUE)
