#Script for the analysis of pheonotypic data

#Starting point is the mean growth of 3 replicates per strain and condition of different parameters:
 # - Growth rate
 # - Lag phase duration
 # - Growth at 24 h (Eff24)
 # - Growth at 40 h (Eff40)
 # - Growth at 48 h (Eff48)
 # - Growth at 65 h (Eff65-end of cultures)
#Eff48 is the parameter chosen for analysis as a proxy for the growth rate and the final biomass produced.

# -----------------------------------------------
# 1. Directories
# -----------------------------------------------

getwd()
setwd("F:/Sergio I/Data_analysis/repglobal")

# -----------------------------------------------
# 2. Loading data and converting to z-score
# -----------------------------------------------

data <- read_excel("F:/Sergio I/Data_analysis/repglobal/20240224_meanparameters_resumen.xlsx")

# Subset the parameter to be represented
# Using bracket indexing to select columns that start with "Eff_" and the "strains" column
Eff_df <- data[, c(grep("^Eff48_", names(data), value = TRUE), "strain...1")] # Dataframe of final time efficiency means

# Rename the column Eff_df$strain...1 to Eff_df$strain
colnames(Eff_df)[colnames(Eff_df) == "strain...1"] <- "strain"

# Convert strain column to numeric
Eff_df$strain <- as.numeric(Eff_df$strain)

allyarro <- Eff_df

# Subset to remove control strains, non-Yarrowia lipolytica strains, and non-existing ones
Eff_df <- subset(Eff_df, !Eff_df$strain >= 131)
print(names(Eff_df))

# Subset to remove the unsuccessful medium, the last Galactose
# Remove media related to inoculum density
datasub <- Eff_df[, -which(names(Eff_df) == "Eff48_YNB_1%_Glucose_Lowinoculum")]
datasub <- Eff_df[, -which(names(Eff_df) ==  "Eff48_YNB_1%_Galactose")]
datasub <- datasub[, -which(names(datasub) ==  "Eff48_YNB_1%_Glucose_0.5gL_FerulicAcid")]

# Check range to determine replacement value for NAs
range(datasub[, 1:42], na.rm = TRUE)

# Several transformations are required for the heatmap
# First compute base 10 log, then z-score
logdata <- log(datasub[, 1:41], base = 10) # Do not include the strain column

range(logdata[,1:41 ], na.rm = TRUE)

## Treatment of inhibitor and nutrient source categories will be handled differently
# Nutrient source <- direct z-score
# Inhibitor <- normalization by growth in glucose followed by z-score

# We have to separate the datasets
mediosliq <- c("CN Intermediate", "Lysine", "Urea", "Glutamine", "Cinnamic Acid", "Nitrates", "YLD",
               "Galactose", "Ribose", "Ethanol", "Fructose", "Glucose", "SDS",
               "Ferulic Acid", "Ethylene Glycol", "Cycloheximide", "Furfural", "Vanillin", "13ºC",
               "Zinc", "Hydroxymethylfurfural", "Cadmium", "37ºC", "4% NaCl", "Copper", "Nickel",
               "Cobalt", "pH 3", "pH7", "Glycerol", "Hexadecane", "Maltose", "Mannose", "Melibiose",
               "Methanol", "Acetate", "Citrate", "Lactate","Propionate", "Tributyrin", "CNHigh")

# Rename the first 41 columns of datasub using the names in mediosliq
colnames(logdata) <- mediosliq
colnames(datasub)[1:41] <- mediosliq

nutrient <- c("Lysine", "Urea", "Glutamine", "Nitrates", "YLD",
              "Galactose", "Ribose", "Ethanol", "Fructose", "Glucose", "Glycerol", "Hexadecane", "Maltose", "Mannose", "Melibiose",
              "Methanol", "Acetate", "Citrate", "Lactate","Propionate", "Tributyrin")

# Create a new vector excluding 'nutrient' elements from 'mediosliq'
mediosinh <- setdiff(mediosliq, nutrient)

# Select logdata columns whose names are in the nutrient vector
datanutrient <- logdata[, colnames(logdata) %in% nutrient]

datainh <- logdata[, colnames(logdata) %in% mediosinh]

# Process the carbon and nitrogen sources
# z-score transformation
zdata <- as.data.frame(scale(datanutrient))

# try <- scale(logdata$`Eff_YNB_1%_Glucose_37ºC`) It seems that z-score is applied column-wise
range(zdata, na.rm = TRUE)

# The lowest value is -7.79, we assign it to NAs
# Replace NA values with -5.592080
zdata <- replace(zdata, is.na(zdata), -5.592080)

# # Remove "Eff_" from column names
# names(zdata) <- gsub("^Eff48_", "", names(zdata))
zdata <- zdata[, 1:41]

# Add strain column
zdata$strains <- datasub$strain

datanutrient <- zdata

#Now the inhibitors
datainh <- datasub[, colnames(datasub) %in% mediosinh]

# z-score transformation

# Divide all columns except 'Glucose' by the 'Glucose' column
inh_divided <- datainh
inh_divided[, colnames(inh_divided) != "Glucose"] <- 
  datainh[, colnames(datainh) != "Glucose"] / datasub$Glucose

zdata <- as.data.frame(scale(inh_divided))

# try <- scale(logdata$`Eff_YNB_1%_Glucose_37ºC`) It seems that z-score is applied column-wise
range(zdata, na.rm = TRUE)

# The lowest value is -7.79, we assign -2.53 to NAs
# Replace NA values with -2.53
zdata <- replace(zdata, is.na(zdata), -2.53)

# # Remove "Eff_" from column names
# names(zdata) <- gsub("^Eff48_", "", names(zdata))
# zdata <- zdata[, 1:20]

# Add strain column
zdata$strains <- datasub$strain
datainh <- zdata

# We merge both datasets
liquido <- cbind(datanutrient[1:21], datainh)

# -----------------------------------------------
# 3. MAtrix building
# -----------------------------------------------
# Generate the matrix
row_names <- liquido[, ncol(liquido)]

# Remove the strains column from the dataframe
# z <- zdata[, c(1:41, 43:56)]
z <- liquido[, c(1:41)]

# Convert the dataframe to a matrix and set row names
mat <- as.matrix(z)
rownames(mat) <- liquido$strains

# Transpose the matrix
mat <- t(mat)

# Create a vector with the strains (columns) to exclude (non-growers, non-Yarrowia)
excluir_columnas <- c("28", "29", "30", "102", "124", "9", "79", "115", "116", "44", "46", "48", "128")

# Get the column names we want to keep
columnas_a_mantener <- setdiff(colnames(mat), excluir_columnas)

# Subset the matrix based on the column names to keep
sub <- mat[, columnas_a_mantener]
mat <- sub
# We now have the working matrix

# We can apply another exclusion for media where none grow or there are no differences,
# as these only add noise
# Create a vector with the row names to exclude
excluir_filas <- c("Furandicarboxylic Acid", "Melibiose", "Maltose", "Methanol", "Galactose", "Piomelanin 240 h", "Ribose")

# Get the row names we want to keep
filas_a_mantener <- setdiff(rownames(mat), excluir_filas)

# Subset the matrix based on the row names to keep
sub <- mat[filas_a_mantener, ]
mat <- sub

#Conditions vectors
# Biomass <- c("Cinnamic Acid", "Ferulic Acid", "Furfural", "Hydroxymethylfurfural", 
#              "Vanillin")
# Inh <- c( "Cadmium", "Cobalt", "Copper", "Cycloheximide",
#           "Ethylene Glycol", "Nickel", "SDS", "Zinc")
# IndusEnv <- c("4% NaCl","13ºC", "37ºC", "pH 3", "pH7", "CN Intermediate", "CN High" )
# phobic <- c("Acetate", "Citrate", "Ethanol", "Glycerol", "Propionate",
#             "Hexadecane", "Lactate", "Tributyrin")
# lipid <- c("lipase CN High", "lipase CN Intermediate", "lipase Glucose",
#            "lipase Glycerol", "lipase pH7", "lipase Urea", "lipid CN High 45 h",
#            "lipid CN High 63h", "lipid Glucose")
# carbon <- c("Fructose", "Glucose", "Mannose", "Ribose", "YLD")
# nitro <- c("Glutamine", "Lysine", "Nitrates", "Urea")
# other <- c("Piomelanin 48 h", "Piomelanin 65 h", "Piomelanin 72 h", "Proteolytic activity")

# -----------------------------------------------
# 4. Prepare metadata
# -----------------------------------------------

#Metadata
popgen <- read_excel("update_manfinetreeposition.xlsx")

#Phylogenetic order
trorder <- read.table("paperlabels.txt", 
                      header = FALSE, 
                      sep = "", 
                      stringsAsFactors = FALSE, 
                      col.names = "order")
trorder$order <- as.factor(trorder$order)
order <- trorder$order
# Subset matrix according to presence in genetic analysis

popgensub <- subset(popgen, popgen$treesamp %in% order)

popgensub <- popgensub[order(match(popgensub$treesamp, trorder$order)), ]

popgensub <- subset(popgensub, popgensub$number %in% colnames(mat))

submat <- mat[, colnames(mat) %in% popgensub$number] ##Subsetting

#Final matrix has 105 strains.

#Order matrix in tree order:
common_cols <- intersect(popgensub$number, colnames(submat))

submat<- submat[, common_cols]

#Matrices of work
inhmat <- submat[rownames(submat) %in% mediosinh, ]
nutmat <- submat[rownames(submat) %in% nutrient,]

#Changing names by collection names:
#COllection names
current_names <- colnames(submat)

# REplace
new_names <- popgensub$Name[match(current_names, popgensub$number)]

# New names
colnames(submat) <- new_names

#nutmat changing
#COllection names
current_names <- colnames(nutmat)

# Replace
new_names <- popgensub$Name[match(current_names, popgensub$number)]

# Newnames
colnames(nutmat) <- new_names

#inhmat
#COllection names
current_names <- colnames(inhmat)

# Replace
new_names <- popgensub$Name[match(current_names, popgensub$number)]

# New names
colnames(inhmat) <- new_names

# -----------------------------------------------
# 5. Solid culture assays
# -----------------------------------------------

solid <- read_excel("F:/Sergio I/solido/solidphenotypying data.xlsx")

#Subsetting strains
solid <- subset(solid, solid$Strain %in% popgensub$number)
rownames(solid) <- solid$Strain

solid <- solid[order(match(solid$Strain, popgensub$number)),]
rownames(solid) <- solid$Strain

#COllection names
current_names <- rownames(solid)

# Reemplazar los nombres según coincidencia con popgensub$number
new_names <- popgensub$Name[match(current_names, popgensub$number)]

# Aplicar los nuevos nombres a nutmat
rownames(solid) <- new_names

# -----------------------------------------------
# 6. Lipase activity
# -----------------------------------------------

# Lipase activity data
lip <- data[,271:276]
rownames(lip) <- data$strain...1
colnames(lip) <- c("CNIntermediate", "Glucose", "Glycerol", "Urea", "pH 7", "CN High")

# Add strain ID as a numeric column
lip$strain <- rownames(lip)
lip$strain <- as.numeric(lip$strain)

# Filter out control strains or strains not of interest (ID >= 131)
lipase <- subset(lip, lip$strain < 131) 
rownames(lipase) <- lipase$strain

# Subset to include only strains of interest (those present in popgensub$number)
lippop <- subset(lipase, lipase$strain %in% popgensub$number)

# Order the rows of lippop to match the order in popgensub$number
lippop <- lippop[order(match(lippop$strain, popgensub$number)),]

# Calculate units of lipase activity
# Define constants
epsilon <- 18000  # Molar extinction coefficient in L·mol^-1·cm^-1
L <- 0.1          # Path length of the well in cm
V <- 0.0002       # Reaction volume in liters (e.g., 200 µL = 0.0002 L)

# Apply formula to calculate enzyme activity units in mol/min
lippop_unidades <- lippop[,1:6] * V / (epsilon * L)
rownames(lippop_unidades) <- lippop$strain

# Convert activity units from mol/min to µmol/min
lippop_unidades <- lippop_unidades * 1e6 * 60

# Check resulting values
print(lippop_unidades)

# Apply z-score normalization across all conditions
lippop_unidades <- as.data.frame(scale(lippop_unidades))

# Collect current row names
current_names <- rownames(lippop_unidades)

# Replace numeric strain IDs with actual strain names from popgensub
new_names <- popgensub$Name[match(current_names, popgensub$number)]

# Apply new row names to the matrix
rownames(lippop_unidades) <- new_names

# -----------------------------------------------
# 7. Final datasets
# -----------------------------------------------

#Lipase activity
colnames(lippop_unidades) <- paste0("lipase_", colnames(lippop_unidades))

#Carbon and nitrogen sources
nut <- t(nutmat) 

#Inhibitors and environmental stres
inhib <- t(inhmat)

#Solid data
sol <- solid[,3:31]

#Piomelanin production
piomel <- data[,266:270]
rownames(piomel) <- data$strain...1
piomel$strain <- rownames(piomel)
piomel <- subset(piomel, piomel$strain %in% popgensub$number)
piomel$strain <- as.numeric(piomel$strain)
piomel <- piomel[order(match(piomel$strain, popgensub$number)),]

#Solid data
sol <- cbind(sol, piomel[,1:5])#Nuestro dataframe de datos discretos
rownames(sol) <- popgensub$Name



#Jump to next pipeline without restarting the R session, objects in the environments will be used