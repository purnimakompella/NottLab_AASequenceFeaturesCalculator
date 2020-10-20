#Also Calculate Fraction of Aromatic Residues
library(stringr)
#install.packages('Peptides') #https://cran.r-project.org/web/packages/Peptides/Peptides.pdf
library(Peptides)

rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

#filepath <- args[1]
#pH <- as.numeric(args[2])

VERBOSE = FALSE

filepath <- '/Users/bioc1660/Documents/Bioinformatics/filesfromTJN/CSVfiles/201014_ddx4vddx3x_fulllength'
pH <- 7.4

print(paste("File path provided: ", filepath))
print(paste("pH provided: ", pH))

files <- list.files(filepath, full.names=TRUE, pattern = "*.csv")

NetChargeTable <- list()
j=1
NetChargeColName <- paste("NetCharge_", pH, sep="")
for (file in files) {
  table <-  read.csv(file)
  df <- data.frame(table)
  if (VERBOSE) print(paste("Length of df: ", nrow(df)))
  pKa_Nterm <- 7.8
  pKa_Cterm <- 3.6
  pKa_Asp <- 4 #D
  pKa_Arg <- 12.5 #R
  pKa_Cys <- 8.14 #C
  pKa_Glu <- 4.5 #E
  pKa_His <- 6.4 #H
  pKa_Tyr <- 9.6 #Y
  pKa_Lys <- 10.4 #K
  
  ASP_pH <- -1*(10**(-(pKa_Asp-pH)))/((10**(-(pKa_Asp-pH)))+1) #D (-1)
  GLU_pH <- (-1)*(10**(-(pKa_Glu-pH)))/((10**(-(pKa_Glu-pH)))+1) #E (-1)
  LYS_pH <- (1)*(10**((pKa_Lys-pH)))/((10**((pKa_Lys-pH)))+1) #K (+1)
  ARG_pH <- (1)*(10**((pKa_Arg-pH)))/((10**((pKa_Arg-pH)))+1) #R (+1)
  HIS_pH <- (1)*(10**((pKa_His-pH)))/((10**((pKa_His-pH)))+1) #H
  CYS_pH <- (-1)*(10**(-(pKa_Cys-pH)))/((10**(-(pKa_Cys-pH)))+1) #C
  TYR_pH <- (-1)*(10**(-(pKa_Tyr-pH)))/((10**(-(pKa_Tyr-pH)))+1) #Y
  Nterm_pH <- (1)*(10**((pKa_Nterm-pH)))/((10**((pKa_Nterm-pH)))+1)
  Cterm_pH <- (-1)*(10**(-(pKa_Cterm-pH)))/((10**(-(pKa_Cterm-pH)))+1)
  
  ASP_simplecharge <- -1 #D
  GLU_simplecharge <- -1 #E
  LYS_simplecharge <- +1 #K
  ARG_simplecharge <- +1 #R
  
  numberchargedresidues = 0
  numbersimplechargedresidues = 0
  numberaromaticresidues = 0
  numA = 0
  numC = 0
  numD = 0
  numE = 0
  numF = 0
  numG = 0
  numH = 0
  numI = 0
  numK = 0
  numL = 0
  numM = 0
  numN = 0
  numP = 0
  numQ = 0
  numR = 0
  numS = 0
  numT = 0
  numV = 0
  numW = 0
  numY = 0
  sequence = ""
  #print(paste("Sequence type:", typeof(sequence)))
  for (i in 1:nrow(df)) 
  {
    #print(paste("Residue", df$Residue[i]))
    #print(paste("Typeof", (typeof(as.character(df$Residue[i])))))
    #sequence = c(sequence, as.character(df$Residue[i]), collapse="")
    sequence = paste(sequence, (as.character(df$Residue[i])), sep="")
    if (df$Residue[i]=='D')
    {
      df$charge[i] <- ASP_pH
      df$simplecharge[i] <- ASP_simplecharge
      numberchargedresidues = numberchargedresidues+1
      numbersimplechargedresidues = numbersimplechargedresidues + 1
      numD = numD + 1
    }
    else if (df$Residue[i]=='R')
    {
      df$charge[i] <- ARG_pH
      df$simplecharge[i] <- ARG_simplecharge
      numberchargedresidues = numberchargedresidues+1
      numbersimplechargedresidues = numbersimplechargedresidues + 1
      numR = numR + 1
    }
    else if (df$Residue[i]=='C')
    {
      df$charge[i] <- CYS_pH
      df$simplecharge[i] = 0
      numberchargedresidues = numberchargedresidues+1
      numC = numC + 1
    }
    else if (df$Residue[i]=='E')
    {
      df$charge[i] <- GLU_pH
      df$simplecharge[i] <- GLU_simplecharge
      numberchargedresidues = numberchargedresidues+1
      numbersimplechargedresidues = numbersimplechargedresidues + 1
      numE = numE + 1
    }
    else if (df$Residue[i]=='H')
    {
      df$charge[i] <- HIS_pH
      df$simplecharge[i] = 0
      numberchargedresidues = numberchargedresidues+1
      numH = numH + 1
    }
    else if (df$Residue[i]=='Y')
    {
      df$charge[i] <- TYR_pH
      df$simplecharge[i] = 0
      numberchargedresidues = numberchargedresidues+1
      numberaromaticresidues = numberaromaticresidues + 1
      numY = numY + 1
    }
    else if (df$Residue[i]=='K')
    {
      df$charge[i] <- LYS_pH
      df$simplecharge[i] <- LYS_simplecharge
      numberchargedresidues = numberchargedresidues+1
      numbersimplechargedresidues = numbersimplechargedresidues + 1
      numK = numK + 1
    }
    else if (df$Residue[i]=='F')
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
      numberaromaticresidues = numberaromaticresidues + 1
      numF = numF + 1
    }
    else if (df$Residue[i]=='W')
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
      numberaromaticresidues = numberaromaticresidues + 1
      numW = numW + 1
    }
    else if (df$Residue[i]=='A')
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
      numA = numA + 1
    }
    else if (df$Residue[i]=='G')
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
      numG = numG + 1
    }
    else if (df$Residue[i]=='I')
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
      numI = numI + 1
    }
    else if (df$Residue[i]=='L')
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
      numL = numL + 1
    }
    else if (df$Residue[i]=='M')
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
      numM = numM + 1
    }
    else if (df$Residue[i]=='N')
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
      numN = numN + 1
    }
    else if (df$Residue[i]=='P')
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
      numP = numP + 1
    }
    else if (df$Residue[i]=='Q')
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
      numQ = numQ + 1
    }
    else if (df$Residue[i]=='S')
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
      numS = numS + 1
    }
    else if (df$Residue[i]=='T')
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
      numT = numT + 1
    }
    else if (df$Residue[i]=='V')
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
      numV = numV + 1
    }
    else
    {
      df$charge[i] = 0
      df$simplecharge[i] = 0
    }
  }
  #print("Sequence:", sequence)
  # print(paste("Outside 2nd for numChargedRes, numAromaticRes: ", numberchargedresidues, numberaromaticresidues))
  currentfile = str_split(basename(file), ".csv")
  proteinname = unlist(currentfile)[1]
  split_conditions =  str_split(unlist(currentfile)[2], "_")
  regionnum = paste(unlist(split_conditions)[6], unlist(split_conditions)[7], sep="_")
  protein_region_name = paste(proteinname, regionnum, sep="_")
  #clean_proteinanme <- gsub(pattern = ".csv", replacement="", proteinname, perl=T)
  
  NetChargeTable$protein[j] = proteinname
  NetChargeTable$protein_region[j] = protein_region_name
  NetChargeTable$Region[j] = gsub(pattern = ".csv", replacement="", unlist(split_conditions)[7], perl=T) 
  NumberofResidues <- nrow(df)
  NetChargeTable$NumberofResidues[j] <- NumberofResidues
  
  NetChargeTable$Sequence[j] <- sequence
  NetChargeTable$pI[j] <- pI(sequence,pKscale="Bjellqvist")
  NetChargeTable$numA[j] <- numA
  NetChargeTable$numC[j] <- numC
  NetChargeTable$numD[j] <- numD
  NetChargeTable$numE[j] <- numE
  NetChargeTable$numF[j] <- numF
  NetChargeTable$numG[j] <- numG
  NetChargeTable$numH[j] <- numH
  NetChargeTable$numI[j] <- numI
  NetChargeTable$numK[j] <- numK
  NetChargeTable$numL[j] <- numL
  NetChargeTable$numM[j] <- numM
  NetChargeTable$numN[j] <- numN
  NetChargeTable$numP[j] <- numP
  NetChargeTable$numQ[j] <- numQ
  NetChargeTable$numR[j] <- numR
  NetChargeTable$numS[j] <- numS
  NetChargeTable$numT[j] <- numT
  NetChargeTable$numV[j] <- numV
  NetChargeTable$numW[j] <- numW
  NetChargeTable$numY[j] <- numY
  
  NetChargeTable$fracA[j] <- numA/NumberofResidues
  NetChargeTable$fracC[j] <- numC/NumberofResidues
  NetChargeTable$fracD[j] <- numD/NumberofResidues
  NetChargeTable$fracE[j] <- numE/NumberofResidues
  NetChargeTable$fracF[j] <- numF/NumberofResidues
  NetChargeTable$fracG[j] <- numG/NumberofResidues
  NetChargeTable$fracH[j] <- numH/NumberofResidues
  NetChargeTable$fracI[j] <- numI/NumberofResidues
  NetChargeTable$fracK[j] <- numK/NumberofResidues
  NetChargeTable$fracL[j] <- numL/NumberofResidues
  NetChargeTable$fracM[j] <- numM/NumberofResidues
  NetChargeTable$fracN[j] <- numN/NumberofResidues
  NetChargeTable$fracP[j] <- numP/NumberofResidues
  NetChargeTable$fracQ[j] <- numQ/NumberofResidues
  NetChargeTable$fracR[j] <- numR/NumberofResidues
  NetChargeTable$fracS[j] <- numS/NumberofResidues
  NetChargeTable$fracT[j] <- numT/NumberofResidues
  NetChargeTable$fracV[j] <- numV/NumberofResidues
  NetChargeTable$fracW[j] <- numW/NumberofResidues
  NetChargeTable$fracY[j] <- numY/NumberofResidues
  
  NetCharge <- Nterm_pH + Cterm_pH + sum(df$charge)
  NetChargeTable$NetCharge[j] <- NetCharge
  NetChargeTable$ChargeDensity[j] <- numberchargedresidues/NumberofResidues
  NetChargeTable$FractionAromatic[j] <- numberaromaticresidues/NumberofResidues
  
  SimpleNetCharge <- sum(df$simplecharge)
  NetChargeTable$SimpleNetCharge[j] = SimpleNetCharge
  NetChargeTable$SimpleChargeDensity[j] <- numbersimplechargedresidues/NumberofResidues
  if (VERBOSE) print(paste("j: ", j))
  j <- j+1
  if (VERBOSE) print(paste("j after addition: ", j))
  if (VERBOSE) print(paste("Size of NC table: ", length(NetChargeTable)))
}

NetChargeTable_df <- as.data.frame(NetChargeTable)
if (VERBOSE) print(paste("Size of NC table outside 2nd for loop: ", nrow(NetChargeTable_df)))
names(NetChargeTable_df)[names(NetChargeTable_df) == "NetCharge"] <- NetChargeColName

outputdirectory = paste(dirname(file),"/Analysis/", sep="")
if (!dir.exists(outputdirectory)){
  dir.create(outputdirectory)
}

filebasename = str_split(basename(file), ".csv")
split_conditions2 =  str_split(unlist(filebasename)[2], "_")
imp_conditions = paste(unlist(split_conditions2)[2], unlist(split_conditions2)[3], unlist(split_conditions2)[4], unlist(split_conditions2)[5], sep="_")
outputfilename = paste(outputdirectory, format(Sys.Date(), format="%y%m%d"), "_NetCharge_pI_pH",pH, "_", imp_conditions, ".csv", sep="")
write.csv(NetChargeTable_df, file = outputfilename, row.names = FALSE)
