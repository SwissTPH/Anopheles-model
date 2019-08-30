library(xlsx)
library(gdata)
library(reshape)
library(dplyr)
#### DATA ####

# read rhythms from Sherrard-Smith(2019) supplementary information 
rhythms_pnas <- 'C:/git_repos/vector-models/LLINparameterisation/pnas.1820646116.sd01.xlsx'

# read rhythms from original Swiss TPH compilation 
tph_rhythms <- read.csv('C:/git_repos/Anopheles-model/downloads/biting_rhythm.csv',sep=';')
# remove rhythm_id and also Huho et al data (which is also in the Sherrard-Smith dataset)
tph_rhythms <- tph_rhythms[ ,(!(names(tph_rhythms)=='rhythm_id'))]
tph_rhythms <- tph_rhythms[tph_rhythms$source_id != 5781,]
varnames <- c("X17:00_18:00", "X18:00_19:00", "X19:00_20:00", "X20:00_21:00", "X21:00_22:00", "X22:00_23:00",
              "X23:00_00:00", "X00:00_01:00", "X01:00_02:00", "X02:00_03:00", "X03:00_04:00", "X04:00_05:00", 
              "X05:00_06:00", "X06:00_07:00") 
names(tph_rhythms) <- c('species','country','sampling','site', varnames, 'comments','source_id','pubmed_id', 'citation')
tph_rhythms$'X16:00_17:00' <- NA
tph_rhythms$'X07:00_08:00' <- NA

# specify the time points from Sherrard-Smith to be included in the analysis files
full_set_of_times <- c("X16.00", "X17.00", "X18.00", "X19.00", 
                             "X20.00", "X21.00", "X22.00", "X23.00", "X00.00", "X01.00", "X02.00", 
                             "X03.00", "X04.00", "X05.00", "X06.00", "X07.00", "X08.00", "X09.00",   
                             "X10.00", "X11.00", "X12.00", "X13.00", "X14.00", "X15.00")
required_subset <- c("X17.00", "X18.00", "X19.00", "X20.00", "X21.00", "X22.00", "X23.00", "X00.00", "X01.00", "X02.00", 
                                                       "X03.00", "X04.00", "X05.00", "X06.00", "X07.00", "X08.00")
fullvarnames <- c("X16:00_17:00",varnames,"X07:00_08:00")


readSherrard_Smith <- function(sheetno, records,sampling) {
  df0 <- t(as.data.frame(read.xls(rhythms_pnas,sheet=sheetno))) 
  colnames(df0) <- full_set_of_times
  required_times <- which(colnames(df0) %in% required_subset)
  required_rows <- seq(from=2,to=(records+1))
  df0 <- df0[required_rows,required_times]
  # writing to csv and reading in again removes unwanted type assignments 
  write.csv(df0,'C:/git_repos/vector-models/LLINparameterisation/temp.csv', row.names = FALSE)
  df <- read.csv('C:/git_repos/vector-models/LLINparameterisation/temp.csv')
  names(df)<-fullvarnames
  df$sampling <- sampling
  if(sheetno < 7){
    df_references <- as.data.frame(read.xls(rhythms_pnas,sheet=sheetno))$X.1
    write.csv(df_references,'C:/git_repos/vector-models/LLINparameterisation/temp.csv', row.names = FALSE)
    df_references0 <- read.csv('C:/git_repos/vector-models/LLINparameterisation/temp.csv')
    df$Reference2 <- as.character(df_references0[!is.na(df_references0)])
    df1 <- left_join(df, temp2, by='Reference2')
    df1$Reference1 <- "-999"
    df <- df1
  }
  return(df)
}

temp2 <- read.csv(file="C:/git_repos/Anopheles-model/raw_data/temporary_reference_information.csv",sep=";")
temp2$Reference1 <- as.character(temp2$Reference1)
temp2$Reference2 <- as.character(temp2$Reference2) 
drops <- c("sampling")
temp2 <- temp2[ , !(names(temp2) %in% drops)]

# read each set of times into dataframes and concatenate
human_indoor <- readSherrard_Smith(sheetno=5, records=22, sampling='IND')
human_bed <- readSherrard_Smith(sheetno=6, records=7, sampling='BED')
# Indoor mosquitoes
df <- readSherrard_Smith(sheetno=7, records=131, sampling='HBI')
df0 <- as.data.frame(read.xls(rhythms_pnas,sheet=7))
df$Reference1 <- sub('X', '', names(df0))[2:132]
mosquito_indoor <- left_join(df, temp2, by='Reference1')
mosquito_indoor$Reference2 <- "-999"
df <- readSherrard_Smith(sheetno=9, records=131, sampling='HBO') 
df0 <- as.data.frame(read.xls(rhythms_pnas,sheet=9))
df$Reference1 <- sub('X', '', names(df0))[2:132]
mosquito_outdoor <- left_join(df, temp2, by='Reference1')
mosquito_outdoor$Reference2 <- "-999"
PNAS_rhythms <- do.call(rbind, list(human_indoor,human_bed,mosquito_indoor,mosquito_outdoor))
drops <- c("Reference1","Reference2","Study")
PNAS_rhythms <- PNAS_rhythms[ , !(names(PNAS_rhythms) %in% drops)]
all_rhythms <- do.call(rbind, list(PNAS_rhythms,tph_rhythms))

# Extract the citations into a separate file, remove duplicates, and make the source_id the unique key
keeps <- c("source_id","pubmed_id","citation")
citations <- unique(all_rhythms[ , (names(all_rhythms) %in% keeps)])
citations <- citations[order(citations$source_id),]
# citations that do not have source_ids from MAP or from previous TPH assignments are assigned values that increment
# starting from the existing maximum value + 1
max_source_id <- max(citations$source_id,na.rm=TRUE)
total_citations <- length(citations$source_id)
valid_source_ids <- length(citations$source_id[!is.na(citations$source_id)])
missing_source_ids <- total_citations - valid_source_ids
citations$source_id <- c(citations$source_id[1:valid_source_ids],(max_source_id+seq(1:missing_source_ids)))
# Merge the source_id back into the all_rhythms data frame and normalise
drops <- c("source_id","pubmed_id")
all_rhythms <- all_rhythms[ , !(names(all_rhythms) %in% drops)]
all_rhythms <- left_join(all_rhythms, citations, by='citation')
 
write.csv(all_rhythms,'C:/git_repos/Anopheles-model/downloads/all_rhythms.csv')

