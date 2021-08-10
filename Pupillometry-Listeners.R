### Borealis String Quartet Experiment: Analysis of Listeners' Pupillometry Data
### Laura Bishop, University of Oslo, email: l.e.bishop@imv.uio.no
### June 29, 2021

# Load libraries ---------------------------------------------------------------------
library(glmmTMB) # Needed for running linear mixed effects models

# Set filepaths ----------------------------------------------------------------------
path.smi <- "Data_to_publish/Pupillometry/Listeners/" # Folder containing eye tracking data
path.difficulty <- "Data_to_publish/Musical_features/Difficulty_ratings/"

# Functions --------------------------------------------------------------------------
bar.labelling <- function(a) {
  b <- ratings[which(ratings$onset <= a & ratings$offset >= a), c("bar")]
  b[1]
}

# Label bars and aggregate data per bar ----------------------------------------------
pupil.list <- data.frame(list.files(path = path.smi))
difficulty <- read.table(paste(path.difficulty, "difficulty_ratings.tsv", sep = ""), header = T)
vars <- read.table("Performer_data_complete.txt")
all.pupil <- data.frame(matrix(rep(0, 11), nrow = 1))

for (i in pupil.list[,1]) {

  ### Set inum = piece name for this file
  iname <- matrix(unlist(strsplit(i, "-")), nrow = 1)
  inum <- gsub(".txt", "", iname[2])
  
  ### Interpolate bar onset times using the bar numbers from "difficulty" and pupil timestamps.
  ### Then use the bar.labelling function to assign bar numbers to pupil data.
  pupil.data <- read.table(paste(path.smi, i, sep = ""))
  ratings <- difficulty[difficulty$piece == inum, ]
  bar.times <- approx(pupil.data$newtime, n = max(ratings$bar)) # Approx bar onsets
  ratings["onset"] <- bar.times$y # Add bar onset times to rating file
  ratings["offset"] <- c(bar.times$y[-1], max(bar.times$y)) # Add bar offset times
  lab.var <- apply(as.matrix(pupil.data$newtime), 1, bar.labelling) # Label bars in filter.dat
  pupil.data["bars"] <- unlist(lab.var) # Add bar labels to data frame
  
  # Aggregate pupil diameters per bar
  pupil.bar <- aggregate(nxmm ~ bars, pupil.data, mean)
  pupil.bar["bars"] <- unique(pupil.data$bars)
  pupil.bar["id"] <- pupil.data$id[1]
  pupil.bar["piece"] <- inum
  
  # Add RMS, cloud diameter, motion, note number/bar
  if (inum == "haydn1") {
    condnum <- 6
  } else if (inum == "haydn2") {
    condnum <- 7
  } else if (inum == "debussy1") {
    condnum <- 8
  } else {
    condnum <- 9
  }
  
  add.vars <- vars[vars$condition == condnum, ]
  pupil.bar["rms"] <- add.vars[add.vars$ins == "Violin1", c("rms")]
  pupil.bar["cloud.diameter"] <- add.vars[add.vars$ins == "Violin1", c("cloud.diameter")]
  head <- aggregate(head ~ bars, add.vars, mean)
  pupil.bar["head"] <- head$head
  arms <- aggregate(arms ~ bars, add.vars, mean, na.action = na.pass, na.rm = T)
  pupil.bar["arms"] <- arms$arms
  technical <- aggregate(technical ~ bars, add.vars, mean, na.action = na.pass, na.rm = T)
  pupil.bar["technical"] <- technical$technical
  express <- aggregate(express ~ bars, add.vars, mean, na.action = na.pass, na.rm = T)
  pupil.bar["express"] <- express$express
  harmonic <- aggregate(harmonic ~ bars, add.vars, mean, na.action = na.pass, na.rm = T)
  pupil.bar["harmonic"] <- harmonic$harmonic
  
  # Add to master frame
  colnames(all.pupil) <- colnames(pupil.bar)
  all.pupil <- rbind(all.pupil, pupil.bar)
}

# Linear mixed effects models -------------------------------------------------------
## This section uses the dataframe "all.pupil" that was created above. We run two models,
## one that includes performers' arm motion (arm.model) and one that includes performers'
## head motion (head.model). The variables bars, piece, and participant ID should be 
## factorized before running the models. It's also necessary to add a constant to the
## pupil diameter variable (nxmm). 
all.pupil$nxmm2 <- all.pupil$nxmm + 2
all.pupil$bars <- as.factor(as.character(all.pupil$bars))
all.pupil$piece <- as.factor(as.character(all.pupil$piece))
all.pupil$id <- as.factor(as.character(all.pupil$id))

## Model containing arm data
arm.model <- glmmTMB(nxmm2 ~ arms + technical + rms + harmonic + cloud.diameter + express
                   + (1|piece) + (1|id) + ar1(bars + 0|piece:id),
                   data = all.pupil, 
                   na.action = na.exclude)
summary(arm.model)

null <- glmmTMB(nxmm2 ~ 1
                + (1|piece) + (1|id) + ar1(bars + 0|piece:id),
                data = all.pupil, 
                na.action = na.exclude)
m1 <- glmmTMB(nxmm2 ~ harmonic
              + (1|piece) + (1|id) + ar1(bars + 0|piece:id),
              data = all.pupil, 
              na.action = na.exclude)
m2 <- glmmTMB(nxmm2 ~ harmonic + express
              + (1|piece) + (1|id) + ar1(bars + 0|piece:id),
              data = all.pupil, 
              na.action = na.exclude)
m3 <- glmmTMB(nxmm2 ~ harmonic + express + technical
              + (1|piece) + (1|id) + ar1(bars + 0|piece:id),
              data = all.pupil, 
              na.action = na.exclude)
m4 <- glmmTMB(nxmm2 ~ harmonic + express + technical + arms
              + (1|piece) + (1|id) + ar1(bars + 0|piece:id),
              data = all.pupil, 
              na.action = na.exclude)
anova(null, m1, m2, m3, m4)

## Model containing head data
head.model <- glmmTMB(nxmm2 ~ head + technical + rms + harmonic + cloud.diameter + express
                   + (1|piece) + (1|id) + ar1(bars + 0|piece:id),
                   data = all.pupil, 
                   na.action = na.exclude)
summary(head.model)

null <- glmmTMB(nxmm2 ~ 1
                + (1|piece) + (1|id) + ar1(bars + 0|piece:id),
                data = all.pupil, 
                na.action = na.exclude)
m1 <- glmmTMB(nxmm2 ~ harmonic
              + (1|piece) + (1|id) + ar1(bars + 0|piece:id),
              data = all.pupil, 
              na.action = na.exclude)
m2 <- glmmTMB(nxmm2 ~ harmonic + express
              + (1|piece) + (1|id) + ar1(bars + 0|piece:id),
              data = all.pupil, 
              na.action = na.exclude)
m3 <- glmmTMB(nxmm2 ~ harmonic + express + technical
              + (1|piece) + (1|id) + ar1(bars + 0|piece:id),
              data = all.pupil, 
              na.action = na.exclude)
m4 <- glmmTMB(nxmm2 ~ harmonic + express + technical + head
              + (1|piece) + (1|id) + ar1(bars + 0|piece:id),
              data = all.pupil, 
              na.action = na.exclude)
anova(null, m1, m2, m3, m4)




































