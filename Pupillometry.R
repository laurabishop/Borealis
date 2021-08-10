### Analysis of pupillometry data for Borealis String Quartet Experiment
### Laura Bishop, University of Oslo, email: l.e.bishop@imv.uio.no
### June 29, 2021

# Load libraries ---------------------------------------------------------------------
library(glmmTMB) # Needed for running linear mixed effects models

# Set filepaths ----------------------------------------------------------------------
path.smi <- "Data_to_publish/Pupillometry/Performers/"
path.motion <- "Data_to_publish/Motion/"
path.rms <- "Data_to_publish/Musical_features/RMS/"
path.tonal <- "Data_to_publish/Musical_features/Cloud_diameter/"
path.difficulty <- "Data_to_publish/Musical_features/Difficulty_ratings/"

# Functions --------------------------------------------------------------------------
bar.labelling <- function(a) {
  b <- ratings[which(ratings$onset <= a & ratings$offset >= a), c("bar")]
  b[1]
}

load_rms_short <- function(trial) {
  rms.frame <- read.csv(paste(path.rms, "trial", trial, "-rms.csv", sep = ""), header = F)
  colnames(rms.frame) <- c("timestamp", "rms", "smrms")
  rms.frame
}

label_rms <- function(rmsdat, ratingsdat) {
  bar.labelling.aud <- function(a) {
    b <- ratingsdat[which(ratingsdat$onset <= a & ratingsdat$offset >= a), c("bar")]
    b[1]
  }
  
  bar.times <- approx(rmsdat$timestamp, n = max(ratingsdat$bar))
  ratingsdat["onset"] <- bar.times$y
  ratingsdat["offset"] <- c(bar.times$y[-1], max(bar.times$y))
  lab.var <- apply(as.matrix(rmsdat$timestamp), 1, bar.labelling.aud)
  rmsdat["bars"] <- unlist(lab.var)
  rmsdat
}

# Assign bar numbers to pupil data, motion data & RMS data --------------------------
## Pupil data
pupil.list <- data.frame(list.files(path = path.smi, pattern = ".txt"))
difficulty <- read.table(paste(path.difficulty, "difficulty_ratings.tsv", sep = ""), header = T) # Load difficulty ratings 
alltrials <- data.frame(matrix(rep(0, 16), nrow = 1)) # Create (empty) master data frame

for (i in pupil.list[,1]) {
  
  ### Set ins = performer's instrument name & inum = performance number for this datafile
  iname <- matrix(unlist(strsplit(i, "-")), nrow = 1)
  ins <- iname[1]
  inum <- gsub("trial", "", iname[2])
  inum <- as.numeric(gsub(".txt", "", inum))
  
  ### Load pupil data & loop through movements (NB. Trial 6 (the concert performance)
  ### has 4 movements; all other trials have 1.)
  trial.full <- read.table(paste(path.smi, i, sep = "/"))
  
  for (m in unique(trial.full$piece)) {
    
    trial.sh <- trial.full[trial.full$piece == m, ]
    
    #### Choose the correct subset of difficulty ratings for this performer and piece.
    if (inum %in% 1:5) {
      ratings <- difficulty[difficulty$ins == ins & difficulty$piece == "haydn1" & 
                              difficulty$bar < 69, ]
    } else {
      if (m == "Haydn1") {
        ratings <- difficulty[difficulty$ins == ins & difficulty$piece == "haydn1", ]
      } else if (m == "Haydn2") {
        ratings <- difficulty[difficulty$ins == ins & difficulty$piece == "haydn2", ]
      } else if (m == "Debussy1") {
        ratings <- difficulty[difficulty$ins == ins & difficulty$piece == "debussy1", ]
      } else {
        ratings <- difficulty[difficulty$ins == ins & difficulty$piece == "debussy2", ]
      }
    }
    
    #### Interpolate bar onset times using a linear interpolation of timestamps from the pupil
    #### data. Add bar onset and offset times to the ratings dataframe, the run the 
    #### "bar.labelling" function (see above under Functions) to add bar number labels to the 
    #### pupil data. 
    bar.times <- approx(trial.sh$newtime, n = max(ratings$bar))
    ratings["onset"] <- bar.times$y
    ratings["offset"] <- c(bar.times$y[-1], max(bar.times$y))
    lab.var <- apply(as.matrix(trial.sh$newtime), 1, bar.labelling)
    trial.full[trial.full$piece == m, c("bars")] <- unlist(lab.var)
  }
  
  ### Add current trial to the parent dataframe "alltrials", which will contain data from
  ### all trials once the loop is run.
  colnames(alltrials) <- colnames(trial.full)
  alltrials <- rbind(alltrials, trial.full)
}
alltrials <- alltrials[-1,]
alltrials$condition <- with(alltrials, ifelse(condition == 6 & piece == "Debussy1", 8,
                                              ifelse(condition == 6 & piece == "Debussy2", 9,
                                                     ifelse(condition == 6 & piece == "Haydn2", 7, 
                                                            condition))))

## RMS data: Data is loaded using the "load_rms_short" function (see above), and bar 
## numbers are assigned using the "label_rms" function (see also above). 
trms1 <- load_rms_short(1)
trms2 <- load_rms_short(2)
trms3 <- load_rms_short(3)
trms4 <- load_rms_short(4)
trms5 <- load_rms_short(5)
trms6 <- load_rms_short("6full")
trms7 <- load_rms_short(7)
trms8 <- load_rms_short(8)
trms9 <- load_rms_short(9)

ratings.h1 <- difficulty[difficulty$ins == "Violin1" & difficulty$piece == "haydn1" & 
                           difficulty$bar < 69, ]
brms1 <- label_rms(trms1, ratings.h1)
brms2 <- label_rms(trms2, ratings.h1)
brms3 <- label_rms(trms3, ratings.h1)
brms4 <- label_rms(trms4, ratings.h1)
brms5 <- label_rms(trms5, ratings.h1)

ratings.h1f <- difficulty[difficulty$ins == "Violin1" & difficulty$piece == "haydn1", ]
brms6 <- label_rms(trms6, ratings.h1f)

ratings.h2 <- difficulty[difficulty$ins == "Violin1" & difficulty$piece == "haydn2", ]
brms7 <- label_rms(trms7, ratings.h2)

ratings.d1 <- difficulty[difficulty$ins == "Violin1" & difficulty$piece == "debussy1", ]
brms8 <- label_rms(trms8, ratings.d1)

ratings.d2 <- difficulty[difficulty$ins == "Violin1" & difficulty$piece == "debussy2", ]
brms9 <- label_rms(trms9, ratings.d2)

## Motion data: Head and arm motion data are loaded separately. Interpolate bar onset times 
## for the ratings dataframe and use those interpolated values to label bars in motion data.
heads <- read.table("headmotion-processed.txt") # Load motion data
arms <- read.table("armmotion-processed.txt") # Load motion data

allheads <- data.frame(matrix(rep(NA, 8), nrow = 1)) # Create (empty) data frame 
colnames(allheads) <- c(colnames(heads), "bars")
allarms <- data.frame(matrix(rep(NA, 5), nrow = 1))
for (cc in unique(heads$condition)) {
  for (ii in unique(heads[heads$condition == cc, c("ins")])) {
    
    cdata <- heads[heads$condition == cc & heads$ins == ii, ]
    adata <- arms[arms$condition == cc & arms$ins == ii, ]
    
    if (cc %in% 1:5) {
      ratings <- difficulty[difficulty$ins == ii & difficulty$piece == "haydn1" & 
                              difficulty$bar < 69, ]
    } else {
      if (cc == 6) {
        ratings <- difficulty[difficulty$ins == ii & difficulty$piece == "haydn1", ]
      } else if (cc == 7) {
        ratings <- difficulty[difficulty$ins == ii & difficulty$piece == "haydn2", ]
      } else if (cc == 8) {
        ratings <- difficulty[difficulty$ins == ii & difficulty$piece == "debussy1", ]
      } else {
        ratings <- difficulty[difficulty$ins == ii & difficulty$piece == "debussy2", ]
      }
    }
    
    bar.times <- approx(cdata$timestamp, n = max(ratings$bar))
    ratings["onset"] <- bar.times$y
    ratings["offset"] <- c(bar.times$y[-1], max(bar.times$y))
    lab.var <- apply(as.matrix(cdata$timestamp), 1, bar.labelling)
    cdata["bars"] <- unlist(lab.var)
    allheads <- rbind(allheads, cdata)
    
    if (nrow(adata) != 0) {
      adata.ag <- aggregate(nvel3d ~ condition + timestamp + ins, adata, sum,
                            na.action = na.pass, na.rm = T)
      adata.ag$nvel3d <- with(adata.ag, ifelse(nvel3d == 0, NA, adata.ag$nvel3d))
      bar.times <- approx(adata.ag$timestamp, n = max(ratings$bar))
      ratings["onset"] <- bar.times$y
      ratings["offset"] <- c(bar.times$y[-1], max(bar.times$y))
      lab.var <- apply(as.matrix(adata.ag$timestamp), 1, bar.labelling)
      adata.ag["bars"] <- unlist(lab.var)
      colnames(allarms) <- colnames(adata.ag)
      allarms <- rbind(allarms, adata.ag)
    }
  }
}
allheads <- allheads[-1,]
allarms <- allarms[-1,]

# Aggregate data per bar -----------------------------------------------------------
## Using the bar-labelled data frames from the previous section, aggregate (average) 
## pupil data, (sum) motion data, and (average) RMS data per bar. 
pupils <- aggregate(nxmm ~ condition + ins + bars, alltrials, mean, 
                    na.action = na.pass, na.rm = T)
pupils$ins <- gsub("vio", "Vio", pupils$ins) # Fix discrepancy in instrument labelling
pupils$ins <- gsub("ce", "Ce", pupils$ins)
heads <- aggregate(nvel3d ~ condition + ins + bars, allheads, sum,
                   na.action = na.pass, na.rm = T) 
armmarks <- aggregate(nvel3d ~ condition + ins + bars, allarms, sum, 
                      na.action = na.pass, na.rm = T)
rms1 <- aggregate(rms ~ bars, brms1, mean)
rms2 <- aggregate(rms ~ bars, brms2, mean)
rms3 <- aggregate(rms ~ bars, brms3, mean)
rms4 <- aggregate(rms ~ bars, brms4, mean)
rms5 <- aggregate(rms ~ bars, brms5, mean)
rms6 <- aggregate(rms ~ bars, brms6, mean)
rms7 <- aggregate(rms ~ bars, brms7, mean)
rms8 <- aggregate(rms ~ bars, brms8, mean)
rms9 <- aggregate(rms ~ bars, brms9, mean)

# Combine all aggregated data into a single data frame -----------------------------
## In this section, data for all predictors are combined into a single data frame for
## use in the mixed effects models (see below). The bar-labelled pupil, motion, and
## RMS data are used as well as the difficulty ratings and tonal tension data. 
difficulty$ins <- gsub("vi", "Vi", difficulty$ins) # Fix instrument labels
difficulty$ins <- gsub("ce", "Ce", difficulty$ins)
fullmaster <- data.frame(matrix(rep(NA, 11), nrow = 1))
colnames(fullmaster) <- c("condition", "ins", "bars", 
                          "pupil", "head", "arms", "rms",
                          "technical", "express", "harmonic",
                          "cloud.diameter") 
tonal.tension.h1 <- read.csv(paste(path.tonal, "haydn1.csv", sep = ""))
tonal.tension.h2 <- read.csv(paste(path.tonal, "haydn2.csv", sep = ""))
tonal.tension.d1 <- read.csv(paste(path.tonal, "debussy1.csv", sep = ""))
tonal.tension.d2 <- read.csv(paste(path.tonal, "debussy2.csv", sep = ""))
correction <- 6.447667028833598

## Here we loop through the performances (9 total since each movement counts now as its
## own performance), and loop through the performers within performances. For each 
## combination of performer/performance, take all data into separate columns, and add
## the dataframe to "fullmaster", which will contain all data for the experiment. 
for (cond in unique(pupils$condition)) {
  for (iin in unique(pupils[pupils$condition == cond, c("ins")])) {
    p.dat <- pupils[pupils$condition == cond & pupils$ins == iin, ]
    h.dat <- heads[heads$condition == cond & heads$ins == iin, ]
    a.dat <- armmarks[armmarks$condition == cond & armmarks$ins == iin, ]
    
    if (cond %in% 1:5) {
      if (cond == 1) {
        m.dat <- rms1
      } else if (cond == 2) {
        m.dat <- rms2
      } else if (cond == 3) {
        m.dat <- rms3
      } else if (cond == 4) {
        m.dat <- rms4
      } else {
        m.dat <- rms5
      }
      ratings.h1 <- difficulty[difficulty$ins == iin & difficulty$piece == "haydn1", ]
      difr <- ratings.h1[ratings.h1$bar < 68, c("technical", "express", "harmonic")]
      tonal.tension <- tonal.tension.h1[1:67, ]
    } else if (cond == 6) {
      m.dat <- rms6
      ratings.h1 <- difficulty[difficulty$ins == iin & difficulty$piece == "haydn1", ]
      difr <- ratings.h1[ratings.h1$bar < 188, c("technical", "express", "harmonic")]
      tonal.tension <- tonal.tension.h1[1:187, ]
    } else if (cond == 7) {
      m.dat <- rms7
      ratings.h2 <- difficulty[difficulty$ins == iin & difficulty$piece == "haydn2", ]
      difr <- ratings.h2[ratings.h2$bar < 74, c("technical", "express", "harmonic")]
      tonal.tension <- tonal.tension.h2[1:73, ]
    } else if (cond == 8) {
      m.dat <- rms8
      ratings.d1 <- difficulty[difficulty$ins == iin & difficulty$piece == "debussy1", ]
      difr <- ratings.d1[ratings.d1$bar < 194, c("technical", "express", "harmonic")]
      tonal.tension <- tonal.tension.d1[1:193, ]
    } else {
      m.dat <- rms9
      ratings.d2 <- difficulty[difficulty$ins == iin & difficulty$piece == "debussy2", ]
      difr <- ratings.d2[ratings.d2$bar < 177, c("technical", "express", "harmonic")]
      tonal.tension <- tonal.tension.d2[1:176, ]
    }
    
    to.add <- cbind(p.dat, h.dat$nvel3d, a.dat$nvel3d, m.dat$rms, difr,
                    tonal.tension[, c("cloud_diameter")])
    colnames(to.add) <- colnames(fullmaster)
    fullmaster <- rbind(fullmaster, to.add)
  }
}
fullmaster <- fullmaster[-1, ]
fullmaster$arms <- with(fullmaster, ifelse(arms == 0, NA, fullmaster$arms))
fullmaster$ins <- factor(fullmaster$ins, levels(factor(fullmaster$ins))[c(3, 4, 2, 1)])
fullmaster$cloud.diameter <- fullmaster$cloud.diameter * correction
write.table(fullmaster, "Performer_data_complete.txt")

# Linear mixed effects models -------------------------------------------------------
## This section uses the dataframe 'fullmaster' that was created above. It is important
## to set condition (i.e., performance number), ins (i.e., performer instrument),
## and bars (musical time) to factors before running the models. The model that includes
## data for performers' arms (and not heads) is called 'arm.model', and the model that
## includes data for performers' heads (and not arms) is called 'head.model'. 

## After running each model, we do a hierarchical analysis to test how each predictor 
## improves the model. Note that tests against g.null must be run on a dataset in which
## NA values are omitted (otherwise, different data will be used at each of g.1, g.2, etc.).
## To solve this issue, we use dataframes called 'armsonly' and 'headsonly' that exclude
## NAs. The models arm.model and head.model could also be run on fullmaster (you should
## get the same results).

## Model containing arm data
fullmaster$condition <- as.factor(as.character(fullmaster$condition))
fullmaster$ins <- as.factor(as.character(fullmaster$ins))
fullmaster$bars <- as.factor(as.character(fullmaster$bars))

armsonly <- fullmaster[!is.na(fullmaster$arms) & 
                         !is.na(fullmaster$technical) & 
                         !is.na(fullmaster$harmonic), ]
armsonly$bars <- as.factor(as.character(armsonly$bars))
arm.model <- glmmTMB(pupil ~ arms + rms + technical + harmonic + cloud.diameter + express
                  + (1|condition) + (1|ins) + ar1(bars + 0|condition:ins),
                  data = armsonly,
                  na.action = na.exclude) 
summary(arm.model)

g.null <- glmmTMB(pupil ~ 1
                  + (1|condition) + (1|ins) + ar1(bars + 0|condition:ins),
                  data = armsonly,
                  na.action = na.exclude) 
g.1 <- glmmTMB(pupil ~ technical
               + (1|condition) + (1|ins) + ar1(bars + 0|condition:ins),
               data = armsonly,
               na.action = na.exclude) 
g.2 <- glmmTMB(pupil ~ technical + harmonic
               + (1|condition) + (1|ins) + ar1(bars + 0|condition:ins),
               data = armsonly,
               na.action = na.exclude) 
g.3 <- glmmTMB(pupil ~ technical + harmonic + cloud.diameter
               + (1|condition) + (1|ins) + ar1(bars + 0|condition:ins),
               data = armsonly,
               na.action = na.exclude) 
g.4 <- glmmTMB(pupil ~ technical + harmonic + cloud.diameter + arms
               + (1|condition) + (1|ins) + ar1(bars + 0|condition:ins),
               data = armsonly,
               na.action = na.exclude) 
anova(g.null, g.1, g.2, g.3, g.4)

## Model containing head data
headsonly <- fullmaster[!is.na(fullmaster$technical) & !is.na(fullmaster$harmonic), ]
head.model <- glmmTMB(pupil ~ head + technical + rms + harmonic + cloud.diameter + express
                  + (1|condition) + (1|ins) + ar1(bars + 0|condition:ins),
                  data = headsonly,
                  na.action = na.exclude)
summary(head.model)

effect0 <- glmmTMB(pupil ~  1
                   + (1|condition) + (1|ins) + ar1(bars + 0|condition:ins), data = headsonly, 
                   na.action = na.exclude)
effect1 <- glmmTMB(pupil ~  technical
                   + (1|condition) + (1|ins) + ar1(bars + 0|condition:ins), data = headsonly, 
                   na.action = na.exclude)
effect2 <- glmmTMB(pupil ~  technical + cloud.diameter
                   + (1|condition) + (1|ins) + ar1(bars + 0|condition:ins), data = headsonly, 
                   na.action = na.exclude)
effect3 <- glmmTMB(pupil ~  technical + cloud.diameter + harmonic
                   + (1|condition) + (1|ins) + ar1(bars + 0|condition:ins), data = headsonly, 
                   na.action = na.exclude)
effect4 <- glmmTMB(pupil ~  technical + cloud.diameter + harmonic + express
                   + (1|condition) + (1|ins) + ar1(bars + 0|condition:ins), data = headsonly, 
                   na.action = na.exclude)
anova(effect0, effect1, effect2, effect3, effect4)

# Test for effects of musical piece and concert time -----------------------------------
## Here, we run a mixed effects model for each performer to test for effects of 
## condition (i.e., performed piece). Conditions are coded as follows: 6 = Haydn1,
## 7 = Haydn2, 8 = Debussy1, 9 = Debussy2. Condition must be factorized before
## running the models. We run three models for each performer in order to test all 
## pairwise combinations of performances. With models II and III, we use the 'relevel'
## function to change the base condition. 
concert <- fullmaster[fullmaster$condition %in% c("6", "7", "8", "9"), ]
concert$condition <- as.factor(as.character(as.numeric(concert$condition)))

v1sub <- concert[concert$ins == "Violin1", ] # Change to test other performers
v1.model <- glmmTMB(pupil ~ condition + (1|condition) + ar1(bars + 0|condition),
                    data = v1sub,
                    na.action = na.exclude) 
summary(v1.model)

v1sub.II  <- within(v1sub, condition <- relevel(condition, ref = 2))
v1.model.II <- glmmTMB(pupil ~ condition + (1|condition) + ar1(bars + 0|condition),
                       data = v1sub.II,
                       na.action = na.exclude) 
summary(v1.model.II)

v1sub.III  <- within(v1sub, condition <- relevel(condition, ref = 3))
v1.model.III <- glmmTMB(pupil ~ condition + (1|condition) + ar1(bars + 0|condition),
                        data = v1sub.III,
                        na.action = na.exclude) 
summary(v1.model.III)