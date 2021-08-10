### Processing of body motion data for Borealis String Quartet Experiment
### Laura Bishop, University of Oslo, email: l.e.bishop@imv.uio.no
### June 29, 2021

### This script should be run to process head and arm motion data before running the 
### script "Pupillometry". This script will write two files: headmotion-processed.txt
### and armmotion-processed.txt, which the Pupillometry script will make use of. 
### This part of the analysis has been placed in a separate script in order to make
### the Pupillometry scrpt easier to follow. 

# Set filepaths ----------------------------------------------------------------------
path.motion <- "Data_to_publish/Motion/" # Folder containing motion data

# Prepare head motion data -----------------------------------------------------------
## Here, we will create a file that contains head motion data from all performers in all 
## performances. This comprehensive file will include the norms of 3D position and velocity data. 
head.list <- data.frame(list.files(path = path.motion, pattern = "heads"))
allheads <- data.frame(matrix(rep(0, 47), nrow = 1))
for (l in head.list[,1]) {
  cond.data <- read.table(paste(path.motion, l, sep = ""))
  
  for (c in c("Violin1", "Violin2", "Viola", "Cello")) {
    pos3d <- cond.data[, grepl("pos", colnames(cond.data)) & grepl(c, colnames(cond.data))]
    pos3D <- sqrt(pos3d[,1]^2 + pos3d[,2]^2 + pos3d[,3]^2)
    cond.data[paste(c, ".pos.3d", sep = "")] <- pos3D
    
    vel3d <- cond.data[, grepl("vel", colnames(cond.data)) & grepl(c, colnames(cond.data))]
    vel3D <- sqrt(vel3d[,1]^2 + vel3d[,2]^2 + vel3d[,3]^2)
    cond.data[paste(c, ".vel.3d", sep = "")] <- vel3D
  }
  
  l.trial <- gsub("condition", "", l)
  l.trial <- gsub("-heads.txt", "", l.trial)
  l.trial1 <- unlist(strsplit(as.character(l.trial), "-"))
  cond.data["condition"] <- l.trial1[1]
  cond.data["timestamp"] <- seq(0, nrow(cond.data)*1/120-1/120, 1/120)
  colnames(allheads) <- colnames(cond.data)
  allheads <- rbind(allheads, cond.data)
}
allheads <- allheads[-1,]

## Next, we reorganize the data into a long format (i.e., so that all performers are listed by
## instrument in one column), and retain only specific columns for the modelling analysis.
heads3d <- allheads[, grepl("3d", colnames(allheads)) | grepl("condition", colnames(allheads))
                    | grepl("movement", colnames(allheads)) | grepl("timestamp", colnames(allheads))]
heads3d.c <- heads3d[, grepl("Cello", colnames(heads3d)) | grepl("condition", colnames(heads3d)) 
                     | grepl("movement", colnames(heads3d)) | grepl("timestamp", colnames(heads3d))]
heads3d.v <- heads3d[, grepl("Viola", colnames(heads3d)) | grepl("condition", colnames(heads3d))
                     | grepl("movement", colnames(heads3d)) | grepl("timestamp", colnames(heads3d))]
heads3d.v1 <- heads3d[, grepl("Violin1", colnames(heads3d)) | grepl("condition", colnames(heads3d))
                      | grepl("movement", colnames(heads3d)) | grepl("timestamp", colnames(heads3d))]
heads3d.v2 <- heads3d[, grepl("Violin2", colnames(heads3d)) | grepl("condition", colnames(heads3d))
                      | grepl("movement", colnames(heads3d)) | grepl("timestamp", colnames(heads3d))]
heads3d.c["ins"] <- "Cello"
heads3d.v["ins"] <- "Viola"
heads3d.v1["ins"] <- "Violin1"
heads3d.v2["ins"] <- "Violin2"
colnames(heads3d.c) <- c("movement", "pos3d", "vel3d", "condition", "timestamp", "ins")
colnames(heads3d.v) <- c("movement", "pos3d", "vel3d", "condition", "timestamp", "ins")
colnames(heads3d.v1) <- c("movement", "pos3d", "vel3d", "condition", "timestamp", "ins")
colnames(heads3d.v2) <- c("movement", "pos3d", "vel3d", "condition", "timestamp", "ins")
allheads3d <- rbind(heads3d.c, heads3d.v, heads3d.v1, heads3d.v2)

## Conditions 6-7 (concert performances) can now be spit up into separate movements (6-9).
allheads3d$condition <- with(allheads3d, ifelse(condition == 7 & movement == 3, 8,
                                                ifelse(condition == 7 & movement == 4, 9,
                                                       ifelse(condition == 6 & movement == 2, 7, 
                                                              condition))))

## Loop through the performers and normalize velocity data. 
wdata2 <- data.frame(matrix(rep(0, 7), nrow = 1))
for (i in unique(allheads3d$ins)) {
  insdata <- allheads3d[allheads3d$ins == i, ]
  insdata["nvel3d"] <- with(insdata, (vel3d - min(vel3d, na.rm = T))/(max(vel3d, na.rm = T) - min(vel3d, na.rm = T)))
  colnames(wdata2) <- colnames(insdata)
  wdata2 <- rbind(wdata2, insdata)
}
wdata2 <- wdata2[-1,]
write.table(wdata2, "headmotion-processed.txt")

# Prepare arm motion data -------------------------------------------------------------
## As for head motion, we will create a dataframe that contains arm motion from all 
## performers in all performances. 
arm.list <- data.frame(list.files(path = path.motion, pattern = "arms", ))
allarms <- data.frame(matrix(rep(0, 11), nrow = 1))
for (l in arm.list[,1]) {
  
  ### Set l.trial1 contains movement number
  l.trial <- gsub("condition", "", l)
  l.trial <- gsub("-arms.txt", "", l.trial)
  l.trial1 <- unlist(strsplit(as.character(l.trial), "-"))
  
  ### Read in arm motion data and loop through the performers by instrument name. For
  ### each performer, calculate the norm of 3d velocity data. Note that the 2nd 
  ### violinist and violist had poorly tracked arms. In these cases, arm marker data 
  ### will be missing from the files located in path.motion. 
  cond.data <- read.table(paste(path.motion, l, sep = ""))
  vel.data <- cond.data[, grepl("vel", colnames(cond.data)) & 
                          grepl("lower", colnames(cond.data))]
  
  for (c in c("Violin1", "Violin2", "Viola", "Cello")) {
    vel.left <- vel.data[, grepl(c, colnames(vel.data)) & grepl("RArm", colnames(vel.data))]
    if (length(vel.left) != 0) {
      vel.left3D <- sqrt(vel.left[,1]^2 + vel.left[,2]^2 + vel.left[,3]^2)
      vel.data[paste("vel", c, "LArm_lower.3d", sep = "")] <- vel.left3D
    } else {
      vel.data[paste("vel", c, "LArm_lower.3d", sep = "")] <- NA
    }
  }
  
  for (c in c("Violin1", "Violin2", "Viola", "Cello")) {
    vel.right <- vel.data[, grepl(c, colnames(vel.data)) & grepl("RArm", colnames(vel.data))]
    if (length(vel.right) != 0) {
      vel.right3D <- sqrt(vel.right[,1]^2 + vel.right[,2]^2 + vel.right[,3]^2)
      vel.data[paste("vel", c, "RArm_lower.3d", sep = "")] <- vel.right3D
    } else {
      vel.data[paste("vel", c, "RArm_lower.3d", sep = "")] <- NA
    }
  }
  
  ### Take only the norm of 3D velocity data and add a timestamp and performance
  ###  number information.
  only3d <- vel.data[, grepl("3d", colnames(vel.data))]
  only3d["condition"] <- l.trial1[1]
  only3d["timestamp"] <- seq(0, nrow(only3d)*1/120-1/120, 1/120)
  only3d["movement"] <- unique(cond.data$movement)
  
  colnames(allarms) <- colnames(only3d)
  allarms <- rbind(allarms, only3d)
}
allarms <- allarms[-1,]

## Reorganize data so that it's in a form that's compatible with the pupil data, then calculate
## normalized velocities. 
velonly <- allarms[, grepl("vel", colnames(allarms)) | grepl("condition", colnames(allarms)) |
                     grepl("timestamp", colnames(allarms)) | grepl("movement", colnames(allarms))]
v1.left <- velonly[, (grepl("Violin1", colnames(velonly)) & grepl("LArm", colnames(velonly))) |
                     grepl("condition", colnames(velonly)) |
                     grepl("timestamp", colnames(velonly)) | grepl("movement", colnames(velonly))]
v2.left <- velonly[, (grepl("Violin2", colnames(velonly)) & grepl("LArm", colnames(velonly)))  |
                     grepl("condition", colnames(velonly)) |
                     grepl("timestamp", colnames(velonly)) | grepl("movement", colnames(velonly))]
v.left <- velonly[, (grepl("Viola", colnames(velonly)) & grepl("LArm", colnames(velonly))) |
                    grepl("condition", colnames(velonly)) |
                    grepl("timestamp", colnames(velonly)) | grepl("movement", colnames(velonly))]
c.left <- velonly[, (grepl("Cello", colnames(velonly)) & grepl("LArm", colnames(velonly))) |
                    grepl("condition", colnames(velonly)) |
                    grepl("timestamp", colnames(velonly)) | grepl("movement", colnames(velonly))]
v1.right <- velonly[, (grepl("Violin1", colnames(velonly)) & grepl("RArm", colnames(velonly))) |
                      grepl("condition", colnames(velonly)) |
                      grepl("timestamp", colnames(velonly)) | grepl("movement", colnames(velonly))]
v2.right <- velonly[, (grepl("Violin2", colnames(velonly)) & grepl("RArm", colnames(velonly))) |
                      grepl("condition", colnames(velonly)) |
                      grepl("timestamp", colnames(velonly)) | grepl("movement", colnames(velonly))]
v.right <- velonly[, (grepl("Viola", colnames(velonly)) & grepl("RArm", colnames(velonly))) |
                     grepl("condition", colnames(velonly)) |
                     grepl("timestamp", colnames(velonly)) | grepl("movement", colnames(velonly))]
c.right <- velonly[, (grepl("Cello", colnames(velonly)) & grepl("RArm", colnames(velonly))) |
                     grepl("condition", colnames(velonly)) |
                     grepl("timestamp", colnames(velonly)) | grepl("movement", colnames(velonly))]
v1.left["ins"] <- "Violin1"
v2.left["ins"] <- "Violin2"
v.left["ins"] <- "Viola"
c.left["ins"] <- "Cello"
v1.right["ins"] <- "Violin1"
v2.right["ins"] <- "Violin2"
v.right["ins"] <- "Viola"
c.right["ins"] <- "Cello"
colnames(v1.left)[1] <- "vel"
colnames(v2.left)[1] <- "vel"
colnames(v.left)[1] <- "vel"
colnames(c.left)[1] <- "vel"
colnames(v1.right)[1] <- "vel"
colnames(v2.right)[1] <- "vel"
colnames(v.right)[1] <- "vel"
colnames(c.right)[1] <- "vel"
v1.left["side"] <- "left"
v2.left["side"] <- "left"
v.left["side"] <- "left"
c.left["side"] <- "left"
v1.right["side"] <- "right"
v2.right["side"] <- "right"
v.right["side"] <- "right"
c.right["side"] <- "right"
arms <- rbind(v1.left, v1.right, v2.left, v2.right, v.left, v.right, c.left, c.right)
arms$condition <- with(arms, ifelse(condition == 7 & movement == 3, 8,
                                    ifelse(condition == 7 & movement == 4, 9,
                                           ifelse(condition == 6 & movement == 2, 7, 
                                                  condition))))
adata2 <- data.frame(matrix(rep(0, 7), nrow = 1))
for (i in unique(arms$ins)) {
  insdata <- arms[arms$ins == i, ]
  insdata["nvel3d"] <- with(insdata, (vel - min(vel, na.rm = T))/
                              (max(vel, na.rm = T) - min(vel, na.rm = T)))
  colnames(adata2) <- colnames(insdata)
  adata2 <- rbind(adata2, insdata)
}
adata2 <- adata2[-1,]
write.table(adata2, "armmotion-processed.txt")
