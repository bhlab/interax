##Dependencies
library(lubridate)
library(dplyr)
library(tidyr)
library(overlap)
library(circular)

##function for generating the activity density of species2
sp2time <- function(data,
                    datetime_col,
                    species_col,
                    species_name,
                    datetime_format,
                    tz = "UTC",
                    bw = "default",
                    export_csv = FALSE,
                    output_file = "activity_density.csv",
                    plot = TRUE) {
  df <- data
  if (!is.null(species_col) && !is.null(species_name)) {
    df <- subset(df, df[[species_col]] == species_name)
  }
  if (is.null(datetime_format)) {
    df[[datetime_col]] <- as.POSIXct(df[[datetime_col]], tz = tz)
    if (all(is.na(df[[datetime_col]]))) {
      stop("Date-time conversion failed. Please mention datetime_format.")
    }
  } else {
    df[[datetime_col]] <- as.POSIXct(df[[datetime_col]],
                                     format = datetime_format,
                                     tz = tz)
  }
  
  hrs  <- as.numeric(format(df[[datetime_col]], "%H"))
  mins <- as.numeric(format(df[[datetime_col]], "%M"))
  secs <- as.numeric(format(df[[datetime_col]], "%S"))
  time_rad <- hrs  * 2 * pi / 24 +
    mins * 2 * pi / (24 * 60)+
    secs * 2 * pi / (24 * 60 * 60)
  time_rad <- circular(time_rad, units = "radians")
  time_rad <- as.numeric(time_rad)
  dens <- densityFit(time_rad, bw = bw)
  if (plot) {
    densityPlot(time_rad, rug = FALSE, main = "Activity Density")
  }
  dens_df <- data.frame(
    radians = dens$x,
    density = dens$y
  )
  if (export_csv) {
    write.csv(dens_df, output_file, row.names = FALSE)
  }
  return(dens_df)
}

#function for randomizing Species 2 sighting within sighting range
rand.data <- function(data.sp1, data.sp2, sp2time){
  randdat <- sample(as.Date(data.sp2$dt), nrow(data.sp2), replace = TRUE)
  x <- sp2time[sp2time$x>=0 & sp2time$x<=24,]
  timelist <- sample(x$x,nrow(data.sp2),replace=TRUE,prob=x$y)
  H <- floor(timelist)
  M.x <- (timelist-H)*60
  M <- floor((timelist-H)*60)
  S <- round((M.x-M)*60)
  randtime <- paste(H,":",M,":",S, sep = "")
  data.sp2$dt <- strptime(paste(as.character(randdat), randtime), "%Y-%m-%d %H:%M:%S")
  data.n <- rbind(data.sp1, data.sp2)
  return(data.n)
}

#function to make matrix from data
make.mat <- function(data, days=10, Species1) {
  data <- data[order(data$dt),]
  sp1.list <- which(data$Species == Species1)
  o.mat <- NULL
  n.sp1 <- NULL
  for (i in 1:length(sp1.list)) {
    sp1cid <- data[sp1.list[i],'Camera']
    dat1 <- data[sp1.list[i]:nrow(data),]
    dat1 <- dat1 %>% filter(Camera == sp1cid)
    sp1seq <- seq(from = dat1$dt[1], length.out = days+1, by = "day")
    dat1 <- dat1 %>% filter(dt >= sp1seq[1] & dt <= sp1seq[days+1]) 
    if (dat1 %>% filter(Species == Species1) %>% nrow == 1) {
      l.sight <- NULL
      for (j in 1:days) {
        c <- dat1[-1,] %>% filter(dt >= sp1seq[j] & dt <= sp1seq[j+1]) %>% nrow
        l.sight <- c(l.sight,c)
      }
      o.mat <- rbind(o.mat,l.sight)
      n.sp1 <- c(n.sp1, sp1.list[i])
    } else {
      #print(paste("T",sp1.list[i]))
    }
  }
  colnames(o.mat) <- paste("Day", c(1:days))
  rownames(o.mat) <- n.sp1
  return(o.mat)
}

#function to make matrix from data for hours
make.mat.hr <- function(data, hours=10, Species1) {
  data <- data[order(data$dt),]
  sp1.list <- which(data$Species == Species1)
  o.mat <- NULL
  n.sp1 <- NULL
  for (i in 1:length(sp1.list)) {
    sp1cid <- data[sp1.list[i],'Camera']
    dat1 <- data[sp1.list[i]:nrow(data),]
    dat1 <- dat1 %>% filter(Camera == sp1cid)
    sp1seq <- seq(from = dat1$dt[1], length.out = hours+1, by = "hour")
    dat1 <- dat1 %>% filter(dt >= sp1seq[1] & dt <= sp1seq[hours+1]) 
    if (dat1 %>% filter(Species == Species1) %>% nrow == 1) {
      l.sight <- NULL
      for (j in 1:hours) {
        c <- dat1[-1,] %>% filter(dt >= sp1seq[j] & dt <= sp1seq[j+1]) %>% nrow
        l.sight <- c(l.sight,c)
      }
      o.mat <- rbind(o.mat,l.sight)
      n.sp1 <- c(n.sp1, sp1.list[i])
    } else {
      #print(paste("T",sp1.list[i]))
    }
  }
  colnames(o.mat) <- paste("Hour", c(1:hours))
  rownames(o.mat) <- n.sp1
  return(o.mat)
}

#function to find out before vs after sightings
bf.af <- function(data, Species1, Species2, hrs=12) {
  #data <- data %>% group_by(CT_ID) %>% arrange(dt, .by_group = TRUE)
  data <- data[order(data$dt),]
  sp1.list <- which(data$Species==Species1)
  before <- 0
  after <- 0
  for (i in 1:length(sp1.list)) {
    det.t <- data$dt[sp1.list[i]]
    det.pt <- det.t + (hrs*60*60)
    det.2pt <- det.t + (2*hrs*60*60)
    det.nt <- det.t - (hrs*60*60)
    det.2nt <- det.t - (2*hrs*60*60)
    sp1cid <- data[sp1.list[i],"Camera"]
    be.t <- data %>% filter(Camera == sp1cid, dt >=  det.2nt & dt <= det.t)
    if (be.t %>% filter(Species == Species1) %>% nrow == 1){
      be.dt <- be.t %>% filter(dt >=  det.nt & dt <= det.t)
      bef <- be.dt %>% filter(Species == Species2) %>% nrow
      before <- before+bef
    }
    af.t <- data %>% filter(Camera == sp1cid, dt >=  det.t & dt <= det.2pt)
    if (af.t %>% filter(Species == Species1) %>% nrow == 1){
      af.dt <- af.t %>% filter(dt >=  det.t & dt <= det.pt)
      aft <- af.dt %>% filter(Species == Species2) %>% nrow
      after <- after+aft
    }
  }
  result <- c("before"= before, "after" = after)
  return(result)
}

#function to list sightings with the time
post.cap <- function(data, Species1, Species2) {
  data <- data[order(data$dt),]
  sp1.list <- which(data$Species == Species1)
  t2en <- NULL
  for (i in 1:length(sp1.list)) {
    sp1cid <- data[sp1.list[i],"Camera"]
    dat1 <- data[sp1.list[i]:nrow(data),]
    dat1 <- dat1[dat1$Camera == sp1cid,]
    if (nrow(dat1) <= 1 | dat1$Species[2] != Species2) next
    en <- difftime(dat1$dt[2], dat1$dt[1], units="hours")
    t2en <- c(t2en, en)
  }
  return(t2en)
}

#function to list sightings with the time and Camera 
# ! TimeDiff in cbind is coerced as character
post.capCT <- function(data, Species1, Species2) {
  data <- data[order(data$dt),]
  sp1.list <- which(data$Species==Species1)
  t2en <- NULL
  for (i in 1:length(sp1.list)) {
    sp1cid <- data[sp1.list[i],"Camera"]
    dat1 <- data[sp1.list[i]:nrow(data),]
    dat1 <- dat1[dat1$Camera == sp1cid,]
    if (nrow(dat1) <= 1 | dat1$Species[2] != Species2) next
    TimeDiff <- difftime(dat1$dt[2], dat1$dt[1], units="hours")
    Camera <- sp1cid
    DateTime <- dat1$DateTime[2]
    en <- cbind(Camera, cbind(DateTime, TimeDiff))
    t2en <- rbind(t2en, en)
  }
  return(t2en)
}

#function for permutation t-test to compare detection probabilites for first 5 days and last 5 days
perm.test <- function (det.p, data.sp1, data.sp2, sp2time, sample1.lab = NULL, sample2.lab = NULL, 
                       B = 999) {
  options(scipen = 999)
  sample1 <- det.p[1:5]
  sample2 <- det.p[6:10]
  if (is.null(sample1.lab) == TRUE) {
    sample1.lab <- "smpl 1"
    sample2.lab <- "smpl 2"
  }
  else {
  }
  n1 <- length(sample1)
  n2 <- length(sample2)
  mean1 <- round(mean(sample1), 2)
  mean2 <- round(mean(sample2), 2)
  error1 <- qnorm(0.975) * sd(sample1)/sqrt(n1)
  error2 <- qnorm(0.975) * sd(sample2)/sqrt(n2)
  sample1_lci <- round(mean1 - error1, 2)
  sample1_uci <- round(mean1 + error1, 2)
  sample2_lci <- round(mean2 - error2, 2)
  sample2_uci <- round(mean2 + error2, 2)
  p.equal.var <- round(t.test(sample1, sample2, var.equal = TRUE)$p.value, 
                       4)
  p.unequal.var <- round(t.test(sample1, sample2, var.equal = FALSE)$p.value, 
                         4)
  pooledData <- c(sample1, sample2)
  nIter <- B
  meanDiff <- numeric(nIter + 1)
  meanDiff[1] <- round(mean1 - mean2, digits = 2)
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (i in 2:B) {
    data.n <- rand.data(data.sp1, data.sp2, sp2time)
    o.mat <- make.mat(data.n)
    det.p2 <- colSums(o.mat)/sum(o.mat)
    sample1.perm <- det.p2[1:5]
    sample2.perm <- det.p2[6:10]
    meanDiff[i] <- mean(sample1.perm) - mean(sample2.perm)
    setTxtProgressBar(pb, i)
  }
  p.lowertail <- (1 + sum(meanDiff[-1] < meanDiff[1]))/(1 + B)
  p.uppertail <- (1 + sum(meanDiff[-1] > meanDiff[1]))/(1 + B)
  two.sided.p <- 2 * min(p.lowertail, p.uppertail)
  graphics::hist(meanDiff, main = paste0("Distribution of permuted mean differences\n(number of permutations: ", 
                                         B, ")"), xlab = "", sub = paste0(sample1.lab, " (n: ", 
                                                                          n1, ") (95% CI lower bound., mean, 95% CI upper bound.): ", 
                                                                          sample1_lci, ", ", mean1, ", ", sample1_uci, "\n", sample2.lab, 
                                                                          " (n: ", n2, ") (95% CI lower bound., mean, 95% CI upper bound.): ", 
                                                                          sample2_lci, ", ", mean2, ", ", sample2_uci, "\nobs. mean diff. (solid dot): ", 
                                                                          meanDiff[1], "; perm. p.value mean ", sample1.lab, " < ", 
                                                                          sample2.lab, ": ", round(p.lowertail, 4), "; perm. p.value mean ", 
                                                                          sample1.lab, " > ", sample2.lab, ": ", round(p.uppertail, 
                                                                                                                       4), "; perm. p.value (2-sided): ", round(two.sided.p, 
                                                                                                                                                                4), "\nregular t-test p-values (2-sided): ", round(p.equal.var, 
                                                                                                                                                                                                                   4), " (equal variance); ", round(p.unequal.var, 4), 
                                                                          " (unequal variance)"), cex.main = 0.85, cex.sub = 0.7)
  rug(meanDiff, col = "#0000FF")
  abline(v = stats::quantile(meanDiff, 0.025), lty = 2, col = "blue")
  abline(v = stats::quantile(meanDiff, 0.975), lty = 2, col = "blue")
  points(x = meanDiff[1], y = 0, pch = 20, col = "black")
  points(x = mean(meanDiff[-1]), y = 0, pch = 1, col = "black")
}

#function for permutation t-test to compare detection probabilities
perm.avoid <- function (data, sp2time, days= 10, Species1, Species2, B = 1000) {
  options(scipen = 1000)
  
  o.mat <- make.mat(data, days, Species1) #make a matrix of sightings
  det.p <- colSums(o.mat)/sum(o.mat) #detection probabilities per day after Species capture
  
  data.sp1 <- data[row.names(o.mat),]
  data.sp2 <- data %>% filter(Species==Species2)
  
  sample.time <- post.cap(data, Species1, Species2)
  
  sample1 <- det.p[1:(days/2)]
  sample2 <- det.p[(1+days/2):days]
  sample.mDiff <- abs(mean(sample1)-mean(sample2))
  
  meanDiff <- numeric(B)
  null.time <- numeric(B)
  null.det.p <- NULL
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (i in 1:B) {
    data.n <- rand.data(data.sp1, data.sp2, sp2time)
    null.time[i] <- mean(post.cap(data.n, Species1, Species2))
    n.mat <- make.mat(data.n, days, Species1)
    det.p2 <- colSums(n.mat)/sum(n.mat)
    sample1.perm <- det.p2[1:(days/2)]
    sample2.perm <- det.p2[(1+days/2):days]
    meanDiff[i] <- abs(mean(sample1.perm) - mean(sample2.perm))
    null.det.p <- rbind(null.det.p, det.p2)
    setTxtProgressBar(pb, i)
  }
  ucomp <- ifelse(meanDiff >= sample.mDiff, 1, 0)
  two.sided.p <- mean(ucomp)
  graphics::boxplot(null.det.p, main = paste0("Distribution of permuted leopard detections\n(number of permutations: ", 
                                              B, ")"), ylab= "Detection probabilities", sub =  paste0("two-sided-p:  ",
                                                                                                      round(two.sided.p,3), "\n",
                                                                                                      "Observed values in magenta dot"))
  points(det.p, pch = 20, col = "magenta")
  
  values <- c("two.sided.p" = two.sided.p, "Observed meantime" = mean(sample.time),
              "Observed mintime" = min(sample.time), "Observed maxtime" = max(sample.time),
              "Null meantime" = mean(null.time), "Null mintime"= min(null.time), "Null maxtime"= max(null.time))
  print(values)
  invisible(list(values,"Observed detection probability"= det.p, 
              "Randomized detection probabilities"= null.det.p))
  
}


#function for permutation t-test to compare detection probabilities at hours interval
perm.avoid.hr <- function (data, sp2time, hours= 10, Species1, Species2, B = 1000) {
  options(scipen = 1000)
  
  o.mat <- make.mat.hr(data, hours, Species1) #make a matrix of sightings
  det.p <- colSums(o.mat)/sum(o.mat) #detection probabilities per day after Species capture
  
  data.sp1 <- data[row.names(o.mat),]
  data.sp2 <- data %>% filter(Species==Species2)
  
  sample.time <- post.cap(data, Species1, Species2)
  
  sample1 <- det.p[1:(hours/2)]
  sample2 <- det.p[(1+hours/2):hours]
  sample.mDiff <- abs(mean(sample1)-mean(sample2))
  
  meanDiff <- numeric(B)
  null.time <- numeric(B)
  null.det.p <- NULL
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  for (i in 1:B) {
    data.n <- rand.data(data.sp1, data.sp2, sp2time)
    null.time[i] <- mean(post.cap(data.n, Species1, Species2))
    n.mat <- make.mat.hr(data.n, hours, Species1)
    det.p2 <- colSums(n.mat)/sum(n.mat)
    sample1.perm <- det.p2[1:(hours/2)]
    sample2.perm <- det.p2[(1+hours/2):hours]
    meanDiff[i] <- abs(mean(sample1.perm) - mean(sample2.perm))
    null.det.p <- rbind(null.det.p, det.p2)
    setTxtProgressBar(pb, i)
  }
  ucomp <- ifelse(meanDiff >= sample.mDiff, 1, 0)
  two.sided.p <- mean(ucomp)
  graphics::boxplot(null.det.p, main = paste0("Distribution of permuted ", Species2, " detections\n(number of permutations: ", 
                                              B, ")"), ylab= "Detection probabilities", sub =  paste0("two-sided-p:  ",
                                                                                                      round(two.sided.p,3), "\n",
                                                                                                      "Observed values in magenta dot"))
  points(det.p, pch = 20, col = "magenta")
  
  values <- c("two.sided.p" = two.sided.p, "Observed meantime" = mean(sample.time),
              "Observed mintime" = min(sample.time), "Observed maxtime" = max(sample.time),
              "Null meantime" = mean(null.time), "Null mintime"= min(null.time), "Null maxtime"= max(null.time))
  print(values)
  invisible(list(values,"Observed detection probability"= det.p, 
                 "Randomized detection probabilities"= null.det.p))
  
}
# Function to list sightings with the time and calculate event duration ## I set default break time 5 mins
time.spent <- function(data, species, breaktime = 5) {
  data <- data %>%
    arrange(dt)
  t.list <- which(data$Species == species)
  results <- data.frame(Station = character(), 
                        event.start = as.POSIXct(character()), 
                        event.end = as.POSIXct(character()), 
                        event.duration = numeric(),
                        stringsAsFactors = FALSE)
  
  if (length(t.list) == 0) {
    return(results)
  }
  start_time <- NULL
  end_time <- NULL
  
  for (i in 1:length(t.list)) {
    current_station <- data$Station[t.list[i]]
    
    if (is.null(start_time)) {
      start_time <- data$dt[t.list[i]]
    }
    if (i == length(t.list) || difftime(data$dt[t.list[i + 1]], data$dt[t.list[i]], units = "mins") > breaktime) {
      end_time <- data$dt[t.list[i]]
      time_spent <- as.numeric(difftime(end_time, start_time, units = "mins"))
      results <- rbind(results, data.frame(Station = current_station, 
                                           event.start = start_time, 
                                           event.end = end_time, 
                                           event.duration = time_spent))
      start_time <- NULL
    }
  }
  return(results)
}

# Function to list sightings with the time and calculate event duration ## default break time 5 mins
time.spent <- function(data, species, breaktime = 5) {
  data <- data %>%
    arrange(dt)
  t.list <- which(data$Species == species)
  results <- data.frame(Station = character(), 
                        event.start = as.POSIXct(character()), 
                        event.end = as.POSIXct(character()), 
                        event.duration = numeric(),
                        stringsAsFactors = FALSE)
  
  if (length(t.list) == 0) {
    return(results)
  }
  start_time <- NULL
  end_time <- NULL
  
  for (i in 1:length(t.list)) {
    current_station <- data$Station[t.list[i]]
    
    if (is.null(start_time)) {
      start_time <- data$dt[t.list[i]]
    }
    if (i == length(t.list) || difftime(data$dt[t.list[i + 1]], data$dt[t.list[i]], units = "mins") > breaktime) {
      end_time <- data$dt[t.list[i]]
      time_spent <- as.numeric(difftime(end_time, start_time, units = "mins"))
      results <- rbind(results, data.frame(Station = current_station, 
                                           event.start = start_time, 
                                           event.end = end_time, 
                                           event.duration = time_spent))
      start_time <- NULL
    }
  }
  return(results)
}

#### Function for making coocurrence matrices using diffeent time intervals. Default time interval is "5 mins"
make_site_sp_mat <- function(data,
                             datetime_col,
                             camera_col,
                             species_col,
                             time_intervals = c("5 mins"),
                             export_csv = FALSE,
                             output_prefix = "site_species_matrix") {
  df <- data
  df$datetime <- suppressWarnings(lubridate::dmy_hms(df[[datetime_col]]))
  if (all(is.na(df$datetime))) {
    stop("Datetime parsing failed. Check datetime format.")
  }
  results <- list()
  for (interval in time_intervals) {
    site_time <- round_date(df$datetime, unit = interval)
    site_time <- paste(df[[camera_col]], site_time, sep = "_")
    site_sp_mat <- df %>%
      mutate(site_time = site_time) %>%
      distinct(site_time, .data[[species_col]]) %>%
      mutate(Present = 1L) %>%
      pivot_wider(
        names_from = .data[[species_col]],
        values_from = Present,
        values_fill = 0L
      )
    results[[interval]] <- site_sp_mat
    if (export_csv) {
      safe_name <- gsub(" ", "_", interval)
      write.csv(site_sp_mat,
                paste0(output_prefix, "_", safe_name, ".csv"),
                row.names = FALSE)
    }
  }
  return(results)
}
