library(tidyverse)
HistoricalRootZoneSoil <- read.csv ("~/Uni/Doctorate/Ch2 Stressor/AWO_2020-2024_HistoricalRootZoneSoilMoisture_NSW.csv")

colnames(HistoricalRootZoneSoil) <- HistoricalRootZoneSoil[17,]
HistoricalRootZoneSoil <- HistoricalRootZoneSoil[-c(1:17),]

colnames(HistoricalRootZoneSoil) <- c("Date","RZSMPerc", "Month")
HistoricalRootZoneSoil$Date <- as.Date(HistoricalRootZoneSoil$Date, format = "%d/%m/%y")
HistoricalRootZoneSoil$RZSMPerc <- as.numeric(HistoricalRootZoneSoil$RZSMPerc)
HistoricalRootZoneSoil$Month <- as.numeric(HistoricalRootZoneSoil$Month)

season_mapping <- c('1' = "Summer", '2' = "Summer", '3' = "Autumn", '4' = "Autumn", '5' = "Autumn", '6' = "Winter", '7' = "Winter", '8' = "Winter", '9' = "Spring", '10' = "Spring", '11' = "Spring", '12' = "Summer")
HistoricalRootZoneSoil$Season <- season_mapping[HistoricalRootZoneSoil$Month] 



ggplot(HistoricalRootZoneSoil, aes(x=Date, y=RZSMPerc, colour = Season)) +
  geom_hline(yintercept=c(5,10,25,50), colour = "black", linetype="dashed") +
  geom_point() +
  #geom_path() +
  theme_minimal()

# < 50% RZSM
count_track_50 <- NULL
count = 0
for (i in 1:nrow(HistoricalRootZoneSoil)){
  RZSM_row <- HistoricalRootZoneSoil[i,2]
  if (RZSM_row <= 50){
    count = count + 1
  } else {
    count_track_50 <- rbind(count_track_50, count)
    count = 0
  }
}

max(count_track_50) # 138 days

# < 25% RZSM
count_track_25 <- NULL
count = 0
for (i in 1:nrow(HistoricalRootZoneSoil)){
  RZSM_row <- HistoricalRootZoneSoil[i,2]
  if (RZSM_row <= 25){
    count = count + 1
  } else {
    count_track_25 <- rbind(count_track_25, count)
    count = 0
  }
}

max(count_track_25) # 122 days

# < 10% RZSM
count_track_10 <- NULL
count = 0
for (i in 1:nrow(HistoricalRootZoneSoil)){
  RZSM_row <- HistoricalRootZoneSoil[i,2]
  if (RZSM_row <= 10){
    count = count + 1
  } else {
    count_track_10 <- rbind(count_track_10, count)
    count = 0
  }
}

max(count_track_10) # 57 days

# < 5% RZSM
count_track_5 <- NULL
count = 0
for (i in 1:nrow(HistoricalRootZoneSoil)){
  RZSM_row <- HistoricalRootZoneSoil[i,2]
  if (RZSM_row <= 5){
    count = count + 1
  } else {
    count_track_5 <- rbind(count_track_5, count)
    count = 0
  }
}

max(count_track_5) # 33 days

# Plot
count_track_5 <- as.data.frame(count_track_5) %>% mutate(perc = "1_count_track_5")
count_track_10 <- as.data.frame(count_track_10) %>% mutate(perc = "2_count_track_10")
count_track_25 <- as.data.frame(count_track_25) %>% mutate(perc = "3_count_track_25")
count_track_50 <- as.data.frame(count_track_50) %>% mutate(perc = "4_count_track_50")

days_below <- rbind(count_track_5, count_track_10, count_track_25, count_track_50)

ggplot(days_below, aes(x=perc, y=V1)) +
  geom_point() + 
  theme_bw() +
  labs(x="Percentage of RZSM threshold", y="Days in a row below percentage")
  #geom_violin() +

# What season do we drought