#Angler Control Index
#wavelet analysis to determine when the fight between fish and angler is controlled by the fish (low ACI) vs the angler (high ACI)
#based on identification of "pump-and-wind" type movements by angler

#code to accompany Willis et al. "Thermocline-associated habitat use and fight dynamics of swordfish during recreational angling events"
#Ciara Willis 2026


#setup----
dat = read.csv("TotalSwordfishData.csv") #your file path here

#remove Phase 0 and Phase 7 (tag in air)
dat = dat %>% filter(!(Phase %in% c(0,7)))

#format time
dat$DateTimeStamp<-as.POSIXct(dat$DateTimeStamp,format="%d/%m/%Y %H:%M")
dat$Timemin <- strftime(dat$DateTimeStamp, format = "%H:%M:%S")

dat$DateTimeStampSec2 = paste(dat$Date, dat$DateTimeStampSec, sep = " ")
dat$Timesec <-as.POSIXct(dat$DateTimeStampSec2,format="%d/%m/%Y %H:%M:%S")

head(dat)


#time elapsed since bait deployment
dat$TimeElapsed = NA
for (i in unique(dat$FishID)) {
  xrowmin = min(which(dat$FishID==i))
  xrowmax = max(which(dat$FishID==i))
  
  xseq = seq(xrowmin:xrowmax)
  
  dat$TimeElapsed[xrowmin:xrowmax] = xseq
  
}

#center time on when fish was hooked
#negative before, positve after
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

dat$TimeElapsed_Hook = NA
dat$TimeElapsed_Hook_01 = NA
for (i in unique(dat$FishID)) {
  
  x = dat %>% filter(FishID == i)
  
  xhookrow = min(which(x$Phase == 3))
  # xhookrow = min(which(dat$FishID == i & dat$Phase == 3))
  
  xtime = c(-((xhookrow-1):1), 0, 1:( dim(x)[1] - xhookrow ))
  length(xtime)
  dim(x)
  
  #standardize TimeElapsed_Hook from 0-1
  xtime_01 = range01(xtime)
  
  xrowmin = min(which(dat$FishID==i))
  xrowmax = max(which(dat$FishID==i))
  
  dat$TimeElapsed_Hook[xrowmin:xrowmax] = xtime
  dat$TimeElapsed_Hook_01[xrowmin:xrowmax] = xtime_01
  
}

dfight = dat %>% dplyr::filter(Phase %in% c(5))


#run wavelet analysis-----

dfight$PWI20_70 = NA
dfight$PWI20_60 = NA

for (i in unique(dat$FishID)[-c(7)]) { 
  print(paste(i, "start"))
  
  dwave = dfight %>% dplyr::filter(FishID == i & Phase %in% c(5))
  
  xrowmin = min(which(dfight$FishID==i))
  xrowmax = max(which(dfight$FishID==i))
  
  #fill in any data gaps so 1 depth reading per s
  dwave_alldates = seq(min(dwave$Timesec), max(dwave$Timesec), by = 1)
  interpdata = approx(dwave$Timesec, dwave$Depth, xout = dwave_alldates)
  
  dwave_interp = data.frame(Timesec = dwave_alldates, Depth = interpdata$y)
  
  #run wavelet analysis
  my.w = analyze.wavelet(dwave_interp, "Depth",
                         loess.span = 0.1,
                         dt = 1, #observation of x is made per time unit. (This defines the time unit.)
                         dj = 1/250, # The period range ... Graphically, the argument dj thus determines the resolution along the period axis
                         lowerPeriod = 8,
                         upperPeriod = 1024,
                         #: This defines the range of periods, expressed in time units, used in the wavelet transformation. Only periods of x within this range will be detected
                         make.pval = TRUE,
                         n.sim = 5, method = "white.noise"
  )
  
  #standardize across each second
  my.w$Power <- apply(my.w$Power, 2, \(x)x/max(x))  
  
  #visualize
  p = ggplot(data.frame(reshape2::melt(my.w$Power), Period = my.w$Period),
             aes(x = Var2, y = Period, fill = value)) +
    theme(legend.text = element_blank())+
    geom_tile() +
    scale_fill_viridis_b(option = "H", n.breaks = 20) +
    scale_y_log10() + geom_hline(yintercept = c(20,70), colour = "white", linetype = "dashed")+ #CW power band limits
    labs(
      title = i,
      x = "Time elapsed (s)", fill = "Power")
  # p
  ggsave(paste0("PowerbyPeriod_",i,".png"), height = 6, width = 10, units = "in")


  #pull out wavelets with periods between 20 and 70 seconds
  rowp20_70 = which(my.w$Period > 20 & my.w$Period < 70) #row index for periods bt 20-70s
  power20_70 = my.w$Power[rowp20_70,]
  
  rowp20_60 = which(my.w$Period > 20 & my.w$Period < 60) #row index 
  power20_60 = my.w$Power[rowp20_60,]
  
  
  #sum each col (sum of power at period 20-70 per time step)
  power20_70_colsums = colSums(power20_70, na.rm = T)
  summary(na.omit(power20_70_colsums))
  length(power20_70_colsums)
  
  power20_60_colsums = colSums(power20_60, na.rm = T)
  summary(power20_60_colsums)
  length(power20_60_colsums)
  
  dwave_interp$period_sumpower20_70 = power20_70_colsums
  dwave_interp$period_sumpower20_60 = power20_60_colsums
  
  head(dwave_interp)
  
  #bring it back to the original dwave df
  dwave2 = merge(dwave, dwave_interp, all.x = T, all.y = F)

  #put the timestamps in order
  dwave3 = dwave2 %>% arrange(Timesec)
  head(dwave3$Timesec) 
  
  #standardize the power20_70_colsums from 0-1
  stand_PWI20_70 = rescale(dwave3$period_sumpower20_70, to = c(0,1))
  stand_PWI20_60 = rescale(dwave3$period_sumpower20_60, to = c(0,1))
  # summary(stand_PWI)
  
  dfight$PWI20_70[xrowmin:xrowmax] = stand_PWI20_70
  dfight$PWI20_60[xrowmin:xrowmax] = stand_PWI20_60
  
  print(paste(i, "stop"))
  
}

#save dfight