#############################################
#                                           #
#           Shaw, D'Aloia, Buston           #
#          Comparison of simulated          #
#              & empirical data             #
#                                           #
#############################################


##Shaw et al. 201X - The evolution of marine larval dispersal kernels in spatially structured habitats:
#                    analytical models, individual-based simulations, and comparisons with empirical estimates
#                    The American Naturalist.

# R code for comparison of simulated data to empirical data (D'Aloia et al. 2015, PNAS)
# https://doi.org/10.1073/pnas.1513754112

#########################################################
#0 Setup
#########################################################

# Set your working directory
# Load these packages

require(reshape2)
require(ggplot2)
require(gridExtra)
require(fitdistrplus)

#########################################################
#1 Generate appdendix Figure E.3
#########################################################

# Here, we will compare the empirical data to simulated data under different levels of 
# breeding resource heterogeneity (none, low, medium, high). For simulated data, we use
# post-survival, precompetition frequencies to most closely match empirical data.

dispersal.freq <-read.csv("dispersal_frequencies.csv", header=TRUE)

# E.3A: Empirical and no patchiness

E3A.data <- dispersal.freq[dispersal.freq$category %in% c("simulated_none", "empirical"),]
  
E3A.data$distance <- factor(E3A.data$distance, levels=c("0-0.5", "0.5-1.5","1.5- 2.5", "2.5-3.5", "3.5-4.5", "4.5-5.5",
                                              "5.5-6.5", "6.5-7.5", "7.5-8.5", "8.5-9.5", "9.5-10.5",
                                              "10.5-11.5", "11.5-12.5", "12.5 - 13.5", "13.5 - 14.5",
                                              "14.5-15.5", "15.5 - 16.5", "16.5 - 17.5", "17.5 - 18.5",
                                              "18.5-19.5", "19.5 - 20.5", "20.5- 21.5", "21.5-22.5",
                                              "22.5 - 23.5", "23.5-24.5", "24.5-25.5", "25.5-26.5",
                                              "26.5-27.5", "27.5-28.5", "28.5-29.5", "29.5 - 30"))

E3A.data$category <- factor(E3A.data$category, levels=c("simulated_none", "empirical"))

E.3A <- ggplot(E3A.data, aes(x = distance,y = rel_freq)) + 
  geom_bar(aes(fill = category), stat="identity", position = "dodge") +
  scale_y_continuous(limits=c(0,0.85),expand = c(0,0))+
  scale_fill_discrete(name=NULL, breaks=c("simulated_none", "empirical"),
                    labels=c("No patchiness", "Empirical")) +
  theme(axis.title.x = element_text(size=18, vjust=0),
        axis.title.y = element_text(size=10, vjust=0.2), axis.text.x = element_text(angle=90),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = c(0.85, 0.8)) +
  labs(title = "(a) No patchiness vs. empirical", y="Relative Frequency of Dispersal Events", x=NULL)
E.3A

# E.3B: Empirical and low patchiness

E3B.data <- dispersal.freq[dispersal.freq$category %in% c("simulated_low", "empirical"),]

E3B.data$distance <- factor(E3B.data$distance, levels=c("0-0.5", "0.5-1.5","1.5- 2.5", "2.5-3.5", "3.5-4.5", "4.5-5.5",
                                              "5.5-6.5", "6.5-7.5", "7.5-8.5", "8.5-9.5", "9.5-10.5",
                                              "10.5-11.5", "11.5-12.5", "12.5 - 13.5", "13.5 - 14.5",
                                              "14.5-15.5", "15.5 - 16.5", "16.5 - 17.5", "17.5 - 18.5",
                                              "18.5-19.5", "19.5 - 20.5", "20.5- 21.5", "21.5-22.5",
                                              "22.5 - 23.5", "23.5-24.5", "24.5-25.5", "25.5-26.5",
                                              "26.5-27.5", "27.5-28.5", "28.5-29.5", "29.5 - 30"))

E3B.data$category <- factor(E3B.data$category, levels=c("simulated_low", "empirical"))

E.3B <- ggplot(E3B.data, aes(x = distance, y = rel_freq)) + 
  geom_bar(aes(fill = category), stat="identity", position = "dodge") +
  scale_y_continuous(limits=c(0,0.85),expand = c(0,0))+
  scale_fill_discrete(name=NULL, breaks=c("simulated_low", "empirical"),
                      labels=c("Low patchiness", "Empirical")) +
  theme(axis.title.x = element_text(size=18, vjust=0),
        axis.title.y = element_text(size=15, vjust=0.2), axis.text.x = element_text(angle=90),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = c(0.85, 0.8)) +
  labs(title = "(b) Low patchiness vs. empirical", y=NULL, x=NULL)
E.3B 

# E.3C: Empirical and medium patchiness

E3C.data <- dispersal.freq[dispersal.freq$category %in% c("simulated_med", "empirical"),]

E3C.data$distance <- factor(E3C.data$distance, levels=c("0-0.5", "0.5-1.5","1.5- 2.5", "2.5-3.5", "3.5-4.5", "4.5-5.5",
                                              "5.5-6.5", "6.5-7.5", "7.5-8.5", "8.5-9.5", "9.5-10.5",
                                              "10.5-11.5", "11.5-12.5", "12.5 - 13.5", "13.5 - 14.5",
                                              "14.5-15.5", "15.5 - 16.5", "16.5 - 17.5", "17.5 - 18.5",
                                              "18.5-19.5", "19.5 - 20.5", "20.5- 21.5", "21.5-22.5",
                                              "22.5 - 23.5", "23.5-24.5", "24.5-25.5", "25.5-26.5",
                                              "26.5-27.5", "27.5-28.5", "28.5-29.5", "29.5 - 30"))

E3C.data$category <- factor(E3C.data$category, levels=c("simulated_med", "empirical"))

E.3C <- ggplot(E3C.data, aes(x = distance, y = rel_freq)) + 
  geom_bar(aes(fill = category), stat="identity", position = "dodge") +
  scale_y_continuous(limits=c(0,0.85),expand = c(0,0))+
  scale_fill_discrete(name=NULL, breaks=c("simulated_med", "empirical"),
                      labels=c("Medium patchiness", "Empirical")) +
  theme(axis.title.x = element_text(size=18, vjust=0),
        axis.title.y = element_text(size=10, vjust=0.2), axis.text.x = element_text(angle=90),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = c(0.85, 0.8)) +
  labs(title = "(c) Medium patchiness vs. empirical", y="Relative Frequency of Dispersal Events", x="Distance Class (km)")
E.3C

# E.3D: Empirical and high patchiness

E3D.data <- dispersal.freq[dispersal.freq$category %in% c("simulated_high", "empirical"),]

E3D.data$distance <- factor(E3D.data$distance, levels=c("0-0.5", "0.5-1.5","1.5- 2.5", "2.5-3.5", "3.5-4.5", "4.5-5.5",
                                              "5.5-6.5", "6.5-7.5", "7.5-8.5", "8.5-9.5", "9.5-10.5",
                                              "10.5-11.5", "11.5-12.5", "12.5 - 13.5", "13.5 - 14.5",
                                              "14.5-15.5", "15.5 - 16.5", "16.5 - 17.5", "17.5 - 18.5",
                                              "18.5-19.5", "19.5 - 20.5", "20.5- 21.5", "21.5-22.5",
                                              "22.5 - 23.5", "23.5-24.5", "24.5-25.5", "25.5-26.5",
                                              "26.5-27.5", "27.5-28.5", "28.5-29.5", "29.5 - 30"))

E3D.data$category <- factor(E3D.data$category, levels=c("simulated_high", "empirical"))

E.3D <- ggplot(E3D.data, aes(x = distance,y = rel_freq)) + 
  geom_bar(aes(fill = category), stat="identity", position = "dodge") +
  scale_y_continuous(limits=c(0,0.85),expand = c(0,0))+
  scale_fill_discrete(name=NULL, breaks=c("simulated_high", "empirical"),
                      labels=c("High patchiness", "Empirical")) +
  theme(axis.title.x = element_text(size=18, vjust=0),
        axis.title.y = element_text(size=18, vjust=0.2), axis.text.x = element_text(angle=90),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = c(0.85, 0.8)) +
  labs(title = "(d) High patchiness vs. empirical", y=NULL, x="Distance Class (km)")
E.3D


grid.arrange(E.3A, E.3B, E.3C, E.3D, ncol = 2, nrow = 2)

###############################################
#2 Fit kernels to simulated and empirical data
##############################################

#Use high-patchy scenario as an example, read in simulated dispersal distances
distances <- read.csv("distances.csv", header=TRUE)
ds1 <- distances$dist_hi
ds1 <- ds1[!is.na(ds1)]

#Fit exponential distribution to simulated data
exp <- fitdist(ds1, "exp", method = "mle")
summary(exp) #rate  2.09793; se = 0.04018904

#Empirical distances
emp.dist <- distances$dist_empirical
emp.dist <- emp.dist[!is.na(emp.dist)]

# Calulate quantiles for simulated pre-survival and empirical data
simulated.quant <- quantile(ds1, c(0.5, 0.95, 0.99))
empirical.quant <- quantile(emp.dist, c(0.5, 0.95, 0.99))

###############################################
#3 Plot Figure 5
##############################################

# Fig. 5 Main plot: Empirical and simulated dispersal kernels
plot(1,xlim=c(0,30),ylim=c(0,2.5),xlab="Distance (km)",ylab="Probability density",col="white", xaxs="i", yaxs="i")
  x<-seq(0,30,0.1)
  lines(x, dexp(x,0.3566428 ), lwd=2, col="darkgray")
  lines(x, dexp(x, 2.09793), lwd=2, col="black")
  legend("topleft", c("Empirical","Simulated"), lty=c(1,1), lwd=c(2,2),col=c("darkgray","black"))

  
#Fig. 5b Inset plot: effect of carrying capacity per site on evolved dispersal pattern

inset.data <- dispersal.freq[dispersal.freq$category %in% c("capacity_1", "capacity_2", "capacity_4"),]

#Truncate to 10 km for better visibility
inset.truncate <- inset.data[! inset.data$distance %in% c("10.5-11.5", "11.5-12.5", "12.5 - 13.5","13.5 - 14.5",
                                                            "14.5-15.5","15.5 - 16.5", "16.5 - 17.5", "17.5 - 18.5", "18.5-19.5",
                                                            "19.5 - 20.5", "20.5- 21.5", "21.5-22.5", "22.5 - 23.5", "23.5-24.5",
                                                            "24.5-25.5", "25.5-26.5","26.5-27.5","27.5-28.5","28.5-29.5","29.5 - 30"), ]

fig_5_inset <- ggplot(inset.truncate, aes(x = distance, y = rel_freq, fill = category)) + 
  geom_bar(colour="black", stat="identity",  width=.75, position = "dodge") +
  scale_fill_manual(values=c("black", "gray", "white"), labels = c("1", "2", "4"))+
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = c(0.8, 0.8)) +
  xlab("Distance (km)") + ylab("Relative Frequency") + labs(fill = "Individuals/Site")
fig_5_inset
