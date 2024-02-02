graphics.off() 
rm(list=ls(all=TRUE)) 

temperature_100 = c(
  58.5,
  57.3,
  59.7,
  59.1,
  58.3,
  56.3,
  57.8,
  58,
  58,
  57.3,
  59.2,
  59.4,
  57.5,
  56.3,
  55.6,
  58.3,
  58.7,
  57.3,
  57.9,
  58,
  55,
  57.6,
  55,
  56.1,
  56.9,
  58.3,
  58,
  58
)

temperature_110 = c(
  57.2,
  59,
  61,
  60.2,
  55.5,
  56.1,
  57.2,
  58.9,
  60.7,
  63,
  60.4,
  62.1,
  61.8,
  55.7,
  57.7,
  57.3,
  62.3,
  63.5,
  53.5,
  53.6,
  55.5,
  56.5,
  54.7,
  56.3,
  49.4,
  52,
  47.5,
  57
)

temperature_120 = c(
  56.7,
  54.5,
  32.1,
  36.2,
  55.6,
  57.2,
  39.4,
  53.2,
  56.2,
  54.4,
  58.3,
  57.1,
  52.7,
  58.1,
  56.4,
  60.2,
  55.6,
  55,
  57.7,
  61.2,
  59.6,
  55.7,
  48.3,
  58.1,
  57.1,
  58.2,
  48.4,
  57
)

# 2. Histogram 

tiff(paste("temperature_histo.tiff", sep="" ), width = 2000, height = 2000,
     units = "px", res = 1000, pointsize = 7)

library(ggplot2)
library(dplyr)
library(hrbrthemes)

# Build dataset with different distributions
data <- data.frame(
  Site = c( rep("Switzerland 100", 28), rep("Switzerland 110", 28), rep("Switzerland 120", 28)),
  value = c( temperature_100, temperature_110, temperature_120)
)

# Represent it
p <- data %>%
  ggplot( aes(x=value, color=Site)) +
  geom_density(lwd = 1.2) +
  scale_color_manual(values=c("deeppink","green", "deepskyblue")) +
  xlab("Temperature (\u00B0C)") + ylab("Density") +
  xlim(30, 70) +
  theme_bw() + theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


 

p

  
dev.off()
