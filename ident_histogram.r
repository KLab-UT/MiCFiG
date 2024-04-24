#install.packages("ggplot2",repos = "http://cran.us.r-project.org")
library("ggplot2")

stats <- read.csv("identities.csv",header = FALSE)

genes <- stats[1]

pident
pident <- as.numeric(unlist(stats[2]))

hist_plot <-ggplot(stats,aes(x=pident)) +
  geom_histogram(binwidth = 1, show.legend = F) +
		 theme_classic() + 
		 scale_color_manual(values = c("gray", "gray")) +  # Change right values to black to bring 
		 scale_fill_manual(values = c("gray", "gray")) +   # shading to values greater than BamA concervancy back
      	         labs(x = "Percent Identity", y = "Number of Hits") +
		 theme(axis.text = element_text(size = 18),  
                 axis.title = element_text(size = 20))
  
  # Adding Vertical lines to histogram
  #geom_vline(aes(xintercept=bamA_proportion), color = "black") +
  #geom_vline(aes(xintercept=tamA_proportion), color = "black")+
  #geom_vline(aes(xintercept=secY_proportion),color = "black") +http://127.0.0.1:43203/graphics/d34d13c5-dfed-4b50-af4f-0ba9aa7fb0e9.png
  
  
  # Labeling vertical lines
  
  #geom_text(aes(x=bamA_proportion-label_h_offset,y = 100, angle = label_angle, label = "BamA: 0.9815")) +
  #geom_text(aes(x=tamA_proportion-label_h_offset,y = 110, angle = label_angle, label = "TamA: 0.9889")) +
  #geom_text(aes(x=secY_proportion-label_h_offset,y = 120, angle = label_angle, label = "SecY: 0.9991")) 

  
hist_plot

