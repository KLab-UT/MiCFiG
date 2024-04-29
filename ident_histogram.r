#install.packages("ggplot2",repos = "http://cran.us.r-project.org")
library("ggplot2")

stats <- read.csv("identities.csv",header = FALSE)
genes <- stats[1]
pident <- as.numeric(unlist(stats[2]))

mean_ident <- mean(pident)
mean_label = paste("mean ",round(mean_ident),"%",sep = "")

max_ident <- max(pident)
max_label = paste("max ",round(max_ident),"%",sep="")
x_label_offset = -30
y_label_offset = -0.0005
label_angle = 0 #90


df <- data.frame(genes,pident)


#####################################
# Density Plot For Percent Identity #
#####################################

densityplot <- ggplot(df) +
  theme_classic() +
  geom_density(
    aes(x=pident,y=after_stat(density)),
    color = "darkblue",
    fill="lightblue",
    linewidth = 1)+
  
  geom_vline(aes(xintercept=mean_ident), color = "black",linetype="dotted") +
  geom_vline(aes(xintercept=max_ident), color = "red",linetype="dotted")+

  geom_text(aes(x=mean_ident - x_label_offset,y = y_label_offset, angle=label_angle,label= mean_label))+
  geom_text(aes(x=max_ident + x_label_offset,y = y_label_offset, angle=label_angle,label= max_label))

  
densityplot

##################################
# Histogram For Percent Identity #
##################################

histplot <- ggplot(df) +
  theme_classic() +
  geom_histogram(
    aes(x=pident),
    binwidth=8,
    color="darkblue",
    fill="lightblue")+
  
  geom_vline(aes(xintercept=mean_ident), color = "black",linetype="dotted") +
  geom_vline(aes(xintercept=max_ident), color = "red",linetype="dotted")+
  
  geom_text(aes(x=mean_ident - x_label_offset,y = -20, angle=label_angle,label= mean_label))+
  geom_text(aes(x=max_ident + x_label_offset,y = -20, angle=label_angle,label= max_label))+
  
  xlab("Percent Identity")+
  ylab("Number of Hits")


histplot

