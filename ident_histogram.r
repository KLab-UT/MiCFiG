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

make_density <- function(){
  densityplot <- ggplot(df) +
    theme_classic() +
    geom_density(
      aes(x=pident,y=after_stat(density)),
      color = "darkblue",
      fill="lightblue",
      linewidth = 1)+
    
    geom_boxplot(
      aes(x=pident),
      width=0.0005,
      fill="#f8fbae",
      position=position_nudge(y=-0.00025))+
    
    geom_vline(aes(xintercept=mean_ident), color = "black",linetype="dotted") +
    geom_vline(aes(xintercept=max_ident), color = "red",linetype="dotted")+
    
    geom_text(aes(x=mean_ident - x_label_offset,y = y_label_offset, angle=label_angle,label= mean_label))+
    geom_text(aes(x=max_ident + x_label_offset,y = y_label_offset, angle=label_angle,label= max_label))
  
  densityplot
}


##################################
# Histogram For Percent Identity #
##################################

make_hist <- function(){
  
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
}



make_scatter <- function(){
  
  scatter <- ggplot(df,aes(x=pident,y=0))+

    geom_point(position = position_jitter(width = 10),
               aes(colour = cut(pident,c(-Inf,100,Inf))),
               shape=1,
               )+

    xlab("Percent Identity")+
    ylab("")+
    
    
    scale_color_manual(name="Spicyness",
                       values=c("(-Inf,100]"="cornflowerblue",
                                "(100,Inf]" = "indianred1"),
                       
                       labels=c("0-100",">100"))+
    
    theme(axis.ticks.y=element_blank(),
            axis.text.y=element_blank()
          )
  scatter
    }

#make_hist()
#make_density()
make_scatter()


