library(data.table)
library(dplyr)
library(ggplot2)



ghist <- function(vec) {
    vec <- mean_res
    df <- data.frame(
         x = vec
    )

    p <- ggplot(df, aes(x = x)) + 
        geom_histogram(aes(y=..density..), colour="black", fill="white")+
        geom_density(alpha=.2, fill="#FF6666") +
        geom_vline(aes(xintercept=mean(x)), color="blue", linetype="dashed", size=1)
    p


}