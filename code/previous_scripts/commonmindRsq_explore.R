library(ggplot2)
library(gridExtra)
library(scales)
library(gplots)
library(grid)
library(ggrepel)
library(viridis)  

load("~/Desktop/glasso2_ssize_test_p.RData")
x = ssize_test_p[,2]
y1 = -log10(ssize_test_p[,3])
y2 = -log10(ssize_test_p[,4])
plot(x,y1)
plot(x,y2)
cor(x,y1); cor(x,y2)
plot(x,ssize_test_p[,2])
plot(log(x),-log10(ssize_test_p[,3]))
cor(x^2,-log10(ssize_test_p[,3]))
plot(sqrt(x),-log10(ssize_test_p[,3]))
plot(x,-log10(ssize_test_p[,3]))

forgg = ssize_test_p[,-3]
forgg[,3] = -log10(forgg[,3])
is_sig = rep("", 44)
for(i in 1:44){
  if(forgg[i,3]>-log10(0.05)){
    is_sig[i] = "Significant"
  }else{
    is_sig[i] = "Non-significant"
  }
}
forgg["sig"] = is_sig
forgg$sig = ordered(forgg$sig, levels = c("Significant", 'Non-significant'))
names(forgg) = c("tissues", "ssize", "logp", "sig")
summary(lm(forgg[,3]~forgg[,2]))$coefficients[8]
pgg = ggplot(forgg, aes(x=ssize, y=logp)) 
pgg = pgg + theme_minimal() 
pgg = pgg + geom_point()
pgg = pgg + ylab('-log10 p-value') + xlab('sample size')
pgg = pgg + geom_label_repel(aes(ssize, logp, fill=sig, label = tissues),fontface = 'bold', color = 'white',box.padding = unit(0.25, "lines")) 
pgg = pgg + scale_fill_viridis(option="magma", discrete=TRUE, begin=0.2, end=0.9)
pgg = pgg + theme(legend.position = "bottom") 
pgg = pgg + theme(legend.text=element_text(size=6), legend.title=element_blank())
pgg = pgg + geom_smooth(method=lm)
pgg = pgg + annotate("text", x=100, y=-0.3, label="P = 4.59e-06", size=6)
pgg = pgg + ggtitle('Improvement in prediction accuracy VS sample size')
pgg
pdf("~/Desktop/MultiGExpr/tmp1.pdf", height = 6, width = 12)
pgg
dev.off()

forgg1 = data.frame(forgg[forgg[,2]<200,3],"Small")
forgg2 = data.frame(forgg[forgg[,2]>=200,3],"Large")
names(forgg1) = c("logp","group")
names(forgg2) = c("logp","group")
gg3 = rbind(forgg1, forgg2)
pgg1 = ggplot(gg3, aes(x=group, y=logp, color=group)) +
  theme_minimal() +
  geom_boxplot(outlier.shape=NA) + ylab('-log10 p-value')+ggtitle('Improvement in prediction accuracy VS sample size')+
  #theme(legend.position = "right",  axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  scale_colour_discrete(name  ="Group",breaks=c("Small", "Large"), labels=c("n < 200", "n >= 200"))
pgg1
pdf("~/Desktop/MultiGExpr/tmp2.pdf", height = 8, width = 6)
pgg1
dev.off()


load("~/Desktop/Rsq_stack_not_na.RData")
for(i in 1:4){
  tmp1 = data.frame(Rsq_stack_not_na[[i]][,1],"single")
  tmp2 = data.frame(Rsq_stack_not_na[[i]][,2],"multi")
  names(tmp1) = c("Rsq","group")
  names(tmp2) = c("Rsq","group")
  tmp0 = rbind(tmp1, tmp2)
  pgg1 = ggplot(tmp0, aes(x=group, y=Rsq, color=group)) +
    theme_minimal() +
    geom_boxplot(outlier.shape=NA) + ylab('Rsq')+ggtitle('Prediction accuracy in CommonMind')+
    scale_colour_discrete(name  ="Group",breaks=c("single", "multi"), labels=c("Single-tissue", "Cross-tissue"))
  pgg1
  pdf("~/Desktop/MultiGExpr/tmp3.pdf", height = 6, width = 4)
  pgg1
  dev.off()
  
}