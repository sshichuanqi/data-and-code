library(openxlsx)
library(corrplot)
library(RColorBrewer)
library(brms)       
library(bayestestR)
library(performance)
library(sjPlot)
library(tidybayes)   
library(ggplot2)
library(ggdist)
library(tidyverse)
library(dplyr)
library(cowplot)

p=read.xlsx('points.xlsx')

#only keeping the rows with all values >= 0
p_selected=data.frame()
for (i in 1:nrow(p)){
  select<-p[i,]>=0
  if ('FALSE' %in% select==F){
    p_selected=rbind(p_selected,p[i,])
  }
}

#AGBD and influencing factors
f=p_selected[,7:20] 
f_scale=as.data.frame(scale(f))

#view pearson correlations among AGBD and its influencing factors
cor_mat=cor(f_scale,method='pearson') 
p_mat=cor.mtest(f_scale) 
mycol=colorRampPalette(colors=c("#d1495b","white","#00798c"))
corrplot.mixed(cor_mat,
               upper='color',upper.col=mycol(20),
               lower='number',lower.col='black',
               p.mat=p_mat$p,sig.level=0.05,insig="blank") #只显示显著的相关系数

corrplot(cor_mat,method='color',type='upper',col=mycol(20),
         addCoef.col="black", 
         #outline='gray40',
         p.mat=p_mat$p,sig.level=0.05,insig="blank",
         tl.pos='d',tl.col='black',tl.cex=1)

#pairs(f_scale,panel=panel.smooth)


###univariate linear regression

dat=data.frame(LST=f$LST,RND=f$RND,MAP=f$MAP,SOM=f$SOM,PAI=f$PAI,H=f$H,AGBD=f$AGBD)

plot_lst <- list()
point_color=c('#9DC4E7','#9DC4E7','#FFD966','#FFD966','#99B078','#99B078')
line_color=c('#2E75B6','#2E75B6','#F2B800','#F2B800','#385723','#385723')

for (i in 1:6) {
  plot<-ggplot(dat, aes_string(x=dat[,i], y='AGBD'))+ 
    geom_point(alpha=0.5,size=1.7,color=point_color[i])+
    labs(x=paste(colnames(dat)[i]),y=paste("AGBD"))+
    theme_classic()+
    theme(plot.margin=margin(t=0.9,r=0.9,b=0.9,l=0.9, unit="cm"),  
          panel.grid=element_blank(),
          axis.line=element_line(linewidth=1),
          axis.ticks=element_line(linewidth=1),
          axis.ticks.length=unit(-0.2,"cm"),
          axis.title.x=element_text(size=24,vjust=-0.5), 
          axis.text.x=element_text(size=22,vjust=-0.5),
          axis.title.y=element_text(size=24,angle=90,hjust=0.5),
          axis.text.y=element_text(size=22,hjust=0.5))
  
  #add fitted line and confidence interval if the regression coefficient is significant
  summary=summary(lm(AGBD~dat[,i],data=dat))
  if (summary$coefficients[2,4]<0.05){
    plot<-plot+stat_smooth(method="lm",formula=y~x,
                           se=T,color=line_color[i],linewidth=1.8) 
  }
  print(plot)
  plot_lst[[i]] <- plot
}

point=plot_grid(plotlist=plot_lst, ncol=3, nrow=2)
point
ggsave("point.png",width=18,height=12,dpi=600)

for (i in 1:6) {
  model.lm<-lm(AGBD~dat[,i],data=dat)
  l <- list(a = as.numeric(format(coef(model.lm)[1], digits = 4)),
            b = as.numeric(format(coef(model.lm)[2], digits = 4)),
            r2 = format(summary(model.lm)$r.squared, digits = 4),
            p = format(summary(model.lm)$coefficients[2,4], digits = 4))
  print(paste0("r2=",l$r2," ","p=",l$p))
}


###Bayesian multivariate regression 

brm_f=brm(AGBD~PAI+H+LST+RND+MAP+SOM,data=f_scale,
          family='gaussian',prior=prior(normal(0,1),class=b)) 

brm_inter=brm(AGBD~PAI*LST+H*LST+PAI*RND+H*RND+MAP+SOM,data=f_scale)

summary(brm_f)
bayes_R2(brm_f)

tab_model(brm_f)
tab_model(brm_inter)

plot(brm_f)
prior_summary(brm_f)
posterior_interval(brm_f,prob=0.95)
ci=posterior_interval(brm_inter,prob=0.80)
round(ci,2)

#draw posterior distributions

posterior_f=spread_draws(brm_f,Intercept,b_PAI,b_H,b_LST,b_RND,b_MAP,b_SOM,sigma)
posterior_f1=pivot_longer(posterior_f,cols=c(b_PAI,b_H,b_LST,b_RND,b_MAP,b_SOM),names_to="predictor") 
posterior_f1$predictor=factor(posterior_f1$predictor,levels=c('b_H','b_PAI','b_SOM','b_MAP','b_RND','b_LST'))

halfeye_f <-ggplot(data=posterior_f1,
                   aes(y=predictor,x=value,fill=predictor))+
  stat_halfeye(point_interval="median_qi", 
               #point_color=NA,            
               position=position_dodge(0.50),
               .width =c(0.95,0.95), #a vector of probabilities to use that determine the widths of the resulting intervals. 
               interval_size_range=c(0.7,0.7),
               normalize = "groups",trim=T
  )+
  #spikes for annotating stat_slabinterval() geometries.
  stat_spike(at=function(x)hdci(x,.width=0.80), #hdci:highest density credible interval
             size=0,linewidth=0.5)+
  # need shared thickness scale so that stat_slab and geom_spike line up
  scale_thickness_shared()

halfeye_f<-halfeye_f+ 
  expand_limits(x=c(-0.5,0.3))+
  scale_x_continuous(breaks=seq(-0.5, 0.3, 0.1))+ 
  labs(x='Standardized effects',y='Influencing factors')+  
  scale_y_discrete(name='Influencing factors',
                   labels=c('b_PAI'='PAI','b_H'='H','b_LST'='LST',
                            'b_RND'='RND','b_MAP'='MAP','b_SOM'='SOM'))+  
  geom_vline(xintercept=0,linetype ="dashed")+
  scale_fill_manual(values =c("b_LST"=alpha('#2E75B6',0.8),                              
                              "b_RND"=alpha('#2E75B6',0.8),                              
                              "b_MAP"=alpha('#F2B800',0.4),
                              "b_SOM"=alpha('#F2B800',0.8),
                              "b_PAI"=alpha('#385723',0.8),
                              "b_H"=alpha('#385723',0.4)))+
  theme_bw() +   
  theme(    
    legend.position="none",    
    panel.grid.major=element_blank(),     
    panel.grid.minor=element_blank(),
    panel.border=element_blank(),
    axis.line.x=element_line(),
    axis.ticks.x=element_line(),
    axis.line.y=element_line(),
    axis.ticks.y=element_line(),
    title=element_text(size=14,family="sans"),    
    text=element_text(size=14,family="sans",color='black'))

halfeye_f

#stacked barplot of the standardized effects of different groups of factors

coef_sum=0.16+0.03+0.4+0.04+0.03+0.08
disturb=0.44/coef_sum
envir=0.11/coef_sum
str=0.19/coef_sum
type=c('Disturbance','Environment','Canopy Structure')
porp=c(disturb,envir,str)
group=rep('A',3)
bar=data.frame(type,porp,group)
bar$type=factor(bar$type,levels=c('Disturbance','Environment','Canopy Structure'))

barplot<-ggplot(data=bar,aes(x=group,y=porp,fill=type))+
  geom_bar(stat='identity',position='fill')+
  geom_text(aes(label=type),
            position=position_fill(vjust=0.5), 
            size=3.6)+
  annotate(geom='text', 
           label='Adj.R^2=0.277',
           x='A',y=1.05,hjust=0.5,size=4)+
  scale_fill_manual(values=c('#9DC4E7','#FFD966','#99B078'))+
  theme_light()+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.y=element_line(),
        axis.text.y=element_text(size=12),
        axis.title.y=element_text(size=12,face='bold',vjust=1.8),
        axis.ticks.y=element_line(colour='black',linewidth=0.6),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.border=element_blank(),
        panel.spacing.y=unit(1,'cm'),
        legend.position='none',
        plot.margin=unit(c(1.5,0.5,1,-0.1),'cm'))+
  ylab('Relative effect of estimates (%)')+
  xlab('')+
  scale_y_continuous(expand=c(0.025,0),position="right")+
  scale_x_discrete(expand=c(0.2,0.3))

barplot

halfeye_bar=cowplot::plot_grid(halfeye_f,barplot,nrow=1,align='h',rel_widths=c(1,0.42))
halfeye_bar

ggsave("halfeye_bar.png",width=7,height=5,dpi=600)


###Bayesian SEM

#construct a disturbance-related composite variable
compositeD=bf(AGBD~LST+RND)
lmD=brm(compositeD, data=f_scale,control=list(adapt_delta=0.98),chains=4)
summary(lmD)
library(bayestestR)
bayestestR::hdi(lmD, ci=c(0.95, .99, .999))
f_scale$compositeD=scale(fixef(lmD)[2,1]*f_scale$LST + fixef(lmD)[3,1]*f_scale$RND)

#construct SEM model
agbd_mod=bf(AGBD~H*compositeD+PAI*compositeD+SOM)
h_mod=bf(H~SOM+compositeD) 
pai_mod=bf(PAI~SOM+compositeD) 
bsem_f<-brm(
  agbd_mod+h_mod+pai_mod+set_rescor(FALSE), 
  data=f_scale,
  #family='gaussian',prior=prior(normal(0,1),class=b),
  cores=parallel::detectCores()-1, 
  chains=4)       

summary(bsem_f)
LOO(bsem_f)
bayestestR::hdi(bsem_f, ci=c(0.95, .99, .999))
bayes_R2(bsem_f)

#draw direct and indirect effects
effects=data.frame(predictor=c('Disturbance','SOM','PAI','H','Disturbance','SOM'),
                   group1=c('Disturbance','Soil','CSC','CSC','Disturbance','Soil'),
                   group2=c(rep('Direct effect',4),rep('Indirect effect',2)),
                   coef=c(0.42,0.08,0.17,0.03,0.15*0.03+0.38*0.17,0.03*0.03+0.12*0.17))
direct=effects[effects$group2=='Direct effect',]
indirect=effects[effects$group2=='Indirect effect',]
direct=direct%>%mutate(predictor=factor(predictor,levels=c('H','PAI','SOM','Disturbance')))
indirect=indirect%>%mutate(predictor=factor(predictor,levels=c('SOM','Disturbance')))

direct_col=ggplot(direct,aes(predictor,coef))+
  coord_flip()+ 
  geom_col(fill=c('#9DC4E7','#FFD966','#99B078','#99B078'),width=0.9)+
  geom_text(aes(y=c(0.42-0.03,0.08-0.03,0.17-0.03,0.03+0.03),label=coef,size=10))+
  ggtitle('Direct effects')+
  ylim(0,0.42)+
  theme_gray()+
  theme(  
    #plot.margin=margin(0.5,0.5,0.5,0.5,unit='cm'),
    legend.position='none',
    panel.border=element_rect(linewidth=0.8), 
    plot.title=element_text(hjust=0.5,size=14,face="bold"),  
    axis.ticks.y=element_line(linewidth=0.8),
    axis.ticks.x=element_blank(),
    axis.text.y=element_text(size=12),
    axis.text.x=element_blank(),
    axis.title=element_blank())
direct_col
#facet_wrap(~group2,ncol=1, nrow=2)

indirect_col=ggplot(indirect,aes(predictor,coef))+
  coord_flip()+
  geom_col(fill=c('#9DC4E7','#FFD966'),width=0.9)+
  geom_text(aes(y=coef+0.03,label=round(coef,2),size=10))+
  ggtitle('Indirect effects')+
  ylab('Standardized effects')+
  ylim(0,0.42)+
  theme_gray()+
  theme( 
    #plot.margin=margin(0.5,0.5,0.5,0.5,unit='cm'),
    legend.position='none',
    panel.border=element_rect(linewidth=0.8),
    plot.title=element_text(hjust=0.5,size=14,face="bold"),  
    axis.ticks.y=element_line(linewidth=0.8),
    axis.ticks.x=element_line(linewidth=0.8),
    axis.text.y=element_text(size=12),
    axis.text.x=element_text(size=12),
    axis.title.y=element_blank(),
    axis.title.x=element_text(hjust=0.5,size=14))
indirect_col

col=cowplot::plot_grid(direct_col,indirect_col,ncol=1,rel_heights=c(1,0.7))
col

ggsave("col.png",width=4,height=5,dpi=600)
