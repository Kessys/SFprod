#Boxplot da eficiencia por regiao
n <- length(dados$ly)
efh1 <- rank(mod.h1$efficiency)/n
efh2 <- rank(mod.h2$efficiency)/n

EFH   <- c("LIML" = efh1, "LIML" = efh2)
MVIC  <- rep("FIML",n)
MVIL  <- rep("LIML",n)
ef    <- c(MVIC, MVIL)
dataH <- data.frame("region" = rep(reg,2), EFH, ef)

setwd("C:/Users/kessy/OneDrive/Área de Trabalho/Mestrado UnB/Dissertação/artigo")
library(ggplot2); library(extrafont); library(ggpubr); library(farver)
#extrafont::loadfonts(device="win")
#font_import(paths = NULL, recursive = TRUE, prompt = TRUE,pattern = NULL)
loadfonts(device = "win")
windowsFonts(Times=windowsFont("TT Times New Roman"))

ggplot(dataH, aes(region, EFH, fill = region)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(fill = "white", color = "black") +
  scale_y_continuous(name = 'Technical efficiency normalized rank') +
  scale_x_discrete(name = 'Region', labels = c("Center-West", "North", "Northeast", "Southeast", "South")) + 
  scale_fill_discrete(name = 'Region') +
  facet_wrap(~ef) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, family = "Times New Roman"),
        text = element_text(size = 20, family = "Times New Roman", colour="black"),
        axis.text.x=element_text(size = 18, family="Times New Roman", colour="black"),
        axis.text.y=element_text(size = 18, family="Times New Roman", colour="black"))
ggsave("efH.png", units = 'in')


ggplot(dataH, aes(region, EFH, fill = region)) + 
  stat_boxplot(geom ='errorbar') +
  geom_boxplot() +
  scale_y_continuous(name = 'Technical efficiency normalized rank') +
  scale_x_discrete(name = 'Region', labels = c("Center-West", "North", "Northeast", "Southeast", "South")) + 
  scale_fill_discrete(name = 'Region') +
  facet_wrap(~ef) +
  theme_bw() +
  theme(plot.title = element_text(size = 16, family = "Times New Roman"),
        text = element_text(size = 18, family = "Times New Roman", colour="black"),
        axis.text.x=element_text(size = 16, family="Times New Roman", colour="black"),
        axis.text.y=element_text(size = 16, family="Times New Roman", colour="black"))
ggsave("efH.png", units = 'in')


