# https://github.com/GuangchuangYu/hexSticker
if (!require('hexSticker')) install.packages('hexSticker')

library(hexSticker)
dire<- "C:/Users/moviedo/OneDrive - Universidade de Santiago de Compostela/GitHub/fda.usc.beta_2019_11_20/inst/script/"
imgurl <-paste0(dire,"DD5.png")
data(tecator)
ldat<-ldata(tecator$y,"x"=fdata.deriv(tecator$absorp.fdata))
class(ldat)


fil2<- "logo09.png"
fil <- fil2
pp<-paste0(dire,fil)
png(pp)
plot(ldat,var.name="Fat")
plot(ldat,col=rep("skyblue",len=nrow(ldat[[1]])))
dev.off()
imgurl <- pp
s2<-sticker(imgurl
        ,package=" fda.usc",p_size = 16,
        p_x=1, p_y=1.67,#, p_width=.6, p_height=.5,
        #s_x=1, s_y=1, s_width=.8, s_height=.8,
        s_x=1, s_y=.94, s_width=.67, s_height=.6,
        #h_color="#4E94B5",, h_fill="white", p_color = "#4E94B5",
        #h_color = "blue", h_fill="white", p_color = "blue",
        h_color = "skyblue", h_fill="white", p_color = "skyblue",
        filename=fil)
s2



sysfonts::font_add_google("Baloo Tammudu", "baloo")
sticker(ut,
        package = "slc rug", 
        p_size = 12.5,
        p_y = 1,
        p_family = "baloo",
        p_color = "#ffffff",
        h_fill = "#453781FF",
        h_color = "#ffffff",
        h_size = 2,
        s_x = 1, s_y = 1, 
        s_width = 1.6, s_height = 1.6,
        url = "github.com/slc-rug",
        u_color	= "#ffffff",
        u_size = 1.2,
        filename="slcrug_hex.png")


#  Functional Data Analysis and Utilities for Statistical Computing

imgurl <- "https://upload.wikimedia.org/wikipedia/commons/thumb/8/8c/Emojione_1F30A.svg/512px-Emojione_1F30A.svg.png"
s2<-sticker(imgurl, package="fda.usc",
            s_x=1, s_y=.8, s_width=.5, s_height=.5,
            h_color="#ff9966", h_fill="white", p_color = "#4E94B5",
            filename="zinbwave.png")
s2


library(ggplot2)
library(joineRML)
library(reshape2)
library(hexSticker)

data(pbc2)
pbc2 <- subset(pbc2, id != "228")

p <- ggplot(aes(x = year, y = log(serBilir)), data = pbc2) +
  geom_line(aes(group = id), alpha = 0.5) +
  geom_line(data = subset(pbc2, id == "11"), colour = "red") +
  geom_line(data = subset(pbc2, id == "96"), colour = "blue") +
  labs(x = "", y = "") +
  theme_void() +
  theme(strip.background = element_blank(),
        strip.text = element_blank())

s<-sticker(p, package = "fda.usc",
        p_size = 15,
        s_x = 1, s_y = .8, s_width = 1.55, s_height = 0.9,
        #p_color = "#FFFFFFDD",
        h_color = "#C74646",
        h_fill = "#46C7C7")
s






library(fda.usc)
data(tecator)
p<-plot(tecator[[1]])
x0<-tecator[[1]]
x<-x0$argvals
y<-x0$data
p <- ggplot(aes(x = x, y = x0), data = pbc2) +
  geom_line(aes(group = id), alpha = 0.5) +
  geom_line(data = subset(pbc2, id == "11"), colour = "red") +
  facet_wrap( ~ variable, scales = "free_y") +
  labs(x = "", y = "") +
  theme_void() +
  theme(strip.background = element_blank(),
        strip.text = element_blank())

library(reshape2)

#ggplot needs a dataframe
data<-y
data <- as.data.frame(data)
#id variable for position in matrix 
data$id <- 1:nrow(data) 
#reshape to long format
plot_data <- melt(data,id.var="id")

#plotx
ggplot(       geom_line(aes(x=x,y=y)))

b<-sticker(p, package = "fda.usc",
        p_size = 8,
        s_x = 1, s_y = .8, s_width = 1.5, s_height = 0.9,
        #p_color = "#FFFFFFDD",
        h_color = "#C76730",
        h_fill = "#3090C7")
b




library(ggplot2)
library(zoo)


