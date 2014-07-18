library(ggplot2)
library(reshape)

setwd("from_axon/tetanus_dat")
soma_v=read.table("soma_v.dat", skip=2)
names(soma_v)<-c("time","soma_v")
gsoma_v<-ggplot(data=soma_v,aes(x=time,y=soma_v))+geom_line()+theme_bw()+
    xlab("time (ms)")+ylab("")

vsyn0<-read.table("v_syn[0].dat", skip=5000, nrow=18000, col.names=c("time", "v_syn"))
vsyn1<-read.table("v_syn[1].dat", skip=5000, nrow=18000, col.names=c("time", "v_syn"))
vsyn2<-read.table("v_syn[2].dat", skip=5000, nrow=18000, col.names=c("time", "v_syn"))
vsyn3<-read.table("v_syn[3].dat", skip=5000, nrow=18000, col.names=c("time", "v_syn"))
vsyn4<-read.table("v_syn[4].dat", skip=5000, nrow=18000, col.names=c("time", "v_syn"))
vsyn5<-read.table("v_syn[5].dat", skip=5000, nrow=18000, col.names=c("time", "v_syn"))
vsyn6<-read.table("v_syn[6].dat", skip=5000, nrow=18000, col.names=c("time", "v_syn"))
vsyn7<-read.table("v_syn[7].dat", skip=5000, nrow=18000, col.names=c("time", "v_syn"))

#v_syn<-rbind(vsyn0, vsyn1, vsyn2, vsyn3, vsyn4, vsyn5, vsyn6, vsyn7)
#flags<- rep(0:7, each=nrow(vsyn0), len=nrow(v_syn))
#v_syn$flags<-flags
somav_tmp<-soma_v[4999:22998,]
names(somav_tmp) <- c("time", "v_syn")
zz<-melt(list(syn1=vsyn1, syn5=vsyn5, syn7=vsyn7, soma=somav_tmp), id.vars="v_syn")
d<- ggplot(zz, aes(value, v_syn, color=L1)) + theme_bw() +geom_line(size=0.7)
colorscal<-c("#a6611a",
    "#dfc27d",
    "#80cdc1",
    "#018571")
d<-d+scale_color_manual(values=colorscal)

casyn0<-read.table("ca_acc_syn[0].dat", skip=2, col.names=c("time", "ca_syn"))
casyn1<-read.table("ca_acc_syn[1].dat", skip=2, col.names=c("time", "ca_syn"))
casyn2<-read.table("ca_acc_syn[2].dat", skip=2, col.names=c("time", "ca_syn"))
casyn3<-read.table("ca_acc_syn[3].dat", skip=2, col.names=c("time", "ca_syn"))
casyn4<-read.table("ca_acc_syn[4].dat", skip=2, col.names=c("time", "ca_syn"))
casyn5<-read.table("ca_acc_syn[5].dat", skip=2, col.names=c("time", "ca_syn"))
casyn6<-read.table("ca_acc_syn[6].dat", skip=2, col.names=c("time", "ca_syn"))
casyn7<-read.table("ca_acc_syn[7].dat", skip=2, col.names=c("time", "ca_syn"))
