Setwd ("D:/data")
Getwd () 
list.files ("GTEx_Analysis_v8_ASE_WASP_counts_by_subject")
x <- list.files ("GTEx_Analysis_v8_ASE_WASP_counts_by_subject")
dir <- paste ("./GTEx_Analysis_v8_ASE_WASP_counts_by_subject/",x,sep="")
n <- length (dir) 
获得文件长度 838
r.data <- read.table(file=dir[1],header=TRUE,sep=",") 
r.data<-r.data[which(r.data$tissue%in% c('ENSG00000106258','ENSG00000160868','ENSG00000160870')),]
r.data <-r.data[which9result.data$gene_id %in% c('LIVER','SRTTRM')),]
for(i in 2:n){new.data <- read.table(file=dir[i], header=TRUE, sep=",")
e.data<-e.data[which(r.data$tissue%in% c('ENSG00000106258','ENSG00000160868','ENSG00000160870')),]
e.data <-e.data[which(r.data$gene_id %in% c('LIVER','SNTTRM')),]
r.data <-rbind(r.data,e.data) }
list1<-result.data [which(result.data$tissue=="LIVER"),]
list1<-result.data [which(result.data$tissue=="SNTTRM"),]
list11<-list1 [which(list1$GENE_ID=="ENSG00000106258"),]
list12<-list1 [which(list1$GENE_ID=="ENSG00000160868"),]
list13<-list1 [which(list1$GENE_ID=="ENSG00000160870"),]
list21<-list2 [which(list2$GENE_ID=="ENSG00000106258"),]
list22<-list2 [which(list2$GENE_ID=="ENSG00000160868"),]
list23<-list2 [which(list2$GENE_ID=="ENSG00000160870"),]
list11$BINOMP[list11$BINOM_P>=0.05]=">=0.05"
list11$BINOMP[list11$BINOM_P<0.05]="<0.05"
list12$BINOMP [list12$BINOM_P>=0.05]=">=0.05"
list11$BINOMP [list11$BINOM_P<0.05]="<0.05"
list13$BINOMP [list13$BINOM_P>=0.05]=">=0.05"
list13$BINOMP [list13$BINOM_P<0.05]="<0.05"
list21$BINOMP [list21$BINOM_P>=0.05]=">=0.05"
list21$BINOMP [list21$BINOM_P<0.05]="<0.05"
list22$BINOMP [list22$BINOM_P>=0.05]=">=0.05"
list22$BINOMP [list22$BINOM_P<0.05]="<0.05"
list23$BINOMP [list23$BINOM_P>=0.05]=">=0.05"
list23$BINOMP [list23$BINOM_P<0.05]="<0.05"
这里对于大于0.05的和小于0.05的位点进行划分。
17
ggplot(list21,aes(x=list21$REF_COUNT,y=list21$ALT_COUNT,color=BINOMPSCALE))+geom_point(alpha=.3)+labs(x="SI3 REF",
y="SI4 ALT")+theme_bw()+ geom_smooth(method="lm")
ggplot(list12,aes(x=list12$REF_COUNT,y=list12$ALT_COUNT,color=factor(BINOMPSCALE)))+geom_point(alpha=.3)+labs(x="
CYP3A4 in liver:REF",y="CYP3A4 in liver:ALT")+=theme_bw()
（重复代码，更换其中的子列表即可）我们可以作出肝脏和小肠中的CYP3A4\5\7的reference和alternative的count数的相对关系的
散点图。 
代码2
我们从所有的肝脏和小肠中筛选出那些属于stoplost synoymoos和missense varient 3’ 和5’的突变位点，这些位点都是潜在的
marker。
list1 < -read.table("D:/data/list1.csv",header=TRUE,sep=",") 
Ls <- list1[which( list1$VARIANT_ANNOTATION %in% c('stop_lost',' synonymous_variant' , 
'3_prime_UTR_variant','5_prime_UTR_variant','missense_variant')),]
library(ggplot2)
ls1 <- ls[which(ls$GENE_ID=="ENSG00000106258"),]
ls2 <- ls[which(ls$GENE_ID=="ENSG00000160868"),]
ls3 <- ls[which(ls$GENE_ID=="ENSG00000160870"),]
ls1$BINOMP[ls1$BINOM_P>=0.05]=">=0.05"
ls1$BINOMP[ls1$BINOM_P<0.05]="<0.05"
ls2$BINOMP[ls2$BINOM_P>=0.05]=">=0.05"
ls2$BINOMP[ls2$BINOM_P<0.05]="<0.05"
ls3$BINOMP[ls3$BINOM_P>=0.05]=">=0.05"
ls3$BINOMP[ls3$BINOM_P<0.05]="<0.05"
BINOMPSCALE1<-factor(ls1$BINOMP)
p1 <- ggplot(ls1,aes(x=ls1$REF_COUNT,y=ls1$ALT_COUNT,color=BINOMPSCALE1))+geom_point(alpha=.7)+labs(x="cyp3a5 
ref_c",y="cyp3a5 alt_c")+theme_bw()+scale_color_manual(values = c("red", "blue"))
BINOMPSCALE2 <- factor(ls2$BINOMP)
p2 <- ggplot(ls2,aes(x=ls2$REF_COUNT,y=ls2$ALT_COUNT,color=BINOMPSCALE2))+geom_point(alpha=.7)+labs(x="cyp3a4 
ref_c",y="cyp3a4 alt_c")+theme_bw()+scale_color_manual(values = c("red", "blue"))
p2
BINOMPSCALE3<-factor(ls3$BINOMP)
p3<-ggplot(ls3,aes(x=ls3$REF_COUNT,y=ls3$ALT_COUNT,color= BINOMPSCALE3))+geom_point(alpha=.7)+labs(x="cyp3a7 ref_c",
y="cyp3a7 alt_c")+theme_bw()+ scale_color_manual(values = c("red", "blue"))
p3
重复之前的作图步骤。
代码3
storea<-c(rep(0,18))
for (i in 1:18)
{for (j in 1:length(names(ud3a4)))
{if(aa4[1,i]==names(ud3a4)[j]){storea[i]<-j}}}
18
ud3a4s<-ud3a4[storea,]
storeb<-c(rep(0,109))
for (i in 1:109)
{for (j in 1:length(names(ud3a5)))
{if(aa5[1,i]==names(ud3a4)[j]){storeb[i]<-j}}}
ud3a5s<-ud3a5[storeb,]
storec<-c(rep(0,117))
for (i in 1:117)
{for (j in 1:length(names(ud3a7)))
{if(aa7[1,i]==names(ud3a7)[j]){storec[i]<-j}}}
ud3a7s<-ud3a7[storea,]
ud3a5s<-as.matrix(ud3a5s)
ud3a4s<-as.matrix(ud3a4s)
ud3a7s<-as.matrix(ud3a7s)
Aqpvalue1<-c(rep(0,6747))
Aqpvalue2<-c(rep(0,6747))
Aqpvalue3<-c(rep(0,6726))
for(i in 1:6747)
{line1<-data.frame(aein=factor(ud3a4s[6748,],levels=c("AEI","NON")),homn=factor(ud3a4s[i,],levels=c("HET","HOM")))
aqpvalue<-chisq.test(table(line1)
Aqpvalue1[i]<-aqpvalue$p.value}
Aqpvalueframe1<-data.frame(Aqpvalue1,ud3a4$POS)
for(i in 1:6747)
{line1<-data.frame(aein=factor(ud3a5s[6748,],levels=c("AEI","NON")),homn=factor(ud3a5s[i,],levels=c("HET","HOM")))
aqpvalue<-chisq.test(table(line1)
Aqpvalue2[i]<-aqpvalue$p.value}
Aqpvalueframe1<-data.frame(Aqpvalue1,ud3a4$POS)
Aqpvalueframe2<-data.frame(Aqpvalue2,ud3a4$POS)
for(i in 1:6726)
{line1<-data.frame(aein=factor(ud3a4s[6727,],levels=c("AEI","NON")),homn=factor(ud3a7s[i,],levels=c("HET","HOM")))
aqpvalue<-chisq.test(table(line1)
Aqpvalue3[i]<-aqpvalue$p.value}
Aqpvalueframe3<-data.frame(Aqpvalue3,ud3a7$POS)
用excel整理为如正文中的形式
第一栏为p值，第二栏为基因hg38坐标下的位置，第三栏为位于的组织。然后我们画出p值的负对数和位置的点图，这样我们方便看到
顺式调控因子相对于CYP3A4\5\7基因的位置。
graph4<-read.table("D:/data/4.csv",header=TRUE,sep=",")
graph5<-read.table("D:/data/5.csv",header=TRUE,sep=",")
graph7<-read.table("D:/data/7.csv",header=TRUE,sep=",")
19
graph4r<-graph4[!is.na(graph4$pval),]
graph5r<-graph5[!is.na(graph5$pval),]
graph7r<-graph7[!is.na(graph7$pval),]
graph4r$logp=(-log10(graph4r$pval))
graph5r$logp=(-log10(graph5r$pval))
graph7r$logp=(-log10(graph7r$pval))
CYP3A5LOG<-ggplot(graph5r,aes(x=graph5r$POS,y=graph5r$logp,color=class5))+labs(x="pos",y="-logpvalue")+theme(axis.
line = element_line(arrow = arrow(length = unit(0.5, 'cm'))))+geom_segment(aes(x=99679998,xend=99648194,y=-0.1,
yend=-0.1),color="black", size=1.1, arrow = arrow(type="closed",length=unit(0.1, "cm") ) ) + annotate("text", label 
= "CYP3A5", x = 99664096, y = -log10(5)+0.5, size=2.5, colour = "black",fontface="bold")+geom_point(alpha=.7)
+scale_color_manual(values = c("red", "blue"))
CYP3A4LOG=ggplot(graph4r,aes(x=graph4r$POS,y=graph4r$logp,color=class4))+labs(x="pos",y="-logpvalue")+theme(axis.
line = element_line(arrow = arrow(length = unit(0.5, 'cm'))))+geom_segment(aes(x=99784248,xend=99756960,y=-0.05,
yend=-0.05),color="black", size=1.1, arrow = arrow(type="closed",length=unit(0.1, "cm") ) ) + annotate("text", 
label = "CYP3A4", x =99770604, y = -log10(4)+0.5, size=2.5, colour = "black",fontface="bold")+geom_point(alpha=.7)
+scale_color_manual(values = c("red", "blue"))
CYP3A7LOG=ggplot(graph7r,aes(x=graph7r$POS,y=graph7r$logp,color=class7))+labs(x="pos",y="-logpvalue")+theme(axis.
line = element_line(arrow = arrow(length = unit(0.5, 'cm'))))+geom_segment(aes(x=99735196,xend=99684957,y=-0.05,
yend=-0.05),color="black", size=1.1, arrow = arrow(type="closed",length=unit(0.1, "cm") ) ) + annotate("text", 
label = "CYP3A7", x = 99710077, y = -log10(4)+0.5, size=2.5, colour = "black",fontface="bold")+geom_point(alpha=.7)
+scale_color_manual(values = c("red", "blue"))
代码4（线性拟合）
fit1<-lm(x11$ALT_COUNT~x11$REF_COUNT)
summary(fit1)
代码5（网页数据爬取）
首先是根据位置搜索SNP的网址，然后是根据SNP搜索组织
加载xml2，rvest,dplyr
web<-read_html(DATA5N$SNP_HTML[i],encoding="utf-8")
w1<-web %>% html_nodes("#main_content") %>% html_nodes("main") %>% html_nodes("div.summary-box.usa-grid-full") %>% 
html_nodes("dl:nth-child(2)") %>% html_nodes("dd:nth-child(4)") %>% html_nodes("div") %>% html_text()
w2<-web %>% html_nodes("#main_content") %>% html_nodes("main") %>% html_nodes("div.summary-box.usa-grid-full") %>% 
html_nodes("dl:nth-child(2)") %>% html_nodes("dd:nth-child(6)") %>% html_text()
unlist然后合并为data frame输出。
在excel表中，按照cit从大到小排列，引用数在2以上的SNP可以进行重点关注。
