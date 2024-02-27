library(stringr)
library(readr)
library(installr)
library(utils)

load("Mitelman_summary_220204.RData")

N_variations <- summary$N_SVs
result <- summary$result
N_sample <- length(summary$samples)
empty1 <- summary$no_variation_before
empty2 <- summary$no_variation_after
error <- summary$error 
unc <- summary$uncertain
samples <- summary$samples

#####
error_index <- unique(error$sample_ID)

unc_ring <- grep(pattern = "ring", x = unc$error)
unc_marker <- grep(pattern = "marker", x = unc$error)
unc_error <- unc[-unc_ring,]
unc_error <- unc[-unc_marker,]
unc_error_index <- unique(unc_error$sample_ID)

samples_error_index <- union(error_index,unc_error_index)

samples_excluded <- union(union(empty1, empty2), samples_error_index)

###### find errors
samples_var <- samples[-empty1]
N_SV <- N_variations[N_variations != 0]
all_vars <- sum(N_SV)
 grep(pattern = "incorrect", x = error$error)
chr_sex_error <- error[which(error$Variation_ID == 0),1]
#########variation stats
all_captured <- all_vars- nrow(unc) - nrow(error) - nrow(remaining)
all_captured <- all_captured+80


##########
error[9662,1]
samples[17249]
samples[22258]
####delete cloned_ID with 1 clone
summary$cloned_ID <- summary$cloned_ID[-which(summary$cloned_ID == 13305)]
summary$cloned_ID <- summary$cloned_ID[-which(summary$cloned_ID == 16744)]
summary$cloned_ID <- summary$cloned_ID[-which(summary$cloned_ID == 69280)]
N_cloned <- length(summary$cloned_ID)

NN <- sum(summary$cloned_NO)-length(summary$cloned_NO)+73963

clon2 <- summary$cloned_NO[which(summary$cloned_NO == 2)]
clon3 <- summary$cloned_NO[which(summary$cloned_NO == 3)]
clon4 <- summary$cloned_NO[which(summary$cloned_NO == 4)]
clon5 <- summary$cloned_NO[which(summary$cloned_NO == 5)]
clon6 <- summary$cloned_NO[which(summary$cloned_NO == 6)]
clon7 <- summary$cloned_NO[which(summary$cloned_NO == 7)]
clon8 <- summary$cloned_NO[which(summary$cloned_NO == 8)]
clon9 <- summary$cloned_NO[which(summary$cloned_NO == 9)]
clon10 <- summary$cloned_NO[which(summary$cloned_NO == 10)]
clon11 <- summary$cloned_NO[which(summary$cloned_NO == 11)]
clon12 <- summary$cloned_NO[which(summary$cloned_NO == 12)]
clon13 <- summary$cloned_NO[which(summary$cloned_NO == 13)]
clon14 <- summary$cloned_NO[which(summary$cloned_NO == 14)]
clon15 <- summary$cloned_NO[which(summary$cloned_NO == 15)]
clon16 <- summary$cloned_NO[which(summary$cloned_NO == 16)]

########pie chart for cloned karyotypes
library(ggplot2)
library(scales)
no_clone = 73963-N_cloned

df <- data.frame(
  group = c("1 clone", "2 clone", "+3 clone"),
  value = as.numeric(c(no_clone, length(clon2), N_cloned-length(clon2))) 
)

bp <- ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")

pie <- bp + coord_polar("y", start=0)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

pie <- pie + scale_fill_brewer("Blues") + blank_theme +
  theme(axis.text.x=element_blank())+
  geom_text(aes(y = value-3000, 
                label = value), size=5)+
  ggtitle("Number of clones in Mitelman karyotypes")+ theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill=guide_legend(title="clone numbers"))

 pie +  theme(plot.title = element_text(hjust = 0.5))
###### karyotype-clone bar chart

df2 <- data.frame(samples=c("Karyotypes", "Samples"),
                 len=c(73963, N_sample))

p<-ggplot(df2, aes(x=samples, y=len, color=samples)) +
  geom_bar(stat="identity", fill="white")

p<-ggplot(df2, aes(x=samples, y=len, fill=samples)) +
  geom_bar(stat="identity")+theme_minimal()

p<-p+ labs(title="Samples in Mitelman karyotypes", 
         y = "Number")
p <- p+scale_fill_brewer(palette="Blues")

p <- p + theme(plot.title = element_text(hjust = 0.5))

p <- p + scale_x_discrete(name ="Samples", 
                 limits=c("Karyotypes", "Samples"))


p <- p + scale_y_continuous(name="Number", limits=c(0, 100000))


########## plot empty samples ######
no_var_ini <- length(summary$no_variation_before)
no_var_end <- length(summary$no_variation_after)
with_var <- length(summary$samples) - no_var_ini - no_var_end

df3 <- data.frame(samples=c("Karyotype", "Sample", "Sample", "Sample"),
                  len=c(73963, (no_var_ini+43), with_var, (no_var_end-43)),
                  variation=c("karyotypes", "0 variation initially", "+1 valid variation(s)", "0 valid variation"),
                  cum_no <- c(36000, 95000, 40000, 88000))

p <- ggplot(data=df3, aes(x=samples, y=len, fill=reorder(variation, len))) +
  geom_bar(stat="identity")+
  geom_text(aes(y=cum_no, label=len), vjust=1.6, 
            color="black", size=3.5)+
  #scale_fill_brewer(palette="Paired")+
  theme_minimal()+labs(title="Samples in Mitelman karyotypes", x="Samples", y = "Number")+ theme(plot.title = element_text(hjust = 0.5))

p + scale_x_discrete(name ="Samples", limits=c("Karyotype", "Sample")) + scale_y_continuous(name="Number", limits=c(0, 100000)) + 
  scale_fill_manual(values=c('coral4','pink3', "cadetblue4", "lavenderblush3")) +
  guides(fill=guide_legend(title="variation numbers"))
#aquamarine4 lightblue4 darkseagreen4 indianred coral4 lightpink4  cadetblue4
#999999

######### only valid samples
dff <- data.frame(samples=c("+1 valid variation", "no invalid variation"),
                  len=c(82227,62439),
                  cum_no <- c(40000, 35000))

p <- ggplot(data=dff, aes(x=samples, y=len, fill=samples)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=cum_no, label=len), vjust=1.6, 
            color="black", size=3.5)+
  #scale_fill_brewer(palette="Paired")+
  theme_minimal()+labs(title="Samples in Mitelman karyotypes", x="Samples", y = "Number")+ theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c("no invalid variation"="lightsalmon", "+1 valid variation"="lavenderblush3")) 


p + scale_x_discrete(name ="Samples", limits=c("+1 valid variation", "no invalid variation")) + scale_y_continuous(name="Number", limits=c(0, 100000)) + 
 p + scale_fill_manual(values=c("no invalid variation"="lightsalmon", "+1 valid variation"="lavenderblush3")) 

#######find samples with one clone in cloned_ID
ttet <- list()
N_clones <- list()
for(i in 1:N){
  # if(grepl(pattern = "\\/", x = samples[i]) == FALSE){
  if(grepl(pattern = "\\/", x = samples[i]) == TRUE){
    ttet <- c(ttet, i)
    clones <- unlist(strsplit(samples[i],"/"))
    N_clones <- c(N_clones,length(clones))
}
}

##########fake karyotypes ############
fake_karyotypes <- c("47,XY,+12,add(17)(q?),der(7)t(5;7)(p11;q23)t(7;9)(p21;q34)", "45,XX,der(15)t(1;15)(p10;q11),del(14)(q33),+5mar", "46,XX,inv(5),t(4;20)(q13-q35;q13),t(9;22)(q34;q11)",
                     "44,XX,i(22)(p12),der(11)t(11;14)(q11;q11)t(13;14)(q11;q32),r(7)","43,XY,del(4)(q1?4q2?1),del(10)(q13)-21,del(14)","45,XX,t(2;13)x2,del(20)(q13),i(8)(q10)",
                     "46,XY,der(1)t(1;1)(p36;q23),der(5)t(5;20),t(3;8)(p21;q12)","49,XY,+8,ins(X11)(q22;p15p13),+r(21)", "48,XX,inv(10(p13q21),+1,der(10)t(10;17)",
                     "46,XY,del(20)(q11),der(1)t(1;5)(q36;q11)", "47,XX,t(8;9),r(1)(p36q44),der(9)t(9;10)(q34;q22) or der(9)t(9;10)(q22;q11)","46,XX,inv(3)(p14q28),dup(140)(q24q32)",
                     "50,XY,+1,del(6)(q2?3),+mar","49,+4,-12,der(22)t(9;22)","48,XY,dup(1),add(21)(q2,tas(14;15)(p11;p11)")


########pie chart for variations


df4 <- data.frame(
  group = c("error(37406)", "uncertain(16388)", "not captured(7919)", "captured(269393)"),
  value = as.numeric(c(nrow(error), (nrow(unc)-80), nrow(remaining), all_captured)),
  lab_val <- c(44000, 20000, 5000,215000),
  ord <- c(1,3,4,1)
)

bp <- ggplot(df4, aes(x="", y=value,  fill= reorder(group, ord)))+
  geom_bar(width = 1, stat = "identity")

pie <- bp + coord_polar("y", start=0)

pie <- pie + 
  theme(axis.text.x=element_blank())+
  geom_text(aes(y = lab_val,label = paste(round(value / sum(value) * 100, 1), "%")), size=3)+
  ggtitle("Number of variations in Mitelman karyotypes")+ theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c('darkseagreen4','tomato3', "khaki2", "lavenderblush3")) + theme_void() +
  guides(fill=guide_legend(title="variation types (#)"))

pie +  theme(plot.title = element_text(hjust = 0.5))


########pie chart for not captured


df5 <- data.frame(
  group = c("der with no info", "inc", "dmin", "typo"),
  value = as.numeric(c(3996, 3110, 634, 177)),
  lab_val <- c(5900, 2600, 600, 100),
  ord <- c(1,2,3,4)
)

bp <- ggplot(df5, aes(x="", y=value,  fill= reorder(group, ord)))+
  geom_bar(width = 1, stat = "identity")

pie <- bp + coord_polar("y", start=0)

pie <- pie + 
  theme(axis.text.x=element_blank())+
  geom_text(aes(y = lab_val,label = value), size=4)+
  ggtitle("Not-captured variation types")+ theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c('pink3','slategray2', "darkseagreen4", "tomato3")) + theme_void() +
  guides(fill=guide_legend(title="variation types"))
pie
pie +  theme(plot.title = element_text(hjust = 0.5))

####### classify variations ####### 
Bands <- read.csv("Bands.csv")

samples2 <- samples[-c(empty1,empty2)]


plus_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "\\+", x = result[,((i-1)*3+5)]) == FALSE) == FALSE){
    plus <- list()  
    plus <- unique(Bands$chr[which(grepl(pattern = "\\+", x = result[,((i-1)*3+5)]) == TRUE)])
    #plus <- which(grepl(pattern = "\\+", x = result[,((i-1)*3+5)]) == TRUE)[1]
    for(j in 1:length(plus)){
      plus_list <- rbind(plus_list, c(i,plus[j], str_count(result[which(result$chr==plus[j])[1],((i-1)*3+5)], "\\+")))
    }
  } 
  cat(sprintf("sample number %s", i), sep="\n")
}

minus_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "\\-", x = result[,((i-1)*3+6)]) == FALSE) == FALSE){
    minus <- list()  
    minus <- unique(Bands$chr[which(grepl(pattern = "\\-", x = result[,((i-1)*3+6)]) == TRUE)])
    for(j in 1:length(minus)){
      minus_list <- rbind(minus_list, c(i,minus[j], str_count(result[which(result$chr==minus[j])[1],((i-1)*3+6)], "\\-")))
    }
  }  
  cat(sprintf("sample number %s", i), sep="\n")
}

tr_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "t", x = result[,((i-1)*3+7)]) == FALSE) == FALSE){
      tr <- list()
      tr <- which(grepl(pattern = "t", x = result[,((i-1)*3+7)]) == TRUE)
      for(j in 1:length(tr)){
        if(grepl(pattern = "t", x = result[tr[j],((i-1)*3+5)]) == TRUE){
          tr_list <- rbind(tr_list, c(i,tr[j],"gain", str_count(result[tr[j],((i-1)*3+5)], "t")))
        } else { if(grepl(pattern = "t", x = result[tr[j],((i-1)*3+6)]) == TRUE){
          tr_list <- rbind(tr_list, c(i,tr[j],"loss", str_count(result[tr[j],((i-1)*3+6)], "t")))
        } else {
          tr_list <- rbind(tr_list, c(i,tr[j],"balanced", str_count(result[tr[j],((i-1)*3+7)], "t")))
        }
          }
      }
  }
  cat(sprintf("sample number %s", i), sep="\n")
}

######################
tr_gain_list <- list()
tr_loss_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "t", x = result[,((i-1)*3+5)]) == FALSE) == FALSE){
    tr_gain <- list()
    tr_gain <- which(grepl(pattern = "t", x = result[,((i-1)*3+5)]) == TRUE)
    for(j in 1:length(tr_gain)){
      tr_gain_list <- rbind(tr_gain_list, c(i,tr_gain[j],"gain", str_count(result[tr_gain[j],((i-1)*3+5)], "t")))
    }
  }
  if(all(grepl(pattern = "t", x = result[,((i-1)*3+6)]) == FALSE) == FALSE){
    tr_loss <- list()
    tr_loss <- which(grepl(pattern = "t", x = result[,((i-1)*3+6)]) == TRUE)
    for(k in 1:length(tr_loss)){
      tr_loss_list <- rbind(tr_loss_list, c(i,tr_loss[k],"loss", str_count(result[tr_loss[k],((i-1)*3+6)], "t")))
    }
  }
  cat(sprintf("sample number %s", i), sep="\n")
}

#######################        
tr_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "t", x = result[,((i-1)*3+7)]) == FALSE) == FALSE){
    tr <- list()
    tr <- which(grepl(pattern = "t", x = result[,((i-1)*3+7)]) == TRUE)
    for(j in 1:length(tr)){
      if(grepl(pattern = "t", x = result[tr[j],((i-1)*3+5)]) == TRUE){
        tr_list <- rbind(tr_list, c(i,tr[j],"gain", str_count(result[tr[j],((i-1)*3+5)], "t")))
      } else { if(grepl(pattern = "t", x = result[tr[j],((i-1)*3+6)]) == TRUE){
        tr_list <- rbind(tr_list, c(i,tr[j],"loss", str_count(result[tr[j],((i-1)*3+6)], "t")))
      } else {
        tr_list <- rbind(tr_list, c(i,tr[j],"balanced", str_count(result[tr[j],((i-1)*3+7)], "t")))
      }
      }
    }
  }
  cat(sprintf("sample number %s", i), sep="\n")
}
############## trans interactions
tr_int_list <- list()
for (i in 1:length(samples2)){
#for (i in 1:100){
  if(all(grepl(pattern = "t", x = result[,((i-1)*3+7)]) == FALSE) == FALSE){
    tr <- list()
    tr <- which(grepl(pattern = "t", x = result[,((i-1)*3+7)]) == TRUE) 
    for(j in 1:length(tr)){  
      sep_terms <- unlist(strsplit( gsub(",","~",result[tr[j],((i-1)*3+7)]), "~" ))
      tr_terms <- sep_terms[which(grepl(pattern = "t", x = sep_terms))] 
      for(k in 1:length(tr_terms)){
        tr_int_list <- rbind(tr_int_list, c(i,tr[j],tr_terms[k]))
      }
    }
  }
  cat(sprintf("sample number %s", i), sep="\n")
}

ex_tr <- list()
for(j in 1:nrow(tr_int_list)){
  if(all(unlist(tr_int_list[j,1]) != samples_excluded) == FALSE){
    ex_tr <- c(unlist(ex_tr),sprintf("%s",j))  
  }
}
tr_int_list2 <- tr_int_list[-as.numeric(ex_tr),]
#tr_int_list2 <- cbind(tr_int_list2, matrix(1,ncol = 1,nrow = nrow(tr_int_list2)))

#tr_mult <- tr_list2[which(tr_list2[,4]!=1),c(1,2,4)]

#for(i in 1:nrow(tr_mult)){
#ss <- list()
#dd <- list()
#ss <- unlist(which(unlist(tr_int_list2[,1]) == unlist(tr_mult[i,1])))
#dd <- unlist(which(unlist(tr_int_list2[ss,2]) == unlist(tr_mult[i,2])))
#tr_int_list2[ss[dd],4] <- unlist(tr_mult[i,3])  
#}
tr_int_list2 < tr_int_list2[,c(1,2,3)]
#######################
add_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "add", x = result[,((i-1)*3+7)]) == FALSE) == FALSE){
    add <- list()
    add <- which(grepl(pattern = "add", x = result[,((i-1)*3+7)]) == TRUE)
      for(j in 1:length(add)){  
        add_list <- rbind(add_list, c(i,add[j], str_count(result[add[j],((i-1)*3+7)], "add")))
    }
  }
  cat(sprintf("sample number %s", i), sep="\n")
}

iso_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "i", x = result[,((i-1)*3+5)]) == FALSE) == FALSE){
    iso <- list()
    iso <- which(grepl(pattern = "i", x = result[,((i-1)*3+5)]) == TRUE)[1]
    if((grepl(pattern = "inv", x = result[iso,((i-1)*3+5)]) == FALSE) &&
       (grepl(pattern = "ins", x = result[iso,((i-1)*3+5)]) == FALSE) &&
       (grepl(pattern = "idic", x = result[iso,((i-1)*3+5)]) == FALSE) &&
       (grepl(pattern = "dic", x = result[iso,((i-1)*3+5)]) == FALSE)){
    if(Bands$band[iso] == "p11"){
      iso_list <- rbind(iso_list, c(i, Bands$chr[iso], "p10", str_count(result[iso,((i-1)*3+5)], "i")))
    } else {
      iso_list <- rbind(iso_list, c(i, Bands$chr[iso],  "q10", str_count(result[iso,((i-1)*3+5)], "i")))
    }
  }
  cat(sprintf("sample number %s", i), sep="\n")
  }
}
    
dup_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "dup", x = result[,((i-1)*3+7)]) == FALSE) == FALSE){
    dup <- list()
    dup <- which(grepl(pattern = "dup", x = result[,((i-1)*3+7)]) == TRUE)
    #unique( dup)
    #if(length(dup))
    for(j in 1:length(dup)){  
      dup_list <- rbind(dup_list, c(i,dup[j], str_count(result[dup[j],((i-1)*3+7)], "dup")))
    }
  }
  cat(sprintf("sample number %s", i), sep="\n")
} 

#del_list <- list()
#for (i in 1:length(samples2)){
#  if(all(grepl(pattern = "del", x = result[,((i-1)*3+7)]) == FALSE) == FALSE){
#    del <- list()
#    del <- unlist(which(grepl(pattern = "del", x = result[,((i-1)*3+7)]) == TRUE))
#    dels <- rbind(del, rep((i-1)*3+7, length(del)))
#    for(j in 1:ncol(del)){
#      sep_terms <- list()
#      sep_terms <- unlist(strsplit(gsub("\\,","~",result[dels[1,j],dels[2,j]]), "~" )) 
#      sep_terms <- sep_terms[which(grepl(pattern = "del", x = sep_terms))]
      
#    for(j in 1:length(del)){  
#      del_list <- rbind(del_list, c(i,del[j], "2bands", str_count(result[del[j],((i-1)*3+7)], "del")))
#    }
#  }
#  cat(sprintf("sample number %s", i), sep="\n")
#}

for (i in 1:length(samples2)){
  if(all(grepl(pattern = "del", x = result[,((i-1)*3+7)]) == FALSE) == FALSE){
    del <- list()
    del <- which(grepl(pattern = "del", x = result[,((i-1)*3+7)]) == TRUE)
    for(j in 1:length(del)){  
      del_list <- rbind(del_list, c(i,del[j], "1bands", str_count(result[del[j],((i-1)*3+7)], "del")))
    }
  }
  cat(sprintf("sample number %s", i), sep="\n")
}

########################

inv_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "inv", x = result[,((i-1)*3+7)]) == FALSE) == FALSE){
    inv <- list()
    inv <- which(grepl(pattern = "inv", x = result[,((i-1)*3+7)]) == TRUE)
    for(j in 1:length(inv)){  
      inv_list <- rbind(inv_list, c(i,inv[j], str_count(result[inv[j],((i-1)*3+7)], "inv")))
    }
  }
  cat(sprintf("sample number %s", i), sep="\n")
} 

insT_list <- list()
insF_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "ins", x = result[,((i-1)*3+7)]) == FALSE) == FALSE){
    insT <- list()
    insF <- list()
    
    insT <- which(grepl(pattern = "insT", x = result[,((i-1)*3+7)]) == TRUE)
    insF <- which(grepl(pattern = "insF", x = result[,((i-1)*3+7)]) == TRUE)
    
    for(j in 1:length(insT)){  
      insT_list <- rbind(insT_list, c(i, "insT", insT[j], str_count(result[insT[j],((i-1)*3+7)], "insT")))
    }
    for(j in 1:length(insF)){  
      insF_list <- rbind(insF_list, c(i, "insF", insF[j], str_count(result[insF[j],((i-1)*3+7)], "insF")))
    }
  }
  cat(sprintf("sample number %s", i), sep="\n")
} 
ins_list <- rbind(insT_list,insF_list)


idic_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "idic", x = result[,((i-1)*3+7)]) == FALSE) == FALSE){
    idic <- list()
    idic <- which(grepl(pattern = "idic", x = result[,((i-1)*3+7)]) == TRUE)
    for(j in 1:length(idic)){
    if(grepl(pattern = "p", x = Bands$band[idic[j]] )){
      idic_list <- rbind(idic_list, c(i, idic[j], "p", str_count(result[idic[j],((i-1)*3+7)], "idic")))
    } else {
      idic_list <- rbind(idic_list, c(i, idic[j], "q", str_count(result[idic[j],((i-1)*3+7)], "idic")))
    }  
    }  
  }    
    cat(sprintf("sample number %s", i), sep="\n")
}


dic_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "dic", x = result[,((i-1)*3+7)]) == FALSE) == FALSE){
    dic <- list()
    dic <- which(grepl(pattern = "dic", x = result[,((i-1)*3+7)]) == TRUE)
    for(j in 1:length(dic)){ 
      if((grepl(pattern = "idic", x = result[dic[j],((i-1)*3+7)]) == FALSE)){
        dic_list <- rbind(dic_list, c(i,dic[j], str_count(result[dic[j],((i-1)*3+7)], "dic")))
    }
    }
  }
  cat(sprintf("sample number %s", i), sep="\n")
}


tas_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "tas", x = result[,((i-1)*3+7)]) == FALSE) == FALSE){
    tas <- list()
    tas <- which(grepl(pattern = "tas", x = result[,((i-1)*3+7)]) == TRUE)
    for(j in 1:length(tas)){ 
      tas_list <- rbind(tas_list, c(i,tas[j], str_count(result[tas[j],((i-1)*3+7)], "tas")))
      
    }
  }
  cat(sprintf("sample number %s", i), sep="\n")
}


hsr_list <- list()
for (i in 1:length(samples2)){
  if(all(grepl(pattern = "hsr", x = result[,((i-1)*3+5)]) == FALSE) == FALSE){
    hsr <- list()
    hsr <- which(grepl(pattern = "hsr", x = result[,((i-1)*3+5)]) == TRUE)
    for(j in 1:length(hsr)){ 
      hsr_list <- rbind(hsr_list, c(i,hsr[j], str_count(result[hsr[j],((i-1)*3+5)], "hsr")))
      
    }
  }
  cat(sprintf("sample number %s", i), sep="\n")
}



classified <- list()
classified$plus <- plus_list
classified$minus <- minus_list
classified$tr <- tr_list
classified$add <- add_list
classified$iso <- iso_list
classified$dup <- dup_list
classified$del <- del_list
classified$inv <- inv_list
classified$ins <- ins_list
classified$idic <- idic_list
classified$dic <- dic_list
classified$tas <- tas_list
classified$hsr <- hsr_list


save(classified, file = sprintf("classified_SVs_v1.RData"))      


###### variation type bar chart
load("classified_SVs_v1.RData")

###########filter problematic samples #####
ex_plus <- list()
for(j in 1:nrow(plus_list)){
  if(all(unlist(plus_list[j,1]) != samples_excluded) == FALSE){
  ex_plus <- c(unlist(ex_plus),sprintf("%s",j))  
  }
}
plus_list2 <- plus_list[-as.numeric(ex_plus),]


ex_minus <- list()
for(j in 1:nrow(minus_list)){
  if(all(unlist(minus_list[j,1]) != samples_excluded) == FALSE){
    ex_minus <- c(unlist(ex_minus),sprintf("%s",j))  
  }
}
minus_list2 <- minus_list[-as.numeric(ex_minus),]


ex_tr <- list()
for(j in 1:nrow(tr_list)){
  if(all(unlist(tr_list[j,1]) != samples_excluded) == FALSE){
    ex_tr <- c(unlist(ex_tr),sprintf("%s",j))  
  }
}
tr_list2 <- tr_list[-as.numeric(ex_tr),]
###################
ex_tr_gain <- list()
for(j in 1:nrow(tr_gain_list)){
  if(all(unlist(tr_gain_list[j,1]) != samples_excluded) == FALSE){
    ex_tr_gain <- c(unlist(ex_tr_gain),sprintf("%s",j))
  }
}
tr_gain_list2 <- tr_gain_list[-as.numeric(ex_tr_gain),]


ex_tr_loss <- list()
for(j in 1:nrow(tr_loss_list)){
  if(all(unlist(tr_loss_list[j,1]) != samples_excluded) == FALSE){
    ex_tr_loss <- c(unlist(ex_tr_loss),sprintf("%s",j))
  }
}
tr_loss_list2 <- tr_loss_list[-as.numeric(ex_tr_loss),]

####################
ex_add <- list()
for(j in 1:nrow(add_list)){
  if(all(unlist(add_list[j,1]) != samples_excluded) == FALSE){
    ex_add <- c(unlist(ex_add),sprintf("%s",j))  
  }
}
add_list2 <- add_list[-as.numeric(ex_add),]


ex_iso <- list()
for(j in 1:nrow(iso_list)){
  if(all(unlist(iso_list[j,1]) != samples_excluded) == FALSE){
    ex_iso <- c(unlist(ex_iso),sprintf("%s",j))  
  }
}
iso_list2 <- iso_list[-as.numeric(ex_iso),]


ex_dup <- list()
for(j in 1:nrow(dup_list)){
  if(all(unlist(dup_list[j,1]) != samples_excluded) == FALSE){
    ex_dup <- c(unlist(ex_dup),sprintf("%s",j))  
  }
}
dup_list2 <- dup_list[-as.numeric(ex_dup),]


ex_del <- list()
for(j in 1:nrow(del_list)){
  if(all(unlist(del_list[j,1]) != samples_excluded) == FALSE){
    ex_del <- c(unlist(ex_del),sprintf("%s",j))  
  }
}
del_list2 <- del_list[-as.numeric(ex_del),]


ex_inv <- list()
for(j in 1:nrow(inv_list)){
  if(all(unlist(inv_list[j,1]) != samples_excluded) == FALSE){
    ex_inv <- c(unlist(ex_inv),sprintf("%s",j))  
  }
}
inv_list2 <- inv_list[-as.numeric(ex_inv),]


ex_ins <- list()
for(j in 1:nrow(ins_list)){
  if(all(unlist(ins_list[j,1]) != samples_excluded) == FALSE){
    ex_ins <- c(unlist(ex_ins),sprintf("%s",j))  
  }
}
ins_list2 <- ins_list[-as.numeric(ex_ins),]


ex_idic <- list()
for(j in 1:nrow(idic_list)){
  if(all(unlist(idic_list[j,1]) != samples_excluded) == FALSE){
    ex_idic <- c(unlist(ex_idic),sprintf("%s",j))  
  }
}
idic_list2 <- idic_list[-as.numeric(ex_idic),]


ex_dic <- list()
for(j in 1:nrow(dic_list)){
  if(all(unlist(dic_list[j,1]) != samples_excluded) == FALSE){
    ex_dic <- c(unlist(ex_dic),sprintf("%s",j))  
  }
}
dic_list2 <- dic_list[-as.numeric(ex_dic),]


ex_tas <- list()
for(j in 1:nrow(tas_list)){
  if(all(unlist(tas_list[j,1]) != samples_excluded) == FALSE){
    ex_tas <- c(unlist(ex_tas),sprintf("%s",j))  
  }
}
tas_list2 <- tas_list[-as.numeric(ex_tas),]


ex_hsr <- list()
for(j in 1:nrow(hsr_list)){
  if(all(unlist(hsr_list[j,1]) != samples_excluded) == FALSE){
    ex_hsr <- c(unlist(ex_hsr),sprintf("%s",j))  
  }
}
hsr_list2 <- hsr_list[-as.numeric(ex_hsr),]

####################
N_plus <- sum(as.numeric(plus_list2[,3]))
N_minus <- sum(as.numeric(minus_list2[,3]))

tr_bal1 <- tr_list[which(tr_list[,3] == "balanced"),]
N_tr_bal1 <- ceiling(sum(as.numeric(tr_bal[,4]))/2)
tr_unbal1 <- tr_list[which(tr_list[,3] != "balanced"),]
N_tr_unbal1 <- ceiling(sum(as.numeric(tr_unbal[,4]))/2)

tr_bal2 <- tr_list2[which(tr_list2[,3] == "balanced"),]
N_tr_bal2 <- ceiling(sum(as.numeric(tr_bal[,4]))/2)
tr_unbal2 <- tr_list2[which(tr_list2[,3] != "balanced"),]
N_tr_unbal2 <- ceiling(sum(as.numeric(tr_unbal[,4]))/2)

#####################
N_add <- sum(as.numeric(add_list2[,3]))

iso_q <- iso_list2[which(iso_list2[,3] == "q10"),]
N_iso_q <- sum(as.numeric(iso_q[,4]))
iso_p <- iso_list2[which(iso_list2[,3] != "q10"),]
N_iso_p <- sum(as.numeric(iso_p[,4]))

N_dup <- sum(as.numeric(dup_list2[,3]))/2  

#del_1 <- del_list[which(del_list[,3] == "2bands"),]
#N_del_1 <- sum(as.numeric(del_1[,4]))/2
#del_2 <- del_list[which(del_list[,3] != "2bands"),]
#N_del_2 <- sum(as.numeric(del_2[,4]))

#del_1 <- del_list[which(del_list[,4] == 1),]
#N_del_1 <- sum(as.numeric(del_1[,4]))/2
#del_2 <- del_list[which(del_list[,3] != "2bands"),]
#N_del_2 <- sum(as.numeric(del_2[,4]))


#########run carefully############


del_sam <- unlist(unique(del_list2[,1]))
del_sta <- vector(length = length(del_sam))
for(k in 1:length(del_sam)){
  del_sta[k] <- length(which(del_list2[,1] == del_sam[k]))
}
unique(del_sta)

N_del_1 <- length(del_sam[which(del_sta == 1)])+length(del_sam[which(del_sta == 3)])+length(del_sam[which(del_sta == 5)])+
  length(del_sam[which(del_sta == 7)])+length(del_sam[which(del_sta == 11)])

N_del_2 <- (sum(as.numeric(del_list2[,4])) - N_del_1)/2

##################

N_inv <- sum(as.numeric(inv_list2[,3]))/2  

N_ins <- sum(as.numeric(ins_list2[which(ins_list2[,2] == "insT"),4]))

N_idic <- sum(as.numeric(idic_list2[,4])) 

N_dic <- ceiling(sum(as.numeric(dic_list2[,3]))/2)

#####
#N_tas <- (sum(as.numeric(tas_list[,3]))/2)-23
N_tas <- ceiling((sum(as.numeric(tas_list2[,3]))/2)-17)

N_hsr <- sum(as.numeric(hsr_list2[,3])) 


df7 <- data.frame(samples=c("+","-", "tr", "tr", "add", 
                            "del", "del"),
                  len=c(N_plus, N_minus, N_tr_bal, N_tr_unbal, N_add,
                        N_del_1+400, N_del_2-400),
#                  len=c(87897, 81284, 29149, 15953, 23536,
#                        N_del_1+300, N_del_2-300),
                  variation=c("+","-", "balanced tr", "unbalanced tr", "add", 
                              "del(1 band)", "del(2 bands)")
                  ,cum_no <- c(43000, 40000, 16000, 40000, 15000, 14000, 7500)
                  )

p <- ggplot(data=df7, aes(x=samples, y=len, fill=reorder(variation, len))) +
  geom_bar(stat="identity")+
  geom_text(aes(y=cum_no, label=len), vjust=1.6, 
            color="black", size=3.5)+
  #scale_fill_brewer(palette="Paired")+
  theme_minimal()+labs(title="Variation frequency in Mitelman")+ theme(plot.title = element_text(hjust = 0.5))

p + scale_x_discrete(name ="Variation type", limits=c("+", "-", "tr", "add", "del")) + scale_y_continuous(name="Frequency", limits=c(0, 90000)) + 
  scale_fill_manual(values=c('coral4','pink3','darkseagreen4', "deepskyblue1", "cadetblue4",  "wheat3", "chocolate3")) +
  guides(fill=guide_legend(title="variation types"))

#aquamarine4 lightblue4 darkseagreen4 indianred coral4 lightpink4  cadetblue4


df8 <- data.frame(samples=c("iso","dup", "inv","ins", "idic","dic", "tas", "hsr"),
                  len=c( N_iso_q+N_iso_p,N_dup, N_inv, N_ins, N_idic, N_dic, N_tas, N_hsr),
                  variation=c("unbalanced","unbalanced", "balanced","balanced", 
                              "unbalanced","unbalanced", "balanced", "unbalanced"),
                  cum_no <- c(3600, 1400, 2200, 600, 375, 900, 340, 420)
                  )

p <- ggplot(data=df8, aes(x=samples, y=len, fill=reorder(variation, len))) +
  geom_bar(stat="identity")+
  geom_text(aes(y=cum_no, label=len), vjust=1.6, 
            color="black", size=3.5) + theme_minimal()+labs(title="Variation frequency in Mitelman")

p <- p + scale_x_discrete(name ="Variation type", limits=c("iso", "dup", "dic", "hsr", "idic", "inv", "ins", "tas")) + 
  scale_y_continuous(name="Frequency", limits=c(0, 7000)) + 
  scale_fill_manual(values=c("cadetblue4","pink3")) +
  guides(fill=guide_legend(title="variation types"))+
#  theme_minimal()+labs(title="Variation frequency in Mitelman", x="Variation type", y = "Frequency")+
  theme(plot.title = element_text(hjust = 0.5))
#aquamarine4 lightblue4 darkseagreen4 indianred coral4 lightpink4  cadetblue4


################bal_unbal

df82 <- data.frame(samples=c("+","-","balanced_tr & tas", "del & dic", "add", 
                            "unbalanced_tr","iso & idic", "inv","dup & hsr","ins"),
                  len=c(N_plus, N_minus , N_tr_bal+N_tas, N_del_1+N_del_2+N_dic, N_add, 
                        N_tr_unbal, N_iso_q+N_iso_p+N_idic, N_inv, N_dup+N_hsr, N_ins),
 #                 len=c(87897, 81284 , 29474, N_del_1+N_del_2+N_dic, 25676, 
#                        15953, 6906, N_inv, 2666, N_ins),
                  variation=c("gain", "loss", "balanced", "loss", "loss", 
                              "loss & gain", "loss & gain", "balanced", "gain", "balanced"),
                  cum_no <- c(33000, 30000, 12000, 5900, 11000,
                              7500, 4000, 3000, 2500, 2300)
)

p <- ggplot(data=df82, aes(x=samples, y=len, fill=reorder(variation, len))) +
  geom_bar(stat="identity")+
  geom_text(aes(y=cum_no, label=len), vjust=1.6, 
            color="black", size=3.5)+
  #scale_fill_brewer(palette="Paired")+
  theme_minimal()+labs(title="Variation frequency in Mitelman without problematic samples")+ theme(plot.title = element_text(hjust = 0.5))

p + scale_x_discrete(name ="Variation type", limits=c("+","-", "balanced_tr & tas", "add", "unbalanced_tr", "del & dic","iso & idic", "inv","dup & hsr","ins")) + 
  scale_y_continuous(name="Frequency", limits=c(0, 90000)) + 
  scale_fill_manual(values=c("pink3","aquamarine3", 'wheat3', "deepskyblue3")) +
  guides(fill=guide_legend(title="variation types"))   

######percent diff

dec_percent <- c(31.7, 32, 32.2, 31.9, 32.8, 32.1, 32.9, 31.2, 32.3, 31.9)


df82 <- data.frame(samples=c("+","-","balanced_tr & tas", "del & dic", "add", 
                             "unbalanced_tr","iso & idic", "inv","dup & hsr","ins"),
                   len= dec_percent,
                   variation=c("gain", "loss", "balanced", "loss", "loss", 
                               "loss & gain", "loss & gain", "balanced", "gain", "balanced"),
                   cum_no <- c(33000, 30000, 12000, 5900, 11000,
                               7500, 4000, 3000, 2500, 2300)
)

p <- ggplot(data=df82, aes(x=samples, y=len, fill=reorder(variation, len))) +
  geom_bar(stat="identity")+
  geom_text(aes(y=len*0.55, label=paste(dec_percent, "%")), vjust=1.6, 
            color="black", size=3.5)+
  #scale_fill_brewer(palette="Paired")+
  theme_minimal()+labs(title="Decrease in the number of variations by deleting problematic samples")+ theme(plot.title = element_text(hjust = 0.5))

p + scale_x_discrete(name ="Variation type", limits=c("+","-", "balanced_tr & tas", "add", "unbalanced_tr", "del & dic","iso & idic", "inv","dup & hsr","ins")) + 
  scale_y_continuous(name="Percentage", limits=c(0, 40)) + 
  scale_fill_manual(values=c("aquamarine3", "deepskyblue3", 'wheat3',"pink3")) +
  guides(fill=guide_legend(title="variation types"))   


################# TOP variations

plus_chr <- unlist(unique(plus_list[,2]))
ppp <- matrix(0,nrow =length(plus_chr), ncol = 2 )
for(i in 1:(length(plus_chr))) {
ppp[i,]  <- c(sum(as.numeric(plus_list[which(plus_list[,2] == plus_chr[i]),3])), plus_chr[i])
}
ppp[,1] <- as.numeric(ppp[,1])
ppp <- ppp[order(as.numeric(ppp[,1]), decreasing = TRUE),]


plus_chr <- unlist(unique(plus_list2[,2]))
ppp2 <- matrix(0,nrow =length(plus_chr), ncol = 2 )
for(i in 1:(length(plus_chr))) {
  ppp2[i,]  <- c(sum(as.numeric(plus_list2[which(plus_list2[,2] == plus_chr[i]),3])), plus_chr[i])
}
ppp2[,1] <- as.numeric(ppp2[,1])
ppp2 <- ppp2[order(as.numeric(ppp2[,1]), decreasing = TRUE),]

chrs <- c("chr1" , "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
          "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
                      
            


df6 <- data.frame(samples=ppp2[,2],
                  len=as.numeric(ppp2[,1])
                  #variation=c("+","-", "balanced tr", "unbalanced tr", "add", 
                   #           "del(1 band)", "del(2 bands)")
                 # ,cum_no <- c(6000, 4200, 3350)
)

p <- ggplot(data=df6, aes(x=samples, y=len)) +
  geom_bar(stat="identity", fill= "deepskyblue3")+
  geom_text(aes(y=as.numeric(len)*0.55, label=len), vjust=1.6, 
            color="black", size=3)+
  #scale_fill_brewer(palette="Paired")+
  theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))

p + scale_x_discrete(name ="Chromosomes", limits=chrs) + scale_y_continuous(name="Frequency", limits=c(0, 7000))+
 theme(text = element_text(size = 12))+labs(title="+ samples (without problematic samples)") 
  #scale_fill_manual(values=c("chocolate3")) 
  #guides(fill=guide_legend(title="variation types"))

#+labs(title="Chromosomes with most +")
#theme(axis.text = element_text(size = 10))

minus_chr <- unlist(unique(minus_list[,2]))
ppp <- matrix(0,nrow =length(minus_chr), ncol = 2 )
for(i in 1:(length(minus_chr))) {
  ppp[i,]  <- c(sum(as.numeric(minus_list[which(minus_list[,2] == minus_chr[i]),3])), minus_chr[i])
}
ppp[,1] <- as.numeric(ppp[,1])
ppp <- ppp[order(as.numeric(ppp[,1]), decreasing = TRUE),]

minus_chr <- unlist(unique(minus_list2[,2]))
ppp2 <- matrix(0,nrow =length(minus_chr), ncol = 2 )
for(i in 1:(length(minus_chr))) {
  ppp2[i,]  <- c(sum(as.numeric(minus_list2[which(minus_list2[,2] == minus_chr[i]),3])), minus_chr[i])
}
ppp2[,1] <- as.numeric(ppp2[,1])
ppp2 <- ppp2[order(as.numeric(ppp2[,1]), decreasing = TRUE),]

df6 <- data.frame(samples=ppp2[,2],
                  len=as.numeric(ppp2[,1])
                  #variation=c("+","-", "balanced tr", "unbalanced tr", "add", 
                  #           "del(1 band)", "del(2 bands)")
                 # ,cum_no <- c(3000, 2800, 2550)
)

p <- ggplot(data=df6, aes(x=samples, y=len)) +
  geom_bar(stat="identity", fill= "wheat3")+
  geom_text(aes(y=as.numeric(len)*0.55, label=len), vjust=1.6, 
            color="black", size=3)+
  #scale_fill_brewer(palette="Paired")+
  theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))

p + scale_x_discrete(name ="Chromosomes", limits=chrs) + scale_y_continuous(name="Frequency", limits=c(0, 4000))+
  theme(text = element_text(size = 12))+labs(title="- samples (without problematic samples)") 



############tr bal

tr_bal_band <- unlist(unique(tr_bal2[,2]))
ppp <- matrix(0,nrow =length(tr_bal_band), ncol = 2 )
for(i in 1:(length(tr_bal_band))) {
  ppp[i,]  <- c(sum(as.numeric(tr_bal2[which(tr_bal2[,2] == tr_bal_band[i]),4])), tr_bal_band[i])
}
ppp[,1] <- as.numeric(ppp[,1])

tas_band <- unlist(unique(tas_list2[,2]))
ttas <- matrix(0,nrow =length(tas_band), ncol = 2 )
for(l in 1:(length(tas_band))) {
  ttas[l,]  <- c(sum(as.numeric(tas_list2[which(tas_list2[,2] == tas_band[l]),3])), tas_band[l])
}
ttas[,1] <- as.numeric(ttas[,1])

for(i in 1:(length(tas_band))){
ind <- which(ppp[,2] == ttas[i,2])
ppp[ind,1] <- as.numeric(ttas[i,1]) + as.numeric(ppp[ind,1])
}
#ppp <- ppp[order(as.numeric(ppp[,1]), decreasing = TRUE),]


ppp <- cbind(ppp, matrix(0, ncol = 1, nrow = 320))
for(j in 1:nrow(ppp)){
 row_ind <- as.numeric(ppp[j,2]) 
 ppp[j,3] <- gsub("chr","",paste(sprintf("%s",Bands$chr[row_ind]),sprintf("%s", Bands$band[row_ind]),sep=""))  
}
ppp[,2] <- as.numeric(ppp[,2])

################################
bandd <- matrix(0,ncol = 1, nrow = 320)
for(j in 1:nrow(Bands)){
  #row_ind <- as.numeric(ppp[j,2]) 
  bandd[j,1] <- gsub("chr","",paste(sprintf("%s",Bands$chr[j]),sprintf("%s", Bands$band[j]),sep=""))  
}
bandd <- bandd[c(1:24,134:159,178:303,25:133,160:177,304:320),] 
#################################

df6 <- data.frame(samples=ppp[,3],
                  len=as.numeric(ppp[,1])
                  #variation=c("+","-", "balanced tr", "unbalanced tr", "add", 
                  #           "del(1 band)", "del(2 bands)")
                  # ,cum_no <- c(3000, 2800, 2550)
)

p <- ggplot(data=df6, aes(x=samples, y=len)) +
  geom_bar(stat="identity", fill= "aquamarine3")+
 # geom_text(aes(y=as.numeric(len)*0.55, label=len), vjust=1.6, 
  #          color="black", size=3)+
  #scale_fill_brewer(palette="Paired")+
  theme_minimal()+ theme(plot.title = element_text(hjust = 0.4))

p + scale_x_discrete(name ="Band", limits=as.factor(bandd) ) + scale_y_continuous(name="Frequency", limits=c(0, 3500))+
  theme(text = element_text(size = 7))+labs(title="Balanced_tr and tas") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 

 
# "aquamarine3", "deepskyblue3"
##############tr unbal
tr_unbal_band <- unlist(unique(tr_unbal2[,2]))
ppp <- matrix(0,nrow =length(tr_unbal_band), ncol = 2 )
for(i in 1:(length(tr_unbal_band))) {
  ppp[i,]  <- c(sum(as.numeric(tr_unbal2[which(tr_unbal2[,2] == tr_unbal_band[i]),4])), tr_unbal_band[i])
}

ppp[,1] <- as.numeric(ppp[,1])

ppp <- cbind(ppp, matrix(0, ncol = 1, nrow = 320))
for(j in 1:nrow(ppp)){
  row_ind <- as.numeric(ppp[j,2]) 
  ppp[j,3] <- gsub("chr","",paste(sprintf("%s",Bands$chr[row_ind]),sprintf("%s", Bands$band[row_ind]),sep=""))  
}


df6 <- data.frame(samples=ppp[,3],
                  len=as.numeric(ppp[,1])
)

p <- ggplot(data=df6, aes(x=samples, y=len)) +
  geom_bar(stat="identity", fill= "pink3")+
  theme_minimal()+ theme(plot.title = element_text(hjust = 0.4))

p + scale_x_discrete(name ="Band", limits=as.factor(bandd) ) + scale_y_continuous(name="Frequency", limits=c(0, 800))+
  theme(text = element_text(size = 7))+labs(title="Unbalanced_tr") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 

#ppp <- ppp[order(as.numeric(ppp[,1]), decreasing = TRUE),]

###################tr unbal gain/loss

tr_gain_band <- tr_unbal2[which(tr_unbal2[,3] == "gain"),]
tr_gain <- unlist(unique(tr_gain_band[,2]))
ppp <- matrix(0,nrow =length(tr_gain), ncol = 2 )
for(i in 1:(length(tr_gain))) {
  ppp[i,]  <- c(sum(as.numeric(tr_gain_band[which(tr_gain_band[,2] == tr_gain[i]),4])), tr_gain[i])
}
ppp[,1] <- as.numeric(ppp[,1])

ppp <- rbind(ppp, matrix(0,ncol = 2, nrow = 5))

ppp[316:320,2] <- c(88,156,167,274,316)

#ppp <- ppp[order(as.numeric(ppp[,1]), decreasing = TRUE),]
ppp <- cbind(ppp, matrix(0, ncol = 1, nrow = 320))
for(j in 1:nrow(ppp)){
  row_ind <- as.numeric(ppp[j,2]) 
  ppp[j,3] <- gsub("chr","",paste(sprintf("%s",Bands$chr[row_ind]),sprintf("%s", Bands$band[row_ind]),sep=""))  
}

ppp <- cbind(ppp, rep("gain", 320))


tr_loss_band <- tr_unbal2[which(tr_unbal2[,3] == "loss"),]
tr_loss <- unlist(unique(tr_loss_band[,2]))
ppp2 <- matrix(0,nrow =length(tr_loss), ncol = 2 )
for(i in 1:(length(tr_loss))) {
  ppp2[i,]  <- c(sum(as.numeric(tr_loss_band[which(tr_loss_band[,2] == tr_loss[i]),4])), tr_loss[i])
}
ppp2[,1] <- as.numeric(ppp2[,1])
ppp2 <- rbind(ppp2, c(0,80))

ppp2 <- cbind(ppp2, matrix(0, ncol = 1, nrow = 320))
for(j in 1:nrow(ppp2)){
  row_ind <- as.numeric(ppp2[j,2]) 
  ppp2[j,3] <- gsub("chr","",paste(sprintf("%s",Bands$chr[row_ind]),sprintf("%s", Bands$band[row_ind]),sep=""))  
}

ppp2[,1] <- -as.numeric(ppp2[,1])
ppp2 <- cbind(ppp2, rep("loss", 320))

unbal_tr_pp <- rbind(ppp,ppp2)


df82 <- data.frame(samples=unbal_tr_pp[,3],
                   len= as.numeric(unbal_tr_pp[,1]),
                   variation=c(rep("gain",320), rep("loss",320))
)

p <- ggplot(data=df82, aes(x=samples, y=len, fill=reorder(variation, len))) +
  geom_bar(stat="identity")+
  theme_minimal()+labs(title="Gain & loss translocation")+ theme(plot.title = element_text(hjust = 0.4))

p + scale_x_discrete(name ="Band", limits=as.factor(bandd)) + 
  scale_y_continuous(name="Frequency", limits=c(-500, 600)) + 
  scale_fill_manual(values=c('wheat3',"deepskyblue3")) +
  guides(fill=guide_legend(title="variation types"))  +
  theme(text = element_text(size = 7))+ theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))


###########tr unbal gain/loss affected regions

tr_gain <- unlist(unique(tr_gain_list2[,2]))
ppp <- matrix(0,nrow =length(tr_gain), ncol = 2 )
for(i in 1:(length(tr_gain))) {
  ppp[i,]  <- c(sum(as.numeric(tr_gain_list2[which(tr_gain_list2[,2] == tr_gain[i]),4])), tr_gain[i])
}
ppp[,1] <- as.numeric(ppp[,1])

#ppp <- ppp[order(as.numeric(ppp[,1]), decreasing = TRUE),]
ppp <- cbind(ppp, matrix(0, ncol = 1, nrow = 320))
for(j in 1:nrow(ppp)){
  row_ind <- as.numeric(ppp[j,2]) 
  ppp[j,3] <- gsub("chr","",paste(sprintf("%s",Bands$chr[row_ind]),sprintf("%s", Bands$band[row_ind]),sep=""))  
}

ppp <- cbind(ppp, rep("gain", 320))


tr_loss <- unlist(unique(tr_loss_list2[,2]))
ppp2 <- matrix(0,nrow =length(tr_loss), ncol = 2 )
for(i in 1:(length(tr_loss))) {
  ppp2[i,]  <- c(sum(as.numeric(tr_loss_list2[which(tr_loss_list2[,2] == tr_loss[i]),4])), tr_loss[i])
}
ppp2[,1] <- as.numeric(ppp2[,1])

ppp2 <- cbind(ppp2, matrix(0, ncol = 1, nrow = 320))
for(j in 1:nrow(ppp2)){
  row_ind <- as.numeric(ppp2[j,2]) 
  ppp2[j,3] <- gsub("chr","",paste(sprintf("%s",Bands$chr[row_ind]),sprintf("%s", Bands$band[row_ind]),sep=""))  
}

ppp2[,1] <- -as.numeric(ppp2[,1])
ppp2 <- cbind(ppp2, rep("loss", 320))

unbal_tr_pp <- rbind(ppp,ppp2)


df82 <- data.frame(samples=unbal_tr_pp[,3],
                   len= as.numeric(unbal_tr_pp[,1]),
                   variation=c(rep("gain",320), rep("loss",320))
)

p <- ggplot(data=df82, aes(x=samples, y=len, fill=reorder(variation, len))) +
  geom_bar(stat="identity")+
  theme_minimal()+labs(title="Gain & loss translocation (affected regions)")+ theme(plot.title = element_text(hjust = 0.4))

p + scale_x_discrete(name ="Band", limits=as.factor(bandd)) + 
  scale_y_continuous(name="Quantity", limits=c(-800, 1800)) + 
  scale_fill_manual(values=c('wheat3',"deepskyblue3")) +
  guides(fill=guide_legend(title="variation types"))  +
  theme(text = element_text(size = 7))+ theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

  ######same fig for bands


bb <- unbal_tr_pp[,3]
bbq <- grep(pattern = "q", x = bb)
bbp <- grep(pattern = "p", x = bb)

for(i in bbq){
  bb[i] <- unlist(strsplit( gsub("q","q~",unbal_tr_pp[i,3]), "~" ))[1]
  
}

for(i in bbp){
  bb[i] <- unlist(strsplit( gsub("p","p~",unbal_tr_pp[i,3]), "~" ))[1]
  
}

arm_band <- matrix(0, ncol = 2, nrow = length(un_b)*2)
un_b <- unique(bb)
arm_band[,1] <- rep(un_b,2)

for(j in 1:length(un_b)){
  pop <- list()
  nen <- list()
  pop <- unlist(which(bb[1:320] == un_b[j]))
  nen <- unlist(which(bb[321:640] == un_b[j]))+320
  arm_band[j,2] <- sum(as.numeric(unbal_tr_pp[pop,1]))
  arm_band[j+length(un_b),2] <- sum(as.numeric(unbal_tr_pp[nen,1]))
  
}

bandarm <- matrix(0, ncol = 1, nrow = 48)
for(i in 1:22){
   
  bandarm[c((2*i)-1,(2*i)),] <- c(sprintf("%sp",i), sprintf("%sq",i))
}
bandarm[45:48,1] <- c("Xp","Xq","Yp","Yq")

df82 <- data.frame(samples=arm_band[,1],
                   len= as.numeric(arm_band[,2]),
                   variation=c(rep("gain",48), rep("loss",48))
)

p <- ggplot(data=df82, aes(x=samples, y=len, fill=reorder(variation, len))) +
  geom_bar(stat="identity")+
  theme_minimal()+labs(title="unbalanced translocation Mitelman")+ theme(plot.title = element_text(hjust = 0.4))

p + scale_x_discrete(name ="chr-arm", limits=as.factor(bandarm)) + 
  scale_y_continuous(name="Quantity", limits=c(-4300, 19000)) + 
  scale_fill_manual(values=c('wheat3',"deepskyblue3")) +
  guides(fill=guide_legend(title="variation types"))  +
  theme(text = element_text(size = 14))+ theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))


###############translocation interaction
tr_int_ind <- as.numeric(unlist(unique(tr_int_list2[,1])))
tr_int_band <- as.numeric(unlist(unique(tr_int_list2[,2])))
tr_int_mat <- matrix(0, nrow = length(tr_int_band), ncol = length(tr_int_band))

for(i in 1:length(tr_int_ind)){
#for(i in 1:200){
qq <- which(tr_int_list2[,1] == tr_int_ind[i])  
uni_tr <- unlist(unique(tr_int_list2[qq,3]))
for(j in 1:length(uni_tr)){
 yy <- qq[which(tr_int_list2[qq,3] == uni_tr[j])]
 zz <- as.numeric(unlist(unique(tr_int_list2[yy,2])))
 tr_int_mat[zz[1],zz[2]] <- ceiling(tr_int_mat[zz[1],zz[2]] + (length(yy)/length(zz)))
 }
}

rownames(tr_int_mat) <- as.factor(bandd)
colnames(tr_int_mat) <- as.factor(bandd)
####################
trr <- unlist(unique(tr_list[,2]))
ppp <- matrix(0,nrow =length(trr), ncol = 2 )
for(i in 1:(length(trr))) {
  ppp[i,]  <- c(sum(as.numeric(tr_list[which(tr_list[,2] == trr[i]),4])), trr[i])
}
ppp[,1] <- as.numeric(ppp[,1])
ppp <- ppp[order(as.numeric(ppp[,1]), decreasing = TRUE),]


df6 <- data.frame(samples=c("balanced","unbalanced", "balanced","unbalanced", "balanced", "unbalanced"),
                  len=c(4609,840, 3983,1127, 3599, 1097),
                  variation=c("14q34","14q34", "22q11","22q11", "9q34", "9q34")
                  ,cum_no <- c(2900,5500,2700,5000, 2500, 4400)
)

p <- ggplot(data=df6, aes(x=variation, y=len, fill=reorder(samples, len))) +
  geom_bar(stat="identity")+
  geom_text(aes(y=cum_no, label=len), vjust=1.6, 
            color="black", size=9)+
  #scale_fill_brewer(palette="Paired")+
  theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))

p + scale_x_discrete(name ="Chromosomes", limits=c("14q34", "22q11","9q34")) + scale_y_continuous(name="Frequency", limits=c(0, 6000))+
  scale_fill_manual(values=c("darkseagreen4","cadetblue4")) + theme(text = element_text(size = 30)) 

#, 'pink3', "deepskyblue1",  





##################
dell <- unlist(unique(del_list[,2]))
ppp <- matrix(0,nrow =length(dell), ncol = 2 )
for(i in 1:(length(dell))) {
  ppp[i,]  <- c(sum(as.numeric(del_list[which(del_list[,2] == dell[i]),4])), dell[i])
}
ppp[,1] <- as.numeric(ppp[,1])
ppp <- ppp[order(as.numeric(ppp[,1]), decreasing = TRUE),]

dell <- unlist(unique(del_1[,2]))
ppp <- matrix(0,nrow =length(dell), ncol = 2 )
for(i in 1:(length(dell))) {
  ppp[i,]  <- c(sum(as.numeric(del_1[which(del_1[,2] == dell[i]),4])), dell[i])
}
ppp[,1] <- as.numeric(ppp[,1])
ppp <- ppp[order(as.numeric(ppp[,1]), decreasing = TRUE),]


df6 <- data.frame(samples=c("1_band","2-band", "1_band","2-band","1_band","2-band"),
                  len=c(137,1188, 34,1232, 681, 509),
                  variation=c("5q13","5q13", "5q33","5q33", "6q21", "6q21")
                  ,cum_no <- c(1330,700,1320,730, 350, 870)
)

p <- ggplot(data=df6, aes(x=variation, y=len, fill=reorder(samples, len))) +
  geom_bar(stat="identity")+
  geom_text(aes(y=cum_no, label=len), vjust=1.6, 
            color="black", size=9)+
  #scale_fill_brewer(palette="Paired")+
  theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))

p + scale_x_discrete(name ="Chromosomes", limits=c("5q13", "5q33","6q21")) + scale_y_continuous(name="Frequency", limits=c(0, 1400))+
  scale_fill_manual(values=c('coral4', 'pink3')) + theme(text = element_text(size = 30)) 



add_chr <- unlist(unique(add_list[,2]))
ppp <- matrix(0,nrow =length(add_chr), ncol = 2 )
for(i in 1:(length(add_chr))) {
  ppp[i,]  <- c(sum(as.numeric(add_list[which(add_list[,2] == add_chr[i]),3])), add_chr[i])
}
ppp[,1] <- as.numeric(ppp[,1])
ppp <- ppp[order(as.numeric(ppp[,1]), decreasing = TRUE),]


df6 <- data.frame(samples=c("14q32", "1p36", "19q13"),
                  len=c(798,725,659)
                  #variation=c("+","-", "balanced tr", "unbalanced tr", "add", 
                  #           "del(1 band)", "del(2 bands)")
                  ,cum_no <- c(400, 360, 330)
)

p <- ggplot(data=df6, aes(x=samples, y=len)) +
  geom_bar(stat="identity", fill="deepskyblue1")+
  geom_text(aes(y=cum_no, label=len), vjust=1.6, 
            color="black", size=9)+
  #scale_fill_brewer(palette="Paired")+
  theme_minimal()+ theme(plot.title = element_text(hjust = 0.5))

p + scale_x_discrete(name ="Chromosomes", limits=c("14q32", "1p36", "19q13")) + scale_y_continuous(name="Frequency", limits=c(0, 800))+
  theme(text = element_text(size = 30)) 


###########sample heatmap
library(ggplot2)
library(pheatmap)
library(gplots)
library(circlize)
library(RColorBrewer)

mat <- c()

hall_heatmap <- pheatmap(hallmark_mat, cluster_rows = TRUE, cluster_cols = FALSE, main = "Hallmark", gaps_col = c(3,13),legend_labels = c("1e-8", "0.2", "0.4", "0.6", "0.8", "1", "FDR value\n\n"),
                         border_color = "grey60", fontsize_row=5, fontsize_col=6, color=colorRampPalette(c("#003300", "forestgreen", "white"))(50))



