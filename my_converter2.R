options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
DATA_PATH <- args[1]
RESULT_PATH <- args[2]

#####test######
#tttt <- "der(5)t(5;20),ins(X;11)(q22;p15p13),+inv(13)(p13q21),ins(X;11)(q22;p15p13),t(9;8),del(12)(p11p13),del(14)(q13),dup(14)(q24q32),der(4)t(3;4)(p2?;q12),-Y,t(9;22)(q34;q12),+i(17)(q10),+der(22)t(9;22),-der(10t(1;2)(q21;q11),+X,+add(6)(p11)"
#tttt <- "der(5)t(5;20),ins(X11)(q22;p15p13),+inv(13(p13q21),ins(X;11)(q22;p15p13),t(9;8),del(12)(p11p13q11),del(14)(q13),dup(140(q24q32),der(4)t(3;4)(p2?;q12),-Y,t(9;22)(q34;q12),+i917)(q10),+der(22)t(9;22),-der(10t(1;2)(q21;q11),+X,+add(6(p11)"
#tttt <-"der(4)t(3;4)(p2?;q12),t(9;22)(q34;q12),+der(22)t(9;22),-der(10t(1;2)(q21;q11),+X"
#tttt <- "45,XX,-8,t(9;22)(q34;q11),der(8)t(8;14)(q11;q32)"
#tttt <- "	47,XX,+12/48,idem,+3/46,idem,-11/45,XX,-11"
#tttt <- "43,XX,-1,-6,-7,-11,-22,+2mar"
#tttt <- "idic(17)(p11),idic(1)(q21)"
#tttt <- "+der(22)t(9;22)(q34;q11)x2,add(18)(p11)x2"
#tttt <- "46,XY,tas(1;7)(q44;p22),add(17)(p?),add(19)(p?),tas(1;19;22)(p36;p13;p13),tas(19;15;21)(q13;q26p13;p12)"
#tttt <- "der(4;11)t(4;11)(p13;q13)hsr(11)(q13)"
#tttt <- "+tas(1;19)(p36;p13)"
#tttt <- "der(2)t(Y;2)(q11;p21)"
#SVs[j] <- "+der(22)t(9;22)(q34;q11)x2,add(18)(p11)x2"
#"der(X)t(X;11)(q22;q23)ins(X;11)(q22;p15p13)"
#"der(13)t(13;15)(q12;q22)t(15;22)(q26;q11)"
#"der(9)t(9;10)(q34;q22) or der(9)t(9;10)(q22;q11)"
#SVs <- c("der(7)t(5;7)(p11;p11)t(7;9)t(5;7)(q23;q34)","t(7;9)(q34;p11)")

#SVs <- unlist(strsplit(tttt, ","))
#N_SVs <- length(SVs)

##########packages###########
library(stringr)
library(readr)
library(installr)
library(utils)
#################functions#################
numbers_only <- function(x) !grepl("\\D", x)

paste_chr <- function(row, col, sv, data, SVj){
  for(i in row){  
  if(data[i,col] == 0){
    if(grepl(pattern = "x", x = SVj) == TRUE){
      sep <- unlist(strsplit(gsub("x","x~",SVj), "~" ))
      times <- sep[which(grepl(pattern = "x", x = sep))+1]
      data[i,col] <- paste(replicate(times, sv), collapse = ",")
    } else {
    data[i,col] <- sv
    }
  } else {
    if(grepl(pattern = "x", x = SVj) == TRUE){
      seps <- unlist(strsplit( gsub("x","x~",SVj), "~" ))
      times <- seps[which(grepl(pattern = "x", x = seps))+1]
      pr <- paste(replicate(times, sv), collapse = ",")
      data[i,col] <- paste(data[i,col], pr, sep=",") 
    } else {
    data[i,col] <- paste(data[i,col], sv, sep=",") 
    }
  }
 }
  return(data)
} 

#paste_chr <- function(row, col, sv, data){
#  for(i in row){  
#    if(data[i,col] == 0){
#      data[i,col] <- sv
#    } else {
#      data[i,col] <- paste(data[i,col], sv, sep=",") 
#    }
#  }
#  return(data)
#} 
qp_term <- function(ind, data){
  a <- data$chr[ind]  
  
  if(grepl(pattern = "p", x = data$band[ind]) == TRUE){
    indices <- c(which(data$chr == sprintf("%s",a))[1] : ind)  
  } else {
    indices <- c(ind : rev(which(data$chr == sprintf("%s",a)))[1])  
  }
  return(indices)
}

number_chr_parser <- function(seg){
  num <- vector(length = length(seg))
  for(j in 1:length(seg)){
  if(grepl(pattern = "X", x = seg[j]) == TRUE){
    num[j] <- "X"
  } else {if(grepl(pattern = "Y", x = seg[j]) == TRUE){
    num[j] <- "Y"
  } else {
    num[j] <- parse_number(seg[j])
         }
  }
  }
  return(num)
} 

####################################
#Bands <- read.csv("Bands.csv", header = TRUE)
Bands <- read.csv(sprintf("%s/Bands.csv", DATA_PATH), header = TRUE)

#karyotypes <- read.table("mitelman_karyshort.txt")
karyotypes <- read.table(sprintf("%s/mitelman_karyshort.txt", DATA_PATH))

#samples <- karyotypes[1:12,1]
#N <- length(samples)
samples <- karyotypes[,1]
N <- length(samples)
converted_table <- as.data.frame(Bands)
colnames(converted_table) <- c("chr", "band", "st", "end")
#########constructing data table#########
cloned_samples <- list()
cloned_numbers <- list()
######check if it is a clone sample
for(i in 1:N){
 # if(grepl(pattern = "\\/", x = samples[i]) == FALSE){
  if(grepl(pattern = "\\/", x = samples[i]) == FALSE){
    colname <- c(sprintf("gain.%s",i), sprintf("loss.%s",i), sprintf("fusion.%s",i))
    new_kar <- data.frame(matrix(0, ncol = 3, nrow = nrow(Bands)))
    colnames(new_kar) <- colname
    converted_table <- cbind(converted_table, new_kar)
  } else {
    ####separating clones
    cloned_samples <- c(cloned_samples, i)
    clones <- unlist(strsplit(samples[i],"/"))
    N_clones <- length(clones)
    cloned_numbers <- c(cloned_numbers,N_clones)
    ##### check if / comes at the end (/ without cloning)
    if(N_clones < 2){
      colname <- c(sprintf("gain.%s",i), sprintf("loss.%s",i), sprintf("fusion.%s",i))
      new_kar <- data.frame(matrix(0, ncol = 3, nrow = nrow(Bands)))
      colnames(new_kar) <- colname
      converted_table <- cbind(converted_table, new_kar)  
    } else{
    # colname <- c(sprintf("gain.%s.",i), sprintf("loss.%s",i), sprintf("fusion.%s",i))

    ####new clones(as separate samples) name in table
    clones_Sv_type <- c(sprintf("gain.%s.1",i), sprintf("loss.%s.1",i),sprintf("fusion.%s.1",i))
    for(k in 2:N_clones){
      clones_Sv_type <- c(clones_Sv_type, sprintf("gain.%s.%s",i,k), sprintf("loss.%s.%s",i,k),sprintf("fusion.%s.%s",i,k))  
    }
    new_kar <- data.frame(matrix(0, ncol = length(clones_Sv_type), nrow = nrow(Bands)))
    colnames(new_kar) <- clones_Sv_type
    converted_table <- cbind(converted_table, new_kar)
    }
  }
#  i = i+N_clones-1
}
cloned_samples <- unlist(cloned_samples)
cloned_numbers <- unlist(cloned_numbers)

######create new samples
i=0
for (j in cloned_samples) {
  j = j+i
  clones <- unlist(strsplit(samples[j],"/"))
  N_clones <- length(clones)
  #####replacing clones with one sample
  samples <- append(samples, clones, after = j)
  samples <- samples[-j]
  i = i+N_clones-1
}

####################################################
N <- length(samples)

error_table <- data.frame(Sample_ID=integer(),
                          Var_ID=integer(),
                          Karyotype=character(),
                          Error_type=character())

uncertain_table <- data.frame(Sample_ID=integer(),
                              Var_ID=integer(),
                              Karyotype=character(),
                              Uncertainty=character())

sex_chr <- matrix(0, nrow = N, ncol = 2)
colnames(sex_chr) <- c("chr_number", "sex")
empty_samples1 <- list()
N_SVs <- vector(length = N)
N_var <- list()
######## Start 
    
for(i in 1:N){
  ##### separate by "," ######
  seperated_k <- unlist(strsplit(samples[i], ","))
  chr_N <- seperated_k[1] #### Number of chrs
  sex_chr[i,1] <- chr_N
      
  ####determine sex ####
  if(grepl(pattern = "idem", x = samples[i]) == TRUE){
    sex_chr[i,2] <- "idem" 
  } else { if(grepl(pattern = "X", x = samples[i]) == TRUE & grepl(pattern = "Y", x = samples[i]) == TRUE){
    sex_chr[i,2] <- "XY"
  } else {
    sex_chr[i,2] <- "XX"   
    }
  } 
  
  ###### check if second element is sex
  if(grepl(pattern = "idem", x = seperated_k[2]) == FALSE && grepl(pattern = "X", x = seperated_k[2]) == FALSE && grepl(pattern = "Y", x = seperated_k[2]) == FALSE){
    error_table <- rbind(error_table, c(i ,0, samples[i], "incorrect format (chr NO/sex)"))
    next
  }    

  #### keep variations after deleting chr number and sex #### 
  SVs <- seperated_k[c(-1,-2)] ####all variations
  N_SVs[i] <- length(SVs)
  #### check to have at least one SV ####
  if(N_SVs[i] == 0){
    empty_samples1 <- rbind(empty_samples1, i) 
    next
  }
  
for(j in 1:N_SVs[i]){
  ###### delete "/" signs #####

  SVs[j] <- gsub("\\/.*","",SVs[j])
  
  ###### delete SVs with ? ####  
  if(grepl(pattern = "\\?", x = SVs[j]) == TRUE){
    error_table <- rbind(error_table, c(i , j, SVs[j], "?")) 
    next
  }
  
  ###### delete SVs with or ####  
  if(grepl(pattern = "or", x = SVs[j]) == TRUE){
    uncertain_table <- rbind(uncertain_table, c(i , j, SVs[j], "or")) 
    next
  } 
  
  ###### save and delete marker chromosome ####
  if(grepl(pattern = "mar", x = SVs[j]) == TRUE){
    uncertain_table <- rbind(uncertain_table, c(i , j, SVs[j], "marker chromosome")) 
    next
  }
  
  ##### save and delete ring chromosome #####
  if((grepl(pattern = "r", x = SVs[j]) == TRUE) && (grepl(pattern = "mar", x = SVs[j]) == FALSE)
     && (grepl(pattern = "der", x = SVs[j]) == FALSE) && (grepl(pattern = "hsr", x = SVs[j]) == FALSE)){
    uncertain_table <- rbind(uncertain_table, c(i , j, SVs[j], "ring chromosome")) 
    next
  }
  
  ##### check x then nothing #####
  if(grepl(pattern = "x", x = SVs[j]) == TRUE){
    seps <- unlist(strsplit(gsub("x","x~",SVs[j]), "~" ))
    multi <- seps[which(grepl(pattern = "x", x = seps))+1]
    if(is.na(multi) == TRUE){
     SVs[j] <- gsub("x","",SVs[j])
    }  
  }
  
  ##### delete x( ######
  if(grepl(pattern = "x\\(", x = SVs[j]) == TRUE){
    error_table <- rbind(error_table, c(i , j, SVs[j], "typo")) 
    next
  } 
  
##### check - in bands ####
  sep_terms <- unlist(strsplit( gsub(")",")~",SVs[j]), "~" ))
  Nk <- length(sep_terms)
  band_error <- list()
  if(Nk > 1){
    for(k in 2:Nk){
      if(grepl(pattern = "\\-", x = sep_terms[k]) == TRUE){
        uncertain_table <- rbind(uncertain_table, c(i , j, SVs[j], "uncertain band"))
        band_error <- 1
        next
      }
    }
  }
  if(is.empty(band_error) == FALSE){
    next
  }

##### + cases ####  
  pos_Svs <- c()
  if(grepl(pattern = "\\+", x = SVs[j]) == TRUE){
    pos_Svs <- gsub("\\+","",SVs[j]) ###delete +
  }
  if(is.null(pos_Svs) == FALSE){
    if((numbers_only(pos_Svs) == TRUE) || (pos_Svs == "X") || (pos_Svs == "Y")){
      converted_table <- paste_chr(which(converted_table$chr == sprintf("chr%s",pos_Svs))
                                   ,(i-1)*3+5,"+", converted_table, SVs[j])
      
    } else { 
    #### + chromosome  
    pos_chr <- gsub("[\\(\\)]", "", regmatches(pos_Svs, gregexpr("\\(.*?\\)", pos_Svs))[[1]])[1] 
#    converted_table <- paste_chr(which(converted_table$chr == sprintf("chr%s",pos_chr))
#                                   ,sprintf("gain.%s",i),"+", converted_table)
    if(grepl(pattern = ";", x = pos_chr)){
      pos_chrs <- unlist(strsplit(pos_chr, ";"))
      converted_table <- paste_chr(which(converted_table$chr == sprintf("chr%s",pos_chrs[1]))
                                   ,(i-1)*3+5,"+", converted_table, SVs[j])
      converted_table <- paste_chr(which(converted_table$chr == sprintf("chr%s",pos_chrs[2]))
                                   ,(i-1)*3+5,"+", converted_table, SVs[j])
    } else {
    converted_table <- paste_chr(which(converted_table$chr == sprintf("chr%s",pos_chr))
                                 ,(i-1)*3+5,"+", converted_table, SVs[j])
    }
    }  
  }    
  
##### - cases ####  
  neg_Svs <- c()
  if(grepl(pattern = "\\-", x = SVs[j]) == TRUE){
    neg_Svs <- gsub("\\-","",SVs[j]) ###delete -
  }
  if(is.null(neg_Svs) == FALSE){
    if((numbers_only(neg_Svs) == TRUE) || (neg_Svs == "X") || (neg_Svs == "Y")){
      converted_table <- paste_chr(which(converted_table$chr == sprintf("chr%s",neg_Svs))
                                   ,(i-1)*3+6 ,"-", converted_table, SVs[j])
    } else { 
      #### - chromosome  
      neg_chr <- gsub("[\\(\\)]", "", regmatches(neg_Svs, gregexpr("\\(.*?\\)", neg_Svs))[[1]])[1] 
      
      if(grepl(pattern = ";", x = neg_chr)){
        neg_chrs <- unlist(strsplit(neg_chr, ";"))
        converted_table <- paste_chr(which(converted_table$chr == sprintf("chr%s",neg_chrs[1]))
                                     ,(i-1)*3+6,"-", converted_table, SVs[j])
        converted_table <- paste_chr(which(converted_table$chr == sprintf("chr%s",neg_chrs[2]))
                                     ,(i-1)*3+6,"-", converted_table, SVs[j])
      } else {
      converted_table <- paste_chr(which(converted_table$chr == sprintf("chr%s",neg_chr))
                                   ,(i-1)*3+6 ,"-", converted_table, SVs[j]) 
      }
    }
  }

##### delete + and - signs #####
  SVs[j] <- gsub("\\-","",SVs[j])
  SVs[j] <- gsub("\\+","",SVs[j])

##### translocation ####
  if(grepl(pattern = "t\\(", x = SVs[j])){ ### t cases
    sep_terms <- unlist(strsplit( gsub(")",")~",SVs[j]), "~" )) ###separate terms within each variation
    #### translocated chromosomes
    tr <- sep_terms[which(grepl(pattern = "t\\(", x = sep_terms))]
    tr_chr <- unlist(strsplit(tr, ";"))
    ###### check bands validity #####
    if(is.null(tr_chr)){
      error_table <- rbind(error_table, c(i , j, SVs[j], "t typo"))
      next
    }
    #
    tr_chr <- number_chr_parser(tr_chr)
    #### check typo in t chromosome
    if(length(tr_chr) != 2*(length(tr))){
      error_table <- rbind(error_table, c(i , j, SVs[j], "t typo"))
      next
    }
    #
    if(all(tr_chr %in% c(1:22, "X", "Y")) == FALSE){
      error_table <- rbind(error_table, c(i , j, SVs[j], "t typo"))
      next
    }
    ### t includes p or q?
    tr_terms <- sep_terms[which(grepl(pattern = "t\\(", x = sep_terms))+1]
    N_trans <- length(tr)
    
    for(k in 1:N_trans){
      #### check if der exist in t term
      if(grepl(pattern = "der", x = tr[k])){
        error_table <- rbind(error_table, c(i , j, SVs[j], "der/t typo"))
        next
      }
      #
      band_term <- list()
    if((grepl(pattern = "p", x = tr_terms[k])) || grepl(pattern = "q", x = tr_terms[k])){
      #### translocated bands
      band_term <- unlist(strsplit(gsub("[\\(\\)]", "", tr_terms[k]), ";")) 
      #### fusion indices for t
      t1_ind <- list()
      t2_ind <- list()
      t1_ind <- intersect(which(converted_table$chr == sprintf("chr%s",tr_chr[((2*k)-1)]))
                          ,which(converted_table$band == sprintf("%s",band_term[1])))
      t2_ind <- intersect(which(converted_table$chr == sprintf("chr%s",tr_chr[(2*k)]))
                          ,which(converted_table$band == sprintf("%s",band_term[2])))
    
      ###### check bands validity #####
      if(is.empty(t1_ind) == TRUE || is.empty(t2_ind) == TRUE){
        error_table <- rbind(error_table, c(i , j, SVs[j], "t invalid bands"))
        next
      }
      #
      converted_table <- paste_chr(t1_ind,(i-1)*3+7 ,sprintf("t%s",j), converted_table, SVs[j]) 
      converted_table <- paste_chr(t2_ind,(i-1)*3+7 ,sprintf("t%s",j), converted_table, SVs[j]) 
      #
      if(grepl(pattern = "der", x = SVs[j])){ ###der cases
        #### derivative chromosome
        der_chr <- sep_terms[which(grepl(pattern = "der", x = sep_terms))]
        der_chr <- number_chr_parser(der_chr)
        ####check if deriviate chr is one of the translocated ones
        if(der_chr == tr_chr[((2*k)-1)] || der_chr == tr_chr[(2*k)]){
        #### check typo in der chromosome
        if((length(which(der_chr == c(1:22, "X", "Y")))) != 1){
          error_table <- rbind(error_table, c(i , j, SVs[j], "der typo"))
          next
        }
        #### find der chr from trans chrs
        der_tr <- tr_chr[c(((2*k)-1),(2*k))]
        neg_tr_band <- band_term[which(der_tr == der_chr)]
        #### loss tr
        tr_neg_ind <- intersect(which(converted_table$chr == sprintf("chr%s",der_chr))
                                 ,which(converted_table$band == sprintf("%s",neg_tr_band))) ###get indices of loss tr
        ###### check bands validity #####
        if(is.empty(tr_neg_ind) == TRUE){
          error_table <- rbind(error_table, c(i , j, SVs[j], "t invalid bands"))
          next
        }
        #
        tr_neg_ind <- qp_term(tr_neg_ind, converted_table)
        converted_table <- paste_chr(tr_neg_ind,(i-1)*3+6 ,sprintf("t%s",j), converted_table, SVs[j]) 
        #### gain tr
        pos_tr_band <- band_term[which(der_tr != der_chr)]
        pos_tr <- der_tr[which(der_tr != der_chr)]
        tr_pos_ind <- intersect(which(converted_table$chr == sprintf("chr%s",pos_tr))
                                ,which(converted_table$band == sprintf("%s",pos_tr_band))) ###get indices of gain tr
        ###### check bands validity #####
        if(is.empty(tr_pos_ind) == TRUE){
          error_table <- rbind(error_table, c(i , j, SVs[j], "t invalid bands"))
          next
        }
        #
        tr_pos_ind <- qp_term(tr_pos_ind, converted_table)
        converted_table <- paste_chr(tr_pos_ind,(i-1)*3+5 ,sprintf("t%s",j), converted_table, SVs[j]) 
      } 
      }
    } else { if(grepl(pattern = "der", x = SVs[j])){
      #### derivative chromosome
      der_chr <- sep_terms[which(grepl(pattern = "der", x = sep_terms))]
      der_chr <- number_chr_parser(der_chr)
      ####check if deriviate chr is one of the translocated ones
      if(der_chr == tr_chr[((2*k)-1)] || der_chr == tr_chr[(2*k)]){
      
      #### check typo in der chromosome
      if((length(which(der_chr == c(1:22, "X", "Y")))) != 1){
        error_table <- rbind(error_table, c(i , j, SVs[j], "der typo"))
        next
      }
      #
      #### check if the same translocation id in another variation of the same karyotype exist
      trans_ind <- which(grepl(pattern = sprintf("t\\(%s;%s\\)",tr_chr[((2*k)-1)],tr_chr[(2*k)]), x = SVs))
      trans_ind <- SVs[trans_ind[which(trans_ind != j )]]
      #new_j <- which(SVs == trans_ind)
      #### filter terms with der
      trans_ind <- trans_ind[which(grepl(pattern = "der", x = trans_ind) == FALSE)]
      #if(is.empty(trans_ind) == FALSE){
        trans_terms <- list()
        trans_terms <- unlist(strsplit( gsub(")",")~",trans_ind), "~" ))
        #### translocated bands
        trans_bands <- trans_terms[which(trans_terms == tr[k])+1]
        ###### check bands validity #####
        if((is.null(trans_bands) == FALSE) && (is.empty(trans_bands) == FALSE) && 
           (grepl(pattern = "p", x = trans_bands) || grepl(pattern = "q", x = trans_bands))){
          tr_bands <- list()
          tr_bands <- unlist(strsplit(gsub("[\\(\\)]", "", trans_bands), ";"))
          #### find der chr from trans chrs
          der_tr <- tr_chr[c(((2*k)-1),(2*k))]
          neg_tr_band <- tr_bands[which(der_tr == der_chr)]
          #### loss tr
          tr_neg_ind <- intersect(which(converted_table$chr == sprintf("chr%s",der_chr))
                                  ,which(converted_table$band == sprintf("%s",neg_tr_band))) ###get indices of loss tr
          ###### check bands validity #####
          if(is.empty(tr_neg_ind) == TRUE){
            error_table <- rbind(error_table, c(i , j, SVs[j], "t bands missing"))
            next
          }
          #
          tr_neg_ind <- qp_term(tr_neg_ind, converted_table)
          converted_table <- paste_chr(tr_neg_ind,(i-1)*3+6 ,sprintf("t%s",j), converted_table, SVs[j]) 
          #### gain tr
          pos_tr_band <- tr_bands[which(der_tr != der_chr)]
          pos_chr <- der_tr[which(der_tr != der_chr)]
          tr_pos_ind <- intersect(which(converted_table$chr == sprintf("chr%s",pos_chr))
                                  ,which(converted_table$band == sprintf("%s",pos_tr_band))) ###get indices of gain tr
          ###### check bands validity #####
          if(is.empty(tr_pos_ind) == TRUE){
            error_table <- rbind(error_table, c(i , j, SVs[j], "t bands missing"))
            next
          }
          #
          tr_pos_ind <- qp_term(tr_pos_ind, converted_table)
          converted_table <- paste_chr(tr_pos_ind,(i-1)*3+5 ,sprintf("t%s",j), converted_table, SVs[j]) 
          
        } else {
          error_table <- rbind(error_table, c(i , j, SVs[j], "t bands missing"))  
        }
      } else {
        error_table <- rbind(error_table, c(i , j, SVs[j], "t bands missing"))   
      }
        ####
    } else {
      error_table <- rbind(error_table, c(i , j, SVs[j], "t bands missing"))   
      }
    }
  }   
}
    
##### add cases ####
  if(grepl(pattern = "add", x = SVs[j]) == TRUE){ ###add cases
  sep_terms <- unlist(strsplit( gsub(")",")~",SVs[j]), "~" )) ###separate terms within each variation
  #### origin chromosome of additional material
  sep1 <- sep_terms[which(grepl(pattern = "add", x = sep_terms))]
  add_chr <- number_chr_parser(sep1)   
  #### check typo in add chromosome
  if((length(which(add_chr == c(1:22, "X", "Y")))) != 1){
    error_table <- rbind(error_table, c(i , j, SVs[j], "add typo"))
    next
  } else {
  #### origin bands
  add_band <- sep_terms[which(grepl(pattern = "add", x = sep_terms))+1]
  add_band <- gsub("[\\(\\)]", "", add_band)
  #### add indix
  add_ind <- intersect(which(converted_table$chr == sprintf("chr%s",add_chr))
                      ,which(converted_table$band == sprintf("%s",add_band)))
  ###### check bands validity #####
  if(is.empty(add_ind) == TRUE){
    error_table <- rbind(error_table, c(i , j, SVs[j], "add invalid bands"))
    next
  }
  #
  ### add loss
  add_loss_ind <- qp_term(add_ind, converted_table)
  converted_table <- paste_chr(add_loss_ind,(i-1)*3+6 ,sprintf("add%s",j), converted_table, SVs[j])
  ### add fusion
  converted_table <- paste_chr(add_ind,(i-1)*3+7 ,sprintf("add%s",j), converted_table, SVs[j]) 
  }
}
  
##### convert ider to isochromosome ####
  if(grepl(pattern = "ider", x = SVs[j]) == TRUE){
    SVs[j] <- gsub("ider","i",SVs[j])
  }
  
##### isochromosome ####
  if(grepl(pattern = "i\\(", x = SVs[j]) == TRUE){ ### i cases
    
    sep_terms <- unlist(strsplit(gsub(")",")~",SVs[j]), "~" )) ###separate terms within each variation
    #### isochromosome chr
    sep1 <- sep_terms[which(grepl(pattern = "i\\(", x = sep_terms))]
    i_chr <- number_chr_parser(sep1)
    #### check typo in i chromosome
    if((length(which(i_chr == c(1:22, "X", "Y")))) != 1){
      error_table <- rbind(error_table, c(i , j, SVs[j], "iso typo"))
      next
    } else {
    #### origin bands
    i_code <- sep_terms[which(grepl(pattern = "i\\(", x = sep_terms))+1]
    i_code <- gsub("[\\(\\)]", "", i_code)
    
    ###### check if there is any isochromosome code #####
    if(is.na(i_code) == TRUE){
      error_table <- rbind(error_table, c(i , j, SVs[j], "iso typo"))
      next
    }
    #
    
    if(i_code == "p10"){
      i_gain_ind <- intersect(which(converted_table$chr == sprintf("chr%s",i_chr))
                               ,which(converted_table$band == "p11")) ###get indices of loss der
      i_gain_ind <- qp_term(i_gain_ind, converted_table)
      
      i_loss_ind <- intersect(which(converted_table$chr == sprintf("chr%s",i_chr))
                              ,which(converted_table$band == "q11")) ###get indices of loss de                                             
      i_loss_ind <- qp_term(i_loss_ind, converted_table)
      converted_table <- paste_chr(i_gain_ind,(i-1)*3+5 ,sprintf("i%s",j), converted_table, SVs[j])
      converted_table <- paste_chr(i_loss_ind,(i-1)*3+6 ,sprintf("i%s",j), converted_table, SVs[j])
      
    } else { if(i_code == "q10"){
      i_gain_ind <- intersect(which(converted_table$chr == sprintf("chr%s",i_chr))
                              ,which(converted_table$band == "q11")) ###get indices of loss der
      i_gain_ind <- qp_term(i_gain_ind, converted_table)
      
      i_loss_ind <- intersect(which(converted_table$chr == sprintf("chr%s",i_chr))
                              ,which(converted_table$band == "p11")) ###get indices of loss de                                             
      i_loss_ind <- qp_term(i_loss_ind, converted_table)
      converted_table <- paste_chr(i_gain_ind,(i-1)*3+5 ,sprintf("i%s",j), converted_table, SVs[j])
      converted_table <- paste_chr(i_loss_ind,(i-1)*3+6 ,sprintf("i%s",j), converted_table, SVs[j])
    } else{
      error_table <- rbind(error_table, c(i , j, SVs[j], "iso invalid bands"))
      next
      }
    }
  }
}
  
##### dup ####
  if(grepl(pattern = "dup", x = SVs[j]) == TRUE){ ###dup cases
    
    sep_terms <- unlist(strsplit(gsub(")",")~",SVs[j]), "~" )) ###separate terms within each variation
    #### dup chr
    sep1 <- sep_terms[which(grepl(pattern = "dup", x = sep_terms))]
    dup_chr <- number_chr_parser(sep1)  
    #### check typo in dup chromosome
    if((length(which(dup_chr == c(1:22, "X", "Y")))) != 1){
      error_table <- rbind(error_table, c(i , j, SVs[j], "dup typo"))
      next
    } #else {
    #### origin bands
    dup_band <- sep_terms[which(grepl(pattern = "dup", x = sep_terms))+1]
    dup_band <- gsub("[\\(\\)]", "", dup_band)
    dup_band1 <- parse_number(dup_band)
    dup_bands <- unlist(strsplit(gsub(dup_band1,sprintf("%s~",dup_band1),dup_band), "~" ))
    ###### dup gain bands
    dup_ind1 <- intersect(which(converted_table$chr == sprintf("chr%s",dup_chr))
                          ,which(converted_table$band == sprintf("%s",dup_bands[1]))) ###get index of first dup band
    
    dup_ind2 <- intersect(which(converted_table$chr == sprintf("chr%s",dup_chr))
                          ,which(converted_table$band == sprintf("%s",dup_bands[2]))) ###get index of second dup band
    ###### check bands validity #####
    if(is.empty(dup_ind1) == TRUE || is.empty(dup_ind2) == TRUE){
      error_table <- rbind(error_table, c(i , j, SVs[j], "dup invalid bands"))
      next
    }
    #
    ##### gain dup
    converted_table <- paste_chr((dup_ind1:dup_ind2),(i-1)*3+5 ,sprintf("dup%s",j), converted_table, SVs[j])
    #### fusion dup
    converted_table <- paste_chr(dup_ind1,(i-1)*3+7 ,sprintf("dup%s",j), converted_table, SVs[j])
    converted_table <- paste_chr(dup_ind2,(i-1)*3+7 ,sprintf("dup%s",j), converted_table, SVs[j])
 # }    
}
  
#### del ####
  if(grepl(pattern = "del", x = SVs[j]) == TRUE){ ###del cases
    
    sep_terms <- unlist(strsplit(gsub(")",")~",SVs[j]), "~" )) ###separate terms within each variation
    #### del chr
    sep1 <- sep_terms[which(grepl(pattern = "del", x = sep_terms))]
    del_chr <- number_chr_parser(sep1) 
    #### check typo in del chromosome
    if((length(which(del_chr == c(1:22, "X", "Y")))) != 1){
      error_table <- rbind(error_table, c(i , j, SVs[j], "del typo"))
      next
    } #else {
    #### origin bands
    del_band <- sep_terms[which(grepl(pattern = "del", x = sep_terms))+1]
    del_band <- gsub("[\\(\\)]", "", del_band)
    del_band1 <- parse_number(del_band)
    del_bands <- unlist(strsplit(gsub(del_band1,sprintf("%s~",del_band1),del_band), "~" ))
    #### check if bands are empty ####
    if(is.na(del_bands) == TRUE){
      error_table <- rbind(error_table, c(i , j, SVs[j], "del invalid bands"))
      next
    }
    #### check if del includes two bands
    if((str_count(del_band, "p") + str_count(del_band, "q")) == 2){
      ###### del bands indices
      del_ind1 <- intersect(which(converted_table$chr == sprintf("chr%s",del_chr))
                            ,which(converted_table$band == sprintf("%s",del_bands[1]))) ###get index of first del band
      
      del_ind2 <- intersect(which(converted_table$chr == sprintf("chr%s",del_chr))
                            ,which(converted_table$band == sprintf("%s",del_bands[2]))) ###get index of second del band
      ###### check bands validity #####
      if(is.empty(del_ind1) == TRUE || is.empty(del_ind2) == TRUE){
        error_table <- rbind(error_table, c(i , j, SVs[j], "del invalid bands"))
        next
      }
      #
      ##### loss del
      converted_table <- paste_chr((del_ind1:del_ind2),(i-1)*3+6 ,sprintf("del%s",j), converted_table, SVs[j])
      #### fusion del
      converted_table <- paste_chr(del_ind1,(i-1)*3+7 ,sprintf("del%s",j), converted_table, SVs[j])
      converted_table <- paste_chr(del_ind2,(i-1)*3+7 ,sprintf("del%s",j), converted_table, SVs[j]) 
      #### check if del includes one band
    } else {if((str_count(del_band, "p") + str_count(del_band, "q")) == 1){
      ###### del band index
      del_ind <- intersect(which(converted_table$chr == sprintf("chr%s",del_chr))
                            ,which(converted_table$band == sprintf("%s",del_bands))) ###get index of first del band
      ###### check bands validity #####
      if(is.empty(del_ind) == TRUE || is.na(del_ind) == TRUE){
        error_table <- rbind(error_table, c(i , j, SVs[j], "del invalid bands"))
        next
      }
      #
      del_ind <- qp_term(del_ind, converted_table)
      ##### loss del
      converted_table <- paste_chr(del_ind,(i-1)*3+6 ,sprintf("del%s",j), converted_table, SVs[j])
      ##### fusion del
      #converted_table <- paste_chr(del_ind1,(3*(i-1)+7),sprintf("del%s",j), converted_table)
      #converted_table <- paste_chr(del_ind2,(3*(i-1)+7),sprintf("del%s",j), converted_table) 
      } else {
        error_table <- rbind(error_table, c(i , j, SVs[j], "del invalid bands")) 
      }
    }
  #}  
}
  
##### inv ####
  if(grepl(pattern = "inv", x = SVs[j]) == TRUE){ ###inv cases
    
    sep_terms <- unlist(strsplit(gsub(")",")~",SVs[j]), "~" )) ###separate terms within each variation
    #### inv chr
    sep1 <- sep_terms[which(grepl(pattern = "inv", x = sep_terms))]
    inv_chr <- number_chr_parser(sep1) 
    #### check typo in inv chromosome
    if((length(which(inv_chr == c(1:22, "X", "Y")))) != 1){
      error_table <- rbind(error_table, c(i , j, SVs[j], "inv typo"))
      next
    } #else {
    #### origin bands
    inv_band <- sep_terms[which(grepl(pattern = "inv", x = sep_terms))+1]
    inv_band <- gsub("[\\(\\)]", "", inv_band)
    inv_band1 <- parse_number(inv_band)
    inv_bands <- unlist(strsplit(gsub(inv_band1,sprintf("%s~",inv_band1),inv_band), "~" ))
    ###### inv bands
    inv_ind1 <- intersect(which(converted_table$chr == sprintf("chr%s",inv_chr))
                          ,which(converted_table$band == sprintf("%s",inv_bands[1]))) ###get index of first inv band
    
    inv_ind2 <- intersect(which(converted_table$chr == sprintf("chr%s",inv_chr))
                          ,which(converted_table$band == sprintf("%s",inv_bands[2]))) ###get index of second inv band
    #### check typo in inv bands
    if(is.empty(inv_ind1) == TRUE || is.empty(inv_ind2) == TRUE){
      error_table <- rbind(error_table, c(i , j, SVs[j], "inv invalid bands"))
      next
    #  
    } #else {
    #### fusion inv
    converted_table <- paste_chr(inv_ind1,(i-1)*3+7 ,sprintf("inv%s",j), converted_table, SVs[j])
    converted_table <- paste_chr(inv_ind2,(i-1)*3+7 ,sprintf("inv%s",j), converted_table, SVs[j])
    #}
  #}
}
     
##### ins ####   
  if(grepl(pattern = "ins", x = SVs[j]) == TRUE){ ###ins cases
    sep_terms <- unlist(strsplit(gsub(")",")~",SVs[j]), "~" )) ###separate terms within each variation
    ####insertion chromosomes
    ins <- sep_terms[which(grepl(pattern = "ins\\(", x = sep_terms))]
    ins_chr <- unlist(strsplit(ins, ";"))
    ###### check bands validity #####
    if(is.null(ins_chr) == TRUE){
      error_table <- rbind(error_table, c(i , j, SVs[j], "ins typo"))
      next
    }
    #
    ins_chr1 <- number_chr_parser(ins_chr[1])
    ins_chr2 <- number_chr_parser(ins_chr[2])
    #### check typo in ins chromosome
    if((length(which(ins_chr1 == c(1:22, "X", "Y"))) + (length(which(ins_chr2 == c(1:22, "X", "Y"))))) != 2){
      error_table <- rbind(error_table, c(i , j, SVs[j], "ins typo"))
      next
    } else {
    ####origin bands
    ins_band <- sep_terms[which(grepl(pattern = "ins\\(", x = sep_terms))+1]
    ins_bands <- unlist(strsplit(gsub("[\\(\\)]", "", ins_band), ";")) ####
    ######where the inserted segment is inserted
    ins_bands1 <- ins_bands[1]
    #####where the inserted segment is coming from
    ins_bands2 <- number_chr_parser(ins_bands[2])
    ins_bands2 <- unlist(strsplit(gsub(ins_bands2,sprintf("%s~",ins_bands2),ins_bands[2]), "~" ))
    ###### ins "to" bands
    ins_ind_to <- intersect(which(converted_table$chr == sprintf("chr%s",ins_chr1))
                          ,which(converted_table$band == sprintf("%s",ins_bands1))) ###get index of ins_to band
    #### check if there is any invalid bands   
    if(is.empty(ins_ind_to) == TRUE){
      error_table <- rbind(error_table, c(i , j, SVs[j], "ins invalid bands"))
      next
    }
    #
    #### fusion ins_to
    converted_table <- paste_chr(ins_ind_to,(i-1)*3+7 ,sprintf("insT%s",j), converted_table, SVs[j])
    ###### ins "from" bands
    ins_ind_from1 <- intersect(which(converted_table$chr == sprintf("chr%s",ins_chr2))
                            ,which(converted_table$band == sprintf("%s",ins_bands2[1]))) ###get index of first ins_from band
    ins_ind_from2 <- intersect(which(converted_table$chr == sprintf("chr%s",ins_chr2))
                               ,which(converted_table$band == sprintf("%s",ins_bands2[2]))) ###get index of second ins_from band
    #### check if there is any invalid bands   
    if(is.empty(ins_ind_from1) == TRUE || is.empty(ins_ind_from2) == TRUE){
      error_table <- rbind(error_table, c(i , j, SVs[j], "ins invalid bands"))
      next
    }
    #
    #### fusion ins_from
    converted_table <- paste_chr(ins_ind_from1,(i-1)*3+7 ,sprintf("insF%s",j), converted_table, SVs[j])
    converted_table <- paste_chr(ins_ind_from2,(i-1)*3+7 ,sprintf("insF%s",j), converted_table, SVs[j])
    
  }
} 

##### idic ####  
  if(grepl(pattern = "idic", x = SVs[j]) == TRUE){ ###idic cases
    sep_terms <- unlist(strsplit(gsub(")",")~",SVs[j]), "~" )) ###separate terms within each variation
    ####isodicentric chromosomes
    idic <- sep_terms[which(grepl(pattern = "idic\\(", x = sep_terms))]
    idic_chr <- number_chr_parser(idic)
    #### check typo in idic chromosome
    if((length(which(idic_chr == c(1:22, "X", "Y")))) != 1){
      error_table <- rbind(error_table, c(i , j, SVs[j], "idic typo"))
      next
    } else {    
    #### origin bands
    idic_band <- sep_terms[which(grepl(pattern = "idic", x = sep_terms))+1]
    idic_band <- gsub("[\\(\\)]", "", idic_band)
    #### idic indix
    idic_ind <- intersect(which(converted_table$chr == sprintf("chr%s",idic_chr))
                           ,which(converted_table$band == sprintf("%s",idic_band)))
    idic_chr_ind <- which(converted_table$chr == sprintf("chr%s",idic_chr))
    ###### check bands validity #####
    if(is.empty(idic_ind) == TRUE){
      error_table <- rbind(error_table, c(i , j, SVs[j], "idic invalid bands"))
      next
    }
    #
    ### check idic arm
    if(grepl(pattern = "p", x = idic_band) == TRUE){
      idic_loss_ind <- qp_term(idic_ind, converted_table)
      idic_gain_ind <- c(idic_ind:(tail(idic_chr_ind,n=1)))
      converted_table <- paste_chr(idic_gain_ind,(i-1)*3+5 ,sprintf("idic%s",j), converted_table, SVs[j])
      converted_table <- paste_chr(idic_loss_ind,(i-1)*3+6 ,sprintf("idic%s",j), converted_table, SVs[j])
      converted_table <- paste_chr(idic_ind,(i-1)*3+7 ,sprintf("idic%s",j), converted_table, SVs[j]) 
    } else { if(grepl(pattern = "q", x = idic_band) == TRUE){
      idic_loss_ind <- qp_term(idic_ind, converted_table)
      idic_gain_ind <- c(idic_chr_ind[1]:idic_ind)
      converted_table <- paste_chr(idic_gain_ind,(i-1)*3+5 ,sprintf("idic%s",j), converted_table, SVs[j])
      converted_table <- paste_chr(idic_loss_ind,(i-1)*3+6 ,sprintf("idic%s",j), converted_table, SVs[j])
      converted_table <- paste_chr(idic_ind,(i-1)*3+7 ,sprintf("idic%s",j), converted_table, SVs[j]) 
    } else {
      error_table <- rbind(error_table, c(i , j, SVs[j], "idic invalid bands"))
      next
    }
    }
  }
}    

##### dic ####  
  if((grepl(pattern = "dic", x = SVs[j]) == TRUE) && (grepl(pattern = "idic", x = SVs[j]) == FALSE)){ ###dic cases
    sep_terms <- unlist(strsplit(gsub(")",")~",SVs[j]), "~" )) ###separate terms within each variation
    #### dicentric chromosomes
    dic <- sep_terms[which(grepl(pattern = "dic", x = sep_terms))]
    dic_chr <- unlist(strsplit(dic, ";"))
    ###### check bands validity #####
    if(is.null(dic_chr) == TRUE){
      error_table <- rbind(error_table, c(i , j, SVs[j], "dic typo"))
      next
    }
    #
    dic_chr1 <- number_chr_parser(dic_chr[1])
    dic_chr2 <- number_chr_parser(dic_chr[2])    
    #### check typo in dic chromosome
    if((length(which(dic_chr1 == c(1:22, "X", "Y"))) + length(which(dic_chr2 == c(1:22, "X", "Y")))) != 2){
      error_table <- rbind(error_table, c(i , j, SVs[j], "dic typo"))
      next
    } else {
      ### dic includes p or q?
      band_term <- sep_terms[which(grepl(pattern = "dic", x = sep_terms))+1]
      if((grepl(pattern = "p", x = band_term) || grepl(pattern = "q", x = band_term))){ 
        #### dicentric bands
        dic2 <- sep_terms[which(grepl(pattern = "dic", x = sep_terms))+1]
        dic_bands <- unlist(strsplit(gsub("[\\(\\)]", "", dic2), ";"))
        #### fusion indices for dic
        dic1_ind <- intersect(which(converted_table$chr == sprintf("chr%s",dic_chr1))
                            ,which(converted_table$band == sprintf("%s",dic_bands[1])))
        dic2_ind <- intersect(which(converted_table$chr == sprintf("chr%s",dic_chr2))
                            ,which(converted_table$band == sprintf("%s",dic_bands[2])))
        ###### check bands validity #####
        if(is.empty(dic1_ind) == TRUE || is.empty(dic2_ind) == TRUE){
          error_table <- rbind(error_table, c(i , j, SVs[j], "dic invalid bands"))
          next
        }
        #
        ###### dic fusion
        converted_table <- paste_chr(dic1_ind,(i-1)*3+7 ,sprintf("dice%s",j), converted_table, SVs[j]) 
        converted_table <- paste_chr(dic2_ind,(i-1)*3+7 ,sprintf("dice%s",j), converted_table, SVs[j]) 
        ###### dic loss
        dic1_loss_ind <- qp_term(dic1_ind, converted_table)
        dic2_loss_ind <- qp_term(dic2_ind, converted_table)
        converted_table <- paste_chr(dic1_loss_ind,(i-1)*3+6 ,sprintf("dice%s",j), converted_table, SVs[j])
        converted_table <- paste_chr(dic2_loss_ind,(i-1)*3+6 ,sprintf("dice%s",j), converted_table, SVs[j])
      } else { 
        error_table <- rbind(error_table, c(i , j, SVs[j], "dic invalid bands"))
        next
      }
    }        
}

##### tas ####  
  if(grepl(pattern = "tas", x = SVs[j])){ ###tas cases
    sep_terms <- unlist(strsplit(gsub(")",")~",SVs[j]), "~" )) ###separate terms within each variation
    #### telomeric chromosomes
    tas <- sep_terms[which(grepl(pattern = "tas", x = sep_terms))]
    tas_chr <- unlist(strsplit(tas, ";"))
    ###### check bands validity #####
    if(is.null(tas_chr)){
      error_table <- rbind(error_table, c(i , j, SVs[j], "tas typo"))
      next
    }
    # 
    N_tas <- length(tas_chr)
    ##### check tas chromosome numbers
    if(N_tas < 2){
      error_table <- rbind(error_table, c(i , j, SVs[j], "tas typo"))
      next
    }
    ################ tas chrs
    tas_chrs <- vector(length = N_tas)
    for(k in 1:N_tas){
    tas_chrs[k] <- number_chr_parser(tas_chr[k])
    }
    #### tas bands
    tas2 <- sep_terms[which(grepl(pattern = "tas", x = sep_terms))+1]
    tas_bands <- unlist(strsplit(gsub("[\\(\\)]", "", tas2), ";"))
    N_bands <- length(tas_bands)
    #### band numbers should be equal or greater than chrs ####
    if(N_tas != N_bands){
      error_table <- rbind(error_table, c(i , j, SVs[j], "tas invalid bands"))
      next
    }
    #### check chr-band separately to see if they are terminal ####
    terminal <- list()
    tas2_ind1 <- list()
    tas2_ind2 <- list()
    tas_ind_p <- list()
    tas_ind_q <- list()
    for (l in 1:N_bands) {
      if((grepl(pattern = "p", x = tas_bands[l]) == TRUE) && (grepl(pattern = "q", x = tas_bands[l]) == TRUE)){
        tas_bands2 <- number_chr_parser(tas_bands[l])
        tas_bands2 <- unlist(strsplit(gsub(tas_bands2,sprintf("%s~",tas_bands2),tas_bands[l]), "~" ))
        tas_chr_ind <- which(converted_table$chr == sprintf("chr%s",tas_chrs[l]))
        tas2_ind1 <- c(tas2_ind1, intersect(tas_chr_ind,which(converted_table$band == sprintf("%s",tas_bands2[1]))))
        tas2_ind2 <- c(tas2_ind2, intersect(tas_chr_ind,which(converted_table$band == sprintf("%s",tas_bands2[2]))))
        ###### check bands validity #####
        if(is.empty(tas2_ind1) == TRUE || is.empty(tas2_ind2) == TRUE){
          error_table <- rbind(error_table, c(i , j, SVs[j], "tas invalid bands"))
          next
        }
        #### check if the band is terminal
        if(rev(tas2_ind1)[1] != tas_chr_ind[1] && rev(tas2_ind1)[1] != rev(tas_chr_ind)[1]){
          error_table <- rbind(error_table, c(i , j, SVs[j], "tas band is not terminal"))
          terminal <- 1
          next
        }
        if(rev(tas2_ind2)[1] != tas_chr_ind[1] && rev(tas2_ind2)[1] != rev(tas_chr_ind)[1]){
          error_table <- rbind(error_table, c(i , j, SVs[j], "tas band is not terminal"))
          terminal <- 1
          next
        }
        #else {
          ###### tas fusion
         # converted_table <- paste_chr(tas2_ind1,(i-1)*3+7 ,sprintf("tas%s",j), converted_table, SVs[j]) 
          #converted_table <- paste_chr(tas2_ind2,(i-1)*3+7 ,sprintf("tas%s",j), converted_table, SVs[j]) 
        #}

      } else {if(grepl(pattern = "p", x = tas_bands[l])){
        tas_chr_ind <- which(converted_table$chr == sprintf("chr%s",tas_chrs[l]))
        tas_ind_p <- c(tas_ind_p, intersect(tas_chr_ind,which(converted_table$band == sprintf("%s",tas_bands[l]))))
        ###### check bands validity #####
        if(is.empty(tas_ind_p)){
          error_table <- rbind(error_table, c(i , j, SVs[j], "tas invalid bands"))
          next
        }
        #### check if the band is terminal
        if(rev(tas_ind_p)[1] != tas_chr_ind[1]){
          error_table <- rbind(error_table, c(i , j, SVs[j], "tas band is not terminal"))
          terminal <- 1
          next
        } 
        #else {
          ###### tas fusion
         # converted_table <- paste_chr(tas_ind,(i-1)*3+7 ,sprintf("tas%s",j), converted_table, SVs[j]) 
        #}  
        
      } else {if(grepl(pattern = "q", x = tas_bands[l])){
        tas_chr_ind <- which(converted_table$chr == sprintf("chr%s",tas_chrs[l]))
        tas_ind_q <- c(tas_ind_q, intersect(tas_chr_ind,which(converted_table$band == sprintf("%s",tas_bands[l]))))
        ###### check bands validity #####
        if(is.empty(tas_ind_q)){
          error_table <- rbind(error_table, c(i , j, SVs[j], "tas invalid bands"))
          next
        }
        #
        #### check if the band is terminal
        if(rev(tas_ind_q) != rev(tas_chr_ind)[1]){
          error_table <- rbind(error_table, c(i , j, SVs[j], "tas band is not terminal"))
          terminal <- 1
          next
        } 
        #else {
          ###### tas fusion
         # converted_table <- paste_chr(tas_ind,(i-1)*3+7 ,sprintf("tas%s",j), converted_table, SVs[j]) 
        #}
        
      } else {
        error_table <- rbind(error_table, c(i , j, SVs[j], "tas band typo"))
        terminal <- 1
        next
      }
      }
    }
    }
    ####### print to table if there is no error ####
    if(is.empty(terminal)){
      tas2_ind1 <- unlist(tas2_ind1)
      tas2_ind2 <- unlist(tas2_ind2)
      tas_ind_p <- unlist(tas_ind_p)
      tas_ind_q <- unlist(tas_ind_q)
      tas_indices <- c(tas2_ind1, tas2_ind2, tas_ind_p, tas_ind_q)
      converted_table <- paste_chr(tas_indices,(i-1)*3+7 ,sprintf("tas%s",j), converted_table, SVs[j]) 
      
    }
}   
  
##### hsr cases ####
  if(grepl(pattern = "hsr", x = SVs[j])){ ###add cases
    sep_terms <- unlist(strsplit( gsub(")",")~",SVs[j]), "~" )) ###separate terms within each variation  
    #### origin chromosome of homogeneously staining region
    sep1 <- sep_terms[which(grepl(pattern = "hsr", x = sep_terms))]
    hsr_chr <- number_chr_parser(sep1) 
    #### check typo in hsr chromosome
    if((length(which(hsr_chr == c(1:22, "X", "Y")))) != 1){
      error_table <- rbind(error_table, c(i , j, SVs[j], "hsr typo"))
      next
    } else {
      #### origin bands
      hsr_band <- sep_terms[which(grepl(pattern = "hsr", x = sep_terms))+1]
      hsr_band <- gsub("[\\(\\)]", "", hsr_band)
      #### hsr indix
      hsr_ind <- intersect(which(converted_table$chr == sprintf("chr%s",hsr_chr))
                           ,which(converted_table$band == sprintf("%s",hsr_band)))
      ###### check bands validity #####
      if(is.empty(hsr_ind)){
        error_table <- rbind(error_table, c(i , j, SVs[j], "hsr invalid bands"))
        next
      }
      #
      ### hsr fusion
      converted_table <- paste_chr(hsr_ind,(i-1)*3+5 ,sprintf("hsr%s",j), converted_table, SVs[j]) 
    }
  }  
    
#######################     
}      
  cat(sprintf("sample number %s", i), sep="\n")
}      
colnames(error_table) <- c("sample_ID", "Variation_ID", "variation", "error")
colnames(uncertain_table) <- c("sample_ID", "Variation_ID", "variation", "error")

########### colname calculation ######
#sample_nmber <- 2762
#which(cloned_samples == sample_nmber)
#tt <- converted_table[,c("chr", "band", "st", "end", "gain.2762", "loss.2762", "fusion.2762")]
#ttt <- cbind(Bands,converted_table[,176:178]) 


#cloned_index <- tail(which(cloned_samples < sample_nmber),n=1)
#sum_clones <- sum(cloned_numbers[1:cloned_index])
#sample_index <- sample_nmber+(sum_clones-cloned_index)



#samples[61]
# cloned_samples[cloned_index]

###### filter no variation samples ######
empty_samples1 <- unlist(empty_samples1) 

empty_samples2 <- list()
for(i in 1:N){
   if(all(converted_table[,((i-1)*3+5):((i-1)*3+7)] == 0)){
     empty_samples2 <- c(empty_samples2, i) 
   }
}
empty_samples2 <- unlist(empty_samples2)

#### delete empty samples from result table ###
empty_indices <- c(((empty_samples2-1)*3+5),((empty_samples2-1)*3+6),((empty_samples2-1)*3+7))
converted_table <- converted_table[,-empty_indices]

###### differentiate before and after empty samples
if(length(empty_samples1) > 0){
  ind <- vector(length = length(empty_samples1))
  for(i in 1:length(empty_samples1)){
  ind[i] <- which(empty_samples2 == empty_samples1[i])
  }
} else {
  ind <- 0
}
empty_samples2 <- empty_samples2[-ind]
#samples[85]

####### summary ######
summary <- list()
summary$N_SVs <- N_SVs
summary$samples <- samples
summary$no_variation_before <- empty_samples1
summary$no_variation_after <- empty_samples2
summary$result <- converted_table
summary$error <- error_table
summary$uncertain <- uncertain_table
summary$sex_chr <- sex_chr
summary$cloned_ID <- cloned_samples
summary$cloned_NO <- cloned_numbers
########### save ##########
#write(summary,"Mitelman_summary_1_15000.R")
save(summary, file = sprintf("%s/Mitelman_summary_220204.RData", RESULT_PATH))      
########## summary ##########
#ssss <- rowSums(summary$result[,c(-1:-4)] != 0)
#maxi <- max(ssss)
#mini <- min(ssss)

#which(ssss == maxi)  
#which(ssss == mini)  

#sorted <- sort(ssss)
#which(ssss == sorted[-1])

########################        






