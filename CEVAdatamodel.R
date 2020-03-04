#####################################
# data model for CEVA data          #
#####################################

library(openxlsx)
library(ggplot2)
data.dir <- "C:/Surfdrive/Projecten/CEVA/"
#read the data of the vaccination group
cevadata1A <- read.xlsx(paste0(data.dir,
                               "P062VectormuneNDtransmissionstudysheddingdata_setA_readible.xlsx"),
                        sheet = "Group1A")
cevadata2A <- read.xlsx(paste0(data.dir, "P062VectormuneNDtransmissionstudysheddingdata_setA_readible.xlsx"),
                        sheet = "Group2A")
cevadata1B <- read.xlsx(paste0(data.dir,
                               "P062VectormuneNDtransmissionstudysheddingdata_setB_readible.xlsx"),
                        sheet = "Group1B")
cevadata2B <- read.xlsx(paste0(data.dir, "P062VectormuneNDtransmissionstudysheddingdata_setB_readible.xlsx"),
                        sheet = "Group2B")

cevadata <- rbind(cevadata1A,cevadata1B,cevadata2A,cevadata2B)

#rename to standard
cevadata$ci <- sapply(cevadata$Challenge, FUN = function(x){ifelse(grepl("contact", x),"contact","seeder")})


#reform data for glm
#count infections based on each type of sample separately
reform.data<- function(data, group.by = NULL, sampleday.vars,
                       positive.if = c("max","min")[1],#max will use at least 1 positive sample, min all samples positive.
                       cut.off = 0){
  #start with a recoding data to binary yes/no
  binary.data <- data; binary.data[, sampleday.vars] <- lapply(binary.data[, sampleday.vars], FUN = function(x){ifelse(x < cut.off,1,0)});
  #group positive samples based on positive.if statement
  aggregate.data <- aggregate.data.frame(binary.data[,sampleday.vars], 
                                         by = list(binary.data$Group, binary.data$bird.id, binary.data$Vaccinated), FUN = positive.if)
  #
  colnames(aggregate.data)[1:3]<- c("Group","bird.id","Vaccinated")
  #aggregate to number of infections prior to time step
  id<-NULL; inf <- NULL; n <- NULL; cases <- NULL;
  for(day in c(1:length(sampleday.vars)))
  {
    #determine the number of infectious at this day
    inf <- rbind(inf, data.frame(dpch = sampleday.vars[day],
                                 ndpch = day,
                                 aggregate.data.frame(aggregate.data[,sampleday.vars[day]], 
                                                          by = list(aggregate.data$Group, aggregate.data$Vaccinated), 
                                                          FUN = sum,na.rm = TRUE)
                                 ))
    #determine the number of birds at this day
    n <- rbind(n, data.frame(dpch = sampleday.vars[day],
                                 ndpch = day,
                                 aggregate.data.frame(aggregate.data[,sampleday.vars[day]], 
                                                           by = list(aggregate.data$Group, aggregate.data$Vaccinated), 
                                                           FUN = function(x) sum( !is.na(x) )) 
                                                      
                             ))
    #determine the number of new cases at this day
     if(day < length(sampleday.vars)){
     cases <- rbind(cases, data.frame(dpch = sampleday.vars[day],
                                ndpch = day,
                                cases =  aggregate.data.frame((
                                  aggregate.data[,sampleday.vars[day+1]]-aggregate.data[,sampleday.vars[day]])>0, 
                                                            by = list(aggregate.data$Group, aggregate.data$Vaccinated), 
                                                            FUN = sum, na.rm = TRUE)
                                  )
                    )
     }else {
       cases <- rbind(cases, 
                      data.frame(dpch = sampleday.vars[day],
                                            ndpch = day,
                                            cases = aggregate.data.frame(aggregate.data[,sampleday.vars[day]], 
                                                                      by = list(aggregate.data$Group, aggregate.data$Vaccinated), 
                                                                      FUN = function(x){0}))
       )
       }
  }
  susceptibles = n$x - inf$x
  #combine infectious, population and cases
  output = cbind(inf,
                 S = susceptibles, 
                 C = cases$cases.x, 
                 N = n$x
                )
  colnames(output)[c(3:8)]<- c("Group","Vaccinated","I","S","C","N")
  
  return(list(output,aggregate.data))
  
  
}

out <- reform.data(cevadata,sampleday.vars = paste0(c(1:14),".dpch"),cut.off = 36)

#infectious period
out.indiv<- out[[2]]
out.indiv$infper <- apply(out[[2]][,c(4:17)], 1,sum, na.rm = T)
#visualize
ggplot(data = out.indiv) +geom_histogram(aes(x = infper, fill = Vaccinated), position = "dodge",binwidth = 1)
#descriptive statistics
mean(out.indiv$infper)
aggregate(out.indiv$infper, list(out.indiv$Vaccinated),mean)
sd(out.indiv$infper)
aggregate(out.indiv$infper, list(out.indiv$Vaccinated),sd)
#test for equal variances
fligner.test(out.indiv$infper, out.indiv$Vaccinated) #conclude unequal variances
#let's compare the mean infectious period
t.test(out.indiv$infper~out.indiv$Vaccinated,var.equal = FALSE) #conclude no difference in mean
#let's compare the sd infectious period
var.test(out.indiv$infper~out.indiv$Vaccinated)

#transmission rate
out.nona <- out[[1]][(!is.na(out[[1]]$S)&!is.na(out[[1]]$C) &out[[1]]$I>0 &out[[1]]$S>0),]

#test glm
fit <- glm(cbind(out.nona$C, out.nona$S-out.nona$C) ~ out.nona$Vaccinated,offset = log(out.nona$I/out.nona$N), family = binomial(link = "cloglog"), data = out.nona, na.action = na.omit)
drop1(fit)
fit.empty <- glm(cbind(out.nona$C, out.nona$S-out.nona$C) ~1,offset = log(out.nona$I/out.nona$N), family = binomial(link = "cloglog"), data = out.nona, na.action = na.omit)
add1(fit.empty, ~out.nona$Vaccinated)
profile.confint <- confint(fit)
#beta non-vaccinated
exp(fit$coefficients[[1]])
exp(profile.confint[1,])
#beta vaccinated
exp(fit$coefficients[[1]]+fit$coefficients[[2]])
exp(profile.confint[2,]+fit$coefficients[[1]])

#calculate R0
alpha = 0.05 #conf.level
#beta
mean.logbeta <- summary(fit)$coefficients[,1] 
se.logbeta<- summary(fit)$coefficients[,2]
#infectious period
mean.inf.per<- aggregate(out.indiv$infper, list(out.indiv$Vaccinated),mean)
sem.inf.per <- aggregate(out.indiv$infper, list(out.indiv$Vaccinated),sd)$x/sqrt(c(sum(out.indiv$Vaccinated =="No"), sum(out.indiv$Vaccinated =="Yes")))
mean.log.inf.per <- aggregate(log(sapply(out.indiv$infper,function(x){max(x,.5)})), list(out.indiv$Vaccinated),mean)
sem.log.inf.per<- aggregate(log(sapply(out.indiv$infper,function(x){max(x,.5)})), list(out.indiv$Vaccinated),sd)$x/sqrt(c(sum(out.indiv$Vaccinated =="No"), sum(out.indiv$Vaccinated =="Yes")))
#R0
mean.R0 <- data.frame(Vaccinated = c("No","Yes"),
                      beta = exp(c(mean.logbeta[1],sum(mean.logbeta))),
                      llbeta = exp(profile.confint[,1]+c(mean.logbeta[1],sum(mean.logbeta))),
                      ulbeta = exp(profile.confint[,2]+c(mean.logbeta[1],sum(mean.logbeta))),
                      infper = mean.inf.per$x,
                      llinfper =  mean.inf.per$x+qnorm(alpha/2,lower.tail = T)*sem.inf.per,
                      ulinfper =  mean.inf.per$x+qnorm(alpha/2,lower.tail = F)*sem.inf.per,
                      logR0 = c(mean.logbeta[1]  ,sum(mean.logbeta))+mean.log.inf.per$x,
                      lllogR0 = (c(mean.logbeta[1]  ,sum(mean.logbeta))+mean.log.inf.per$x) + qnorm(alpha/2,lower.tail = T)*(se.logbeta + sem.log.inf.per),
                      ullogR0 = (c(mean.logbeta[1]  ,sum(mean.logbeta))+mean.log.inf.per$x) + qnorm(alpha/2,lower.tail = F)*(se.logbeta + sem.log.inf.per)
                      )

mean.R0$R <- exp(mean.R0$logR0)
mean.R0$llR <- exp(mean.R0$lllogR0)
mean.R0$ulR <- exp(mean.R0$ullogR0)
mean.R0
