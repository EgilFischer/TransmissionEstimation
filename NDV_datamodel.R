####################################################################
#       NDV data model                                             #
#       Casts a data set with measurements pre day into a data set #
#       that can be used for parameters estimation                 #
####################################################################

#fill in groups and types of birds when NA
fill.type.na.rename <- function(data){
  new.data <- data;
  records <- length(data[,1]);
  for(i in c(1:records))
  {
      types <- c(ifelse(!is.na(new.data[i,]$Group),new.data[i,]$Group, types[1]),#set group
                 ifelse(!is.na(new.data[i,]$Subgroup),new.data[i,]$Subgroup, types[2]),#set group
                 ifelse(!is.na(new.data[i,]$Role.in.the.study),new.data[i,]$Role.in.the.study, types[3]),#set group
                 ifelse(!is.na(new.data[i,]$bird.id),new.data[i,]$bird.id, types[4]))#set group
      new.data[i,c(1:4)]<- types
        
  }
  colnames(new.data)[c(1:3,5)]<- c("Vaccinated", "Group","ci","Sample")
  new.data$Vaccinated <- ifelse(new.data$Vaccinated =="Vaccinated", "Yes", "No")
  new.data$Challenge <- ifelse(new.data$ci =="seeder", "seeder", "contact")
  return(new.data[!is.na(data$Sample.type),])#remove empty records
}

#uses file with states(positive / negative) already determined manually from Ct value
reform.SIR.data <- function(data, 
                            sampleday.vars,
                            exclude.seeder = T
){
  reform.data(data,
              sampleday.vars = sampleday.vars,
              exclude.seeder = exclude.seeder,
              SIR.state = T)
}

#count infections based on each type of sample separately 
reform.data<- function(data, 
                       model = c("SIS","SIR","SI")[1],#transmission model default SIS model
                       inf.rule =1, #rule with number of consecutive positive samples to be determined infectious
                       rec.rule =1, #rule with number of consecutive negative samples to be determined recovered
                       sampleday.vars,
                       group.by = NULL, 
                       positive.if = c("max","min")[1],#max will use at least 1 positive sample, min all samples positive.
                       cut.off = 0,
                       exclude.seeder = T,
                       SIR.state = F#False means that SIR status needs to be determined
)
{
  #remove manual classification which was added for the publication
  data <-data[(data$Sample=="category"&!is.na(data$Sample)),];
  data[,sampleday.vars]<- sapply(data[,sampleday.vars], function(x){ifelse(x=="S",0,ifelse(x=="I", 1,ifelse(x =="R",2,ifelse(x == "NI", NA, NA))))});
  #
  if(!SIR.state){
    #start with a recoding data to binary yes/no
    binary.data <- data; binary.data[, sampleday.vars] <- lapply(binary.data[, sampleday.vars], 
                                                                 FUN = function(x){ifelse(x < cut.off,1,0)});
    #group positive samples based on positive.if statement
    aggregate.data <- aggregate.data.frame(binary.data[,sampleday.vars], 
                                           by = list(binary.data$Group, 
                                                     binary.data$bird.id, 
                                                     binary.data$Vaccinated,
                                                     binary.data$ci),
                                           FUN = positive.if)
    #recode based on transmission model
    for(ani in c(1:length(aggregate.data[,1])))
    {
      
      aggregate.data[ani,sampleday.vars]<-SIR.state(as.vector(aggregate.data[ani,sampleday.vars][1,]),
                                                    model = model,
                                                    inf.rule = inf.rule,
                                                    rec.rule = rec.rule )
      
    }
    
    #change the colnames
    colnames(aggregate.data)[1:4]<- c("Group","bird.id","Vaccinated","Challenge")
  }
  else{
    aggregate.data <- data #already as input presented.
  }
  #aggregate to number of infections prior to time step TODO is to put this in a separate function
  id<-NULL; inf <- NULL; rec<- NULL;n <- NULL; cases <- NULL;susceptibles <- NULL;
  for(day in c(1:length(sampleday.vars)))
  {
    
    #determine the number of infectious on this day
    inf <- rbind(inf, data.frame(dpch = sampleday.vars[day],
                                 ndpch = day,
                                 aggregate.data.frame(aggregate.data[,sampleday.vars[day]], 
                                                      by = list(aggregate.data$Group, aggregate.data$Vaccinated), 
                                                      FUN = function(x){sum(x==1,na.rm = TRUE)})
    ))
    #determine the number of recoverd on this day
    rec <- rbind(rec, data.frame(dpch = sampleday.vars[day],
                                 ndpch = day,
                                 aggregate.data.frame(aggregate.data[,sampleday.vars[day]], 
                                                      by = list(aggregate.data$Group, aggregate.data$Vaccinated), 
                                                      FUN = function(x){sum(x==2, na.rm = TRUE)})
    ))
    #determine the number of birds on this day
    n <- rbind(n, data.frame(dpch = sampleday.vars[day],
                             ndpch = day,
                             aggregate.data.frame(aggregate.data[,sampleday.vars[day]], 
                                                  by = list(aggregate.data$Group, aggregate.data$Vaccinated), 
                                                  FUN = function(x) sum( !is.na(x) )) 
                             
    ))
    
    #determine the number of susceptibles
    if(exclude.seeder){
      susceptibles <- rbind(susceptibles, data.frame(dpch = sampleday.vars[day],
                                                     ndpch = day,
                                                     aggregate.data.frame(aggregate.data[aggregate.data$Challenge == "contact",sampleday.vars[day]], 
                                                                          by = list(aggregate.data$Group[aggregate.data$Challenge == "contact"], 
                                                                                    aggregate.data$Vaccinated[aggregate.data$Challenge == "contact"]), 
                                                                          FUN = function(x) {sum(x==0, na.rm = TRUE)}))) 
      
    }else
    {
      susceptibles <- rbind(susceptibles, data.frame(dpch = sampleday.vars[day],
                                                     ndpch = day,
                                                     aggregate.data.frame(aggregate.data[,sampleday.vars[day]], 
                                                                          by = list(aggregate.data$Group, aggregate.data$Vaccinated), 
                                                                          FUN = function(x) {sum(x==0, na.rm = TRUE)})))
    }
    
    #determine the number of new cases at this day
    if(day < length(sampleday.vars)){
      if(exclude.seeder){
        cases <- rbind(cases, data.frame(dpch = sampleday.vars[day],
                                         ndpch = day,
                                         cases =  aggregate.data.frame((
                                           aggregate.data[aggregate.data$Challenge == "contact",sampleday.vars[day+1]]-aggregate.data[aggregate.data$Challenge == "contact",sampleday.vars[day]])>0, 
                                           by = list(aggregate.data$Group[aggregate.data$Challenge == "contact"], aggregate.data$Vaccinated[aggregate.data$Challenge == "contact"]), 
                                           FUN = sum, na.rm = TRUE)   )      )
      }
      else
      {
        cases <- rbind(cases, data.frame(dpch = sampleday.vars[day],
                                         ndpch = day,
                                         cases =  aggregate.data.frame((
                                           aggregate.data[,sampleday.vars[day+1]]-aggregate.data[,sampleday.vars[day]])>0, 
                                           by = list(aggregate.data$Group, aggregate.data$Vaccinated), 
                                           FUN = sum, na.rm = TRUE)))
      }
    }else {
      
      #determine cases and potential exclude the seeders
      if(exclude.seeder){
        cases <- rbind(cases, 
                       data.frame(dpch = sampleday.vars[day],
                                  ndpch = day,
                                  cases = aggregate.data.frame(aggregate.data[aggregate.data$Challenge == "contact",sampleday.vars[day]], 
                                                               by = list(aggregate.data$Group[aggregate.data$Challenge == "contact"],
                                                                         aggregate.data$Vaccinated[aggregate.data$Challenge == "contact"]), 
                                                               FUN = function(x){0}))
        )}
      else{
        cases <- rbind(cases, 
                       data.frame(dpch = sampleday.vars[day],
                                  ndpch = day,
                                  cases = aggregate.data.frame(aggregate.data[,sampleday.vars[day]], 
                                                               by = list(aggregate.data$Group, aggregate.data$Vaccinated), 
                                                               FUN = function(x){0}))
        )}
    }
  }
  
  #combine infectious, population and cases
  output = cbind(inf,
                 S = susceptibles$x, 
                 C = cases$cases.x, 
                 N = n$x)
  
  colnames(output)[c(3:8)]<- c("Group","Vaccinated","I","S","C","N")
  
  return(list(output,aggregate.data))
  
  
}

###Determine state based on model ####
SIR.state<- function(in.data,#vector of consecutive samples
                     model,#type of transmission model to use
                     inf.rule =1,
                     rec.rule =1,
                     only.last = T #if only last cell is positive remove
){
  out.data <- in.data;
  #determine potential infection and recovery events
  change <- c(0,in.data[2:length(in.data)]-in.data[1:(length(in.data)-1)])
  #then based on the model recode it to being susceptible (0), infectious (1) or recovered (2)
  
  if(model == "SI")
  {
    #first positive samples and loop such that infection rule is fullfilled
    found = F; index = 1;index.first <- length(change)+1
    while(index <= length(change) & found == F & max(change==1,na.rm = T) == 1 )
    {
      index.first <- c(1:length(out.data))[change==1][index]
      if(!is.na(index.first)){
        if(min(in.data[index.first:min(index.first+inf.rule-1, length(in.data))])==1){found = T}
        else
        {
          #false positive removal
          out.data[index.first]<- 0;
          #add 
          index = index + 1;}
      }else{index = index + 1;}
    }
    #get the NA's
    index.nas <- c(1:length(out.data))[is.na(out.data)]
    #set all to one
    if(!is.na(index.first)){if(index.first <= length(out.data)){out.data[index.first:length(out.data)] <- 1}}
    out.data[index.nas] <- NA
  }
  
  if(model == "SIR")
  {
    #first positive samples and loop such that infection rule is fullfilled
    found = F; index = 1;index.first <- length(change)+1
    while(index <= length(change) & found == F & max(change==1,na.rm = T) == 1)
    {
      index.first <- c(1:length(out.data))[change==1][index]
      if(!is.na(index.first)){
        if(min(in.data[index.first:min(index.first+inf.rule-1, length(in.data))])==1){
          found = T
          #remove changes before the first index 
          change[1:(index.first-1)] <- 0
        }
        else{
          #false positive removal
          out.data[index.first]<- 0;
          #add 
          index = index + 1;}
      }else{index = index + 1;}
    }
    
    #first negative samples and loop such that infection rule is fullfilled
    found = F; index = 1;index.last <- length(change)+1
    while(index <= length(change) & found == F)
    {
      index.last <- c(1:length(out.data))[change== -1 & !is.na(change)][index]
      if(!is.na(index.last))
      {
        if(max(in.data[index.last:min(index.last+rec.rule-1, length(in.data))])==0)
        {
          found = T
        }
        else{
          #false negative removal
          out.data[index.last]<- 1;
          #add 
          index = index + 1;}
      } else{index = index + 1}
    }
    #get the NA's
    index.nas <- c(1:length(out.data))[is.na(out.data)]
    if(is.na(index.first)|index.first> length(out.data))
    {
      index.first <- NA #no infection
      index.last<- NA #no recovery without infection
    }
    
    #set all to one
    if(!is.na(index.first)){
      out.data[index.first:length(out.data)] <- 1
      if(!is.na(index.last)){ out.data[(index.last):length(out.data)] <- 2}}
    out.data[index.nas] <- NA
  }
  
  
  #if(model == "SIS") {stop("not yet implemented")}
  
  if(only.last){
    if(sum(as.numeric(out.data[1,]),na.rm = T)==1){
      if(tail(out.data[!is.na(out.data)], n = 1) == 1)
      {
        out.data[!is.na(out.data)]<- 0}#set all to be 0
    }
    
  }
  
  return(out.data)
}

