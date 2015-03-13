Make more changes, does it get automatically detected?

Make changes to this file and see it it reflects

# Reading Crude data: 
#        -Has missing values
#        -Given as time stamps (time format)
# Transit time between stops (Xi) is difference between arrival times of stops
# Missing value imputation using linear interpolation based on distance
# Predicting future Arrival times with reference to present stop (Yi-Yi') based on transit time of the past using Multiple linear regression - Automatic regression

# rm(list=ls()) #will remove ALL objects previously stored in workspace
# source("D:/transit/progfun.r")

load_library <- function()
{
library(leaps)
library(pls)
library(zoo)
library(grid)
library(lmtest)
library(MASS)
library(nortest)
library(mvtnorm)
library(multcomp)
library(RColorBrewer)
library(latticeExtra)
library(reshape)
library(plyr)
library(colorspace)
library(HH)
library(lmtest)
return()
}

Init_Global_variables <-function()
{
dates <<-vector()               # Global variables
Mon <<-vector()
Tue <<-vector()
Wed <<-vector()
Thu <<-vector()
Fri <<-vector()
Sat <<-vector()

}

# dprep<-read.table("crude45.data",header=T,sep="\t")

# Generating Date matrix & Weekday variable 

dayvariable<-function(dprep)
{
dtuple<<-dprep
dates <- as.Date(paste(dtuple$Date,"-2012",sep=""),"%d-%b-%Y")

for (j in 1:length(dtuple$Date))  
{
if(weekdays(dates[j], abbreviate = TRUE) == "Mon") 
{ 
Mon<<-cbind(Mon,1)
Tue<<-cbind(Tue,0)
Wed<<-cbind(Wed,0)
Thu<<-cbind(Thu,0)
Fri<<-cbind(Fri,0)
Sat<<-cbind(Sat,0)
 }
else if(weekdays(dates[j], abbreviate = TRUE) == "Tue") 
{ 
Mon<<-cbind(Mon,0)
Tue<<-cbind(Tue,1)
Wed<<-cbind(Wed,0)
Thu<<-cbind(Thu,0)
Fri<<-cbind(Fri,0)
Sat<<-cbind(Sat,0)
 }
else if(weekdays(dates[j], abbreviate = TRUE) == "Wed") 
{ 
Mon<<-cbind(Mon,0)
Tue<<-cbind(Tue,0)
Wed<<-cbind(Wed,1)
Thu<<-cbind(Thu,0)
Fri<<-cbind(Fri,0)
Sat<<-cbind(Sat,0)
 }
else if(weekdays(dates[j], abbreviate = TRUE) == "Thu") 
{ 
Mon<<-cbind(Mon,0)
Tue<<-cbind(Tue,0)
Wed<<-cbind(Wed,0)
Thu<<-cbind(Thu,1)
Fri<<-cbind(Fri,0)
Sat<<-cbind(Sat,0)
 }
else if(weekdays(dates[j], abbreviate = TRUE) == "Fri") 
{ 
Mon<<-cbind(Mon,0)
Tue<<-cbind(Tue,0)
Wed<<-cbind(Wed,0)
Thu<<-cbind(Thu,0)
Fri<<-cbind(Fri,1)
Sat<<-cbind(Sat,0)
 }
else if(weekdays(dates[j], abbreviate = TRUE) == "Sat") 
{ 
Mon<<-cbind(Mon,0)
Tue<<-cbind(Tue,0)
Wed<<-cbind(Wed,0)
Thu<<-cbind(Thu,0)
Fri<<-cbind(Fri,0)
Sat<<-cbind(Sat,1)
}
}

Mon<<-as.vector(Mon)
Tue<<-as.vector(Tue)
Wed<<-as.vector(Wed)
Thu<<-as.vector(Thu)
Fri<<-as.vector(Fri)
Sat<<-as.vector(Sat)

# return(cbind(Mon,Tue,Wed,Thu,Fri,Sat))
}

# Converting Arrival time (Yi) to Seconds data
Time_conversion<-function(dtuple)
{
r<-c("7:00:00")                                                 # Reference at 7 am 
tr<-as.POSIXct(r,"%H:%M:%S",tz="")                              # standard time format 

yi<-list()
Yi<-list()

for (j in 2:length(dprep))       
{ 
     ki<-vector()
     gi<-vector()
     ti<-vector()
        gi <- as.character(dprep[,j])
         
     ti<-as.POSIXct(gi,"%H:%M:%S",tz="")
     for (i in 1:length(ti))
        {
        if(is.null(ti[i]))
           {
           }
        else
           {
            ki[i]<-as.numeric(difftime(ti[i],tr,unit=c("secs")))
           }         
        } 
   
      yi<-as.matrix(cbind(yi,ki,deparse.level=0))
      Yi<-as.matrix(cbind(Yi,ki,deparse.level=0))               # Back-up (Unprocessed) Yi

	  
}
Yi<-as.data.frame(Yi)
return(Yi)
}
# Interpolation function:

spl<-function (tprev,tnext,dprev,dnext,dcurr)
{
tcurr <- as.numeric((tnext-tprev)*(dcurr-dprev)/(dnext-dprev) + tprev)
return(tcurr)
}

# dist<-c(1.19,6.23,2.58,2.5,2.28,1.28,2.09,1.52,0.75,0.66,0.4,3.53)

init_distance_vector <-function()
{
dref <<-c(5,6.19,12.42,15,17.5,19.78,21.06,23.15,24.67,25.42,26.08,26.48,30.01)
}
#############################################################################################################################

# Imputation of missing values based on distance based linear interpolation 

# ispy<-yi

brute_imputation_y <-function(ispy)
{
for(i in 1:length(yi[,1]))
{ 
for(j in 1:length(yi[i,])) 
{
if(is.na(yi[i,j]))                                                                  # enters only if missing values
{
    if(j==1)                                                                        # only if first column                               
	{
	  if(!is.na(yi[i,j+1]) & !is.na(yi[i,j+2]))                                     # 2 & 3 column are not NULL
	            {
	   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+2], dprev=dref[j+1],
                         tnext=as.numeric(yi[i,j+2]), tprev<-as.numeric(yi[i,j+1])) 
                }
      else if(is.na(yi[i,j+1]) & !is.na(yi[i,j+2]) & !is.na(yi[i,j+3]))                  # 2 column is null 
                {
      	ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+3], dprev=dref[j+2],
                         tnext=as.numeric(yi[i,j+3]), tprev<-as.numeric(yi[i,j+2]))   
                }
	  else if(!is.na(yi[i,j+1]) & is.na(yi[i,j+2]) & !is.na(yi[i,j+3]))                  # 3 column is null
	            {
        ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+3], dprev=dref[j+1],
                         tnext=as.numeric(yi[i,j+3]), tprev<-as.numeric(yi[i,j+1]))    
                }
	  else if(!is.na(yi[i,j+4]) & !is.na(yi[i,j+3]) )                  # 3 & 4 column is not null
	            {
        ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+4], dprev=dref[j+3],
                         tnext=as.numeric(yi[i,j+4]), tprev<-as.numeric(yi[i,j+3]))    
                }
	   else if(!is.na(yi[i,j+4]) & !is.na(yi[i,j+2]) )                  # 2 & 4 column is not null
	            {
        ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+4], dprev=dref[j+2],
                         tnext=as.numeric(yi[i,j+4]), tprev<-as.numeric(yi[i,j+2]))    
                }
	   else if(!is.na(yi[i,j+5]) & !is.na(yi[i,j+3]) )                  # 3 & 5 column is not null
	            {
        ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+5], dprev=dref[j+3],
                         tnext=as.numeric(yi[i,j+5]), tprev<-as.numeric(yi[i,j+3]))    
                }
		else if(!is.na(yi[i,j+5]) & !is.na(yi[i,j+4]) )                  # 4 & 5 column is not null
	            {
        ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+5], dprev=dref[j+4],
                         tnext=as.numeric(yi[i,j+5]), tprev<-as.numeric(yi[i,j+4]))    
                }
				
				
	}
	else if(j==length(yi[i,]))                                                           # only if last column
	{
	  if(!is.na(yi[i,j-1]) & !is.na(yi[i,j-2]))                                     # 2nd & 3rd last column are not NULL
	            {
	   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j-1], dprev=dref[j-2],
                         tnext=as.numeric(yi[i,j-1]), tprev<-as.numeric(yi[i,j-2])) 
                }
      else if(is.na(yi[i,j-1]) & !is.na(yi[i,j-2]) & !is.na(yi[i,j-3]))                  # 2nd last column is NULL
                {
      	ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j-2], dprev=dref[j-3],
                         tnext=as.numeric(yi[i,j-2]), tprev<-as.numeric(yi[i,j-3]))   
                }
	  else if(!is.na(yi[i,j-1]) & is.na(yi[i,j-2]) & !is.na(yi[i,j-3]))                  # 3rd last column is NULL
	            {
        ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j-1], dprev=dref[j-3],
                         tnext=as.numeric(yi[i,j-1]), tprev<-as.numeric(yi[i,j-3]))    
                }
	  else if(!is.na(yi[i,j-3]) & is.na(yi[i,j-2]) & !is.na(yi[i,j-4]))                  # 3rd last column is NULL
	            {
        ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j-3], dprev=dref[j-4],
                         tnext=as.numeric(yi[i,j-3]), tprev<-as.numeric(yi[i,j-4]))    
                }
	  else if(!is.na(yi[i,j-5]) & !is.na(yi[i,j-4]))                  # 3rd last column is NULL
	            {
        ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j-4], dprev=dref[j-5],
                         tnext=as.numeric(yi[i,j-4]), tprev<-as.numeric(yi[i,j-5]))    
                }
	}
    else
    {
        if(!is.na(yi[i,j-1]) & !is.na(yi[i,j+1]))                                   # nearest neighbour are not NULL
                {
         ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+1], dprev=dref[j-1],
                         tnext=as.numeric(yi[i,j+1]), tprev<-as.numeric(yi[i,j-1]))
                }
        else if(!is.na(yi[i,j-1]) & is.na(yi[i,j+1]))                                    # right side is NULL
		{
            if(j==(length(yi[i,])-1))              	 						        # consider 12th column
		    {
			    if(!is.na(yi[i,j-2]))         
                {
                  ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j-1],dprev=dref[j-2],
                              tnext=as.numeric(yi[i,j-1]), tprev<-as.numeric(yi[i,j-2]))
				}
				else if(!is.na(yi[i,j-3]))         
                {
                  ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j-1],dprev=dref[j-3],
                              tnext=as.numeric(yi[i,j-1]), tprev<-as.numeric(yi[i,j-3]))
			    }
			}
		    else if(!j==(length(yi[i,])-1)) 
			{
			if(!is.na(yi[i,j+2]))
			    {
                ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+2],dprev=dref[j-1],
                              tnext=as.numeric(yi[i,j+2]), tprev<-as.numeric(yi[i,j-1]))
                }
			else if(!is.na(yi[i,j+3]))
			    {
                ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+3],dprev=dref[j-1],
                              tnext=as.numeric(yi[i,j+3]), tprev<-as.numeric(yi[i,j-1]))
                }
			else if(!is.na(yi[i,j+4]))
			    {
                ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+4],dprev=dref[j-1],
                              tnext=as.numeric(yi[i,j+4]), tprev<-as.numeric(yi[i,j-1]))
                }	
		    }
			else if(!j==2)
			{
			if(!is.na(yi[i,j-2]))         
                {
                ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j-1],dprev=dref[j-2],
                              tnext=as.numeric(yi[i,j-1]), tprev<-as.numeric(yi[i,j-2]))
			    }
		    else if(!j==3 & !is.na(yi[i,j-3]))         
                {
                ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j-1],dprev=dref[j-3],
                              tnext=as.numeric(yi[i,j-1]), tprev<-as.numeric(yi[i,j-3]))
			    }
		    else if(j==3 & !is.na(yi[i,j+2]))         
                {
                ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+2],dprev=dref[j-1],
                              tnext=as.numeric(yi[i,j+2]), tprev<-as.numeric(yi[i,j+1]))
			    }
		    else if(j==3 & !is.na(yi[i,j+3]))         
                {
                ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+3],dprev=dref[j-1],
                              tnext=as.numeric(yi[i,j+3]), tprev<-as.numeric(yi[i,j+1]))
			    }
			}
			else if(j==2 & !is.na(yi[i,j+2]))         
                {
                ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+2],dprev=dref[j-1],
                              tnext=as.numeric(yi[i,j+2]), tprev<-as.numeric(yi[i,j+1]))
			    }
			else if(j==2 & !is.na(yi[i,j+3]))         
                {
                ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+3],dprev=dref[j-1],
                              tnext=as.numeric(yi[i,j+3]), tprev<-as.numeric(yi[i,j+1]))
			    }
			
        }
                 		
        else
		if(is.na(yi[i,j-1]) & !is.na(yi[i,j+1]))                                    # left side is also NULL
		{
		    if(j==2)
			{
			    if(!is.na(yi[i,j+2]))
                {
                ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+2],dprev=dref[j+1],
                              tnext=as.numeric(yi[i,j+2]), tprev<-as.numeric(yi[i,j+1]))
			    }
			    else
                if(!is.na(yi[i,j+3]))
                {
                ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+3],dprev=dref[j+1],
                              tnext=as.numeric(yi[i,j+3]), tprev<-as.numeric(yi[i,j+1]))
			    }
            }
			else if(j==3)
			{
			    if(!is.na(yi[i,j+2]))
                {
                ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+2],dprev=dref[j+1],
                              tnext=as.numeric(yi[i,j+2]), tprev<-as.numeric(yi[i,j+1]))
			    }
			    else
                if(!is.na(yi[i,j+3]))
                {
                ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+3],dprev=dref[j+1],
                              tnext=as.numeric(yi[i,j+3]), tprev<-as.numeric(yi[i,j+1]))
			    }
            }
			else if(j==(length(yi[i,])-1))
			{
				if(!is.na(yi[i,j-2])) 
			    {
               ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+1],dprev=dref[j-2],
                              tnext=as.numeric(yi[i,j+1]), tprev<-as.numeric(yi[i,j-2]))
			    }
				else if(!is.na(yi[i,j-3])) 
			    {
               ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+1],dprev=dref[j-3],
                              tnext=as.numeric(yi[i,j+1]), tprev<-as.numeric(yi[i,j-3]))
			    }
	        }
			else if(j==(length(yi[i,])-2))
			{
			    if(!is.na(yi[i,j-2])) 
			    {
               ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+1],dprev=dref[j-2],
                              tnext=as.numeric(yi[i,j+1]), tprev<-as.numeric(yi[i,j-2]))
			    }
				else if(!is.na(yi[i,j-3])) 
			    {
               ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+1],dprev=dref[j-3],
                              tnext=as.numeric(yi[i,j+1]), tprev<-as.numeric(yi[i,j-3]))
			    }
				else if(!is.na(yi[i,j+2])) 
			    {
               ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+2],dprev=dref[j+1],
                              tnext=as.numeric(yi[i,j+2]), tprev<-as.numeric(yi[i,j+1]))
			    }
			}
            else
			{
			if(!is.na(yi[i,j-2])) 
			    {
               ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+1],dprev=dref[j-2],
                              tnext=as.numeric(yi[i,j+1]), tprev<-as.numeric(yi[i,j-2]))
			    }
            else if(!is.na(yi[i,j+2])) 
			    {
               ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+2],dprev=dref[j+1],
                              tnext=as.numeric(yi[i,j+2]), tprev<-as.numeric(yi[i,j+1]))
			    }
            
			}
			
	    }
        
		if(is.na(yi[i,j+1]) & is.na(yi[i,j-1]))                            #Right and left are NULL
		{
		    if(j==2 | j==3)                                                # Taking care of extreme stops (first few and last few) 
			{
			    if(!is.na(yi[i,j+2]) & !is.na(yi[i,j+3]))
                {
                   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+3],dprev=dref[j+2],
                              tnext=as.numeric(yi[i,j+3]), tprev<-as.numeric(yi[i,j+2]))
			    }
				else
				if(!is.na(yi[i,j+4]) & !is.na(yi[i,j+3]))
                {
                   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+4],dprev=dref[j+3],
                              tnext=as.numeric(yi[i,j+4]), tprev<-as.numeric(yi[i,j+3]))
			    }
		        else
				if(!is.na(yi[i,j+2]) & !is.na(yi[i,j+4]))
                {
                   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+4],dprev=dref[j+2],
                              tnext=as.numeric(yi[i,j+4]), tprev<-as.numeric(yi[i,j+2]))
			    }
				
		    }
			else
			if(j==(length(yi[i,])-1) | j==(length(yi[i,])-2))
			{
			    if(!is.na(yi[i,j-2]) & !is.na(yi[i,j-3]))
                {
                   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j-2],dprev=dref[j-3],
                              tnext=as.numeric(yi[i,j-2]), tprev<-as.numeric(yi[i,j-3]))
			    }
				else
				if(!is.na(yi[i,j-2]) & !is.na(yi[i,j-4]))
                {
                   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j-2],dprev=dref[j-4],
                              tnext=as.numeric(yi[i,j-2]), tprev<-as.numeric(yi[i,j-4]))
			    }
				else
				if(!is.na(yi[i,j-4]) & !is.na(yi[i,j-3]))
                {
                   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j-3],dprev=dref[j-4],
                              tnext=as.numeric(yi[i,j-3]), tprev<-as.numeric(yi[i,j-4]))
			    }
				
				
			}
            else
            {
			    if(!is.na(yi[i,j-2]) & !is.na(yi[i,j+2]))
                { 
                   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+2],dprev=dref[j-2],
                              tnext=as.numeric(yi[i,j+2]), tprev<-as.numeric(yi[i,j-2]))
			    }
				else
				if(!is.na(yi[i,j-3]) & !is.na(yi[i,j+3]))
                { 
                   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+3],dprev=dref[j-3],
                              tnext=as.numeric(yi[i,j+3]), tprev<-as.numeric(yi[i,j-3]))
			    }
				else
				if(!is.na(yi[i,j-3]) & !is.na(yi[i,j+2]))
                { 
                   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+2],dprev=dref[j-3],
                              tnext=as.numeric(yi[i,j+2]), tprev<-as.numeric(yi[i,j-3]))
			    }
				else
				if(!is.na(yi[i,j-2]) & !is.na(yi[i,j+3]))
                { 
                   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+3],dprev=dref[j-2],
                              tnext=as.numeric(yi[i,j+3]), tprev<-as.numeric(yi[i,j-2]))
			    }
				else
				if(!is.na(yi[i,j+2]) & !is.na(yi[i,j+3]))
                { 
                   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j+3],dprev=dref[j+2],
                              tnext=as.numeric(yi[i,j+3]), tprev<-as.numeric(yi[i,j+2]))
			    }
				else
				if(!is.na(yi[i,j-2]) & !is.na(yi[i,j-3]))
                { 
                   ispy[i,j]<-spl(dcurr=dref[j], dnext=dref[j-2],dprev=dref[j-3],
                              tnext=as.numeric(yi[i,j-2]), tprev<-as.numeric(yi[i,j-3]))
			    }
						    #        ispy[i,j]<-c("wat")                       # debugging keyword
				else
		        {
		        ispy[i,j]<-mean(ispy[i,])
                }	
            }		
        }			
    }
}
# print(ispy[,j])                                             # Column debugger
}
}
# print(ispy)


counter<-0

for (j in 2:length(ispy[1,]))
{ 
     for (i in 1:length(ispy[,j]))
        {
        if(as.numeric(ispy[i,j]) - as.numeric(ispy[i,j-1])<=0)
           {
		   ispy[i,j]<- as.numeric(ispy[i,j-1])+10
		   counter<-counter + 1 
           }
        } 
}

# print(ispy)
# print("Imputation black spots for Y:")
# print(counter)

# Creating data frame for numeric matrix
ispY<-ispy  
ispY <- as.data.frame(ispY)                                                 # Backup 

for (j in 1:length(ispY[1,]))
{
ispY[,j]<- round(as.numeric(ispY[,j]),0)
}

# print(ispY[,4] - ispY[,3])  

return(ispY)
}
                                        
#########################################################################################################################

distance_vector <-function(ispY)
{
 
dstops<-vector()
for(k in 1:(length(dref)-1))
{
dstops[k]<-(dref[k+1] - dref[k])
}

### Changing starting distance dref[1]
# dstops<-c(1.19,6.23,2.58,2.5,2.28,1.28,2.09,1.52,0.75,0.66,0.4,3.53)

dref[1]<<-round((as.numeric(dref[length(dref)])-as.numeric(dref[2]))/(mean(as.numeric(ispY[,length(ispY)]))-mean(as.numeric(ispY[,2])))*mean(as.numeric(ispY[,1])),2)

for(j in 2:length(dref))
{
dref[j]<<-(dref[j-1]+dstops[j-1])
}

}
# print(dref)

# Generating Xi from ispy

generate_x <-function(ispy)
{
xi<-list()
Xi<-list()
Ki<-vector()
Ki<-ispy[,1]
Xi<-cbind(Xi,Ki,deparse.level = 0)
xi<-cbind(xi,Ki,deparse.level = 0)

for (j in 2:length(ispy[1,]))
{ 
     Ki<-vector() 
    
     for (i in 1:length(ispy[,j]))
        {
        if(is.null(ispy[i,j]) | is.null(ispy[i,j-1]))
           {
           }
        else
           {
            Ki[i]<-(as.numeric(ispy[i,j]) - as.numeric(ispy[i,j-1]))
           }         
        } 
      xi<-cbind(Xi,Ki,deparse.level = 0)
      Xi<-cbind(Xi,Ki,deparse.level = 0)                            # Backup
}

Xi<-as.data.frame(Xi)

## If Imputation fails, fill floor and ceiling values in Xi 

counter<-0

for (j in 2:length(Xi[1,]))
{ 
     for (i in 1:length(Xi[,j]))
        {
        if(Xi[i,j]<20)
           {
		   Xi[i,j]<-mean(as.numeric(Xi[,j]))
		   counter<-counter + 1 
           }
        else if(Xi[i,j]>1200)
           {
            Xi[i,j]<-mean(as.numeric(Xi[,j]))
			 counter<-counter + 1 
           }         
        } 
      
}

# print(Xi)
# print("Imputation black spots:")
# print(counter)

for (j in 1:length(Xi[1,]))
{
Xi[,j]<- round(as.numeric(Xi[,j]),0)
}

# print(Xi[,4] + Xi[,3])                                             # Testing if Xi is converted to numeric

return(Xi)
}

# Predicting future Arrival times with reference to present stop (Yi-Yi') based on 
# transit time of the past using Multiple linear regression - Automatic regression

naming_dataframe_x <-function(mat)
{
colX<-vector()

for(e in 1:length(mat))
{
colX<-c(colX,paste("X",e,sep=""))
}
colnames(mat)<-colX

return(mat)
}

naming_dataframe_y <-function(mat)
{
colY<-vector()

for(e in 1:length(mat))
{
colY<-c(colY,paste("Y",e,sep=""))
}
colnames(mat)<-colY

return(mat)
}

nor_metric<-function(x){                                                      ############ #NN 
print("AD test")
print("CVM test")
print("Shapiro test")
print("Lillie test")
print("SF test")
print("BP test")
print("DW test")
print(ad.test(x$residuals/(summary(x)$sigma))$p.value)
print(cvm.test(x$residuals/(summary(x)$sigma))$p.value)
print(shapiro.test(x$residuals/(summary(x)$sigma))$p.value)
print(lillie.test(x$residuals/(summary(x)$sigma))$p.value)
print(sf.test(x$residuals/(summary(x)$sigma))$p.value)
print(as.numeric(bptest(x)$p.value))
print(as.numeric(dwtest(x)$p.value))
}

model_build <- function(mod_x,mod_y)
{

rsq<<-vector()
Fst<<-vector()
Fdf1<<-vector()
Fdf2<<-vector()
Fpv<<-vector()
coeff <<-list()
nor_all <<-matrix()
stepmat<<-list()
summarymat<<-list()
bpres<<-list()
cooksvec<<-vector()
lambda <<-vector()
lambtak <<-vector()
stepmodel<<-list()
bpresinit<<-list()
resmeaninit<<-vector()
resmean<<-vector()
cookedrows<<-vector()
cooksdist<<-list()
postcooks<<-vector()
cox_t<<-vector()
nor_all_init<<-matrix()
nor_ad<<-vector()
nor_shapiro<<-vector()
nor_cvm<<-vector()
nor_sf<<-vector()
nor_lillie<<-vector()
nor_ad_in <<-vector()
nor_shapiro_in<<-vector()
nor_cvm_in<<-vector()
nor_sf_in<<-vector()
nor_lillie_in<<-vector()
bpres_middle<<-vector()
dw_test<<-vector()
coxed_m<<-list()
dw_test_init<<-vector()
isZ_Cox<<-list()
p=1

 for (j in 1:((length(mod_y))-1))
{
      for (i in (j+1):(length(mod_y)))
      {       
              			
			list <- vector()
			    print(p)
			  for (k in 1:j)
              {
               xstr <- paste("mod_x[,", k,"]",sep="")
			   list <- c(list, xstr)
		      }
	 		 leftside <- paste("mod_y[,",i,"]-mod_y[,",j, "]~", sep="")
			 form <- as.formula(paste(leftside ,paste( "Mon + Tue + Wed + Thu + ", paste(list, collapse = "+"), sep="")))
    #         form <- as.formula(paste(leftside , paste(list, collapse = "+")))
			 coxed<-lm(form)
		    	boxcoxed<-boxcox(coxed,plotit=F,lambda = seq(-5,5, 1/5))
			 coxed_m[[p]]<<-coxed
			 bpresinit[p]<<-as.numeric(bptest(coxed)$p.value)
			 var_unequal<<-as.numeric(bptest(coxed)$p.value)
			 dw_test_init[p]<<-as.numeric(dwtest(coxed)$p.value)
			 err_corr<-as.numeric(dwtest(coxed)$p.value)
			 resmeaninit[p]<<-t.test(resid(coxed),mu=0,alternative ="two.sided")$p.value
			 nor_ad_in[p]<<-ad.test(coxed$residuals/(summary(coxed)$sigma))$p.value
              nor_shapiro_in[p]<<-shapiro.test(coxed$residuals/(summary(coxed)$sigma))$p.value
			  nor_cvm_in[p]<<-cvm.test(coxed$residuals/(summary(coxed)$sigma))$p.value
			  nor_lillie_in[p]<<-lillie.test(coxed$residuals/(summary(coxed)$sigma))$p.value
			  nor_sf_in[p]<<-sf.test(coxed$residuals/(summary(coxed)$sigma))$p.value
			  nor_ad_i<-ad.test(coxed$residuals/(summary(coxed)$sigma))$p.value
              nor_shapiro_i<-shapiro.test(coxed$residuals/(summary(coxed)$sigma))$p.value
			  nor_cvm_i<-cvm.test(coxed$residuals/(summary(coxed)$sigma))$p.value
			  nor_lillie_i<-lillie.test(coxed$residuals/(summary(coxed)$sigma))$p.value
			  nor_sf_i<-sf.test(coxed$residuals/(summary(coxed)$sigma))$p.value
	#	     mean_nonzero <-t.test(resid(coxed),mu=0,alternative ="two.sided")$p.value			  
			count_normal<-0
            if(nor_ad_i < 0.05)
            { count_normal = count_normal + 1 }
  			if(nor_shapiro_i < 0.05)
            { count_normal = count_normal + 1 }
			if(nor_cvm_i < 0.05)
            { count_normal = count_normal + 1 }
			if(nor_lillie_i < 0.05)
            { count_normal = count_normal + 1 }
			if(nor_sf_i < 0.05)
            { count_normal = count_normal + 1 }
			
			  

	#	    outcrude <<- file.path("D:","transit","Crude Out Plots Cox",paste("plot_" ,p,"pl_Y" ,i,"Y",j,".jpg", sep = ""))
   # 		jpeg(file=outcrude)
	#	      par(mfrow=c(2,2))
    #         plot(coxed,which=c(1,4,5))
#		     qqnorm(resid(coxed))
	#	    qqline(resid(coxed))
    #        dev.off()
			
			flagcox<<-1
			
		if(var_unequal <= 0.05 | count_normal >= 3 | err_corr <= 0.05)	
		{
			 flagcox<<-with(boxcoxed, x[which.max(y)])
        }	
		    print("Lambda:")
			 print(flagcox)
    #        coxpath <<- file.path("D:","transit","Cox Plots Cox",paste(p,"pl_Y" ,i,"Y",j, ".jpg", sep = ""))
    #            jpeg(file=coxpath)
    #            iotitle = paste("Cox plot", "Y",i,"-Y",j)
    #          boxcox(coxed,plotit=T,lambda = seq(-5,5, 1/5))
     #           boxcox(coxed,plotit=T)
	#			dev.off()

			 isMon<-Mon
			 isTue<-Tue
			 isWed<-Wed
			 isThu<-Thu
			 isFri<-Fri
			 isY<-mod_y
            isX<-mod_x
			isX <-naming_dataframe_x(mat=isX)
			isY <-naming_dataframe_y(mat=isY)
            isZ<-cbind(isX,isY)		
				attach(isZ)
            r_list <- vector()
				
			  for (k in 1:j)
              {
			   r_xstr <- paste("X", k,sep="")
			   r_list <- c(r_list, r_xstr)
			  		 
	#		      inputspath <<- file.path("D:","transit","IPlots Cox",paste(p,"pl_Y" ,i,"Y",j,"X",k, ".jpg", sep = ""))
    #            jpeg(file=inputspath)
    #            iotitle = paste("Outputs versus inputs", "Y",i,"-Y",j,": X",k)
    #            plot(isX[,k],isY[,i]-isY[,j], main = iotitle)
    #			abline(lm(isY[,i]-isY[,j]~isX[,k]), col="red") # regression line (y~x)
    #            lines(lowess(isY[,i]-isY[,j],isX[,k]), col="blue") # lowess line (x,y) 
    #            dev.off()		
		     	}  
				  
			if(flagcox <= -2)
             {
               leftpart <- paste("I((Y",i,"-Y",j, ")^(flagcox))~", sep="")
			   cox_t[p]<<-flagcox
			 } 
			  else if(flagcox > -2 && flagcox <= -0.5)
			 {
			 leftpart <- paste("I((Y",i,"-Y",j, ")^(-1))~", sep="")
			 cox_t[p]<<- -1
			 }
	#		 else if(flagcox > -1 && flagcox <= -0.5)
	#		 {
	#		# leftpart <- paste("I((Y",i,"-Y",j, ")^(-0.5))~", sep="")
	#		  leftpart <- paste("I((Y",i,"-Y",j, ")^(flagcox))~", sep="")
	#		# leftpart <- paste("I(log(Y",i,"-Y",j, "))~", sep="")
	#		# cox_t[p]<<- -0.5
	#		 cox_t[p]<<- flagcox
	#		  
	#		 }
			 else if(flagcox > -0.5 && flagcox <= 0)
			 {
			 leftpart <- paste("I(log(Y",i,"-Y",j, "))~", sep="")
			 cox_t[p]<<- "log_op"
			 }
			 else if(flagcox > 0 && flagcox <= 0.5)
			 {
			 leftpart <- paste("I((Y",i,"-Y",j, ")^(0.5))~", sep="")
			 cox_t[p]<<- 0.5
			 }
			 else if(flagcox > 0.5 && flagcox < 2)
			 {
			 leftpart <- paste("Y",i,"-Y",j, "~", sep="")
			 cox_t[p]<<- 1
			 }
			 else if(flagcox >= 2 )
			 {
			 leftpart <- paste("I((Y",i,"-Y",j, ")^(flagcox))~", sep="")
			 cox_t[p]<<-flagcox
			 }
			 else
			 {
			 print("Something is wrong in Y")
			 }
			 form <- as.formula(paste(leftpart ,paste( "isMon + isTue + isWed + isThu + ", paste(r_list, collapse = "+"), sep="")))
	#		 form <- as.formula(paste(leftpart , paste(r_list, collapse = "+")))
  #			 print(form)	
			  cooked <- lm(form)
             bpres_middle[p]<<-as.numeric(bptest(cooked)$p.value) 		  
#  print(summary(cooked))
			  cooksvec<-cooks.distance(cooked)
			  rowsremoved<-0
		for(g in 1:length(cooksvec))                       
		{
		if(cooksvec[g]> 4/length(cooksvec))
		{
		isZ<-isZ[-g,]
		isMon<-isMon[-g]
		isTue<-isTue[-g]
		isWed<-isWed[-g]
		isThu<-isThu[-g]
		isFri<-isFri[-g]
		rowsremoved<-rowsremoved + 1
		}
		}
		cooksdist[[p]]<<-cooksvec
		
		detach(isZ)
		attach(isZ)
        isZ_Cox[[p]]<<-cbind(isZ,isMon,isTue,isWed,isThu)		
		
		      model<-lm(form)	
              step <- stepAIC(model, direction="both",AICc=TRUE,trace=FALSE)
          #    stepmat[[p]] <- stepAIC(model, direction="both",AICc=TRUE,trace=FALSE)
              print(summary(step))
			  summarymat[[p]]<<-summary(step)
              rsq<<-c(rsq,summary(step)$r.squared)
              coeff[[p]]<<-step$coefficients
              nor_ad[p]<<-ad.test(step$residuals/(summary(step)$sigma))$p.value
              nor_shapiro[p]<<-shapiro.test(step$residuals/(summary(step)$sigma))$p.value
			  nor_cvm[p]<<-cvm.test(step$residuals/(summary(step)$sigma))$p.value
			  nor_lillie[p]<<-lillie.test(step$residuals/(summary(step)$sigma))$p.value
			  nor_sf[p]<<-sf.test(step$residuals/(summary(step)$sigma))$p.value
			  
			 # stepmodel[[p]]<<-lm(step$model)
			  stepmodel[[p]]<<-step
			  bpres[p]<<-as.numeric(bptest(step)$p.value)
			  dw_test[p]<<-as.numeric(dwtest(step)$p.value)
			  resmean[p]<<-t.test(resid(step),mu=0,alternative ="two.sided")$p.value
			  lambda[p]<<-flagcox
         if (is.null(summary(step)$fstatistic[1])) {
                Fst[p]<<-0
                   Fdf1[p]<<-0
                   Fdf2[p]<<-0
                   Fpv[p]<<-0.99999
              }
         else {
                   Fst[p]<<-summary(step)$fstatistic[1]
                   Fdf1[p]<<-summary(step)$fstatistic[2]
                   Fdf2[p]<<-summary(step)$fstatistic[3]
                   Fpv[p]<<-(1-pf(Fst[p],Fdf1[p],Fdf2[p]))
              }
      
   #      residualpath <<- file.path("D:","transit","Res Plots Cox",paste("plot_" ,p,"pl_Y" ,i,"Y",j,".jpg", sep = ""))
   ##      jpeg(file=residualpath)
#		 par(mfrow=c(2,2))
 #        resititle = paste("Residuals versus fitted", p)
  #       plot(step$fitted.values,step$residuals, main = resititle)
  #       dev.off()
		 
		outpath <<- file.path("D:","transit","Out Plots Cox",paste("plot_" ,p,"pl_Y" ,i,"Y",j, ".jpg", sep = ""))
    		jpeg(file=outpath)
		      par(mfrow=c(2,2))
             plot(step,which=c(1,4,5))
		 qqnorm(resid(step))
 #		  qqline(resid(step))
         dev.off()
        
 detach(isZ)
					             
	p=p+1 
	}

}
}

# pairs(~Xi[,1]+Xi[,2]+Xi[,3]+Xi[,4]+Xi[,5]+Xi[,6]+Xi[,7]+Xi[,8]+Xi[,9]+Xi[,10]+Xi[,11]+Xi[,12],
# main="Scatterplot Matrix")

#  Checking Eigen Values XTX for Input matrix. Insignificant value means collinearity ??

#Mi<-stdize(as.matrix(Xi[,-13]))
#e<-eigen(t(Mi) %*% Mi)
#signif(e$values,4)

# Large Kappa values means nonlinearities there

# kapmodT<-model.matrix(~Ti[,1]+Ti[,2]+Ti[,3]+Ti[,4]+Ti[,5]+Ti[,6]+Ti[,7]+Ti[,8]+Ti[,9]+Ti[,10]+Ti[,11]+Ti[,12])
# kappa(kapmodT)

# One or more Eigen values near 0.00 means colinearities
#eigen(cov(Mi))$values

# Smaller the determinant, more is the collinearity 
# det(cov(Mi))

# VIF for Xi
#vif(Xi[,-13])
#### No multi-collinearity found

# To measure lamba accurate of Box Cox
#with(ff, x[which.max(y)])
#ff<-boxcox(lm(ispY[,4]-ispY[,3]~Xi[,2]+Xi[,3]+Xi[,1]),plotit=T)

#############################      EXECUTION OF CODE   ##########################################################

load_library()
Init_Global_variables()

dprep <<-read.table("crude45.data",header=T,sep=",")
# dprep <<-read.csv("crude50.csv",header=T,sep=",")
dayvariable(dprep)

yi<-Time_conversion(dtuple=dprep)
print("Time converted:")
print(yi)
init_distance_vector()
print(dref)
ispY <-brute_imputation_y(ispy=yi)
print("Imputed Y:")
print(ispY)
distance_vector(ispY=ispY)
print(dref)
Xi<-generate_x(ispy=ispY)
print("Imputed X:")
print(Xi)
#Xi <-naming_dataframe_x(mat=Xi)
#print("Named matrices:")
print(Xi)
#ispY <-naming_dataframe_y(mat=ispY)
print(ispY)
x_in<-Xi
y_in<-ispY
model_build(mod_x=x_in,mod_y=y_in)
#print("Models:")
#print(stepmodel)

nor_all_init<<-round(cbind(nor_ad_in,nor_lillie_in,nor_shapiro_in,nor_sf_in,nor_cvm_in),4)
nor_all<<-round(cbind(nor_ad,nor_lillie,nor_shapiro,nor_sf,nor_cvm),4)

#write.csv(round(nor_all,4),"D:/transit/Results/Normality Cox_F.csv")
#write.csv(round(Fpv,4),"D:/transit/Results/F-test Cox_F.csv")
#write.csv(round(as.numeric(bpres),4),"D:/transit/Results/BP Cox_F.csv")
#write.csv(round(as.numeric(dw_test),4),"D:/transit/Results/DW Cox_F.csv")
#write.csv(round(rsq,4),"D:/transit/Results/R-sq Cox_F.csv")
