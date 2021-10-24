#Parameters
#TH:        Threshold base omega
#windows:   Windows w
#follow:    Following time steps (jumps deleted) tau
#J:         Jump threshold Theta
#slice:     Period Number N
#Qstart:    Start Date, e.g. "2018-01-01"
#Qend:      End Date, e.g. "2018-07-31"

#Install package
PackageInstall()

#Read GHI Data and transfer to Clearness Index Z
data=ReadClearnessIndex()
date=data$date
Z=data$Z
GHI=data$GHI

#Get every-day start and end points 
start_end=DayTime(date,Z)
starttime=start_end$model1.fitted.values
endtime=start_end$model2.fitted.values

#Get jump sign modelling
sign=JumpSign(Qstart,Qend,TH,windows,follow,J,slice,starttime,endtime,date,Z)
model_sign=sign$model_sign
sign$accurate_rate

#Modelling and calculate parameters
parameters=Modelling(Qstart,Qend,TH,windows,follow,J,slice,starttime,endtime,date,Z)

#Clustering
#part:  Clustering numbers obtained by gap statistics method
cluster=Cluster(part,parameters)

#Simulation of a specific day
#Q: the data of a specific day
#Obtain the order of the specific day
n=as.numeric(as.Date(Q)-as.Date("2018-01-01"))+1 
#Obtain start and end points on the specific day
startpoint=starttime[n]
endpoint=endtime[n]
#Obtain GHI, time and Clearness Index data on the specific day
dayZ=Z[(date(date)==Q)&(hour(date)*60+minute(date)>=startpoint)&(hour(date)*60+minute(date)<=endpoint)]
daytime=date[(date(date)==Q)&(hour(date)*60+minute(date)>=startpoint)&(hour(date)*60+minute(date)<=endpoint)]
dayGHI=GHI[(date(date)==Q)&(hour(date)*60+minute(date)>=startpoint)&(hour(date)*60+minute(date)<=endpoint)]
#Obtain the order of the specific day
n2=as.numeric(as.Date(Q)-as.Date(Qstart))+1 
#Obtain the regime switching process on the specific day
ClusterNumber=cluster$cll[((n2-1)*slice+1):(n2*slice)]
#Simulation results of solar irradaince on the specific day
simulation_para=SimulationParameter(ClusterNumber,cluster)
