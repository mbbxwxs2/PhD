PackageInstall=function()
{
  library("magrittr")
  library("lubridate") #Using ms() to convert to minutes
  library("MASS")
  library("freqparcoord")
  library("moments")
  library("MLmetrics")
  library("ggplot2")
  library("RColorBrewer")
  library("hydroGOF")
  library("data.table")
  library("R.utils")
  library("fitdistrplus")
}

ReadClearnessIndex=function()
{
  #Read GHI & Date from 2018-01-01 to 2018-07-31 
  setwd("~/Desktop/Solar/Data")
  data=read.csv("solardata.csv",header=F)
  date=as.POSIXct(paste(data$V1,data$V2), format="%d/%m/%Y %H:%M")
  GHI=data$V4
  date=as.numeric(date)
  t=2
  cond=FALSE
  m=0
  while (cond==FALSE) 
  {
    if (day(as.POSIXct(date[t],origin = "1970-01-01"))==day(as.POSIXct(date[t-1],origin = "1970-01-01")))
      if (date[t]-date[t-1]>60)
      {
        aver=(GHI[t]+GHI[t-1])/2
        GHI=insert(GHI,ats=t,values=aver)
        aver=(Temp[t]+Temp[t-1])/2
        Temp=insert(Temp,ats=t,values=aver)
        date=insert(date,ats=t,values=date[t-1]+60)
      }
    #Go to next day
    t=t+1
    if (t==length(GHI)) cond=TRUE
  }
  date=as.POSIXct(date,origin = "1970-01-01")
  #Calculate Clearness Index
  #Solar time
  #LONG: longitude (57.5522)
  long=57.4684*pi/180
  lst=15*pi/180*round(long/(15*pi/180))
  n=as.numeric(date(date)-as.Date("2018-01-01"))+1 
  B=2*pi*(n-1)/365  
  E=229.2*(0.000075+0.001868*cos(B)-0.032077*sin(B)-0.014615*cos(2*B)-0.04089*sin(2*B))
  localtime=hour(date)*60+minute(date)
  solartime=localtime+4*180/pi*(lst-long)+E
  w=(solartime-720)/4*pi/180
  #Zenith Angle cos(theta(t))
  #lat:  latitude (-20.2230)
  phi=20.2230*pi/180
  sigma=-23.45*pi/180*sin((n+284)/365*2*pi)
  zenith=cos(phi)*cos(sigma)*cos(w)+sin(phi)*sin(sigma)
  #Clearness index Z(t)
  E0=1+0.033*cos(n/365*2*pi)
  G=1367*E0*zenith
  Z=GHI/G
  data=data.frame(GHI,date,Z,Temp)
  return(data)
}

#Function of Clearness Index (To calculate the start and end time)
ClearnessIndex=function(date,localtime,root)
{
  long=57.4684*pi/180
  lst=15*pi/180*round(long/(15*pi/180))
  n=as.numeric(date(date)-as.Date("2018-01-01"))+1 
  B=2*pi*(n-1)/365  
  E=229.2*(0.000075+0.001868*cos(B)-0.032077*sin(B)-0.014615*cos(2*B)-0.04089*sin(2*B))
  solartime=localtime+4*180/pi*(lst-long)+E
  w=(solartime-720)/4*pi/180
  phi=-20.2230*pi/180
  sigma=23.45*pi/180*sin((n+284)/365*2*pi)
  zenith=cos(phi)*cos(sigma)*cos(w)+sin(phi)*sin(sigma)
  E0=1+0.033*cos(n/365*2*pi)
  G=1367*E0*zenith
  return (G-root)
}

#Calculate Start and End times
DayTime=function(date,Z)
{
  Q=seq(as.Date("2018-01-01"), as.Date("2018-07-31"), "days")
  sunrise_time=rep(0,length(Q))
  sunset_time=rep(0,length(Q))
  for (i in 1:length(Q))
  {
    if (length(Z[(date(date))==Q[i]])>0)
    {
      maximum=optimize(ClearnessIndex, interval=c(0, 1440), date=Q[i],root=0,maximum=TRUE)$objective
      percentage=maximum*0.1
      sunrise_time[i]=uniroot(ClearnessIndex,c(100,600),date=Q[i],root=percentage)$root
      sunset_time[i]=uniroot(ClearnessIndex,c(700,1440),date=Q[i],root=percentage)$root
    }
  }
  ##Sunrise Time
  x=1:length(sunrise_time)
  x=x[sunrise_time!=0]
  y1=sunrise_time[sunrise_time!=0]
  model1=lm(y1 ~ poly(x,4))
  #Sunset
  x=1:length(sunset_time)
  x=x[sunset_time!=0]
  y2=sunset_time[sunset_time!=0]
  model2=lm(y2 ~ poly(x,4))
  output=data.frame(model1$fitted.values,model2$fitted.values)
  return(output)
}

#Jump Sign Modelling
JumpSign=function(start,end,TH,windows,follow,J,slice,starttime,endtime,date,Z)
{
  #Initialize
  posnum=negnum=pos=neg=0
  neg=0
  pos=0
  #Data set according to Estimated Start & end point
  QD=seq(as.Date(start), as.Date(end), "days")
  for (u in 1:length(QD))
  {
    Q=QD[u]
    n=as.numeric(date(QD[u])-as.Date("2018-01-01"))+1 
    if (n>91) n=n-16
    startpoint=starttime[n]
    endpoint=endtime[n]
    dayZ=Z[(date(date)==QD[u])&(hour(date)*60+minute(date)>=startpoint)&(hour(date)*60+minute(date)<=endpoint)]
    daytime=date[(date(date)==QD[u])&(hour(date)*60+minute(date)>=startpoint)&(hour(date)*60+minute(date)<=endpoint)]
    if (length(dayZ)==0) next 
    #Signal Process 
    ###Need: dayZ, daytime, TH, windows
    ###Get: Signal
    mZ=dayZ 
    mmdayZ=dayZ
    mtime=daytime
    mZZ=mZ
    l=length(mZZ)
    signal=rep(0,l)
    id=1:l
    idd=id
    run=0 #Stop times
    old=0
    while (run<=50) 
    {
      mZZ=mZ
      mdayZ=mmdayZ
      idd=id
      l=length(mZZ)
      i=windows+1
      for (tt in (windows+1):l)
      {
        move=TRUE #Operate or not
        windowfix=mdayZ[(idd[tt]-windows):(idd[tt]-1)][mdayZ[(idd[tt]-windows):(idd[tt]-1)]!=-100]#Fix Windows
        if (length(windowfix)>1) 
        {
          MA=mean(windowfix)
          V=sd(windowfix)
        }  else move=FALSE
        dif=abs(mZZ[tt]-MA)
        if ((dif>TH*V)&(move))
        {
          if (mZZ[tt]>MA) signal[idd[tt]]=1 else signal[idd[tt]]=-1
          mZ=mZ[-i]
          mtime=mtime[-i]
          id=id[-i]
          mmdayZ[idd[tt]]=-100
        } else i=i+1
      }
      if (old==sum(signal!=0)) break else old=sum(signal!=0)
      run=run+1
    }
    #Jump Process
    ###Need: dayZ, daytime, signal, J
    ###Get: jump
    jump=rep(0,(length(dayZ)-1)) #jump_k=Z_(k+1)-Z_k
    jtime=daytime[2:length(daytime)] #jtime_k is the time of jump_k
    distance=0
    for (i in 2:length(dayZ)) distance[i-1]=dayZ[i]-dayZ[i-1] #Distance=Z_(k+1)-Z_k
    #If D>J Then Jump filtered
    for (i in 2:length(signal)) if ((signal[i]!=0) &(abs(distance[i-1])>J)) jump[i-1]=distance[i-1]
    #Delate Follow
    ###Need: jump, follow
    ###Get: mjump (Pure jump to estimate jump size distribution)
    ll=length(jump)
    mjump=jump
    tt=1
    while (tt<=ll) 
    {
      if (mjump[tt]!=0)
      {
        for (j in 1:follow) if ((tt+j)<=ll) if (mjump[tt+j]*mjump[tt]<0) mjump[tt+j]=0
      }
      tt=tt+1
    }
    #Jump size distribution
    ###Need: mjump
    #pjump is the pure jump size
    #ptime is the jump time corresponding to the jump, ptime_k is the time of jump_k, which is between (Z_k, Z_(k+1))
    #for (i in 1:length(mjump)) mjump[i]=mjump[i]/dayZ[i]
    pjump=mjump
    ptime=1:length(pjump)
    i=1
    tt=length(pjump)
    while (tt>0) {
      if (pjump[i]==0)
      {
        pjump=pjump[-i]
        ptime=ptime[-i]
      }
      else
        i=i+1
      tt=tt-1
    }
    #Q-Q plot of positive jump
    pnum=length(pjump[pjump>0])
    nnum=length(pjump[pjump<0])
    if (pnum!=0) pos[(posnum+1):(posnum+pnum)]=dayZ[ptime[pjump>0]]
    if (nnum!=0) neg[(negnum+1):(negnum+nnum)]=dayZ[ptime[pjump<0]]
    posnum=posnum+pnum
    negnum=negnum+nnum
  }
  #Logistic
  total=c(pos,neg)
  sign=c(rep(1,length(pos)),rep(0,length(neg)))
  sign=sign[order(total)]
  total=sort(total)
  model_sign=glm(sign~total,family=binomial(link='logit'))
  model.prob=predict(model_sign,type="response")
  model.pred=rep(0,length(sign))
  model.pred[model.prob>.5]=1
  accurate_rate=mean(model.pred==sign)
  output=list(model_sign=model_sign,accurate_rate=accurate_rate)
  return(output)
}

##O-U bridge
fill_data=function(a,b,gap,theta,mu,sigma)
{
  #Wiener process
  #set.seed(1)
  Z=sigma*rnorm(n=gap-1,sd =1)
  Z=c(0, cumsum(Z))
  #Brownian bridge
  BrownianBridge=rep(0,gap)
  T=gap-1
  for (i in 1:gap) BrownianBridge[i]=(b-a)/T*(i-1)+Z[i]-Z[gap]/T*(i-1)
  #O-U bridge
  result=rep(0,gap-1)
  result[1]=a-mu
  for (i in 1:T) result[i+1]=result[i]-theta*result[i]+2*theta*(b*exp(theta*(i-1+T))-result[i]*exp(2*theta*(i-1)))/(exp(2*theta*T)-exp(2*theta*(i-1)))+BrownianBridge[i+1]-BrownianBridge[i]
  result=result+mu
  result[gap]=b
  return(result)
}

#Estimate Parameters in Modelling
Modelling=function(Q1,Q2,TH,windows,follow,J,slice,starttime,endtime,date,Z)
{
  QD=seq(as.Date(Q1), as.Date(Q2), "days")
  mean_dayZ=sd_dayZ=lambda=pospara1=pospara2=negpara1=negpara2=mu=sigma=theta=rep(0,length(QD)*slice)
  for (day in 1:length(QD))
  {
    n=as.numeric(date(QD[day])-as.Date("2018-01-01"))+1 
    if (n>91) n=n-16
    startpoint=starttime[n]
    endpoint=endtime[n]
    dayZ=Z[(date(date)==QD[day])&(hour(date)*60+minute(date)>=startpoint)&(hour(date)*60+minute(date)<=endpoint)]
    daytime=date[(date(date)==QD[day])&(hour(date)*60+minute(date)>=startpoint)&(hour(date)*60+minute(date)<=endpoint)]
    if (length(dayZ)==0) next 
    #Initialize
    l=length(date)
    oo=day #day's number 
    #Signal Process 
    ###Need: dayZ, daytime, TH, windows
    ###Get: Signal
    mZ=dayZ 
    mmdayZ=dayZ
    mtime=daytime
    mZZ=mZ
    l=length(mZZ)
    signal=rep(0,l)
    id=1:l
    idd=id
    run=0 #Stop times
    old=0
    while (run<=50) 
    {
      mZZ=mZ
      mdayZ=mmdayZ
      idd=id
      l=length(mZZ)
      i=windows+1
      for (tt in (windows+1):l)
      {
        move=TRUE #Operate or not
        windowfix=mdayZ[(idd[tt]-windows):(idd[tt]-1)][mdayZ[(idd[tt]-windows):(idd[tt]-1)]!=-100]#Fix Windows
        if (length(windowfix)>1) 
        {
          MA=mean(windowfix)
          V=sd(windowfix)
        }  else move=FALSE
        dif=abs(mZZ[tt]-MA)
        if ((dif>TH*V)&(move))
        {
          if (mZZ[tt]>MA) signal[idd[tt]]=1 else signal[idd[tt]]=-1
          mZ=mZ[-i]
          mtime=mtime[-i]
          id=id[-i]
          mmdayZ[idd[tt]]=-100
        } else i=i+1
      }
      if (old==sum(signal!=0)) break else old=sum(signal!=0)
      run=run+1
    }
    #Jump Process
    ###Need: dayZ, daytime, signal, J
    ###Get: jump
    jump=rep(0,(length(dayZ)-1)) #jump_k=Z_(k+1)-Z_k
    jtime=daytime[2:length(daytime)] #jtime_k is the time of jump_k
    distance=0
    for (i in 2:length(dayZ)) distance[i-1]=dayZ[i]-dayZ[i-1] #Distance=Z_(k+1)-Z_k
    #If D>J Then Jump filtered
    for (i in 2:length(signal)) if ((signal[i]!=0) &(abs(distance[i-1])>J)) jump[i-1]=distance[i-1]
    #Delate Follow
    ###Need: jump, follow
    ###Get: mjump (Pure jump to estimate jump size distribution)
    ll=length(jump)
    mjump=jump
    tt=1
    while (tt<=ll) 
    {
      if (mjump[tt]!=0)
      {
        for (j in 1:follow) if ((tt+j)<=ll) if (mjump[tt+j]*mjump[tt]<0) mjump[tt+j]=0
      }
      tt=tt+1
    }
    #Jump size distribution
    ###Need: mjump
    #pjump is the pure jump size
    #ptime is the jump time corresponding to the jump, ptime_k is the time of jump_k, which is between (Z_k, Z_(k+1))
    #for (i in 1:length(mjump)) mjump[i]=mjump[i]/dayZ[i]
    pjump=mjump
    ptime=1:length(pjump)
    i=1
    tt=length(pjump)
    while (tt>0) {
      if (pjump[i]==0)
      {
        pjump=pjump[-i]
        ptime=ptime[-i]
      }
      else
        i=i+1
      tt=tt-1
    }
    ##Clustering Criteria: Mu, s.d., Lambda
    slice_l=round(length(dayZ)/slice)
    for (q in 1:(slice-1))
    {
      mean_dayZ[q+(oo-1)*slice]=mean(dayZ[((q-1)*slice_l+1):(q*slice_l)])
      sd_dayZ[q+(oo-1)*slice]=var(dayZ[((q-1)*slice_l+1):(q*slice_l)])
      lambda[q+(oo-1)*slice]=sum(mjump[((q-1)*slice_l+1):(q*slice_l)]!=0)/slice_l
    }
    mean_dayZ[slice+(oo-1)*slice]=mean(dayZ[((slice-1)*slice_l+1):length(dayZ)])
    sd_dayZ[slice+(oo-1)*slice]=var(dayZ[((slice-1)*slice_l+1):length(dayZ)])
    lambda[slice+(oo-1)*slice]=sum(mjump[((slice-1)*slice_l+1):(length(mjump))]!=0)/(length(mjump)-(slice-1)*slice_l-1)
    #Jump size
    #pjump is the pure jump size
    #ptime is the jump time corresponding to the jump, ptime_k is the time of jump_k, which is between (Z_k, Z_(k+1))
    for (q in 1:slice)
    {
      if (q<slice)
      {
        pjump=mjump[((q-1)*slice_l+1):(q*slice_l)]
        ptime=((q-1)*slice_l+1):(q*slice_l)
      } else {
        pjump=mjump[((q-1)*slice_l+1):length(mjump)]
        ptime=((q-1)*slice_l+1):length(mjump)
      }
      i=1
      tt=length(pjump)
      while (tt>0) {
        if (pjump[i]==0)
        {
          pjump=pjump[-i]
          ptime=ptime[-i]
        }
        else
          i=i+1
        tt=tt-1
      }
      pos=pjump[pjump>0]
      if (length(pos)>1)
      {
        fit_params=fitdistr(pos,"lognormal")
        pospara=fit_params$estimate
        pospara1[q+(oo-1)*slice]=pospara[1]
        pospara2[q+(oo-1)*slice]=pospara[2]
      } else
      {
        pospara1[q+(oo-1)*slice]=-100
        pospara2[q+(oo-1)*slice]=-100
      }
      #Q-Q plot of negative jump
      neg=-pjump[pjump<0]
      if (length(neg)>1)
      {
        fit_params <- fitdistr(neg,"lognormal")
        negpara=fit_params$estimate
        negpara1[q+(oo-1)*slice]=negpara[1]
        negpara2[q+(oo-1)*slice]=negpara[2]
      } else 
      {
        negpara1[q+(oo-1)*slice]=-100
        negpara2[q+(oo-1)*slice]=-100
      } 
    }
    #After remove jumps
    slice_l=round(length(dayZ)/slice)
    #1-(n-1) slices
    for (q in 1:slice)
    {
      #Get diffusion component
      if (q<slice)
      {
        l=length(dayZ[((q-1)*slice_l+1):(q*slice_l)])
        diffusion=dayZ[((q-1)*slice_l+1):(q*slice_l)]
        mt=((q-1)*slice_l+1):(q*slice_l)
        smjump=mjump[((q-1)*slice_l+1):(q*slice_l)]    
      } else {
        l=length(dayZ[((slice-1)*slice_l+1):length(dayZ)])
        diffusion=dayZ[((slice-1)*slice_l+1):length(dayZ)]
        mt=((slice-1)*slice_l+1):length(dayZ)
        smjump=mjump[((slice-1)*slice_l+1):(length(dayZ)-1)]  
      }
      i=1
      j=1
      while (i<l)
      {
        if (smjump[j]!=0) 
        {
          diffusion[i+1]=-1
          mt[i+1]=-1
          #diftime=diftime[-(i+1)]
          l=l-1
        } else i=i+1
        j=j+1
      }
      len=length(diffusion)
      if (sum(mt==-1)!=0)
      {
        #Fill the unequal data
        fillstep=0
        fillstart=0
        for (i in 2:(len-1))
        {
          if (diffusion[i]==-1 && mt[i]==-1)
          {
            if (fillstep==0) fillstart=i
            fillstep=fillstep+1
          }
          if (fillstep!=0&& diffusion[i]!=-1 && mt[i]!=-1)    
          {
            fillvalue=(diffusion[i]-diffusion[fillstart-1])/(i-fillstart+1)
            for (j in fillstart:(i-1)) diffusion[j]=diffusion[j-1]+fillvalue
            fillstep=0
          }
        } 
      }
      #Mean Reversion Component 
      #Get mu_1,sigma_1,theta
      tmu=ttheta=tsigma=0
      s0=s1=s00=s01=0
      s0=sum(diffusion[1:(len-1)])/(len-1)
      s1=sum(diffusion[2:len])/(len-1)
      for (i in 2:len)
      {
        s00=s00+diffusion[i-1]^2
        s01=s01+diffusion[i-1]*diffusion[i]
      }
      s00=s00/(len-1)
      s01=s01/(len-1)
      tmu=(s1*s00-s0*s01)/(s0*s1-s0^2-s01+s00)
      jundge=FALSE
      if ((s0-tmu)/(s1-tmu)>0) ttheta=log((s0-tmu)/(s1-tmu)) else jundge=TRUE
      beta=(1-exp(-ttheta))/ttheta
      m=0
      for (i in 2:len)
      {
        mm=(1-exp(-ttheta))/ttheta
        m=tmu*ttheta*mm+diffusion[i-1]*(1-ttheta*mm)
        tsigma=tsigma+(diffusion[i]-m)^2
      }
      tsigma=tsigma/((len-1)*beta*(1-0.5*ttheta*beta))
      if ((ttheta>=-0.05)&&(!jundge)&&(tmu>=-0.1)&&(tmu<1.2))
      {
        theta[q+(oo-1)*slice]=ttheta
        mu[q+(oo-1)*slice]=tmu
        sigma[q+(oo-1)*slice]=sqrt(tsigma)
      } else {
        #Random Walk Model (Ordinary least squares procedure)
        #Get mu_2, sigma_2
        y0=x1=0
        for (i in 1:(length(diffusion)-1))
        {
          y0[i]=(diffusion[i+1]-diffusion[i])
          x1[i]=1
        }
        f=lm(y0~-1+x1)
        summary(f)
        # Calculate \mu, \theta, \sigma
        beta1=f$coefficients[1]
        tmu=beta1
        tsigma=sd(f$residuals)
        theta[q+(oo-1)*slice]=-1
        mu[q+(oo-1)*slice]=tmu
        sigma[q+(oo-1)*slice]=tsigma
      }
      if (sum(mt==-1)!=0)
      {
        ##Convergency (re-fill the data)
        repeat
        {
          update_theta=theta[q+(oo-1)*slice]
          update_mu=mu[q+(oo-1)*slice]
          update_sigma=sigma[q+(oo-1)*slice]
          #Fill the unequal data 
          fill_start=1
          for (i in 2:len)
            if (mt[i]!=-1)
            {
              if (i-fill_start>1) 
              {
                fill_value=fill_data(diffusion[fill_start],diffusion[i],i-fill_start+1,update_theta,update_mu,update_sigma)
                diffusion[fill_start:i]=fill_value
              }
              fill_start=i
            }
          #Mean Reversion Component 
          #Get mu_1,sigma_1,theta
          tmu=ttheta=tsigma=0
          s0=s1=s00=s01=0
          s0=sum(diffusion[1:(len-1)])/(len-1)
          s1=sum(diffusion[2:len])/(len-1)
          for (i in 2:len)
          {
            s00=s00+diffusion[i-1]^2
            s01=s01+diffusion[i-1]*diffusion[i]
          }
          s00=s00/(len-1)
          s01=s01/(len-1)
          tmu=(s1*s00-s0*s01)/(s0*s1-s0^2-s01+s00)
          jundge=FALSE
          if ((s0-tmu)/(s1-tmu)>0) ttheta=log((s0-tmu)/(s1-tmu)) else jundge=TRUE
          beta=(1-exp(-ttheta))/ttheta
          m=0
          for (i in 2:len)
          {
            mm=(1-exp(-ttheta))/ttheta
            m=tmu*ttheta*mm+diffusion[i-1]*(1-ttheta*mm)
            tsigma=tsigma+(diffusion[i]-m)^2
          }
          tsigma=tsigma/((len-1)*beta*(1-0.5*ttheta*beta))
          if ((ttheta>=-0.05)&&(!jundge)&&(tmu>=0)&&(tmu<1))
          {
            theta[q+(oo-1)*slice]=ttheta
            mu[q+(oo-1)*slice]=tmu
            sigma[q+(oo-1)*slice]=sqrt(tsigma)
          } else {
            #Random Walk Model (Ordinary least squares procedure)
            #Get mu_2, sigma_2
            y0=x1=0
            for (i in 1:(length(diffusion)-1))
            {
              y0[i]=(diffusion[i+1]-diffusion[i])
              x1[i]=1
            }
            f=lm(y0~-1+x1)
            summary(f)
            # Calculate \mu, \theta, \sigma
            beta1=f$coefficients[1]
            tmu=beta1
            tsigma=sd(f$residuals)
            theta[q+(oo-1)*slice]=-1
            mu[q+(oo-1)*slice]=tmu
            sigma[q+(oo-1)*slice]=tsigma
          }
          epsilon_theta=abs(update_theta-theta[q+(oo-1)*slice])
          epsilon_mu=abs(update_mu-mu[q+(oo-1)*slice])
          epsilon_sigma=abs(update_sigma-sigma[q+(oo-1)*slice])
          if (epsilon_mu<0.01|epsilon_theta<0.01|epsilon_sigma<0.01) break else
          {
            update_theta=theta[q+(oo-1)*slice]
            update_mu=mu[q+(oo-1)*slice]
            update_sigma=sigma[q+(oo-1)*slice]
          }
        }
      }
    }
  }
  output=list(mean_dayZ=mean_dayZ,sd_dayZ=sd_dayZ,lambda=lambda,pospara1=pospara1,pospara2=pospara2,negpara1=negpara1,negpara2=negpara2,theta=theta,mu=mu,sigma=sigma,mjump=mjump)
  return(output)
}

Cluster=function(part,parameters)
{
  #Initialize
  mean_dayZ=parameters$mean_dayZ[parameters$mean_dayZ!=0]
  sd_dayZ=parameters$sd_dayZ[parameters$mean_dayZ!=0]
  lambda=parameters$lambda[parameters$mean_dayZ!=0]
  pospara1=parameters$pospara1[parameters$mean_dayZ!=0]
  pospara2=parameters$pospara2[parameters$mean_dayZ!=0]
  negpara1=parameters$negpara1[parameters$mean_dayZ!=0]
  negpara2=parameters$negpara2[parameters$mean_dayZ!=0]
  theta=parameters$theta[parameters$mean_dayZ!=0]
  mu=parameters$mu[parameters$mean_dayZ!=0]
  sigma=parameters$sigma[parameters$mean_dayZ!=0]
  ##Propotion (s.d.)
  day_inf=data.frame(scale(mean_dayZ),scale(sd_dayZ),scale(lambda))
  cll=kmeans(day_inf,part)$cluster
  l=length(cll)
  for (i in 1:l) if (theta[i]==-1) cll[i]=cll[i]+part
  for (i in 1:(2*part-1))
    if (sum(cll==i)==0) {
      for (j in (i+1):(part*2)) cll[cll==j]=cll[cll==j]-1
    }
  #Clustering for mean-reverting process and SDE dX=mu*dt+sigma*dB
  group_lambda=group_mu=group_theta=group_sigma=group_pos1=group_pos2=group_neg1=group_neg2=rep(0,max(cll))
  for (i in 1:max(cll))
  {
    group_lambda[i]=median(parameters$lambda[cll==i])
    group_theta[i]=median(theta[cll==i])
    group_mu[i]=median(mu[cll==i])
    group_sigma[i]=median(sigma[cll==i])
    pos1=pospara1[cll==i]
    pos2=pospara2[cll==i]
    neg1=negpara1[cll==i]
    neg2=negpara2[cll==i]
    if (length(pos1[pos1!=-100])<1) group_pos1[i]=-100 else group_pos1[i]=mean(pos1[pos1!=-100])
    if (length(pos2[pos2!=-100])<1) group_pos2[i]=-100 else group_pos2[i]=mean(pos2[pos2!=-100])
    if (length(neg1[neg1!=-100])<1) group_neg1[i]=-100 else group_neg1[i]=mean(neg1[neg1!=-100])
    if (length(neg2[neg2!=-100])<1) group_neg2[i]=-100 else group_neg2[i]=mean(neg2[neg2!=-100])
  }
  output=list(cll=cll,group_lambda=group_lambda,group_mu=group_mu,group_theta=group_theta,group_sigma=group_sigma,group_pos1=group_pos1,group_pos2=group_pos2,group_neg1=group_neg1,group_neg2=group_neg2)
  return(output)
}

SimulationParameter=function(ClusterNumber,cluster)
{
  #Initialization
  label=ClusterNumber
  slambda=smu=ssigma=spospara1=spospara2=snegpara1=snegpara2=stheta=0
  for (i in 1:length(label))
  {
    slambda[i]=cluster$group_lambda[label[i]]
    smu[i]=cluster$group_mu[label[i]]
    ssigma[i]=cluster$group_sigma[label[i]]
    stheta[i]=cluster$group_theta[label[i]]
    spospara1[i]=cluster$group_pos1[label[i]]
    spospara2[i]=cluster$group_pos2[label[i]]
    snegpara1[i]=cluster$group_neg1[label[i]]
    snegpara2[i]=cluster$group_neg2[label[i]]
  }
  output=data.frame(slambda,smu,ssigma,spospara1,spospara2,snegpara1,snegpara2,stheta)
  return(output)
}

#Simulation
simulation=function(dayZ,slice,simulation_para,ClusterNumber)
{
  #Initialize
  slambda=simulation_para$slambda
  smu=simulation_para$smu
  ssigma=simulation_para$ssigma
  spospara1=simulation_para$spospara1
  spospara2=simulation_para$spospara2
  snegpara1=simulation_para$snegpara1
  snegpara2=simulation_para$snegpara2
  stheta=simulation_para$stheta
  ###
  sZ=0
  sZ[1]=dayZ[1]
  sjump=rep(0,length(dayZ))
  slice_l=round(length(dayZ)/slice)
  for (i in 1:slice)
  {
    if (i<slice) simnumber=slice_l else simnumber=length(dayZ)-slice_l*(slice-1)-1
    B=rnorm(length(dayZ))
    prob=rpois(simnumber,slambda[i])
    if(ClusterNumber[i]<=part)
    {
      for (j in 1:simnumber)
      {
        if (prob[j]==1)
        {
          #Only positive jump
          if ((spospara1[i]!=-100)&(spospara2[i]!=-100)&(snegpara1[i]==-100)&(snegpara2[i]==-100))
            sjump[j+(i-1)*slice_l]=rlnorm(1,spospara1[i],spospara2[i]) #postive jump
          #Only negative jump
          else if ((spospara1[i]==-100)&(spospara2[i]==-100)&(snegpara1[i]!=-100)&(snegpara2[i]!=-100))
            sjump[j+(i-1)*slice_l]=-rlnorm(1,snegpara1[i],snegpara2[i]) #negative jump}
          else  
            #Both jumps
            if ((spospara1[i]!=-100)&(spospara2[i]!=-100)&(snegpara1[i]!=-100)&(snegpara2[i]!=-100))
            {
              prob1=runif(1,0,1)
              if (prob1<predict(model_sign, list(total = sZ[j+(i-1)*slice_l]),type="response"))
                sjump[j+(i-1)*slice_l]=rlnorm(1,spospara1[i],spospara2[i]) #postive jump
              else 
                sjump[j+(i-1)*slice_l]=-rlnorm(1,snegpara1[i],snegpara2[i]) #negative jump}
            }          
        }
        #Cap at 1.5 and -0.02
        if (sZ[j+(i-1)*slice_l]+stheta[i]*(smu[i]-sZ[j+(i-1)*slice_l])+ssigma[i]*B[j+(i-1)*slice_l]+sjump[j+(i-1)*slice_l]>1) 
        {
          B[j+(i-1)*slice_l]=rnorm(1,0,1)
          if (sZ[j+(i-1)*slice_l]+stheta[i]*(smu[i]-sZ[j+(i-1)*slice_l])+ssigma[i]*B[j+(i-1)*slice_l]+sjump[j+(i-1)*slice_l]>1.5) return(-1)
        }
        if (sZ[j+(i-1)*slice_l]+stheta[i]*(smu[i]-sZ[j+(i-1)*slice_l])+ssigma[i]*B[j+(i-1)*slice_l]+sjump[j+(i-1)*slice_l]<0) 
        {
          B[j+(i-1)*slice_l]=rnorm(1,0,1)
          if (sZ[j+(i-1)*slice_l]+stheta[i]*(smu[i]-sZ[j+(i-1)*slice_l])+ssigma[i]*B[j+(i-1)*slice_l]+sjump[j+(i-1)*slice_l]< (-0.02)) return(-1)
        }
        sZ[j+1+(i-1)*slice_l]=sZ[j+(i-1)*slice_l]+stheta[i]*(smu[i]-sZ[j+(i-1)*slice_l])+ssigma[i]*B[j+(i-1)*slice_l]+sjump[j+(i-1)*slice_l]
      }
    } else {
      for (j in 1:simnumber)
      {
        if (prob[j]==1)
        {
          #Only positive jump
          if ((spospara1[i]!=-100)&(spospara2[i]!=-100)&(snegpara1[i]==-100)&(snegpara2[i]==-100))
            sjump[j+(i-1)*slice_l]=rlnorm(1,spospara1[i],spospara2[i]) #postive jump
          #Only negative jump
          else if ((spospara1[i]==-100)&(spospara2[i]==-100)&(snegpara1[i]!=-100)&(snegpara2[i]!=-100))
            sjump[j+(i-1)*slice_l]=-rlnorm(1,snegpara1[i],snegpara2[i]) #negative jump}
          else  
            #both jumps
            if ((spospara1[i]!=-100)&(spospara2[i]!=-100)&(snegpara1[i]!=-100)&(snegpara2[i]!=-100))
            {
              prob1=runif(1,0,1)
              if (prob1<predict(model_sign, list(total = sZ[j+(i-1)*slice_l]),type="response"))
                sjump[j+(i-1)*slice_l]=rlnorm(1,spospara1[i],spospara2[i]) #postive jump
              else 
                sjump[j+(i-1)*slice_l]=-rlnorm(1,snegpara1[i],snegpara2[i]) #negative jump}
            }          
        }          
        #Cap at 1.5 and -0.02
        if (sZ[j+(i-1)*slice_l]+smu[i]+ssigma[i]*B[j+(i-1)*slice_l]+sjump[j+(i-1)*slice_l]>1) 
        {
          B[j+(i-1)*slice_l]=rnorm(1,0,1)
          if (sZ[j+(i-1)*slice_l]+smu[i]+ssigma[i]*B[j+(i-1)*slice_l]+sjump[j+(i-1)*slice_l]>1.5) return(-1)
        }
        if (sZ[j+(i-1)*slice_l]+smu[i]+ssigma[i]*B[j+(i-1)*slice_l]+sjump[j+(i-1)*slice_l]<0) 
        {
          B[j+(i-1)*slice_l]=rnorm(1,0,1)
          if (sZ[j+(i-1)*slice_l]+smu[i]+ssigma[i]*B[j+(i-1)*slice_l]+sjump[j+(i-1)*slice_l]<(-0.02)) return(-1)
        }
        sZ[j+1+(i-1)*slice_l]=sZ[j+(i-1)*slice_l]+smu[i]+ssigma[i]*B[j+(i-1)*slice_l]+sjump[j+(i-1)*slice_l]
      }
    }
  }  
  output=data.frame(sZ,sjump)
  return(output)
}


InverseGHI=function(sZ,daytime)
{
  date=daytime
  #Solar time
  #LONG=57.5522
  long=57.4684*pi/180
  lst=15*pi/180*round(long/(15*pi/180))
  n=as.numeric(date(date)-as.Date("2018-01-01"))+1 
  #B=2*pi*(n-81)/365
  #E=9.87*sin(2*B)-7.53*cos(B)-1.5*sin(B)
  B=2*pi*(n-1)/365  
  E=229.2*(0.000075+0.001868*cos(B)-0.032077*sin(B)-0.014615*cos(2*B)-0.04089*sin(2*B))
  localtime=hour(date)*60+minute(date)
  solartime=localtime+4*180/pi*(lst-long)+E
  w=(solartime-720)/4*pi/180
  #Zenith Angle cos(theta(t))
  phi=20.2230*pi/180
  sigma=-23.45*pi/180*sin((n+284)/365*2*pi)
  zenith=cos(phi)*cos(sigma)*cos(w)+sin(phi)*sin(sigma)
  #Clearness index Z(t)
  E0=1+0.033*cos(n/365*2*pi)
  G=1367*E0*zenith
  output=sZ*G
  return(output)
}

#Search for the longest matching labels
Mapping=function(u,part,sim,cll) #u is the forecast number now
{
  prop=0
  nprop=0
  leng=0
  l=length(cll)
  for (i in 1:(l-1))
  {
    same=TRUE
    samp=0
    while (same==TRUE) 
    {
      samp=samp+1
      if ((l-i-samp+1<1)|(u-samp)<1) same=FALSE 
      if (same&&(sim[u-samp]!=cll[l-i+1-samp])) same=FALSE
    }
    if (samp>1)
    {
      nprop=nprop+1
      prop[nprop]=cll[l-i+1]
      leng[nprop]=samp-1  
    }
  }
  chose=prop[leng==max(leng)]
  option=rep(0,part)
  for (i in 1:length(chose)) option[chose[i]]=option[chose[i]]+1
  option=option/sum(option)
  return(option)
}  

GHIanalyse=function(Q,QD,run_time,Pplot,l,sim)
{
  checkn=which(QD==Q)
  ClusterNumber=sim[(l+1+(checkn-1)*slice):(l+checkn*slice)]  
  n=as.numeric(as.Date(Q)-as.Date("2018-01-01"))+1 
  if (n>91) n=n-16
  startpoint=starttime[n]
  endpoint=endtime[n]
  dayZ=Z[(date(date)==Q)&(hour(date)*60+minute(date)>=startpoint)&(hour(date)*60+minute(date)<=endpoint)]
  daytime=date[(date(date)==Q)&(hour(date)*60+minute(date)>=startpoint)&(hour(date)*60+minute(date)<=endpoint)]
  dayGHI=GHI[(date(date)==Q)&(hour(date)*60+minute(date)>=startpoint)&(hour(date)*60+minute(date)<=endpoint)]
  n2=as.numeric(as.Date(Q)-as.Date(Qstart))+1 
  datef=data.frame(x=daytime,y=dayGHI)
  gplot=ggplot(datef,aes(x=x,y=y)) +
    #geom_line(color=brewer.pal(9, "Greys")[7],size=1)+
    geom_line(color=brewer.pal(7, "Set1")[2],size=1)+
    theme(text=element_text(size=16),
          legend.text=element_text(size=16), 
          axis.text=element_text(size=18, colour="black"),
          #panel.background = element_rect(fill = 'white', colour = 'black'),
          axis.title.y = element_text(size=18, vjust = 2), 
          axis.title.x = element_text(size=18, vjust = 0),
          legend.key=element_rect(fill="white")) +
    scale_x_datetime(date_labels = "%H:%M",date_breaks = "2 hour")+ 
    xlab("Time") + ylab("GHI")
  ks_value=MaxAE_value=MAPE_value=RMSE_value=MAE_value=NRMSE=skew_value=kurtosis_value=sd_value=mean_value=rep(0,run_time)
  for (i in 1:run_time)
  {
    simulation_para=SimulationParameter(ClusterNumber,cluster)
    k=0
    repeat
    {
      k=k+1
      output=simulation(dayZ,slice,simulation_para,ClusterNumber)
      if (length(output)!=1) break
      if (k>=2000) return(-1)
    }
    sZ=output$sZ
    sGHI=InverseGHI(sZ,daytime)
    gplot=gplot+geom_line(y=sGHI,color=brewer.pal(7, "Set2")[2],size=1)
    #Tests
    ks_value[i]=ks.test(jitter(dayGHI),jitter(sGHI))$p.value
    MaxAE_value[i]=max(abs(sGHI-dayGHI))
    MAPE_value[i]=MAPE(sGHI,dayGHI)
    RMSE_value[i]=RMSE(sGHI,dayGHI)
    MAE_value[i]=MAE(sGHI,dayGHI)
    NRMSE[i]=nrmse(sGHI,dayGHI,norm="maxmin")
    #Property
    mean_value[i]=mean(sGHI)
    sd_value[i]=sd(sGHI)
    skew_value[i]=skewness(sGHI)
    kurtosis_value[i]=kurtosis(sGHI)
  }
  #gplot=gplot+geom_line(y=dayGHI,color=brewer.pal(9, "Greys")[7],size=1)
  # Add a ribbon with the confidence band
  result=data.frame(ks=mean(ks_value),rmse=mean(RMSE_value),nrmse=mean(NRMSE),mae=mean(MAE_value),maxae=mean(MaxAE_value),mape=mean(MAPE_value),mean=mean(mean_value),sd=mean(sd_value),skew=mean(skew_value),kurt=mean(kurtosis_value),histmean=mean(dayGHI),histsd=sd(dayGHI),histskew=skewness(dayGHI),histkurt=kurtosis(dayGHI))
  if (Pplot==TRUE) print(gplot)
  return(result)
}

CIanalyse=function(Q,QD,run_time,Pplot,l,sim)
{
  checkn=which(QD==Q)
  ClusterNumber=sim[(l+1+(checkn-1)*slice):(l+checkn*slice)]  
  n=as.numeric(as.Date(Q)-as.Date("2018-01-01"))+1 
  if (n>91) n=n-16
  startpoint=starttime[n]
  endpoint=endtime[n]
  dayZ=Z[(date(date)==Q)&(hour(date)*60+minute(date)>=startpoint)&(hour(date)*60+minute(date)<=endpoint)]
  daytime=date[(date(date)==Q)&(hour(date)*60+minute(date)>=startpoint)&(hour(date)*60+minute(date)<=endpoint)]
  n2=as.numeric(as.Date(Q)-as.Date(Qstart))+1 
  datef=data.frame(x=daytime,y=dayZ)
  gplot=ggplot(datef,aes(x=x,y=y)) +
    #geom_line(color=brewer.pal(9, "Greys")[7],size=1)+
    geom_line(color=brewer.pal(7, "Set1")[2],size=1)+
    theme(text=element_text(size=16),
          legend.text=element_text(size=16), 
          axis.text=element_text(size=18, colour="black"),
          #panel.background = element_rect(fill = 'white', colour = 'black'),
          axis.title.y = element_text(size=18, vjust = 2), 
          axis.title.x = element_text(size=18, vjust = 0),
          legend.key=element_rect(fill="white")) +
    scale_x_datetime(date_labels = "%H:%M",date_breaks = "2 hour")+ 
    xlab("Time") + ylab("Clearness Index")
  ks_value=MaxAE_value=MAPE_value=RMSE_value=MAE_value=NRMSE=skew_value=kurtosis_value=sd_value=mean_value=rep(0,run_time)
  for (i in 1:run_time)
  {
    simulation_para=SimulationParameter(ClusterNumber,cluster)
    k=0
    repeat
    {
      k=k+1
      output=simulation(dayZ,slice,simulation_para,ClusterNumber)
      if (length(output)!=1) break
      if (k>=2000) return(-1)
    }
    sZ=output$sZ
    #gplot=gplot+geom_line(y=sZ,color="red",size=0.5)
    gplot=gplot+geom_line(y=sZ,color=brewer.pal(7, "Set2")[2],size=1)
    #Tests
    ks_value[i]=ks.test(jitter(dayZ),jitter(sZ))$p.value
    MaxAE_value[i]=max(abs(sZ-dayZ))
    MAPE_value[i]=MAPE(sZ,dayZ)
    RMSE_value[i]=RMSE(sZ,dayZ)
    MAE_value[i]=MAE(sZ,dayZ)
    NRMSE[i]=nrmse(sZ,dayZ,norm="maxmin")
    #Property
    mean_value[i]=mean(sZ)
    sd_value[i]=sd(sZ)
    skew_value[i]=skewness(sZ)
    kurtosis_value[i]=kurtosis(sZ)
  }
  #gplot=gplot+geom_line(y=dayZ,color=brewer.pal(9, "Greys")[7],size=1)
  # Add a ribbon with the confidence band
  result=data.frame(ks=mean(ks_value),rmse=mean(RMSE_value),nrmse=mean(NRMSE),mae=mean(MAE_value),
                    maxae=mean(MaxAE_value),mape=mean(MAPE_value),mean=mean(mean_value),
                    sd=mean(sd_value),skew=mean(skew_value),kurt=mean(kurtosis_value),
                    histmean=mean(dayZ),histsd=sd(dayZ),histskew=skewness(dayZ),histkurt=kurtosis(dayZ))
  if (Pplot==TRUE) print(gplot)
  return(result)
}


