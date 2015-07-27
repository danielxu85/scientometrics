#R code for incoming citations

# We will fit a varying coefficient Spine Poisson regression with radial kernel.
# The belief we have is that the number of incoming citations have a similar structure when
# we parametrize them as a function of |t - i| where t is the outgoing citation year, and i is the incoming
# citation year. A Spline Poisson regression models counts of incoming citations as random fluctuations around a 
# smooth function of |t-i|. Varying coefficients with a radial kernel allows this smooth function to also change 
# smoothly in t. 


#load data, and libraries ####

load("../data/timeseries.Rdata")
# Nomenclature Definition: From A -> To B, which means that document A cites document B
library(splines)
library(reshape)
library(ggplot2)

# Functions ####

fit_pub_vcm <- function(dat,lam=1,d=3,xx=0:10)
{
  #Fits a vcm model to incoming citation data
  #lam is a parameter controlling degree of connection in t (higher lambda is less connection)
  #d is the degrees of freedom of the spline fit for |t-i| (higher d leads to more exact, but less reliable fit)
  #xx describes formatting of the output
  
  dat.melt = melt(dat)
  
  dat.melt$diff = dat.melt$From.Pub.Year - dat.melt$To.Pub.Year
  
  temp.idx = which(dat.melt$diff < 0)
  
  dat.melt = dat.melt[-temp.idx,]
  
  #plot(dat.melt$diff,dat.melt$value,col=dat.melt$From.Pub.Year - 1999)
  
  getrespt <- function(t,xx=0:10)
  {
    wts = exp(-(dat.melt$From.Pub.Year - t)^2 / lam)
    
    ns.temp = ns(dat.melt$diff,df=d)
    
    dat.melt.temp = data.frame(cbind(dat.melt$value,ns.temp))
    
    colnames(dat.melt.temp)[1] = "y"
    
    temp.fit = glm(y ~ .,data=dat.melt.temp,family="poisson",weights=wts)
    
    if(length(xx) < 3)
    {
      xx1 = seq(0,t-xx[1],1)
      if(length(xx) > 1)
      {
        if(length(xx1) < xx[2])
          xx1 = 0:xx[2]
      }
      xx = xx1
    }
    
    temp.pred = predict.glm(temp.fit,newdata = data.frame(predict(ns.temp,newx=xx)),type="response",se.fit=TRUE)
    
    ret = data.frame(tt = t - xx,fit=temp.pred$fit,se.fit=temp.pred$se.fit,year=rep(t,length(xx)))
    return(ret)
  }
  
  years = unique(dat.melt$From.Pub.Year)
  resps <- sapply(years,getrespt,xx=xx)
  
  ret = data.frame(tt = unlist(resps[1,]),fit = unlist(resps[2,]), se.fit = unlist(resps[3,]), year = unlist(resps[4,]))
  
  get_t_vals <- function(t)
  {
    vals.temp = dat.melt[dat.melt$From.Pub.Year == t,]
    to.temp = ret[ret$year == t,"tt"]
    vals.temp = vals.temp$value[match(to.temp,vals.temp$To.Pub.Year)]
    vals.temp
  }
  
  ret$value = c(unlist(sapply(years,get_t_vals)))
  
  return(ret)
}

makeVCMTPlot <- function(dat,name,x.size,width,height,ncol=NULL,every.other=TRUE)
{
  #makes a wrapped plot of VCM fits for a citation database
  #dat is output from fit_pub_vcm
  #name is the plot name
  #x.size is the size of x-axis text
  #ncol is number of columns in faceting
  #every.other will only label every other x-axis label if set to TRUE (to avoid overcrowding)
  
  xlabs2 = paste("'",sapply(dat$tt,function(x){substr(x,3,4)}),"\n(",dat$tt - dat$year,")",sep="")
  dat$xlabs = factor(xlabs2,levels=xlabs2)
  dat$xx = dat$year - dat$tt + 1
  xlab1 = dat$xlabs
  if(every.other)
    xlab1 = xlab1[seq(1,length(levels(dat$xlabs)),by=2)]
  p <- ggplot(dat) + geom_blank(aes(x=xlabs, y=fit)) + scale_x_discrete(breaks=xlab1,labels=xlab1)
  p <- p + geom_line(aes(x=xx,y=fit),col=cbPalette[7]) + geom_ribbon(aes(x=xx,ymin=fit - 2*se.fit,ymax=fit+2*se.fit),col=NA,alpha=.25)
  p <- p + geom_point(aes(x=xlabs,y=value),col=cbPalette[3]) 
  p + labs(x = "year" , y= "Number Citations" , title = "") + theme_bw() + facet_wrap( ~ year , scales = "free", ncol=ncol) + theme(axis.text.x = element_text(size=x.size))
  fname = paste(paste(name,".png",sep=""))
  Sys.sleep(5)
  ggsave(file=fname,width=width,height=height)
}

# Fit VCM models ####

ma.vcm <- fit_pub_vcm(ma_timeseries,xx=1995)
phd.vcm <- fit_pub_vcm(PhD_timeseries,xx=2001)
journal.vcm <- fit_pub_vcm(journal_timeseries[,1:(dim(journal_timeseries)[2]-1)],xx=1980)

# Plots of Fits ####

makeVCMTPlot(phd.vcm,"incoming_phd",10,12,10)
makeVCMTPlot(ma.vcm,"incoming_ma",10,12,10)
makeVCMTPlot(journal.vcm,"incoming_journal",7,14,14,4)

# F tests to compare 4 different incoming citation models ####

get_pub_RSS <- function(dat,lam=1,d=3,type=1,dmax=50)
{
  #Fits a model to incoming citation data
  #Type 1: is the most general, fitting a separate NS with df=d for each source year, number of target citations parametrized by target relative to source year
  #Type 2: incoming citation counts solely depend on target year, and not on source year, e.g. all citations target papers from 1996
  #Type 3: incoming citation counts solely depend on target year relative to source year, e.g. each target year only uses citations from 3 years ago
  #Type 4: is a VCM model with coefficients of an NS for target year relative to source, smoothly varying with source year, somewhere in between type 1 and 2
  #lam is a parameter controlling degree of connection in t (higher lambda is less connection)
  #d is the degrees of freedom of the spline fit for |t-i| (higher d leads to more exact, but less reliable fit)
  #xx describes formatting of the output
  
  dat.melt = melt(dat)
  
  dat.melt$diff = dat.melt$From.Pub.Year - dat.melt$To.Pub.Year
  
  temp.idx = which(dat.melt$diff < 0)
  
  dat.melt = dat.melt[-temp.idx,]
  
  #plot(dat.melt$diff,dat.melt$value,col=dat.melt$From.Pub.Year - 1999)
  if(type == 1)
  {
    getrss.source <- function(t)
    {
      dat.temp = dat.melt[dat.melt$From.Pub.Year == t,]
      
      ns.temp = ns(dat.temp$diff,df=d)
      dat.temp = data.frame(cbind(dat.temp$value,ns.temp))
      colnames(dat.temp)[1] = "y"
      
      temp.fit = glm(y ~ .,data=dat.temp,family="poisson")
    
      ret = c(temp.fit$deviance,d+1)
      names(ret) = c("dev","df")
      return(ret)
    }
    
    rss = sapply(unique(dat.melt$From.Pub.Year),getrss.source)
    ret = rowSums(rss)
    ret = c(ret,type)
    names(ret) = c("dev","df","type")
    return(ret)

  } else if(type == 2)
  {
     getrss.target <- function(d2)
     {
       ns.temp = ns(dat.melt$To.Pub.Year,df=d2)
       dat.temp = data.frame(cbind(dat.melt$value,ns.temp))
       colnames(dat.temp)[1] = "y"
       
       temp.fit = glm(y ~ .,data=dat.temp,family="poisson")
       
       ret = c(temp.fit$deviance,d2+1,temp.fit$aic)
       names(ret) = c("dev","df","aic")
       return(ret)
     }
     
     d.idx = 3:dmax
     rss = sapply(d.idx,getrss.target)
     ret = rss[1:2,which.min(rss[3,])]
     ret = c(ret,type)
     names(ret) = c("dev","df","type")
     return(ret)
     
  } else if(type == 3)
  {
    getrss.relative <- function(d2)
    {
      ns.temp = ns(dat.melt$diff,df=d2)
      dat.temp = data.frame(cbind(dat.melt$value,ns.temp))
      colnames(dat.temp)[1] = "y"
      
      temp.fit = glm(y ~ .,data=dat.temp,family="poisson")
      
      ret = c(temp.fit$deviance,d2+1,temp.fit$aic)
      names(ret) = c("dev","df","aic")
      return(ret)
    }
    
    d.idx = 3:dmax
    rss = sapply(d.idx,getrss.relative)
    ret = rss[1:2,which.min(rss[3,])]
    ret = c(ret,type)
    names(ret) = c("dev","df","type")
    return(ret)
    
  } else 
  {
    
    getrss.vcm <- function(t)
    {
      wts = exp(-(dat.melt$From.Pub.Year - t)^2 / lam)
      
      ns.temp = ns(dat.melt$diff,df=d)
      
      dat.melt.temp = data.frame(cbind(dat.melt$value,ns.temp))
      
      colnames(dat.melt.temp)[1] = "y"
      
      temp.fit = glm(y ~ .,data=dat.melt.temp,family="poisson",weights=wts)
      
      wts = wts / sum(wts)
      df = cbind(rep(1,dim(ns.temp)[1]),ns.temp)
      df = sum(diag(df %*% solve(t(df) %*% diag(wts) %*% df) %*% t(df) %*% diag(wts)))
      
      resid.t = which(dat.melt$From.Pub.Year == t)
      resid.t = sum(resid(temp.fit)[resid.t]^2)
      
      ret = c(resid.t,df)
      names(ret) = c("dev","df")
      return(ret)
      
    }
    
    rss = sapply(unique(dat.melt$From.Pub.Year),getrss.vcm)
    ret = rowSums(rss)
    ret = c(ret,type)
    names(ret) = c("dev","df","type")
    return(ret)
  }
}

ma.rss = sapply(1:4,get_pub_RSS,dat=ma_timeseries,lam=1,d=3)
ma.rss4 = sapply(seq(.01,1,.05),get_pub_RSS,dat=ma_timeseries,type=4,d=3)

phd.rss = sapply(1:4,get_pub_RSS,dat=PhD_timeseries,lam=.5,d=3,dmax=10)
phd.rss4 = sapply(seq(.3,1,.05),get_pub_RSS,dat=PhD_timeseries,type=4,d=3)

journal.rss = sapply(1:4,get_pub_RSS,dat=journal_timeseries,lam=1,d=3)
journal.rss4 = sapply(seq(.2,1,.05),get_pub_RSS,dat=journal_timeseries,type=4,d=3)

ma.p = ma.rss[,2:3] - ma.rss[,4]
ma.p = 1 - pchisq(ma.p[1,],ma.p[2,]*-1)

phd.p = phd.rss[,2:3] - phd.rss[,4]
phd.p = 1 - pchisq(phd.p[1,],phd.p[2,]*-1)

journal.p = journal.rss[,2:3] - journal.rss[,4]
journal.p = 1 - pchisq(journal.p[1,],journal.p[2,]*-1)


#get_pub_RSS(journal_timeseries)

#check new fits with smaller parameters ####

ma.vcm2 <- fit_pub_vcm(ma_timeseries,xx=1995,lam=seq(.01,1,.05)[6])
phd.vcm2 <- fit_pub_vcm(PhD_timeseries,xx=2001,lam=.5)
journal.vcm2 <- fit_pub_vcm(journal_timeseries[,1:(dim(journal_timeseries)[2]-1)],xx=1980,lam=.2)

# Plots of Fits ####

makeVCMTPlot(phd.vcm2,"incoming_phd2",10,12,10)
makeVCMTPlot(ma.vcm2,"incoming_ma2",10,12,10)
makeVCMTPlot(journal.vcm2,"incoming_journal2",7,14,14,4)

# table of deviance p-values ####

ma.vec = c(ma.rss[1:2,2],ma.rss[1:2,3],ma.rss[1:2,4],1,round(ma.p,3))

PhD.vec = c(phd.rss[1:2,2],phd.rss[1:2,3],phd.rss[1:2,4],0.5,round(phd.p,3))

journal.vec = c(journal.rss[1:2,2],journal.rss[1:2,3],journal.rss[1:2,4],1,round(journal.p,3))

inc.tab = rbind(ma.vec,journal.vec,PhD.vec)

colnames(inc.tab) = c("Deviance 1","Df 1","Dev. 2","Df 2","Dev. 3","Df 3","alpha","P(D3<D1)","P(D3<D2)")
rownames(inc.tab) = c("ma","journal","PhD")

library(xtable)

print(xtable(inc.tab,label="tab:tab1",caption="Supplementary Table 1"))

#getting distributions of null generalized log likelihood for both hypothesis tests
#get_pub_RSS <- function(dat,lam=1,d=3,type=1,dmax=50)
get_null_loglik <- function(dat,lam=1,d=3,type=2,B=1000)
{
  #Fits a model to incoming citation data
  #Type 2: incoming citation counts solely depend on target year, and not on source year, e.g. all citations target papers from 1996
  #Type 3: incoming citation counts solely depend on target year relative to source year, e.g. each target year only uses citations from 3 years ago
  #both of these will be null hypotheses compared against the type 4 model:
  #Type 4: is a VCM model with coefficients of an NS for target year relative to source, smoothly varying with source year, somewhere in between type 1 and 2
  #lam is a parameter controlling degree of connection in t (higher lambda is less connection)
  #d is the degrees of freedom of the spline fit for |t-i| (higher d leads to more exact, but less reliable fit)
  #xx describes formatting of the output
  
  dat.melt = melt(dat)
  
  dat.melt$diff = dat.melt$From.Pub.Year - dat.melt$To.Pub.Year
  
  temp.idx = which(dat.melt$diff < 0)
  
  dat.melt = dat.melt[-temp.idx,]
  
  #get null model
  if(type == 2)
  {
    get.nullmodel <- function(d2,dat.melt)
    {
      ns.temp = ns(dat.melt$To.Pub.Year,df=d2)
      dat.temp = data.frame(cbind(dat.melt$value,ns.temp))
      colnames(dat.temp)[1] = "y"
      
      temp.fit = glm(y ~ .,data=dat.temp,family="poisson")
      
      ret = list(fit=temp.fit,df=d2)
      return(ret)
    }
    
  } else if(type == 3)
  {
    get.nullmodel <- function(d2,dat.melt)
    {
      ns.temp = ns(dat.melt$diff,df=d2)
      dat.temp = data.frame(cbind(dat.melt$value,ns.temp))
      colnames(dat.temp)[1] = "y"
      
      temp.fit = glm(y ~ .,data=dat.temp,family="poisson")
      
      ret = list(fit=temp.fit,df=d2)
      return(ret)
    }
  }
  
  model0 = get.nullmodel(d,dat.melt)
  
  getmu.vcm <- function(t,d2,dat.melt)
  {
    wts = exp(-(dat.melt$From.Pub.Year - t)^2 / lam)
    
    ns.temp = ns(dat.melt$diff,df=d2)
    
    dat.melt.temp = data.frame(cbind(dat.melt$value,ns.temp))
    
    colnames(dat.melt.temp)[1] = "y"
    
    temp.fit = try(glm(y ~ .,data=dat.melt.temp,family="poisson",weights=wts),silent=TRUE)
    
    #now need to predict solely for those y with from == t
    idx = which(dat.melt$From.Pub.Year == t)
    
    if(class(temp.fit) == "try-error")
    {
      mu = 0*idx
    } else
    {
      mu = predict(temp.fit,dat.melt.temp[which(dat.melt$From.Pub.Year == t),],type="response")  
    }
    
    ret = cbind(dat.melt$To.Pub.Year[idx],rep(t,length(idx)),mu)
    colnames(ret) = c("to","from","value")
    
    return(ret)
  }
  
  compute.lnb <- function(sim=TRUE)
  {
    #first draw from mode0 to get new data
    if(sim)
    {
      yb = predict(model0$fit,type="response")
      yb = rpois(length(yb),yb)
    } else
    {
      yb = dat.melt$value
    }
    
    dat.meltb = dat.melt
    dat.meltb$value = yb
    
    #get null model mean for new data yb
    model0b = get.nullmodel(d,dat.meltb)
    mu0b = predict(model0b$fit,type="response")
    
    #get VCM model mean for new data yb
    mu.vcmb = do.call("rbind",sapply(unique(dat.meltb$From.Pub.Year),getmu.vcm,dat.melt=dat.meltb,d2=d))
    mu.vcmb = mu.vcmb[,"value"]
    
    ln = log(mu.vcmb / mu0b) * yb + (mu0b - mu.vcmb)
    return(-1 * sum(ln))
  }
  
  lnbs = replicate(B,compute.lnb())
  ln.star = compute.lnb(FALSE)
  return(list(lnbs=lnbs,ln.star=ln.star))
}

fname = "vcm_lnb.Rdata"
lnb.ma.2 = get_null_loglik(ma_timeseries,type=2,B=1000)
lnb.ma.3 = get_null_loglik(ma_timeseries,type=3,B=1000)

save(lnb.ma.2,lnb.ma.3,file=fname)

lnb.phd.2 = get_null_loglik(PhD_timeseries,type=2,B=1000,lam=.5)
lnb.phd.3 = get_null_loglik(ma_timeseries,type=3,B=1000,lam=.5)

save(lnb.ma.2,lnb.ma.3,lnb.phd.2,lnb.phd.3,file=fname)

lnb.journal.2 = get_null_loglik(journal_timeseries[,1:(dim(journal_timeseries)[2]-1)],type=2,B=1000)
lnb.journal.3 = get_null_loglik(journal_timeseries[,1:(dim(journal_timeseries)[2]-1)],type=3,B=1000)

save(lnb.ma.2,lnb.ma.3,lnb.phd.2,lnb.phd.3,lnb.journal.2,lnb.journal.3,file=fname)

#print(xtable(inc.tab,label="tab:relevance",caption="Relevance table for all bins found significant"),size="small")

#look at the lnb distributions

sort(lnb.ma.2$lnbs)
#lnb = -1 * ( ln - l0) = log( L0 / LN )
#So when it's postitive -> L0 is a better fit than LN
#And when it's negative -> LN is a better fit that L0
#I ran simulations under the null that L0 is true, hence sometimes LN does really poorly

png("ln_histograms.png",width=1000,height=1000)
par(mfcol=c(2,3))
xlims = c(-20,50)
breakz = 20
cexz = 2.5
idx = which(lnb.ma.2$lnbs < 1e10)
name = "Histogram of D(H1,H0)b: \n"
name = paste(name,"MA, Stag")
hist(-1*lnb.ma.2$lnbs[idx],breaks=breakz,main=name,xlim=xlims,cex.main=cexz,xlab="",ylab="")

#sort(lnb.ma.3$lnbs)
name = "Histogram of D(H1,H0)b: \n"
name = paste(name,"MA, Flow")
hist(-1*lnb.ma.3$lnbs,main=name,breaks=breakz,xlim=xlims,cex.main=cexz,xlab="",ylab="")

idx = which(lnb.phd.2$lnbs < 1e10)
name = "Histogram of D(H1,H0)b: \n"
name = paste(name,"PhD, Stag")
hist(-1*lnb.phd.2$lnbs[idx],breaks=breakz,main=name,xlim=xlims,cex.main=cexz,xlab="",ylab="")

#sort(lnb.ma.3$lnbs)
name = "Histogram of D(H1,H0)b: \n"
name = paste(name,"PhD, Flow")
hist(-1*lnb.phd.3$lnbs[which(lnb.phd.3$lnbs < 500)],main=name,breaks=breakz,xlim=xlims,cex.main=cexz,xlab="",ylab="")

idx = which(lnb.journal.2$lnbs < 1e10)
name = "Histogram of D(H1,H0)b: \n"
name = paste(name,"Journal, Stag")
hist(-1*lnb.journal.2$lnbs[idx],breaks=breakz,main=name,xlim=xlims,cex.main=cexz,xlab="",ylab="")

#sort(lnb.ma.3$lnbs)
name = "Histogram of D(H1,H0)b: \n"
name = paste(name,"Journal, Flow")
hist(-1*lnb.journal.3$lnbs,main=name,breaks=breakz,xlim=xlims,cex.main=cexz,xlab="",ylab="")
dev.off()

#get p-values for each one
getpvallnb <- function(lnb)
{
  ret = ( sum(lnb$lnbs < lnb$ln.star,na.rm=T) ) / ( length(lnb$lnbs) - sum(is.na(lnb$lnbs)) )
  if(ret == 0)
    return(1 / ( length(lnb$lnbs) - sum(is.na(lnb$lnbs)) ))
  return(ret)
}

getpvallnb(lnb.ma.2)
getpvallnb(lnb.ma.3)

getpvallnb(lnb.phd.2)
getpvallnb(lnb.phd.3)

getpvallnb(lnb.journal.2)
getpvallnb(lnb.journal.3)

#all p-vals are <= 0.001

