require(deSolve)
require(ggplot2)
require(matlab)

dataloader <- read.csv("owid-covid-data (3).csv")


loadcountry <- function(country, code) {
  country <- dataloader[dataloader$iso_code == code,c(5:6, 9,23)]
  x <- 1 : nrow(country)
  rownames(country) <- x
  country <- cbind(country, x)
  country <- cbind(country, log10(country$new_cases))
  return(country)
}
netherlands <- loadcountry(netherlands, "NLD")
denmark <- loadcountry(denmark, "DNK")
sweden <- loadcountry(sweden, "SWE")
norway <- loadcountry(norway, "NOR")
uk <- loadcountry(uk, "GBR")
spain <- loadcountry(spain, "ESP")
spain <- rbind(spain, c(0,0,0,0,0,0))
germany <- loadcountry(germany, "DEU")
italy <- loadcountry(italy, "ITA")
france <- loadcountry(france, "FRA")
canada <- loadcountry(canada, "CAN")
usa <- loadcountry(usa, "USA")


plotcountry <- function(country, name, yn) {
  lim <- 10000
  if(yn == TRUE) {
    lim <- 100000
  }
  plot(country$new_cases~country$x, log = "y" , ylim = c(1,lim), xlim = c(0,350), main = name, xlab = "time starting on 12/31/2019", ylab = "X2", col = "blue")
}


asolver <- function(r,y,n) {
  a = (r * y) / n
  print(a)
  return(a)
}

nsolver <- function(p,r,y) {
  n = p / (1 - ((1 + log(r)) / r))
  print(n)
  return(n)
}

ysolver <- function(r,val,p) {
  y = val / (r - 1)
  print(y)
  return(y)
}

sirmodel <- function(t,state,params) {
  x1 = state[1]
  x2 = state[2]
  alpha = params[1]
  gamma = params[2]
  dx1dt <- -alpha * x1 * x2
  dx2dt <- alpha * x1 * x2 - gamma * x2
  dxdt <- c(dx1dt,dx2dt)
  list(dxdt)
}

solvefunction <- function(vals,p,val,times, name, valstart, valfinal, x2start) {
  ee <- 100
  j <- 1
  name <- name[name$x > valstart,]
  name <- name[name$x <= valfinal,]
  name <- name[name$x > ee,]
  for(i in vals) {
    y <- ysolver(i,val,p)
    n <- nsolver(p, i, y)
    a <- asolver(i,y,n)
    xstart <- c(x1 = n, x2 = x2start)
    params <- c(a,y)
    out <- as.data.frame(
      ode(func=sirmodel, y=xstart, times = times, parms = params)
    )
    if(j == 1) {
      lines(x2~time,data = out, col = "red", lwd = 2)
      out <- out[out$time > ee,]
      out <- out[out$time <= 263,]
      if(nrow(name) == 162) {
        name <- rbind(name, c(0,0,0,0,0,0,0,0,0))
      }
      resid1 <- out$x2
      name <- cbind(name, resid1)
    }
    else if(j == 2) {
      lines(x2~time,data = out, col = "black", lwd = 2)
      out <- out[out$time > ee,]
      out <- out[out$time <= 263,]
      if(nrow(name) == 162) {
        name <- rbind(name, c(0,0,0,0,0,0,0,0,0))
      }
      resid2 <- out$x2
      name <- cbind(name, resid2)
    }
    else if(j==3) {
      lines(x2~time,data = out, col = "red", lwd = 2)
      out <- out[out$time > ee,]
      out <- out[out$time <= 263,]
      if(nrow(name) == 162) {
        name <- rbind(name, c(0,0,0,0,0,0,0,0,0))
      }
      resid3 <- out$x2
      name <- cbind(name, resid3)
    }
    j <- j + 1
  }
  return(name)
}

improvedcasesolver <- function(timeStart,timeEnd,p,slope,country,vals, caseresiduals,x2start) {
  times2 <- seq(timeStart+1,timeEnd,1)
  p2 <- p
  slope2 <- slope
  caseresiduals <- solvefunction(vals,p2,slope2,times2,country,timeStart,timeEnd, x2start)
  return(caseresiduals)
  }
deathsolver <- function(country,percent,scalednum,vals,num,p, name) {
  plot(country$new_deaths~rownames(country), log = "y", ylim = c(1,10000),xlim = c(50,300), main = name, xlab = "time starting on 12/31/2019", ylab = "X4", col = "red")
  scaled <- country$new_cases * percent
  country <- cbind(country, scaled)
  times <- seq(scalednum,300,1)
  name <- solvefunction(vals,p*percent,num,times,country, 1, 300, 1)
  return(name)
}

ploterror <- function(caseresiduals, name, bool, length) {
  caseresiduals <- cbind(caseresiduals, sqrt((((caseresiduals$resid3-caseresiduals$resid2) * (caseresiduals$resid3-caseresiduals$resid2)) + ((caseresiduals$resid1-caseresiduals$resid2) * (caseresiduals$resid1-caseresiduals$resid2)) / 2)))
  caseresiduals <- cbind(caseresiduals, rep(0,length))
  colnames(caseresiduals)[11] <- "pvalues"
  colnames(caseresiduals)[10] <- "stdev"
  caseresiduals <- cbind(caseresiduals, c(1:length))
  if(bool == TRUE) {yn = caseresiduals$new_cases }
  if(bool == FALSE) {yn = caseresiduals$new_deaths}
  for(i in caseresiduals$`c(1:length)`) {
    if(yn[i] > caseresiduals$resid2[i]) {
      caseresiduals$pvalues[i] <- pnorm(yn[i], mean = caseresiduals$resid2[i], sd = caseresiduals$stdev[i], lower.tail = FALSE)
    }
    else if(yn[i] < caseresiduals$resid2[i]) {
      caseresiduals$pvalues[i] <- pnorm(yn[i], mean = caseresiduals$resid2[i], sd = caseresiduals$stdev[i], lower.tail = TRUE)
    }
  }
  caseresiduals <- caseresiduals[yn > 0,]
  plot(caseresiduals$pvalues~caseresiduals$x, main = name, xlab = "time starting on 12/31/2019", ylab = "p value", col = "darkorchid1", cex = 1.25)
  abline(h=0.05, col = "blue", cex = 2)
}

doallplots <- function(country, timeStart, timeEnd, valsStart, valsEnd, p, slope, delta, newStart, countryName, p2,slope2,x2start, bool) {
  plotcountry(country, paste(countryName, "- Cases"), bool)
  times <- seq(timeStart, timeEnd, 1)
  vals <- seq(valsStart, valsEnd, 0.5)
  caseresiduals <- solvefunction(vals, p, slope, times, country, 1, timeEnd,x2start)
  caseresiduals2 <- improvedcasesolver(timeEnd,350,p2,slope2,netherlands, vals, caseresiduals,x2start)
  deathresiduals <- deathsolver(country, delta, newStart, vals, slope, p, paste(countryName, "- Deaths"))
  caseresiduals3 <- merge(caseresiduals, caseresiduals2,all = TRUE)
  print(ggplot(caseresiduals3) + geom_point(aes(x=x, y=new_cases-resid2), stat="identity", fill="skyblue", alpha=0.7) + geom_errorbar(aes(x=x,ymin=new_cases-resid3,ymax=new_cases-resid1), width=0.4) + ggtitle(paste("Residual Plot-Cases:",countryName)) + xlab("time starting on 12/31/2019") + ylab("Residual"))
  print(ggplot(deathresiduals) + geom_point(aes(x=x, y=new_deaths-resid2), stat="identity", fill="skyblue", alpha=0.7) + geom_errorbar(aes(x=x,ymin=new_deaths-resid3,ymax=new_deaths-resid1), width=0.4) + ggtitle(paste("Residual Plot-Deaths", countryName)) + xlab("time starting on 12/31/2019") + ylab("Residual"))
  ploterror(caseresiduals3, paste("Associated p value for Cases-", countryName), TRUE, 163)
  ploterror(deathresiduals, paste("Associated p value for Deaths-", countryName), FALSE, 163)
}

doallplots(netherlands, 68,185,3.9,4.9,1200,0.193,0.15,64, "Netherlands", 2000, 0.09,30, FALSE)
doallplots(denmark, 77,200,3.7,4.7,300,0.192,0.08,77, "Denmark", 800,0.09,20, FALSE)
doallplots(sweden, 52,263,3.6,4.6,600,0.192,0.18,76, "Sweden", 0, 0, 1, FALSE)
#no further model possible for Sweden
doallplots(norway, 65,205,3.6,4.6,250,0.25,0.05,85, "Norway",150,0.11,10, FALSE)
doallplots(uk, 75,210,10,11,5300,0.255,0.21,69, "United Kingdom",6000,0.068,220, FALSE)
doallplots(germany, 65,205,4.2,5.2,6300,0.25,0.05,73, "Germany",2000,0.085,150, FALSE)
doallplots(italy, 60,210,5,6,6000,0.25,0.14,54, "Italy",1800,0.1,100, FALSE)
doallplots(canada, 80,220,8,9,1800,0.2,0.1,80, "Canada",1500,0.06,150, FALSE)
doallplots(france, 66,170,4,5,4600,0.233,0.19,62, "France",12000,0.07,100, FALSE)
doallplots(usa, 87,170,19,20,25000,0.31,0.1,70, "USA",60000,0.18,11000, TRUE)
doallplots(spain, 65,170,5,6,8200,0.31,0.11,63, "Spain",12000,0.08,120, FALSE)

#fin