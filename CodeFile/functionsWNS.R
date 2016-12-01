## functionsWNS.R



########### Bat Paramerters ###########

load.bat <- function(species = speciesOption){
  ### Function that loads the parameter set for the chosen bat species ###
  ## Argument: species - species to use, takes one of the following values:
  ## "M.lucifugus" for little brown bats
  ## "M.myotis" for European myotis
  ## "E.serotinus" for serotine bats
  ## "E.fuscus" for big brown bats
  ## "M.yumanensis" for Yuma myotis
  ## "M.californicus" for California myotis
  ### Note: Can be used to call parameters for anything existing as a named column in species.parms.csv
  pars <- read.csv("parameterFiles/species.parms.csv",sep=",",head=T)
  if(species%in%names(pars)){
    # Extract and return parameter set for the desired species
    pset <- subset(pars, select = c("Parameter",species))
  }else{
    warning("Unknown species option selected.")
  }
  par <- subset(pset,Parameter != "Ta") # Remove Ta!
  batParams <- get(species,par)
  names(batParams) <- par$Parameter
  ### PARAMETERS FOR THE BAT ENERGETICS MODEL ###
  # RMR             Rest met rate
  # TMRmin          Torpor met rate
  # Teu             Temp euthermic
  # Tlc             lower crit temp during euthermia
  # Ttormin         Temp defended to stop freezing
  # Ceu             euthermic thermal conductance
  # Ct              torpor conductance
  # S               spec heat capacity of tissues
  # k               constant
  # ttormax         max torpor time
  # tar             time to arouse
  # teu             time euthermic
  # mass            total bat mass
  
  return(batParams)	
}

########### Fungal Growth Dynamics ############

load.fung <- function(growth=growthOption){
  ### Function that loads the chosen parameter set for fungal growth and scales ###
  ##Argument: growth - rate parameter set to use, takes one of the following values:
    ## "Chaturvedi" (greater growth) from Chaturvedi et al. PLoS One
    ## "Verant" (slower growth) from Verant et al. PLoS One
  
  pars <- read.csv("parameterFiles/rate.parms.csv")
  if(growth%in%names(pars)){
    # Extract and return desired parameter set
    pset <- subset(pars, select = c("pars",growth))
  }else{
    warning("Unknown option selected for growth rate parameters.")
  }
  rate.par <- pset	#load rates from growth option
  scaled <- read.csv("parameterFiles/humid.parms.csv",header=T)
  colnames(scaled) <- colnames(rate.par)
  parms <- rbind(scaled,rate.par)
  growthParams <- get(growth,parms)
  names(growthParams) <- c("mu1","mu2","beta1","beta2","beta3")
  ### PARAMETERS FOR THE FUNGAL GROWTH MODEL ###
  # mu1             scaling parameter for Michaelis-Menton function/Rel. Humidity
  # mu2             scaling parameter for Michaelis-Menton function/Rel. Humidity
  # beta1           temp-dept hourly rate shape par
  # beta2           temp-dept hourly rate shape par
  # beta3           temp-dept hourly rate shape par - max growth
  return(growthParams)
}

##### Model Engines ####
detModel <- function(t,y,parms){
  require(deSolve)
  ### Fucntion for creating the deterministic model of fungal growth
  ##Arguments: t <- temperature
  ##Arguments: y <- dependant varriables
  ##Arguments: parms <- parameters related to fundal growth
  with(c(as.list(y),parms),{
    ttor  <-  calcTorporTime(Ttor,FungalArea,WNS,parms)
    
    dpTdt <- pE/teu - pT/ttor # change in TorporProp (pT) / dt
    dpEdt <- pT/ttor - pE/teu # change in EuthermicProp (pE) / dt
    dJdt <- Eeu*pE + Etor*pT + Ear*pT/ttor # change in EnergyConsumed / dt
    dFdt <- growth*pT # change in FungalArea / dt
    
    list(c(dpTdt,dpEdt,dJdt,dFdt))
  })
}

max_to_current <- function(x) {
  cummax(x)[-1]
}

dynamicEnergyPd <- function(env.df, WNS=TRUE, modelParams = c(growthParams,batParams)){
  require(deSolve); require(data.table)
  ### Function for determining the growth area of Pd, and the total amount of enegry consumed
  ##Arguments: Ta <- range of ambient temperture
  ##Arguments: twinter <- range of time of witner
  ##Arguments: Hd <- range of ambient relative humidity
  ##Arguments: WNS <- is WNS in the model (Logical)
  ##Arguments: modelParams <- model parameters loaded through functionsWNS.R
  blah <- with(as.list(modelParams),{
    Ta <- env.df[[1]]
    Hd <- env.df[[2]]
    if(beta3>=Teu){
      warning("The model assumes fungal growth does not to occur at euthermic temperatures.\nThis assumption is violated in the specified parameter range.")
    }		
    Ttor <- ifelse (Ta > Ttormin,Ta,Ttormin) # determine torpid body temperature, given ambient temperature
    Tb  <-  ifelse(Ttor < Teu,Ttor,Teu) # determine body temperature
    values <- c(Ttor = Ttor, WNS = WNS,
                #- Fungal growth rate -#
                growth = fungalGrowthRate(Tb,modelParams)*scaleFungalGrowthRate(Hd,modelParams), # fungal growth rate as a function of body temperature and humidity
                #- Energy use -#
                Eeu  =  calcEnergyPerTimeEuthermic(Ta, modelParams),
                Etor  =  calcEnergyPerTimeTorpor(Ta, modelParams),
                Ear  =  calcEnergyArousal(Ttor, modelParams),
                #- Predefined parameters -#
                modelParams)
    # run deterministic model with dynamic fungal growth\
    ts <- data.table(lsoda( y = c(pT=1,
                                  pE=0,
                                  EnergyConsumed=0,
                                  FungalArea=0),
                            times = twinter,
                            func = detModel,
                            parms = values))
    
    # EnergyConsumed should be monotone increasing, so we should just be able to grab the energy
    # consumed at each time point, but in case it can go down in some other variation of the model
    # we'll instead take maximum energy consumed up to the given point in time.  The same applies
    # to fungal area.
    Ewinter  <-  max_to_current(ts$EnergyConsumed)    # energy used over full winter duration
    fatConsumed  <-  convertToFat(Ewinter) # convert Ewinter to grams of fat consumed (over full winter duration)
    results <- data.table(Ta = Ta, Humidity = Hd, cbind(GfatConsumed=c(0,fatConsumed), Pdgrowth=c(0,max_to_current(ts$FungalArea)), time=ts$time)) # results - grams of fat consumed, plus area of total Pd growth (Pd)
    return(results)
  })
}

########### Model functions ############

# how Q varies with ambient temperature (ambTemp = Ta)
calcQ <- function(ambTemp,k1=1.6,k2=0.26,k3=0.006){
  k1 + k2*ambTemp - k3*ambTemp^2
}
# Methods equation 2 - how fungal growth rate varies with body temperature (bodyTemp = Tb)
fungalGrowthRate <- function(bodyTemp, params, Tmin=0){
  with(as.list(params),{
    ifelse(bodyTemp>beta3|bodyTemp<=Tmin,0,beta1*(bodyTemp-Tmin)*(1-exp(beta2*(bodyTemp-beta3))))
  })}
# Methods equation 3 - how fungal growth rate scales with relative humidity (pctRH = Hd)
scaleFungalGrowthRate <- function(pctRH, params){
  with(as.list(params),
       mu1*pctRH/(1+(mu2*pctRH))
  )}
# Methods equation 6 - Calculate energy per time at euthermic body temperature, given ambient temperature
calcEnergyPerTimeEuthermic <- function(ambTemp, params){
  with(as.list(params),
       RMR + (Tlc-ambTemp)*Ceu
  )}
# Methods equations 7 and 8 - Calculate energy per time at torpid body temperature, given ambient temperature
calcEnergyPerTimeTorpor <- function(ambTemp, params, Q=calcQ(ambTemp)){
  with(as.list(params),{
    ifelse (ambTemp > Ttormin, 
            TMRmin*Q^((ambTemp-Ttormin)/10),
            TMRmin + (Ttormin-ambTemp)*Ct)
  })}
# Calculate energy needed for arousal, given torpid body temperature
calcEnergyArousal <- function(torpidTemp, params){
  with(as.list(params),
       (Teu - torpidTemp)*S
  )}
# Methods equations 9 and 10 - Caluculate time in torpor as a function of ambient temperature and fungal area
calcTorporTime <- function(ambTemp, areaPd, inf, params, Q=calcQ(ambTemp)){
  with(as.list(params),{
    # how time in torpor varies with ambient temp (Pd absent)
    time <- ifelse(ambTemp > Ttormin, 
                   (ttormax/Q^((ambTemp-Ttormin)/10)), 
                   ttormax/(1+(Ttormin-ambTemp)*Ct/TMRmin))
    if(inf==TRUE){
      areaPd <- ifelse (areaPd < 1, 1, areaPd)
      time  <-  time/areaPd # fungal growth reduces ttor
    }
    return(time)
  })
}
# conversion for energy into grams of fat
convertToFat <- function(energy,k1=20.1,k2=39.3,k3=100){
  energy*k1/(k2*k3)
}


energy.Pd  <-  function(Ta, twinter, Hd, WNS, modelParams = c(growthParams,batParams)){
  ### Function that calculates energy use and colony area as a function of environement and infection status ###
  ## Argument: Ta - ambient temperature (degrees C), takes a numeric value
  ## Argument: twinter - winter duration in hours, takes a numeric value (should be <8767)
  ## Argument: Hd - ambient relative humidity (percentage), takes a numeric value <=100
  ## Argument: WNS - infection status, takes a logical value (TRUE if infected; FALSE if free of fungus)
  ## Argument: modelParams - parameters for fungal growth and bat energetics, as specified 
  with(as.list(modelParams),{
    Ttor <- ifelse (Ta > Ttormin,Ta,Ttormin) # determine torpid body temperature, given ambient temperature
    Tb  <-  ifelse(Ttor < Teu,Ttor,Teu) # determine body temperature
    #- Fungal growth -#
    growthrate  <- fungalGrowthRate(Tb,modelParams) # fungal growth rate as a function of body temperature
    growthrate_H <- scaleFungalGrowthRate(Hd,modelParams) # scaling factor for fungal growth based on humidity
    Pd <- growthrate*twinter*growthrate_H # Methods equation 4 - total Pd growth
    #- Energy use -#
    Eeu  <-  calcEnergyPerTimeEuthermic(Ta, modelParams)
    Etor  <-  calcEnergyPerTimeTorpor(Ta, modelParams)
    Ear  <-  calcEnergyArousal(Ttor, modelParams) 
    ttor  <-  calcTorporTime(Ta,Pd,WNS,modelParams)
    Ebout  <-  Eeu*teu + Etor*ttor + Ear*tar # Methods equation 5 - calculate total energy use for a torpor/arousal cycle
    Ewinter  <-  (twinter/(ttor+teu+tar))*Ebout    # calculate energy used over full winter duration
    fatConsumed  <-  convertToFat(Ewinter) # convert Ewinter to grams of fat consumed (over full winter duration)
    results <- c(grams=fatConsumed, area=Pd) # results - grams of fat consumed, plus area of total Pd growth (Pd area in CM)
    return(results)
  })
}


