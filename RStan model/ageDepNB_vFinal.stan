// Description: Bayesian model of age-dependent transmission when there is
// contact-based transmission and environmental transmission
// The model reads in daily age-distributed hospitalisation data (and non-age distr removal data)
// as well as population pyramid and contact matrix

// note convention throughout is that variables with D_ in front all involve data coming as input
data {
  int <lower=0> D_NDays; // numDays for array lengths
  int <lower=0> D_NAgeGroups; // numAgegroups for array lengths
  int <lower=0> D_cumulHospByAge[D_NAgeGroups,68]; // the daily hosp data is cumulative
  int <lower=0> D_ageDisIn[D_NAgeGroups]; // age distribution
  real <lower=0> D_meanAges[D_NAgeGroups]; // useful to have the mean ages of the cohorts
  real <lower=0> D_removedsByAge_frac[68]; // the daily removed data is cumulative
  real <lower=0> D_contactMPre[D_NAgeGroups,D_NAgeGroups]; // prelockdown CM
  real <lower=0> D_contactMPost[D_NAgeGroups,D_NAgeGroups]; // Post lockdown not used!
  int <lower=0, upper=4> D_runType; // binary for changing assumption, not currently working
  int <lower=0, upper=10> D_negTau; // look forward delay
  int <lower=0, upper=10> D_posTau; // look backward delay
  real <lower=0> D_meancontact; // useful to have means (by national pop) of age cohorts
} 

parameters {
  real<lower=0> transMC; //
  real<lower=0> transSC; //
  real<lower=0> transMU; //
  real<lower=0> transSU; //
  real<lower=0> phi_inv; //
  real<lower=1, upper=100000> relSev_inv[(!(D_runType>0)) ? D_NAgeGroups : 0];// conditionally defining these parameters in an on/off manner as there are multiple run scenarios
  real<lower=1, upper=100> coefExp_inv[(D_runType>0) ? 1 : 0];   
  real<lower=0, upper=0.5> ageExp[(D_runType>0) ? 1 : 0];   
  real<lower=-1, upper=1> deviationFromExp[(D_runType>0) ? D_NAgeGroups : 0];   
  real<lower=-1, upper=1> deviationFromUE[(D_runType==2) ? D_NAgeGroups : 0];   
  real<lower=-1, upper=1> deviationFromUT[(D_runType>2) ? D_NAgeGroups : 0];   
  real<lower=-1, upper=1> deviationFromUT_S[(D_runType==4) ? D_NAgeGroups : 0];  
  //real<lower=0, upper=10> juvFac2[D_NAgeGroups]; KEEP FOR USE IF EXPLORING GRADIENT IN INFECTIOUSNESS
}
transformed parameters {
  real phi=1./phi_inv;
  real<lower=0, upper=1> coefExp[(D_runType>0) ? 1 : 0];  
  real<lower=0, upper=1> ProbSevere[(!(D_runType>0)) ? D_NAgeGroups : 0];
  vector[D_NAgeGroups] ageDis;
  vector[D_NAgeGroups] contacts_vec;
  vector[D_NAgeGroups] SevereInc;
  vector[D_NAgeGroups] MildInc;  
  vector[D_NAgeGroups] transPrevC;
  vector[D_NAgeGroups] transPrevU;
  vector[D_NAgeGroups] multSevereToAll;
  real transCoef;
  real propS;
  real transMain;
  real multInf[D_NAgeGroups];
  real NBintensity[D_NAgeGroups,(D_NDays)];
  
  if (D_runType>0) 
  {
    coefExp[1]=1./coefExp_inv[1];   
  }
  else{
    for (ss in 1:D_NAgeGroups) 
    {
      ProbSevere[ss]=1/relSev_inv[ss];
    }
  }
    for (ss in 1:D_NAgeGroups) 
  {
    if (D_runType>0){
  multSevereToAll[ss]= (coefExp[1]*exp(ageExp[1]*D_meanAges[ss])/exp(ageExp[1]*D_meanAges[8]))+deviationFromExp[ss];
    }else{
      multSevereToAll[ss]= ProbSevere[ss];
    }
    if (D_runType==2){
      multInf[ss]=1+deviationFromUE[ss]; 
    }else{
      multInf[ss]=1;  // modelling of uniform infectiousness (above for gradient through deviation)  
    }
  } 
  ageDis=to_vector(D_ageDisIn);

  for (kk in (1+D_negTau):(D_NDays+D_negTau)) {        //sampling distribution
      // different components:
      SevereInc=to_vector(D_cumulHospByAge[,(kk)])*(1-D_removedsByAge_frac[(kk)]);
      // scale up from severe to mild incidence via prob severe
      MildInc=(((SevereInc)./(to_vector(multSevereToAll))))-SevereInc;
      if (D_runType==3) {
        transPrevC=(((SevereInc*transSC)+((transMC*MildInc).*(to_vector({1,1,1,1,1,1,1,1})+to_vector(deviationFromUT))))./ageDis); 
        transPrevU=((SevereInc*transSU)+((transMU*MildInc).*(to_vector({1,1,1,1,1,1,1,1})+to_vector(deviationFromUT))));        
      }else if (D_runType==4) {  
        transPrevC=(((SevereInc*transSC).*(to_vector({1,1,1,1,1,1,1,1})+to_vector(deviationFromUT_S))+((transMC*MildInc).*(to_vector({1,1,1,1,1,1,1,1})+to_vector(deviationFromUT))))./ageDis); 
        transPrevU=((SevereInc*transSU).*(to_vector({1,1,1,1,1,1,1,1})+to_vector(deviationFromUT_S))+((transMU*MildInc).*(to_vector({1,1,1,1,1,1,1,1})+to_vector(deviationFromUT))));        
        
      }else{
        transPrevC=(((SevereInc*transSC)+(MildInc*transMC))./ageDis); 
        transPrevU=((SevereInc*transSU)+(MildInc*transMU)); 
      }
      
      
      for (jj in 1:D_NAgeGroups) 
      { 
        contacts_vec=to_vector(D_contactMPre[jj,]); // grabbing relevant row of CM for age jj
        propS=fmax(1-(((D_cumulHospByAge[jj,(kk)])/(multSevereToAll[jj]))*(1/ageDis[jj])),0); //prob. jj is suscept.
        transCoef=D_ageDisIn[jj]*multSevereToAll[jj]*multInf[jj];  //various components of tranmission rates
        transMain=sum(((transPrevC).*(contacts_vec))+(transPrevU*(D_meancontact/mean(ageDis)))); // main transm term  
        NBintensity[jj,(kk-D_negTau)]=transCoef*propS*transMain;  // combining all the bits for overall transmission
        
        
        // below block is for investigating exposure patterns (NOT in MS) see also gen. quant. section
        //NBintensity2[jj]=NBintensity2[jj]+sum(((transPrevC).*(contacts_vec)));
        //if ((kk-D_negTau)<2) {
        //  NBintensity3a[jj]=NBintensity3a[jj]+sum(((transPrevC).*(contacts_vec)));
        //}
        //if ((kk-D_negTau)<9) {  
        //  NBintensity3[jj]=NBintensity3[jj]+sum(((transPrevC).*(contacts_vec)));
        //}else if ((kk-D_negTau)<17) {
        //  NBintensity4[jj]=NBintensity4[jj]+sum(((transPrevC).*(contacts_vec)));
        //}else if ((kk-D_negTau)<25) {
        //  NBintensity5[jj]=NBintensity5[jj]+sum(((transPrevC).*(contacts_vec)));
        //}
          
      }
  }
}

model {
  // many of these vectors are only defined because their equivalent arrays (often declared earlier) are the default but some things are easier as vectors 
  int D_SCasesByAge;
  // defining Bayesian priors (default is gamma, which is the conjugate prior of Poisson; we are using NB but it is not the probability based formulation hence gamma appropriate)
  transMC ~ gamma(1,5.0);
  transSC ~ gamma(1,5.0);
  transMU ~ gamma(1,5.0);
  transSU ~ gamma(1,5.0);
  phi_inv ~ normal(2,0.5);
  
  // Cauchy priors are here defined for inverse of parameters that represent probabilities only
  if (D_runType>0){
    ageExp[1] ~ gamma(1,5.0); 
    coefExp_inv[1] ~ cauchy(0., 1); 
  }
  for (ss in 1:D_NAgeGroups) 
  {
    if (D_runType>0){
        deviationFromExp[ss] ~ normal(0,0.01); 
    if (D_runType==2){  
        deviationFromUE[ss] ~ normal(0,0.01); 
    } if (D_runType>2){  
        deviationFromUT[ss] ~ normal(0,0.01);  
    } if (D_runType==4){  
    deviationFromUT_S[ss] ~ normal(0,0.01); 
    }  
    }else{
      relSev_inv[ss] ~ cauchy(0., 1);  // for initial Fig. where each age has its own separately modelled Psevere
    }
  } 

  for (pp in 1:(D_NDays)) {  
      for (jj in 1:D_NAgeGroups) 
      { 
        D_SCasesByAge=D_cumulHospByAge[jj,pp+1+D_posTau]-D_cumulHospByAge[jj,pp+D_posTau]; // getting new cases from the data
        D_SCasesByAge ~ neg_binomial_2(NBintensity[jj,pp],pow(NBintensity[jj,pp],phi)); // fitting as negative binomial with dispersion determined by mean and a fitted dispersion factor
      }
  }



  
}
generated quantities {  // bits in here dont influence the running but stan will add them to stuff that is outputted
real additProtect[(D_runType>0) ? D_NAgeGroups : 0];
real relTransC;
real relTransU;
real relTransCvU_S;
real relTransCvU_M;
int epi_pred[(D_runType==1) ? (D_NDays) : 0]; 
int tmpEpiTot;

relTransC=transSC/transMC;
relTransU=transSU/transMU;
relTransCvU_S=transSC/transSU;
relTransCvU_M=transMC/transMU;

if (D_runType>0){
  for (ss in 1:D_NAgeGroups) 
    additProtect[ss]=deviationFromExp[ss]/(coefExp[1]*exp(ageExp[1]*D_meanAges[ss])/exp(ageExp[1]*D_meanAges[8]));
}

if (D_runType==1) {
  
  for (kk in 1:((D_NDays))) 
  { 
    tmpEpiTot=0;
    for (jj in 1:D_NAgeGroups) 
    {
      if (NBintensity[jj,kk]>0) {
        tmpEpiTot=tmpEpiTot + neg_binomial_2_rng(NBintensity[jj,kk],pow(NBintensity[jj,kk],phi));
      }else{
        tmpEpiTot=tmpEpiTot + 0;
      }
    }
    epi_pred[kk] = tmpEpiTot;
  }
  
}

}


