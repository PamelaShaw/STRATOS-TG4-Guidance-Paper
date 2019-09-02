title1 "OPEN Study: Regressing BMI on Log Potassium Density, Age and Gender";
title2 "Maximum Likelihood Estimation using Proc Calis";

libname lib "h:\documents\freedman\stratos\data";

options center linesize=120 pagesize=42;

proc format;
  value gendfmt
    0 = "Male"
    1 = "Female";
  run;

 /* get open data */

data open;
  set lib.open;

 /* redefine gender: 0 = male, 1 = female */

  gender = gender - 1;

 /* log potassium density */

  ffq   = log(potass_ffq1 / energy_ffq1);
  biom1 = log(potass_ur1 / tee1);
  biom2 = log(potass_ur2 / tee1);

 /* for observations 251-484, set biomarker equal to missing */

  if (_n_ > 250) then do;
    biom1 = .;
    biom2 = .;
    end;
  run;

 /* calcluate means */

proc means data=open n mean std var min max;
  var bmi age gender ffq biom1-biom2;
  run;

 /* fit mle model using proc calis */

ods exclude modelstatement endogenousvar exogenousvar stdmanifesteq stdlatenteq sqmultcorr;

proc calis data=open method=fiml stderr nostand pshort outest=outest;
  nloptions omethod=newrap gconv=1e-12 fconv=1e-12 maxfunc=1000 maxiter=200 history;
  var bmi age gender ffq biom1-biom2;

  lineqs

 /* age */

  age = mean_age * intercept + d_age,

 /* gender */

  gender = mean_gender * intercept + d_gender,

 /* unobserved X */

  f_x = gamma0 * intercept + gamma1 * age + gamma2 * gender + d_phi,

 /* response bmi */

  bmi = beta0 * intercept + beta1 * f_x + beta2 * age + beta3 * gender + d_epsilon,

 /* ffq */

  ffq = alpha0 * intercept + alpha1 * f_x + alpha2 * age + alpha3 * gender + d_u,

 /* unbiased biomarker */

  biom1 = 0 * intercept + f_x + d_delta1,
  biom2 = 0 * intercept + f_x + d_delta2;

  variance
    d_age     = var_age,
    d_gender  = var_gender,
    d_phi     = var_phi,
    d_epsilon = var_epsilon,
    d_u       = var_u,
    d_delta1-d_delta2 = var_delta [...];

  cov
    d_age * d_gender = cov_age_gender;

  parameters beta0-beta3 alpha0-alpha3 gamma0-gamma2
             var_epsilon var_phi var_u
             mean_age mean_gender var_age var_gender cov_age_gender;

 /* ODS statements for output */

  ods output addparms=addparms;
  run;

ods exclude none;
