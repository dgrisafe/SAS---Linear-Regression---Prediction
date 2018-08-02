*==================================================================*
Program: Linear Regression Exploratory Data Analysis
Purpose: Prediction modeling data analysis
Date: 06/15/18
Note1: if comments are in ALL CAPS it is a place to manually enter code
*==================================================================*;

/* Fresh start when running SAS code */

	*clear log and output;
	dm 'log;clear;output;clear;';

	*clear the Results window within your program;
	dm 'odsresults; clear';

	*avoid getting too many labels error;
	ods graphics on / LABELMAX=2000;


*==================================================================*
Initialize dataset and look at all variables
*==================================================================*;

libname datset 'C:\Users\Dom\Documents\Screening Linear\datasets';
libname home 'C:\Users\Dom\Documents\Screening Linear';

*initialize dataset;
data stuff;
	*PATH TO APPROPRIATE DATASET;
	set datset.chol_313;
run;

*look at for all variables dataset;
proc contents;
run;


*<<<<<<<<<<<<<<<==================================================================*
 <<<<<<<<<<<<<<<	DEFINE GLOBAL MACROS:
 <<<<<<<<<<<<<<<		id variables
 <<<<<<<<<<<<<<<		outcome
 <<<<<<<<<<<<<<<		exposure of interest
 <<<<<<<<<<<<<<<		continuous covariates
 <<<<<<<<<<<<<<<		categorical covariates
 <<<<<<<<<<<<<<<		all variables
*<<<<<<<<<<<<<<<==================================================================*;

%let idvar		= ID;
%let outco		= SBP;
*all possible predictors (x variables);
%let contpred	= CHOL AGE BMI DBP HDL HT LDL SKIN TG VLDL WT;
%let catpred	= FEMALE PROBAND;
*all variables, except id variables;
%let allvar 	= &outco &contpred &catpred;


*==================================================================*
Preliminary data formatting:

	Define formats

	Format and label variables

	Remove unrealistic observations:
		Age range 0 to 122 - Guiness world record
		Weight max 635 kg = 1400 lbs - Guiness world record
		Highest recorded blood pressure in young adults preforming heavy, dynamic weight lifting 370/360 mmHg/mmHg - Narloch (1995) et. al.;
*==================================================================*;
title;

*DEFINE CODED FORMATS;
proc format;
	value ynf 0='No' 1='Yes';
	value FMLF 0='MALE' 1='FEMALE';
run;

data stuff;
	set stuff;

	*ASSIGN CODED FORMATS;
	format proband ynf. female FMLF.;

	*ASSIGN VARIABLE NAMES;
	label CHOL = 'Total cholesterol (dg/mL)';
	label AGE = 'Age (years)';
	label BMI = 'Body mass index (kg/m^2)';
	label DBP = 'Diastolic blood pressure (mmHg)';
	label HDL = 'High density lipoprotein (mg/dL)';
	label HT = 'Height (in)';
	label LDL = 'Low density lipoprotein (mg/dL)';
	label SKIN = 'Blood glucose levels (mg/dL)';
	label TG = 'Triglycerides (mg/dL)';
	label VLDL = 'Very low density lipoprotein (mg/dL)';
	label WT = 'Weight (lbs)';
	Label FEMALE = 'Gender (Female)';
	label PROBAND = 'Proband?';
	label SBP = 'Systolic blood pressure (mmHg)';

	*REMOVE UNREALISTIC DATA OBSERVATIONS;
	
		*Guiness book world record oldest person 122 yo;
		if AGE < 0 or AGE > 122 then AGE = .;

		*Guiness book world record heaviest person 635 kg = 1400 lbs;
		if WT > 1400 then WT = .;

		*Highest recorded blood pressure in young adults preforming heavy, dynamic weight lifting: 370/360 Narloch (1995) et. al.;
		if SBP > 370 then SBP = .;
		if DBP > 360 then DBP = .;

	*make temporary id to use in merging proc means output;
	tempid=1;
	
run;
proc print data=stuff (obs=10);
	title "Check to make sure dataset looks okay";
	title2 "Before splitting into Training (~80%) and Holding (~20%) sets)";
run;


*==================================================================*
Check for and delete duplicates
*==================================================================*;
title;

proc sort data = stuff nodupkey dupout = dups;
	*identify variables to check for duplicates;
 	by &idvar &allvar;
run;
proc print data=dups (obs=20);
	title 'First 20 duplicate observations';
run;


*====================================================================================================================================*
Make training and holdout datasets
	80% for training set
	20% for validation set
*====================================================================================================================================*;
title;

*make sure data is sorted by id variable;
data stuff;
	set stuff;
proc sort data=stuff;
	by &idvar;
run;

*choose ~80% of dataset for training set;
data stuffTrain;
	set stuff; if _n_ <=150; *make training set of first 150 (~80%) observations;
proc print data=stuffTrain (obs=10);
	title "Traning set preview";

*choose ~20% of dataset for holdout set;
data stuffHoldout;
	set stuff; if _n_ >150; *make training set of remaining (~20%) observations;
proc print data=stuffHoldout (obs=10);
	title "Holdout set preview";
run;

*summary statistics of training and holdout sets;
proc means data=stuffTrain;
	title "Training set summary statistics";
proc means data=stuffHoldout;
	title "Holdout set summary statistics";
run;


*====================================================================================================================================*
*====================================================================================================================================*
		EXPLORATORY DATA ANALYSIS:
*====================================================================================================================================*
*====================================================================================================================================*;

*==================================================================*
Table 1

Continuous outcome and covariates normality assessment & summary statistics:
	n
	means
	standard deviations
	ranges

Categorical covariates:
	n
	frequencies
	percentages
	cumulative frequencies
	number missing values
*==================================================================*;
title;

*view means and summary statistics of Continuous outcome, exposure, and covariates;
proc means data=stuffTrain maxdec=2;
	var &outco &contpred;
	*class female; *good if there's a group to compare by;
	title 'Continuous outcome, expos and covariates means, std dev, ranges';
	title2 'Use for Table 1';
run;

*view frequencies of Categorical varibles;
proc freq data=stuffTrain;
	tables &catpred;
	title 'Categorical variables frequencies';
	title2 'Use for Table 1';
run;


*==================================================================*
Continuous outcome and covariates normality assessment & summary statistics:
	histogram, sideways
	boxplot
	QQ-plot
	histogram, proper orientation

Categorical covariates:
	boxplots
*==================================================================*;
title;


*ASSESS LINEARITY of outcome to see if need to transform;
proc univariate data=stuffTrain normal plots;
	var &outco;
	histogram;
	title "&outco (outcome) need to be transformed to satisfy linearity?";
	title2 "Look at extreme 5 observations to assess for impossible values";
	title3 "Normality? Does mean = median = mode?";
run;


*ASSESS NORMALITY of outcome:
	*histograms of continuous outcome variable;
	*tests and plots to assess nomality;
%macro univar(var);
	proc univariate data=stuffTrain normal plots;
		var &var;
		histogram;
		title "&var univariate analysis";
		title2 "Look at extreme 5 observations to assess for impossible values";
		title3 "Independent variables do NOT have to satisfy Normality (mean = median = mode)";
	run;
%mend univar;

*Assess potential predictor variables as well;
%univar(AGE);
%univar(BMI);
%univar(DBP);
%univar(HDL);
%univar(HT);
%univar(LDL);
%univar(SKIN);
%univar(TG);
%univar(VLDL);
%univar(WT);
%univar(FEMALE);
%univar(PROBAND);

*boxplots of Categorical variables;
%macro boxplot(outcount, catpred);
	proc sgplot data=stuffTrain;
		vbox &outcount / category=&catpred;
		title "Boxplots of &outcount by &catpred (Categorical variable)";
	proc sgplot;
		histogram &outcount / group=&catpred transparency=0.75;
		density &outcount / type=kernel group=&catpred;
		title "Comparative histograms of &outcount by &catpred (Categorical variable)";
	run;
%mend boxplot;
%boxplot(&outco,female)
%boxplot(&outco,proband)


*====================================================================================================================================*
*====================================================================================================================================*
		UNIVARIATE ANALYSIS:
*====================================================================================================================================*
*====================================================================================================================================*;

*==================================================================*
Center continuous variables on their means in Training set
*==================================================================*;
title;

*make means of all variables;
proc means data = stuffTrain mean std n;
	by tempid; *all subjects have tempid=1;
	var &allvar;
	output out = stuffmeans mean= / autoname;
	title "Means of all variables (in Training set) for centering";
run;

*center all continuous independent variables on their means;
data stuffTrain;

	*merge means with dataset;
	merge stuffTrain stuffmeans;
	by tempid;

	*CENTER VALUES;
	CHOL_Cent 		= CHOL 		- CHOL_mean;
	AGE_Cent 		= AGE 		- AGE_mean;
	BMI_Cent 		= BMI 		- BMI_mean;
	DBP_Cent 		= DBP 		- DBP_mean;
	HDL_Cent 		= HDL 		- HDL_mean;
	HT_Cent 		= HT		- HT_mean;
	LDL_Cent 		= LDL 		- LDL_mean;
	SKIN_Cent 		= SKIN 		- SKIN_mean;
	TG_Cent 		= TG 		- TG_mean;
	VLDL_Cent 		= VLDL 		- VLDL_mean;
	WT_Cent 		= WT 		- WT_mean;

	*ASSIGN CENTERED VARIABLE NAMES;
	label CHOL_Cent = 'Total cholesterol (dg/mL), centered';
	label AGE_Cent = 'Age (years), centered';
	label BMI_Cent = 'Body mass index (kg/m^2), centered';
	label DBP_Cent = 'Diastolic blood pressure (mmHg), centered';
	label HDL_Cent = 'High density lipoprotein (mg/dL), centered';
	label HT_Cent = 'Height (in), centered';
	label LDL_Cent = 'Low density lipoprotein (mg/dL), centered';
	label SKIN_Cent = 'Blood glucose levels (mg/dL, centered)';
	label TG_Cent = 'Triglycerides (mg/dL), centered';
	label VLDL_Cent = 'Very low density lipoprotein (mg/dL), centered';
	label WT_Cent = 'Weight (lbs), centered';

	*TRANSFORMATIONS OF OUTCOME - IF NEEDED AFTER LOOKING AT:
		1) SCATTER OF OUTCOME VS PRIMARY EXPOSURE
		2) RESID VS. PRED PLOT AFTER FINAL MODEL);

		*Tukey's Ladder - outcome^lambda;
		*outco_log = log(&outco); 	*lambda = 0;
		*outco_sqrt = sqrt(&outco); *lambda = 1/2 (sqrt);
		*outco_2 = &outco*&outco;	*lambda = 2;

		*Note: In ESTIMATION, only consider transformations of the outcome after looking at the scatter of outcome vs main exposure 
		(not covariates). Covariates are not transformed in modeling!;

*check to make sure dataset correct;
proc print data=stuffTrain (obs=10);
	title "Check to make sure Training dataset looks okay";
	title2 "After centering on means of variables in Training set";
run;

*look at univariate analysis of transformed variable;
*%univar(&outco);
*%univar(outco_log);
*%univar(outco_sqrt);
*%univar(outco_2);


*==================================================================*
Center continuous variables on (training set) means in Holdout set
*==================================================================*;
title;

*center on exactly the same means for the HOLDOUT DATASET so it is comprable;
data stuffHoldout;

	*merge training set means with holdout dataset;
	merge stuffHoldout stuffmeans;
	by tempid;

	*CENTER VALUES;
	CHOL_Cent 		= CHOL 		- CHOL_mean;
	AGE_Cent 		= AGE 		- AGE_mean;
	BMI_Cent 		= BMI 		- BMI_mean;
	DBP_Cent 		= DBP 		- DBP_mean;
	HDL_Cent 		= HDL 		- HDL_mean;
	HT_Cent 		= HT		- HT_mean;
	LDL_Cent 		= LDL 		- LDL_mean;
	SKIN_Cent 		= SKIN 		- SKIN_mean;
	TG_Cent 		= TG 		- TG_mean;
	VLDL_Cent 		= VLDL 		- VLDL_mean;
	WT_Cent 		= WT 		- WT_mean;

	*ASSIGN CENTERED VARIABLE NAMES;
	label CHOL_Cent = 'Total cholesterol (dg/mL), centered';
	label AGE_Cent = 'Age (years), centered';
	label BMI_Cent = 'Body mass index (kg/m^2), centered';
	label DBP_Cent = 'Diastolic blood pressure (mmHg), centered';
	label HDL_Cent = 'High density lipoprotein (mg/dL), centered';
	label HT_Cent = 'Height (in), centered';
	label LDL_Cent = 'Low density lipoprotein (mg/dL), centered';
	label SKIN_Cent = 'Blood glucose levels (mg/dL, centered)';
	label TG_Cent = 'Triglycerides (mg/dL), centered';
	label VLDL_Cent = 'Very low density lipoprotein (mg/dL), centered';
	label WT_Cent = 'Weight (lbs), centered';

	*REPEAT TRANSFORMATIONS OF OUTCOME IN HOLDOUT SET - IF YOU DID THIS IN THE TRAINING SET:
		1) SCATTER OF OUTCOME VS PRIMARY EXPOSURE
		2) RESID VS. PRED PLOT AFTER FINAL MODEL);

		*Tukey's Ladder - outcome^lambda;
		*outco_log = log(&outco); 	*lambda = 0;
		*outco_sqrt = sqrt(&outco); *lambda = 1/2 (sqrt);
		*outco_2 = &outco*&outco;	*lambda = 2;

		*Note: In PREDICTION, only consider transformations of the outcome 
		PREDICTORS (X vars) are not USUALLY transformed in PREDICTION modeling!;

*check to make sure dataset correct;
proc print data=stuffHoldout (obs=10);
	title "Check to make sure Holdout dataset looks okay";
	title2 "After centering on means of variables in Training set";
run;


*<<<<<<<<<<<<<<<==================================================================*
 <<<<<<<<<<<<<<<	DEFINE GLOBAL MACROS FOR CENTERED VARIABLES:
 <<<<<<<<<<<<<<<		outcome, NOT centered BUT maybe transformed
 <<<<<<<<<<<<<<<		potential continuous predictors, centered
 <<<<<<<<<<<<<<<		potential categorical predictors, NOT centered
 <<<<<<<<<<<<<<<		all variables, centered, plus (transformed?) outcome
*<<<<<<<<<<<<<<<==================================================================*;

%let outTrans	= &OUTCO; *change if transforming outcome above;

*ALL CONTINUOUS VARIABLES AFTER CENTERING;
%let contCent	= CHOL_Cent AGE_Cent BMI_Cent DBP_Cent HDL_Cent HT_Cent LDL_Cent SKIN_Cent TG_Cent VLDL_Cent WT_Cent;

*all variables (w/centered) except id variables;
%let allCent	= &outTrans &contcent &catpred;


*==================================================================*
Check centering of continuous variables on (training set) means worked
	in Training set
	and in Holdout set
*==================================================================*;
title;

*check to make sure centering worked in Training set;
proc means data = stuffTrain mean std n;
	var &allCent;
	title "Check centering of variables in Training set worked";
	title2 "Means of continuous variables should be zero";
run;

*check to make sure centering worked on holdout set;
proc means data = stuffHoldout mean std n;
	var &allCent;
	title "Check centering of variables in Holdout set worked";
	title2 "Means of continuous variables should not necessarily be zero in Holdout set,";
	title3 "bc continuous variables were centered on (similar but different) means of Training set variables";
run;


*==================================================================*
Scatter plots of univariate relationships, w/macro
	*just in case need to add a higher order term of any predictor variable vs the outcome
*==================================================================*;
title;

*scatter plots macro;
%macro scatthat(x,y,dataset);
	*simple scatter plot;
	proc sgplot data=&dataset;
		scatter x=&x y=&y / datalabel = &idvar;
		reg x=&x y=&y;
		title "Scatter plot of &y on &x";
		title2 "In &dataset dataset";
	run;
%mend scatthat;
%scatthat(CHOL_Cent,&outTrans,stuffTrain)
%scatthat(AGE_Cent,&outTrans,stuffTrain)
%scatthat(BMI_Cent,&outTrans,stuffTrain)
%scatthat(DBP_Cent,&outTrans,stuffTrain)
%scatthat(HDL_Cent,&outTrans,stuffTrain)
%scatthat(HT_Cent,&outTrans,stuffTrain)
%scatthat(LDL_Cent,&outTrans,stuffTrain)
%scatthat(SKIN_Cent,&outTrans,stuffTrain)
%scatthat(TG_Cent,&outTrans,stuffTrain)
%scatthat(VLDL_Cent,&outTrans,stuffTrain)
%scatthat(WT_Cent,&outTrans,stuffTrain)
%scatthat(FEMALE,&outTrans,stuffTrain)
%scatthat(PROBAND,&outTrans,stuffTrain)

*IF SEE A PREDICTOR VARIABLE THAT MAY BE BETTER ILLUSTRATED AS HIGHER ORDER TERM, ADD IT TO THE &contCent MACRO ABOVE;
*DONT FORGET: IF HIGHER ORDER TERM IS SELECTED LATER, FIRST ORDER TERM MUST BE INCLUDED AS WELL

*==================================================================*
Correlations between all independent variables and continuous outcome
	COLLINEARITY CAN LEAD TO CONFUSING RESULTS WHEN USING MODEL SELECTION TECHNIQUES
	YOU ALREADY CENTERED CONTINUOUS X VARIABLES, WHICH REDUCES COLLINIARITY AMONG HIGHER ORDER TERMS
*==================================================================*;
title;

*Pearsons correlations of all variables;
%macro corrthat(vars);
	proc corr data=stuffTrain best = 6; *shows the top 6 correlations for each variable;
	  var &vars;
	  title "Pearson correlations: Percent variability shared = R^2 * 100%";
	  title2 "> 0.60-0.70 b/w exposure or covariates suggests collinearity";
	  *By squaring the correlation and then multiplying by 100%, 
	  you can determine what percentage of the variability is shared.  
	  Let’s round 0.59678 to be 0.6, which when squared would be .36, multiplied by 100 would be 36%.  
	  Hence read shares about 36% of its variability with write.;
	run;
%mend corrthat;
%corrthat(&allCent)

*EXPLAIN IF ANY VARIABLES WILL BE REMOVED DUE TO HIGH CORRELATION (I.E. THEY ARE COLLINEAR)
	REMOVE WEIGHT BECAUSE IT IS COLLINEAR WITH BMI (R = 0.866)
	REMOVE LDL BECAUSE IT IS COLLINEAR WITH CHOLESTEROL (R = 0.971)
	REMOVE TG BECAUSE IT IS COLLINEAR WITH VLDL (r = 0.906)


*<<<<<<<<<<<<<<<==================================================================*
 <<<<<<<<<<<<<<<	DEFINE GLOBAL MACROS FOR ALL POTENTIAL PREDICTORS - MAXIMUM MODEL:
 <<<<<<<<<<<<<<<		outcome, NOT changed
 <<<<<<<<<<<<<<<		potential continuous predictors, centered
 <<<<<<<<<<<<<<<			(After dropping continuous covariates with collinearity)
 <<<<<<<<<<<<<<<		potential categorical predictors, NO change
 <<<<<<<<<<<<<<<		select all potential predictor variables to start, including higher order terms
 <<<<<<<<<<<<<<<		all variables, centered, plus outcome
*<<<<<<<<<<<<<<<==================================================================*;

*DECIDE COMPLETE LIST OF X VARIABLES TO CONSIDER AS SET OF POTENTIAL PREDICTORS OF Y;
	*maximum # X variables: 	k = n - 2;
	*reasonable # X variables: 	k = n/5;
	*MAY WANT TO INCLUDE ADDITIONAL VARIABLES: BMI, AGE^2, HEIGHT^2, BMI^2
		*USE A PRIORI KNOWLEDGE TO DECIDE ON HIGHER ORDER TERMS;
		*ALSO LOOK AT SCATTER PLOTS OF OUTCOME ON EACH POTENTIAL PREDICTOR VARIABLE (next section);
%let contMax	= CHOL_Cent AGE_Cent BMI_Cent DBP_Cent HDL_Cent HT_Cent SKIN_Cent VLDL_Cent WT_Cent LDL_Cent TG_Cent;

*all variables (w/centered) except id variables;
%let allMax	= &outTrans &contMax &catpred;

*==================================================================*
Step-wise model selection - NO interaction terms, at first
	Put in all the different selection approachs commented out just in case
		"stepwise was/was not used, similar results observed in forwards and backwards

	Selecting "Best" model: 
		R-squared, adjusted - maximize while having most parimonious model possible (fewest x variables)
		F-test - does model 1 provide a statistically significant improvement over model 2? (extra sums of squares p. 5-25)
		Mallow's Cp - smallest Cp usually best model
			Cp = k + 1 (for max model)
			for smaller models Cp usually larger
*==================================================================*;
title;

/*
*Best subsets elimination selection;
proc reg data=stuffTrain;
*need to make interaction terms and dummy variables for nominal variabels with >1 outcome;
	model &outTrans = &contMax &catpred / selection=rsquare cp; *add &highPred if necessary;
	title "Best subsets elimination for &outco";
	title2 "Best model has lowest Cp values";
	title3 "Beyond the chosen model, the R-squared is changing little";		
run;
quit;
	proc glm data=stuffTrain;
		*MANUALLY ENTER the best model from the best subsets selection method above (w/ lowest Cp);
		model &outTrans =  ;
		title "Best subsets chosen model for &outco";	
	run;
quit;

title;
*Forward elimination selection;
proc glmselect data=stuffTrain;
	*class = XXX; *if theres any nominal or ordinal variables with more than 2 possible outcomes;
	model &outTrans = &contMax &catpred / selection=forward 
	details=all hierarchy=single sle=0.15 showpvalues select=sl stop=sl; *add &highPred if necessary;
	title "Forward elimination for &outco";	
run;
quit;

title;
*Backwards elimination selection;
proc glmselect data=stuffTrain;
	*class = XXX; *if theres any nominal or ordinal variables with more than 2 possible outcomes;
	model &outTrans = &contMax &catpred / selection=backward 
	details=all hierarchy=single sls=0.15 showpvalues select=sl stop=sl; *add &highPred if necessary;
	title "Backward elimination for &outco";	
run;
quit;
*/

title;
*Stepwise elimination selection;
proc glmselect data=stuffTrain;
	*class = XXX; *if theres any nominal or ordinal variables with more than 2 possible outcomes;
	model &outTrans = &contMax &catpred / selection=stepwise 
	details=all hierarchy=single sls=0.15 sle=0.15 showpvalues select=sl stop=sl; *add &highPred if necessary;
	title "Stepwise elimination for &outco";
run;
quit;

*<<<<<<<<<<<<<<<	DEFINE GLOBAL MACRO FOR ALL RELEVANT FIRST AND HIGHER ORDER VARIABLES IN MODEL (NO INTERACTION TERMS YET THO);
%let firstPred = AGE_CENT BMI_CENT DBP_CENT; 
*DONT FORGET TO ADD RELEVANT FIRST ORDER TERM IF HIGHER ORDER TERM SELECTED;


*NOTES - MODEL SELECTION:
	1) In Best Subsets (proc reg), if only some of the indicator variables are included for nominal/ordinal variables,
		must include or exculde all dummy variables (judgement call)
	2) If a higher order term is selected for the model, the lower order term must also be included
	3) If colliniearity issues arise due to above rules (e.g. BMI and weight in model), 
		decide which to include and which to exclude
	4) Want to force 2 variables to remain in the model? 
		Put them first and second in the variable list and say "include=2" in the model options (p.10-32);


*==================================================================*
Step-wise model selection -  INTERACTION terms of all SIGNIFICANT PREDICTORS
	Put in all the different selection approachs commented out just in case
	"stepwise was/was not used, similar results observed in forwards and backwards

	Selecting "Best" model: 
		R-squared, adjusted - maximize while having most parimonious model possible (fewest x variables)
		F-test - does model 1 provide a statistically significant improvement over model 2? (extra sums of squares p. 5-25)
		Mallow's Cp - smallest Cp usually best model
			Cp = k + 1 (for max model)
			for smaller models Cp usually larger
*==================================================================*;
title;

*<<<<<<<<<<<<<<<	DEFINE GLOBAL MACRO FOR ALL POSSIBLE INTERACTIONS AMONG RELEVANT FIRST AND HIGHER ORDER, NON-INTERACTION TERMS;
*<<<<<<<<<<<<<<<	Maximum # possible interaction terms = N(N-1)/2;
%let ixnMax = 
	AGE_CENT*DBP_CENT AGE_CENT*BMI_CENT
	BMI_CENT*DBP_CENT
;

*Stepwise elimination selection;
proc glmselect data=stuffTrain;
	*class = XXX; *if theres any nominal or ordinal variables with more than 2 possible outcomes;
	model &outTrans = &firstPred &ixnMax / selection=stepwise 
	details=all hierarchy=single sls=0.15 sle=0.15 showpvalues select=sl stop=sl;
	title "Assessing INTERACTION terms w/Stepwise elimination for &outco";
	title2 "Possible interaction terms = n(n-1)/2";
	title3 "E.g. if n = 5, there are 5(5-1)/2 = 10 possible interaction terms";
	score data=stuffTrain out=scoreTrain p=predTrain r=residTrain; *training output datasets w/ predicted outcome and residuals;
	score data=stuffHoldout out=scoreHold p=predHold r=residHold; *holdout output data sets w/ predicted outcome and residuals;
run;
quit;

proc print data=scoreTrain (obs=10);
	title "Preliminary: Training set predictions";
	title2 "Model built from Training set used to make Predicted values from Training set predictors";
proc print data=scoreHold (obs=10);
	title "Preliminary: Holding set predictions";
	title2 "Model built from Training set used to make Predicted values from Holding set predictors";
run;

*Calculating R-squared for Preliminary final model on Training data;
proc corr data=scoreTrain;
	var &outTrans predTrain;
	title "Training: Pearson's correlation b/w predicted and observed &outTrans in Training set";
	title2 "Square Pearsons to get R-squared";
	title3 "These R-squareds are preliminary, and should be reported after assessing LINE assumptions and Sensitivity analysis";
run;

*Calculating R-squared for Preliminary final model on Holdout data;
proc corr data=scoreHold;
	var &outTrans predHold;
	title "Holding: Pearson's correlation b/w predicted and observed &outTrans in Holding set";
	title2 "Square Pearsons to get R-squared";
	title3 "These R-squareds are preliminary, and should be reported after assessing LINE assumptions and Sensitivity analysis";
run;


*<<<<<<<<<<<<<<<==================================================================*
 <<<<<<<<<<<<<<<	DEFINE GLOBAL MACROS FOR VARIABLES IN PRELIMINARY FINAL MODEL:
 <<<<<<<<<<<<<<<		outcome, no change
 <<<<<<<<<<<<<<<		first order predictor variables, no change
 <<<<<<<<<<<<<<<		select all potential predictor variables to start, including higher order terms
 <<<<<<<<<<<<<<<		all variables, centered, plus outcome
*<<<<<<<<<<<<<<<==================================================================*;

*all significant interaction terms from final stepwise selection;
%let ixnFinal	= AGE_Cent*DBP_Cent;

*all variables in preliminary final model (w/centered) except id variables;
%let preFinMod	= &outTrans &firstPred &ixnFinal;


*==================================================================*
Preliminary Final Model
	Assess LINE assumptions
	Model Diagnostics: Tolerance, Cooks D, Leverage
*==================================================================*;
title;

*macro to produce the Prelminary Main Effects model;
%macro prelimFinalModelReg(out,var,ixn,ident);
	proc glm data=stuffTrain plot=diagnostics(label);
		id &ident;
		model &out=&var &ixn / clparm tolerance;
		title "Preliminary final model of &out on &var,";
		title2 "w/interaction terms &ixnFinal";
		title3 "Tolerance < 0.1 indicates substantial collinearity";
		output out=prefinalmodel p=Predicted r=Residual student=StudentResid cookd=CooksD 
				h=Leverage rstudent=JackknifeResid dffits=DFFITs covratio=CovRatio;
	run;
	quit;
	proc sgplot data=prefinalmodel;
		scatter x=predicted y=residual;
		loess  x=predicted y=residual / clm smooth=0.4;
		refline 0;
		title "Assess Linearity:"; 
		title2 "Assess Homoscedasticity: Are Variances (of residuals) Equal?";
		title3 "Residuals vs. Predicted &out to look for pattern (Loess) or unequal variance";
	proc univariate data=prefinalmodel normal plots;
		var residual;
		title "Are Residuals Normally Distributed?";
		title2 "Assess Normality of &out on &var,";
		title3 "w/interaction terms &ixnFinal";
	run;
%mend prelimFinalModelReg;

*Preliminary Final Model that contains all interaction terms that maintain significance;
%prelimFinalModelReg(&outTrans,&firstPred,&ixnFinal,&idvar)

*WRITE OUT INTERPRETATIONS OF LINE ASSUMPTIONS AND TOLERANCE (COLLINEARITY) DIAGNOSTIC
	1) LINEARITY IS NOT SATISFIED, LOOKING AT ONE POINT TOWARDS THE EXTREME PREDICTED VALUE
	2) INDEPENDENCE IS ASSUMED BECAUSE EACH SUBJECT IS THEIR OWN. THIS MAY BE AN INAPPROPRIATE ASSUMPTION BECAUSE MANY SUBJECTS ARE RELATED
	3) NORMAILITY IS ASSUMED BY LOOKING AT THE GRAPHS OF THE RESIDUALS AND THE QQ PLOT
	4) EQUAL VARIANCES/HOMOSCEDASTICITY IS ASSUMED BECAUSE THE RESIDUALS SEEM TO BE SPREAD RELATIVELY EVENLY THROUGHOUT PREDICTED VALUES
		EXCEPT FOR THE ONE OUTLIER THAT LOOKS LIKE ITS RESPONSIBLE FOR THE MODEL VIOLATING LINEARITY

	5) NO TYPE 1 TOLERANCE VALUES ARE BELOW 0.1 FOR ANY VARIABLES, SUGGESTING NO STRONG CORRELATION BW ANY PREDICTORS;
;

*look graphically for outliers;
	*diagnostics;
		%scatthat(leverage,predicted,prefinalmodel)
		%scatthat(cooksd,predicted,prefinalmodel)
		%scatthat(DFFITS,predicted,prefinalmodel)
	*scatter plot of residuals vs predicted;
		%scatthat(residual,predicted,prefinalmodel)
		%scatthat(studentResid,predicted,prefinalmodel)
		%scatthat(jackknifeResid,predicted,prefinalmodel)

*assess most extreme diagnostic values;
%macro extremeDiag(diag);
	data diagnosticPFM;
		set prefinalmodel;
		keep id &outTrans predicted residual &diag;
	proc sort data=diagnosticPFM;
		by descending &diag;
	proc print data=diagnosticPFM(obs = 25);
		title "&diag: 25 most extreme observations";
		title2 "n = # observations ; k = # model parameters";
		title3 "Leverage > 2(k+1)/n	; CooksD > 1 ; DFFITs > 2*sqrt(k/n)";
run;
%mend;
*diagnostics;
	%extremeDiag(leverage)
	%extremeDiag(cooksd)
	%extremeDiag(DFFITS)
*residuals;
*	%extremeDiag(residual);
*	%extremeDiag(studentResid);
*	%extremeDiag(jackknifeResid);

*extreme predicted values;
%extremeDiag(predicted)

*look for collinearity (again) in the variables in the preliminary final model;
%corrthat(&outTrans &firstPred)


*==================================================================*
Sensitivity analysis
*==================================================================*;
title;

*macro to produce the Sensitivity analysis;
*PUT THE ID NUMBERS OF THE OUTLIERS IN THE GREEN VALUES FOR 42 AND 115;
%macro SensitivityAnalysis(out,var,ixn,ident);
	proc glm data=stuffTrain(where=(&ident ne 42)) plot=diagnostics(label);
		id &ident;
		model &out=&var &ixn / clparm;
		title "Sensitivity analysis of Final Model (&out on &var), w/adjusting vars and ixns";
		title2 "Eliminating potentially influential points";
		output out=SensitivityAnalysis p=Predicted r=Residual student=StudentResid cookd=CooksD 
				h=Leverage rstudent=JackknifeResid dffits=DFFITs covratio=CovRatio;
	run;
	quit;
	proc sgplot data=SensitivityAnalysis;
		scatter x=predicted y=residual;
		loess  x=predicted y=residual / clm smooth=0.4;
		refline 0;
		title "Assess Linearity:"; 
		title2 "Assess Homoscedasticity: Are Variances (of residuals) Equal?";
		title3 "Residuals vs. Predicted &out to look for pattern (Loess) or unequal variance";
	proc univariate data=SensitivityAnalysis normal plots;
		var residual;
		title "Are Residuals Normally Distributed?";
		title2 "Assess Normality of &out on &var,";
		title3 "w/interaction terms &ixn";
	run;
%mend SensitivityAnalysis;

*Sensitivity Analysis that contains all interaction terms that maintain significance;
%SensitivityAnalysis(&outTrans,&firstPred,&ixnFinal,&idvar)


*==================================================================*
Model validation on holdout dataset
*==================================================================*;
title;

*remove unrealistic outliers from observations;
data stuffTrain;
	set stuffTrain(where=(&idvar ^= 42));
run;

*re-fit final model from stuffTrain dataset;
proc glm data=stuffTrain;
	id &idvar;
	model &outTrans = &firstPred &ixnFinal / clparm;
	store finalModel;
	title "Validating Final model based on Training dataset";
run;

*get predicted values of final model on training dataset;
proc plm restore=finalModel;
	score data=stuffTrain out=scoreTrainFin Predicted=predTrainFin;
run;
*Calculating R-squared for Final model on Training data;
proc corr data=scoreTrainFin;
	var &outTrans predTrainFin;
	title "Training: Pearson's correlation b/w predicted and observed &outTrans in Training set";
	title2 "Square Pearsons to get R-squared";
	title3 "These are final Pearsons R, need to be squared to be reported";
run;

*get predicted values of final model on holdout dataset;
proc plm restore=finalModel;
	score data=stuffHoldout out=scoreHoldFin Predicted=predHoldFin;
run;
*Calculating R-squared for Preliminary final model on Holdout data;
proc corr data=scoreHoldFin;
	var &outTrans predHoldFin;
	title "Holding: Pearson's correlation b/w predicted and observed &outTrans in Holding set";
	title2 "Square Pearsons to get R-squared";
	title3 "These are final Pearsons R, need to be squared to be reported";
run;
