#include "def.h"

double solve_BendersPre(int restype)
{
	int i,j,k;
	int index,index1;  // auxiliar indices to fill in the constraint matrix
	int nodecount;
	//Variables to call cplex
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double    *rhs;    // right hand side of constraints ................................
	char      *sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int       *matbeg; // index of first non-zero element in each row...................
	int       *matind; // associated column of each non-zero element ...................
	double    *matval; // coefficient values for the non-zero elements of constraints....
	int       status;  // optimization status......................... .................
	char probname[16]; // problem name for cplex .......................................
	char Printinst[40];
	double    value;   // objevtive value of solution ..................................
	double *x;
	char *ctype;//The C type for converting into a MIP
	int addcuts, com_check,cutnz,num_cuts;
	double *cutval;
	int    *cutind;
	int		*optact;
	int		*feasact;
	int		Flg_stop=0;
	int iterations;
	int numlocrows;
	double LB;
	int eliminate;
	int exitflag;
	int	acum;
	int acum2;
	double *cut_sets;
	double    mincost;
    int    minindex;
	int		*candidate;
	int    acumcuts;
	int lastLPenter=0;
	int  firsttime=0;
	double checkint;
	int flagint;
	int intloopLP;//counter for the interior loop of the lift and project cuts.
	
	candidate=create_int_vector(M);
	cutind=create_int_vector(M+1);
	cutval=create_double_vector(M+1);
	flag_frac=1;
	First_time=1;
	flag_keepcut=1;
		
	//Initialize CPLEX envPreironment
	envPre = CPXopenCPLEX (&status); 
	if ( envPre == NULL ) {
		char  errmsg[1024];
		printf ("Could not open CPLEX. \n");
		CPXgeterrorstring (envPre, status, errmsg);
		printf ("%s", errmsg);
	}

    // Create the problem in CPLEX 
	strcpy(probname,"MPPre"); //Copy the name placed as the problem name
	lpPre = CPXcreateprob (envPre, &status, probname);  // Initialize the pointer that signals the problem 
	
	if ( envPre == NULL ) {
		char  errmsg[1024];
		printf ("Could not create lpPre. \n");
		CPXgeterrorstring (envPre, status, errmsg);
		printf ("%s", errmsg);
	}
	
	/************We entered into cplex the arrays with the design variable****************/
	status = CPXnewcols (envPre, lpPre, M, yobj, ylb, yub, NULL, NULL);
	if( status )  fprintf (stderr,"CPXnewcols failed.\n");
	/**********************************************************************************/
	/************We entered into cplex the arrays with the routing cost variables****************/
	status = CPXnewcols (envPre, lpPre, K, zobj, zlb, zub, NULL, NULL);
	if( status )  getchar();
	/**********************************************************************************/	
	
	/**********************************Here we can choose to add these*********************/
	/**********************************************************************************/
	numrows = 1;
	numnz = M;
	d_vector(&rhs,numrows,"open_cplex:2");
	c_vector(&sense,numrows,"open_cplex:3");
	i_vector(&matbeg,numrows,"open_cplex:4");
	i_vector(&matind,numnz,"open_cplex:6");
	d_vector(&matval,numnz,"open_cplex:7");

	/************First we write that there must be an arc leaving the origin of each commodity*******/
	for(k=0;k<K;k++){	
		index = 0;//columns
		index1 = 0;//rows 	
		sense[index1]='G';      
		matbeg[index1] = index;
		rhs[index1++]=1;
		for(j=0;j<M;j++){			 
			if(com[k].i==edges[j].i){ //if node is origin
					 matind[index] = pos_y[edges[j].i][edges[j].j];
					 matval[index] = 1;			 
					 index++;
			}
		}
		UpdFeas_pool(1,matval, matind, index);
		status = CPXaddrows (envPre, lpPre, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
		if( status ) fprintf (stderr,"CPXaddrows failed.\n");
	}
	/************Next we write that there must be an arc arriving at each destination*******/
	for(k=0;k<K;k++){	
		index = 0;//columns
		index1 = 0;//rows 	
		sense[index1]='G';      
		matbeg[index1] = index;
		rhs[index1++]=1;
		for(j=0;j<M;j++){			 
			if(com[k].j==edges[j].j){ //if node is origin
					 matind[index] = pos_y[edges[j].i][edges[j].j];
					 matval[index] = 1;			 
					 index++;
			}
		}
		UpdFeas_pool(1,matval, matind, index);
		status = CPXaddrows (envPre, lpPre, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
		if( status ) fprintf (stderr,"CPXaddrows failed.\n");
	}	
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);


	d_vector(&rhs,1,"open_cplex:2");
	c_vector(&sense,1,"open_cplex:3");
	i_vector(&matbeg,1,"open_cplex:4");
	numcols = CPXgetnumcols (envPre, lpPre);
	d_vector(&x,numcols,"open_cplex:0");	
	CPXchgobjsen(envPre, lpPre,CPX_MIN);
	reducedcost=create_double_vector(numcols);
	/*Changing some of the precision*/
	CPXsetintparam(envPre,CPX_PARAM_NUMERICALEMPHASIS,1); //Numerical precision of the time
	CPXsetintparam(envPre,CPX_PARAM_THREADS, 1);
	CPXsetintparam(envPre,CPX_PARAM_LPMETHOD,CPX_ALG_DUAL);//Use the Primal Algorithm
	iterations=0;
	//start_time = clock();
	/************************************************/
	do{
		/**************************Solving whatever cuts we have*************************/
		num_cuts=0;
		addcuts=0;
		violation=0.00001;
		cutsettol=0.9999; //Cutset inequality
		CPXlpopt(envPre,lpPre);  //solve the linear program
		status=CPXgetstat(envPre,lpPre);
		if(status==1){
			printf("Solution found with ");
			CPXgetobjval(envPre,lpPre,&value);
			LB=value;
			printf("LP Relax: %.2f\n",LB);			
			CPXgetx(envPre,lpPre,x,0, numcols-1);  // obtain the values of the decision variables
			for(j=0;j<M;j++){                  //Updating y_frac to obtain the cuts in the dual subproblem		 
				ybar_frac[edges[j].i][edges[j].j] = x[pos_y[edges[j].i][edges[j].j]];					
			}//printf("/******************************/\n\n");
		}
		else{
			//printf("Unknown stopping criterion (%d) \n",status);
			//getchar();
			Flg_stop=1;
		}
		/*************************Check for violated Feasibility Cuts*******************************************************/
		sense[0]='G';      
		matbeg[0] = 0;
		rhs[0]=1;
		if(1==1){
			for(i=0;i<num_origins;i++){
				for(j=0;j<origins[i].dim;j++){
					addcuts=maxFlowAndMinCut(origins[i].s, origins[i].d[j], ybar_frac, cutval,cutind);
					if(addcuts==1){
						//printf("Added a feasiblity cut\n");
						com_check=origins[i].ind[j][0];
						status = CPXaddrows (envPre, lpPre, 0, 1, feasnumnz, rhs, sense, matbeg, cutind, cutval, NULL, NULL);					
						UpdFeas_pool(2,cutval, cutind, feasnumnz);
						cuts_SC[2]++;
						num_cuts=num_cuts+addcuts;
					}
				}
			}
		}
		/*************************Check for violated Optimality Cuts*******************************************************/
		if ( num_cuts == 0) {
			solve_net_sub(x);         //MIN_COST flow M&W
			for(k=optcount-vioptcuts;k<optcount;k++){					//Add the violated cuts
				rhs[0]=optcuts[k].RHS;
				status = CPXaddrows (envPre, lpPre, 0, 1, optcuts[k].numnz, rhs, sense, matbeg, optcuts[k].ind, optcuts[k].coeff, NULL, NULL);
				cuts_SC[3]++;
			}
			if(vioptcuts==0){						//Try again except this time with the actual value
				//printf("Added an Exact Optimality cut\n");
				solve_net_exact(x);
				for(k=optcount-vioptcuts;k<optcount;k++){					//Add the violated cuts
					rhs[0]=optcuts[k].RHS;
					status = CPXaddrows (envPre, lpPre, 0, 1, optcuts[k].numnz, rhs, sense, matbeg, optcuts[k].ind, optcuts[k].coeff, NULL, NULL);
					cuts_SC[3]++;
				}
			}
			if(vioptcuts==0){
				RHSLPcutoff=0;
				for(j=0;j<M;j++){ //Keeping the LP solution to use as corepoint for the root node.
					LPsol[j]=x[j];
					RHSLPcutoff+=x[j];
				}
				if(firsttime==0){//Keeping the LP solution to be used as a corepoint	
					firsttime++;
					/***********Here we will run our heuristic here******************/ 
					//if(useheur==1 && method==2)master_p(ylb, indarcinpath);
					/****************************************************************/
				}
				else{
					split_var=frac_index(x);
					if(LifProunds<LifPlimit/* && lastLPenter<iterations-1*/ && split_var>=0){  //Lift and Project in preprocess
						lastLPenter=iterations;						
						ConstrMatrix_LP(x,10000);				//we allow for extra ones to appear
						for(intloopLP=0;intloopLP<split_var/depthLP;intloopLP++){
							if(LiftandP(splitvars[intloopLP], x)==1){
								status = CPXaddrows (envPre, lpPre, 0, 1,  LPcuts[LifPcount].numnz, &LPcuts[LifPcount].RHS, sense, matbeg, LPcuts[LifPcount].ind, LPcuts[LifPcount].coeff, NULL, NULL);
								printf("> %lf  with %d / %d non zeros\n", LPcuts[LifPcount].RHS,LPcuts[LifPcount].numnz, M+K);
								if( status )  fprintf (stderr,"CPXaddrows failed.\n");
								/********Adding now the lift and project cuts to the set of optimality cuts for it to be used again-- Added to the code on January 20th 2019********/
										optcuts[optcount].RHS =LPcuts[LifPcount].RHS;
										optcuts[optcount].numnz=LPcuts[LifPcount].numnz;
										optcuts[optcount].fromfrac =1;
										for(i=0;i<LPcuts[LifPcount].numnz;i++){
											optcuts[optcount].coeff[i]=LPcuts[LifPcount].coeff[i];
											optcuts[optcount].ind[i]=LPcuts[LifPcount].ind[i];
										}
										if(optcount+1==optlim*M*K-1){							//Now increasing optcount while checking there is stil space
											//getchar();
											optlim++;
											optcuts= (CUTPOOL *) realloc(optcuts,optlim*M*K*sizeof *optcuts); 
										}							
										optcount++;
										added_LP++;
										optcuts[optcount].coeff=create_double_vector(M+K+1);
										optcuts[optcount].ind=create_int_vector(M+K+1);
										free(LifPcut);
							}
							else{
								printf("No violated constraint found\n");
								Flg_stop=1;
							}
						}
						/*************************************************************************************************************/
						//printf("Added a Lift and project Cut\n");
						LifProunds++;
						TotalLPcuts++;
					}					
					else Flg_stop=1;
				}
			}
		}
		iterations++;
	}while(Flg_stop==0);
	free(matbeg);
	free(sense);
	free(rhs);
	LifPlimit=0;
	RHSLPcutoff=0;
	
	/***********Passing the LP solution and checking if integral******************/
	flagint=0;
	for(j=0;j<M;j++){ //Keeping the LP solution to use as corepoint for the root node.
		LPsol[j]=x[j];
		if(x[j]>0.5){
			checkint=1.0-x[j];	
		}
		else{
			checkint=x[j];	
		}
		if(ABS(checkint)<0.0001) flagint++;
		RHSLPcutoff+=x[j];	
	}
	if(flagint==M) bestobj=LB;
	end_time = clock();
	cputimeBD = (double)(end_time - start_time) / CLOCKS_PER_SEC;
	presolve_ind=0;
	if(100*(bestobj-LB)/ABS(LB)>.01){
		ConstrMatrix_LP(x,1.0); //If we haven't found the optimum then continue
		for(j=0;j<M;j++){
			stabilizer[j]=LPsol[j]+.01;
			core_point[edges[j].i][edges[j].j]=stabilizer[j];
			indprior[j]=j;
			priority[j]=(int) (x[j]*edges[j].f);
		}

		/********Now we will do variable elimination*********/
		//status = CPXgetdj (envPre, lpPre, reducedcost, 0, numcols-1); //Get the information of the reduced cost
		//CPXgetobjval(envPre,lpPre,&value);
		//for(j=0;j<M;j++){
		//	forelimintest[j]=value+reducedcost[pos_y[edges[j].i][edges[j].j]];
		//	//printf("Testing arc (%d,%d) for elimination %lf+%lf=%lf >? %lf\n", edges[j].i,edges[j].j, value,reducedcost[pos_y[edges[j].i][edges[j].j]], forelimintest[j], bestobj);
		//	if(forelimintest[j]>bestobj+.00001){
		//		yub[j]=0;
		//		/*printf("Fixed arc (%d,%d) to 0\n", edges[j].i,edges[j].j);
		//		getchar();*/
		//	}		  
		//}
		//getchar();
	/***************************************************************/
	}
	//getchar();
	/****Clean out the memory****/
	free(x);
	if ( lpPre != NULL ) {
			  status = CPXfreeprob (envPre, &lpPre);
			  if ( status ) {
				fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
			  }
	}
	if ( envPre != NULL ) {
			  status = CPXcloseCPLEX (&envPre);
			  if ( status ) {
				char  errmsg[1024];
				fprintf (stderr, "Could not close CPLEX environment.\n");
				CPXgeterrorstring (envPre, status, errmsg);
				fprintf (stderr, "%s", errmsg);
			  }
	}
	return LB;
}

int ConstrMatrix_LP( double *x, double Tol_Act){
	int i,j,k;
	int NbActive;
	int index,bigInd;
	double slack_1;
	/**Initializing so that there are no active cuts***/
	num_feasact=0;
	num_optact=0;
	/* Here we define the largest of all of them*/
	bigInd=MAX(feascount,optcount);
	index=0;
	do{
		/********************************First we'll check for active feasibility cuts**************************************There doesn't seem to be any active feasiblility cuts */
		if(index<feascount){
			slack_1=0;
			j=0;
			while(j<feascuts[index].numcol && (slack_1<=1.1)){
				slack_1+=x[feascuts[index].ind[j]];				
				j++;
			}
			if(j==feascuts[index].numcol && ABS(slack_1-1)<=Tol_Act){    //Checking if it's active
				feasact[num_feasact]=index;          
				num_feasact++;
			}
		}
		if(index<optcount){
			slack_1=-optcuts[index].RHS;
			j=0;
			while(j<optcuts[index].numnz && (slack_1<=0.1)){
				slack_1+=optcuts[index].coeff[j]*x[optcuts[index].ind[j]];
				j++;
			}
			if(j==optcuts[index].numnz && ABS(slack_1)<=Tol_Act){
				optact[num_optact]=index;
				num_optact++;
			}/*else{ printf("rejected optimality cut with slack %lf",slack_1);getchar();}*/
		}
		index++;
	}while(index<bigInd);
	NbActive=num_feasact+num_optact;
	//printf("%d feasibility cuts out of %d are active\n%d optimality cuts out of %d are active\nFor a total of %d active Benders cuts\n",num_feasact,feascount,num_optact,optcount,NbActive);
	return NbActive;
}

int frac_index(double *x){    //Function returns the most first index it finds that is most fractional
	int i,j, most_frac, temp, sizeList;
	double difference,currdiff;
	double *weight;
	weight=create_double_vector(M);
	most_frac=-1;
	difference=0;
	sizeList=0;
	for(i=0;i<M;i++){
		if(x[i]>0.01 && x[i]<0.99){
			splitvars[sizeList]=i;
			weight[sizeList]=(1-ABS(x[i]-0.5))*edges[i].f;
			sizeList++;
		}
	}
	/******************Now we begin sorting******************/
	for (i = 0; i < sizeList; i++)
    {
        for (j = 0; j < (sizeList - i - 1); j++)
        {
            if (weight[j] < weight[j + 1])
            {
                temp = weight[j];
                weight[j] = weight[j + 1];
                weight[j + 1] = temp;
				/*****Now change the index set***/
				temp = splitvars[j];
                splitvars[j] = splitvars[j + 1];
                splitvars[j + 1] = temp;
            }
        }
    }
	/*printf("Sorted list is \n");
	for(i=0;i<sizeList;i++){
		printf("Arc %d with solution value %lf and weight %lf\n", splitvars[i],x[splitvars[i]],weight[i]);
	}*/
	free(weight);
	return sizeList;
}

int LiftandP(int split_var, double *soln){
	int i,j,k;
	int index,index1;  // auxiliar indices to fill in the constraint matrix
	int *countindfeas; //Array with the number of indices explored thus far feasibility cuts
	int *countindopt; //Array with the number of indices explored thus far optimality cuts
	double cut_val;
	int nodecount;
	//Variables to call cplex
	CPXLPptr  lpaux;      // data strucutre to store a problem in cplex ...................
	CPXENVptr envaux;     // cplex environment.............................................
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double    *obj;    // objective function coefficients ..............................
	double    *rhs;    // right hand side of constraints ................................
	char      *sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int       *matbeg; // index of first non-zero element in each row...................
	int       *matind; // associated column of each non-zero element ...................
	double    *matval; // coefficient values for the non-zero elements of constraints....
	double    *lb;     // lower bounds of variables.....................................
	double    *ub;     // upper bounds of variables.....................................
	int       status;  // optimization status......................... .................
	char probname[16]; // problem name for cplex .......................................
	char      *ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	double    *CGLPsol;
	int       num_alpha, num_u,num_v, num_final;
	int       iteration;
	int *pos_alpha;
	int *pos_u;
	int *pos_v;
	int pos_u0;
	int pos_v0;
	int pos_beta;
	double mk;
	int    UPlimit;
	

	/*****************We declare the calculated solutions***********************/	
	double alphacalc1,alphacalc2, betacalc1,betacalc2;
	/*****************************Initializing the new position matrices *******************************/
	
	countBigMatrix=num_feasact+num_optact+2*M+K;  //tHE 2*M+K is because of the bounding constraints on the arcs and the artificial commodity variables
	countindfeas=create_int_vector(num_feasact);
	countindopt=create_int_vector(num_optact);
	pos_alpha=create_int_vector(M+K);
	pos_u=create_int_vector(countBigMatrix);
	pos_v=create_int_vector(countBigMatrix);
	LifPcut=create_double_vector(M+K+1);
	
	/*************************************************************************************************/
	
	/*****************************Going to initialize the cplex environment*************************/
	//Initialize CPLEX environment
	envaux = CPXopenCPLEX (&status); 
	if ( envaux == NULL ) {
		char  errmsg[1024];
		printf ("Could not open CPLEX. \n");
		CPXgeterrorstring (envaux, status, errmsg);
		printf ("%s", errmsg);
	}

    // Create the problem in CPLEX 
	strcpy(probname,"LPMP"); //Copy the name placed as the problem name
	lpaux = CPXcreateprob (envaux, &status, probname);  // Initialize the pointer that signals the problem 
	
	if ( envaux == NULL ) {
		char  errmsg[1024];
		printf ("Could not create LP. \n");
		CPXgeterrorstring (envaux, status, errmsg);
		printf ("%s", errmsg);
	}
	/*************************************************************************************************/

	/*****************************Defining all the variables *******************************/
	/******If we haven't decided which one to perform disjunction on then we will determine it ourselves*****/	                                        
    index1 = 0;  // index of columns		//Define alpha variables 
	numcols = M+K;
	d_vector(&obj,numcols,"open_cplex:1"); //objective function value
	d_vector(&lb,numcols,"open_cplex:8");//lower bound
	d_vector(&ub,numcols,"open_cplex:9");//upper bound
    for(i=0;i<M;i++){
	   pos_alpha[i] = index1;
	   obj[index1] = soln[pos_y[edges[i].i][edges[i].j]];
	   //printf("Value of arc (%d,%d) is %lf\n", edges[i].i,edges[i].j, soln[pos_y[edges[i].i][edges[i].j]]);
	   lb[index1] = -CPX_INFBOUND;
       ub[index1] = CPX_INFBOUND;
       index1++;
	}
	for(i=M;i<M+K;i++){
	   pos_alpha[i] = index1;
	   obj[index1] = soln[pos_art[i-M]];
	   //printf("Value of routing cost (%d,%d) is %lf\n", com[i-M].i,com[i-M].j, soln[pos_art[i-M]]);
	   lb[index1] = -CPX_INFBOUND;
       ub[index1] = CPX_INFBOUND;
	   index1++;	  
	}
	//getchar();
	status = CPXnewcols (envaux, lpaux, index1, obj, lb, ub, NULL, NULL);
	if( status )  fprintf (stderr,"CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	num_alpha = index1;//number of Z variables
	//printf("There are %d alpha variables accumulateing to %d total variables so far\n",index1, num_alpha);

	                                        //Define the u variables
    index1 = 0;  // index of columns
	numcols = countBigMatrix;
	d_vector(&obj,numcols,"open_cplex:1");
	d_vector(&lb,numcols,"open_cplex:8");
	d_vector(&ub,numcols,"open_cplex:9");	
	for(i=0;i<countBigMatrix;i++){
    		pos_u[i] = num_alpha+ index1; // position of the additional variables
			obj[index1] = 0;
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
	}		
    status = CPXnewcols (envaux, lpaux, index1, obj, lb, ub, NULL, NULL);
	if( status ) 
      fprintf (stderr,"CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	num_u = num_alpha+index1;//number of x variables
	//printf("There are %d u variables accumulateing to %d total variables so far\n",index1,num_u);

	                                        //Define the v variables
    index1 = 0;  // index of columns
	numcols = countBigMatrix;
	d_vector(&obj,numcols,"open_cplex:1");
	d_vector(&lb,numcols,"open_cplex:8");
	d_vector(&ub,numcols,"open_cplex:9");
	for(i=0;i<countBigMatrix;i++){
    		pos_v[i] = num_u + index1; // position of the additional variables
			obj[index1] = 0;
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
	}		
    status = CPXnewcols (envaux, lpaux, index1, obj, lb, ub, NULL, NULL);
	if( status ) 
      fprintf (stderr,"CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	num_v = num_u+index1;//number of x variables
	//printf("There are %d alpha v accumulateing to %d total variables so far\n",index1,num_v);
									//Define the other three variables
	index1 = 0;  // index of columns
	numcols = 3;
	d_vector(&obj,numcols,"open_cplex:1");
	d_vector(&lb,numcols,"open_cplex:8");
	d_vector(&ub,numcols,"open_cplex:9");
    pos_u0 = num_v + index1; // position of the u0 variable
	obj[index1] = 0;
	lb[index1] =0;
	ub[index1] = CPX_INFBOUND;
	index1++;
	pos_v0 = num_v + index1; // position of the v0 variable
	obj[index1] = 0;
	lb[index1] =0;
	ub[index1] = CPX_INFBOUND;
	index1++;
	pos_beta = num_v + index1; // position of the beta variable
	obj[index1] = -1;
	lb[index1] =-CPX_INFBOUND;
	ub[index1] = CPX_INFBOUND;
	index1++;
    status = CPXnewcols (envaux, lpaux, index1, obj, lb, ub, NULL, NULL);
	if( status ) 
      fprintf (stderr,"CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	num_final = num_v+index1;//number of x variables
	//printf("There are %d alpha variables accumulateing to %d total variables so far\n",index1,num_final);
	//getchar();
	/*************************************************************************************************/
	/***********************Here, let's chose what's the variable we will use***************************/

	/*************************************************************************************************/
	                                             //First set of constraints with the us
	for(j=0;j<M+K;j++){
		numrows = 1;
		numnz = countBigMatrix;
		d_vector(&rhs,numrows,"open_cplex:2");
		c_vector(&sense,numrows,"open_cplex:3");
		i_vector(&matbeg,numrows,"open_cplex:4");
		i_vector(&matind,numnz,"open_cplex:6");
		d_vector(&matval,numnz,"open_cplex:7");
		index = 0;//columns
		index1 = 0;//rows   
		sense[index1]='G';      //These are the alpha and u constraints
		matbeg[index1] = index;
		rhs[index1++]=0;
		matind[index] = pos_alpha[j];
		matval[index] = 1;			 
		index++;
		if(j==split_var){
			matind[index] = pos_u0;
			matval[index] = 1;
			index++;
		}
		if(j<M){				//We only do this for the variables that are arcs not the continuous artificial ones.
			for(i=0;i<num_feasact;i++){	
				for(k=countindfeas[i];k<feascuts[feasact[i]].numcol;k++){
					if(feascuts[feasact[i]].ind[k]==j){ //if A[i][j] is different from 0
						matind[index] = pos_u[i];
						matval[index] = -feascuts[feasact[i]].coeff[k];
						index++;
						countindfeas[i]=k;
						//printf("Feasibility constraint [%d] will be in restriction %d with coefficient %lf with pos_u[%d]=%d\n", feasact[i], j, matval[index-1],i,matind[index-1]); 
					}
					if(feascuts[feasact[i]].ind[k]>j) break;
				}
			}
			
		}
		for(i=0;i<num_optact;i++){	
			for(k=countindopt[i];k<optcuts[optact[i]].numnz;k++){
				if(optcuts[optact[i]].ind[k]==j){ //if A[i][j] is different from 0
					matind[index] = pos_u[num_feasact+i];
					matval[index] = -optcuts[optact[i]].coeff[k];
					//printf("Optimality constraint [%d] will be in restriction %d with coefficient %lf with pos_u[%d]=%d\n", optact[i], j, matval[index],i,matind[index]); 
					/*if(pos_u[num_feasact+i]<0 || optcuts[optact[i]].coeff[k]<0){
							printf("Something wrong at optcuts i=%d and j=%d\n",i,j);
							getchar();
					}*/
					index++;
					countindopt[i]=k;
				}
				if(optcuts[optact[i]].ind[k]>j) break;
			}
		}
		//if(j>=M)getchar();
/***********Now we will be adding the usual bounding constraints that are necessary**************************/
		if(j<M){		//This is if it falls within the arc variables
			//First it's the more than or equal to 0 of the arc variables
			matind[index] = pos_u[num_feasact+num_optact+j];
			matval[index] = -1;
			index++;
			//Next it's the less than or equal to 1 of the arc variables
			matind[index] = pos_u[num_feasact+num_optact+M+j];
			matval[index] = 1;
			index++;
		}
		else{ //Now it's the more than or equal to 0 of the continuous variables
			matind[index] = pos_u[num_feasact+num_optact+M+j];
			matval[index] = -1;
			index++;
		}
		status = CPXaddrows (envaux, lpaux, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
		if( status ) fprintf (stderr,"CPXaddrows failed.\n");
		free(matbeg);
		free(matind);
		free(matval);
		free(sense);
		free(rhs);
/**********************************************************************************************************/
	}
	numrows = 2;
	numnz = 2*countBigMatrix;
	d_vector(&rhs,numrows,"open_cplex:2");
	c_vector(&sense,numrows,"open_cplex:3");
	i_vector(&matbeg,numrows,"open_cplex:4");
	i_vector(&matind,numnz,"open_cplex:6");
	d_vector(&matval,numnz,"open_cplex:7");
	index = 0;//columns
	index1 = 0;//rows 
	sense[index1]='G';      //These are the Beta and u constraints with >
	matbeg[index1] = index;
	rhs[index1++]=0;
	matind[index] = pos_beta;
	matval[index] = -1;			 
	index++;
	for(i=0;i<num_feasact;i++){		 		 
		matind[index] = pos_u[i];
		matval[index] = 1;
		//printf("RHS %d is %lf\n", i, A[i][M+K]); 
		index++;
	}
	for(i=0;i<num_optact;i++){		 		 
		matind[index] = pos_u[num_feasact+i];
		matval[index] = optcuts[optact[i]].RHS;
		//printf("RHS %d is %lf\n", i, A[i][M+K]); 
		index++;
	}
/****************Adding the right hand sides of the bounding constraints************************************/
	for(i=0;i<M;i++){
		matind[index] = pos_u[num_feasact+num_optact+M+i];
		matval[index] = -1; 
		index++;
	}
	for(i=0;i<K;i++){
		matind[index] = pos_u[num_feasact+num_optact+2*M+i];
		matval[index] = shortest[i].costwdem; 
		index++;
	}
/**********************************************************************************************************/
	sense[index1]='L';      //These are the Beta and u constraints with <
	matbeg[index1] = index;
	rhs[index1++]=0;
	matind[index] = pos_beta;
	matval[index] = -1;			 
	index++;
	for(i=0;i<num_feasact;i++){		 		 
		matind[index] = pos_u[i];
		matval[index] = 1;
		index++;
	}
	for(i=0;i<num_optact;i++){		 		 
		matind[index] = pos_u[num_feasact+i];
		matval[index] = optcuts[optact[i]].RHS; 
		index++;
	}
/****************Adding the right hand sides of the bounding constraints************************************/
	for(i=0;i<M;i++){
		matind[index] = pos_u[num_feasact+num_optact+M+i];
		matval[index] = -1; 
		index++;
	}
	for(i=0;i<K;i++){
		matind[index] = pos_u[num_feasact+num_optact+2*M+i];
		matval[index] = shortest[i].costwdem; 
		index++;
	}
/**********************************************************************************************************/
	status = CPXaddrows (envaux, lpaux, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if( status ) fprintf (stderr,"CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	//printf("Entered Lift and Project with (%d,%d)\n", edges[split_var].i,edges[split_var].j);
/**********************Reset the counters of the countind****************************************/
	UPlimit=num_feasact;
	if(num_optact>UPlimit)UPlimit=num_optact;
	for(i=0;i<UPlimit;i++){
		if (i<num_feasact) countindfeas[i]=0;
		if(i<num_optact)  countindopt[i]=0;
	}

/*************************************************************************************************/
	                                           //Second set of constraints, with the Vs
	numrows = M+K+2;
	numnz = (M+K+2)*countBigMatrix;
	d_vector(&rhs,numrows,"open_cplex:2");
	c_vector(&sense,numrows,"open_cplex:3");
	i_vector(&matbeg,numrows,"open_cplex:4");
	i_vector(&matind,numnz,"open_cplex:6");
	d_vector(&matval,numnz,"open_cplex:7");
	index = 0;//columns
	index1 = 0;//rows   
	for(j=0;j<M+K;j++){
		sense[index1]='G';      //First we write that there must be an arc leaving the origin of each commodity
		matbeg[index1] = index;
		rhs[index1++]=0;
		matind[index] = pos_alpha[j];
		matval[index] = 1;			 
		index++;
		if(j==split_var){
			matind[index] = pos_v0;
			matval[index] = -1;
			index++;
		}
		if(j<M){				//We only do this for the variables that are arcs not the continuous artificial ones.
			for(i=0;i<num_feasact;i++){	
				for(k=countindfeas[i];k<=feascuts[feasact[i]].numcol;k++){
					if(feascuts[feasact[i]].ind[k]==j){ //if A[i][j] is different from 0
						matind[index] = pos_v[i];
						matval[index] = -feascuts[feasact[i]].coeff[k];
						index++;
						countindfeas[i]=k;
					}
					if(feascuts[feasact[i]].ind[k]>j) break;
				}
			}
		}
		for(i=0;i<num_optact;i++){	
			for(k=countindopt[i];k<=optcuts[optact[i]].numnz;k++){
				if(optcuts[optact[i]].ind[k]==j){ //if A[i][j] is different from 0
					matind[index] = pos_v[num_feasact+i];
					matval[index] = -optcuts[optact[i]].coeff[k];
					index++;
					countindopt[i]=k;
				}
				if(optcuts[optact[i]].ind[k]>j) break;
			}
		}
/***********Now we will be adding the usual bounding constraints that are necessary**************************/
		if(j<M){		//This is if it falls within the arc variables
			//First it's the more than or equal to 0 of the arc variables
			matind[index] = pos_v[num_feasact+num_optact+j];
			matval[index] = -1;
			index++;
			//Next it's the less than or equal to 1 of the arc variables
			matind[index] = pos_v[num_feasact+num_optact+M+j];
			matval[index] = 1;
			index++;
		}
		else{ //Now it's the more than or equal to 0 of the continuous variables
			matind[index] = pos_v[num_feasact+num_optact+M+j];
			matval[index] = -1;
			index++;
		}
/**********************************************************************************************************/
	}
	sense[index1]='G';      //These are the Beta and u constraints with >
	matbeg[index1] = index;
	rhs[index1++]=0;
	matind[index] = pos_beta;
	matval[index] = -1;			 
	index++;
	for(i=0;i<num_feasact;i++){		 		 
		matind[index] = pos_v[i];
		matval[index] = 1;
		//printf("RHS %d is %lf\n", i, A[i][M+K]); 
		index++;
	}
	for(i=0;i<num_optact;i++){		 		 
		matind[index] = pos_v[num_feasact+i];
		matval[index] = optcuts[optact[i]].RHS;
		//printf("RHS %d is %lf\n", i, A[i][M+K]); 
		index++;
	}
/****************Adding the right hand sides of the bounding constraints************************************/
	for(i=0;i<M;i++){
		matind[index] = pos_v[num_feasact+num_optact+M+i];
		matval[index] = -1; 
		index++;
	}
	for(i=0;i<K;i++){
		matind[index] = pos_v[num_feasact+num_optact+2*M+i];
		matval[index] =  shortest[i].costwdem; 
		index++;
	}
/**********************************************************************************************************/
	matind[index] = pos_v0;
	matval[index] = 1;
	index++;
	sense[index1]='L';      //These are the Beta and u constraints with <
	matbeg[index1] = index;
	rhs[index1++]=0;
	matind[index] = pos_beta;
	matval[index] = -1;			 
	index++;
	for(i=0;i<num_feasact;i++){		 		 
		matind[index] = pos_v[i];
		matval[index] = 1;
		//printf("RHS %d is %lf\n", i, A[i][M+K]); 
		index++;
	}
	for(i=0;i<num_optact;i++){		 		 
		matind[index] = pos_v[num_feasact+i];
		matval[index] = optcuts[optact[i]].RHS;
		//printf("RHS %d is %lf\n", i, A[i][M+K]); 
		index++;
	}
/****************Adding the right hand sides of the bounding constraints************************************/
	for(i=0;i<M;i++){
		matind[index] = pos_v[num_feasact+num_optact+M+i];
		matval[index] = -1; 
		index++;
	}
	for(i=0;i<K;i++){
		matind[index] = pos_v[num_feasact+num_optact+2*M+i];
		matval[index] = shortest[i].costwdem; 
		index++;
	}
/**********************************************************************************************************/
	matind[index] = pos_v0;
	matval[index] = 1;
	index++;
	status = CPXaddrows (envaux, lpaux, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if( status ) fprintf (stderr,"CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

/*****************We need to now add the normalization constraint**************************************/
	numrows = 2;
	numnz = 2*(2+2*countBigMatrix);
	d_vector(&rhs,numrows,"open_cplex:2");
	c_vector(&sense,numrows,"open_cplex:3");
	i_vector(&matbeg,numrows,"open_cplex:4");
	i_vector(&matind,numnz,"open_cplex:6");
	d_vector(&matval,numnz,"open_cplex:7");
	index = 0;//columns
	index1 = 0;//rows 
	sense[index1]='L';      //First we write that there must be an arc leaving the origin of each commodity
	matbeg[index1] = index;
	rhs[index1++]=1;
	for(j=0;j<countBigMatrix;j++){		
		matind[index] = pos_u[j];
		matval[index] = 1;			 
		index++;
		matind[index] = pos_v[j];
		matval[index] = 1;			 
		index++;
	}
	matind[index] = pos_v0;
	matval[index] = 1;			 
	index++;
	matind[index] = pos_u0;
	matval[index] = 1;			 
	index++;
	sense[index1]='G';      //First we write that there must be an arc leaving the origin of each commodity
	matbeg[index1] = index;
	rhs[index1++]=1;
	for(j=0;j<countBigMatrix;j++){		
		matind[index] = pos_u[j];
		matval[index] = 1;			 
		index++;
		matind[index] = pos_v[j];
		matval[index] = 1;			 
		index++;
	}
	matind[index] = pos_v0;
	matval[index] = 1;			 
	index++;
	matind[index] = pos_u0;
	matval[index] = 1;			 
	index++;
	status = CPXaddrows (envaux, lpaux, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	if( status ) printf ("CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	/*************************************************************************************************/
	/**********************Reset the counters of the countind****************************************/
	UPlimit=num_feasact;
	if(num_optact>UPlimit)UPlimit=num_optact;
	for(i=0;i<UPlimit;i++){
		if (i<num_feasact) countindfeas[i]=0;
		if(i<num_optact)  countindopt[i]=0;
	}

/*************************************************************************************************/
	CPXchgobjsen(envaux, lpaux,CPX_MIN);
	//CPXsetintparam(envaux,CPX_PARAM_LPMETHOD,CPX_ALG_PRIMAL );//Use the Primal Algorithm

	CPXsetintparam(envaux,CPX_PARAM_NUMERICALEMPHASIS,1); //Numerical precision of the time    
    //CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.    
    //CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,1);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfeas
    //CPXsetdblparam(env,CPX_PARAM_TILIM,86400); // time limit
    //CPXsetdblparam(env,CPX_PARAM_TRELIM, 14000); // B&B memory limit
    //CPXsetdblparam(env,CPX_PARAM_EPGAP, 0.0000000001); // e-optimal solution (%gap)
    //CPXsetdblparam(env,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)    
	//CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use	
	//CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
    //CPXsetintparam(env,CPX_PARAM_HEURFREQ, -1); //heuristic frequency and intensisty 	
	//CPXsetdblparam(env,CPX_PARAM_CUTUP,UpperBound+.01); // provide an initial upper bound
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(env,CPX_PARAM_PARALLELMODE, 1); 	
	//CPXsetintparam(env,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_BESTBOUND);  // MIP emphasis: optimality, feasibility, moving best bound
	/*Additional parameter settings for using the branch and cut*/
   /************************************************/
	//CPXwriteprob(envaux,lpaux,"CGLP.lp",NULL);   
	CPXlpopt(envaux,lpaux);  //solve the linear program	
    status=CPXgetstat(envaux,lpaux);
	if(status==1||status==2){
		CPXgetobjval(envaux,lpaux,&value);
		//printf("The optimal LP solution of the active constraints is %lf\n",value);
		numcols = CPXgetnumcols (envaux,lpaux);
		/*printf("The here are %d columns\n",numcols);
		getchar();*/
		d_vector(&CGLPsol,numcols,"open_cplex:0");			
		CPXgetx(envaux,lpaux,CGLPsol,0, numcols-1);
		/*for(j=0;j<numcols;j++){
			printf("CGLPsol[%d]=%lf\n",j,CGLPsol[j]);
			printf("LifPcut[%d]=%lf\n",j,LifPcut[j]);
			getchar();
		}*/
		
		alphacalc=create_double_vector(M+K);
		if(value<-0.0001){
			value=1;
			LPcuts[LifPcount].numnz=0;
			LPcuts[LifPcount].coeff=create_double_vector(M+K);
			LPcuts[LifPcount].ind=create_int_vector(M+K);
			for(j=0;j<M+K;j++){
				LifPcut[j]=CGLPsol[pos_alpha[j]];
				//printf("CGLPsol[%d]=%lf\n",pos_alpha[j],CGLPsol[pos_alpha[j]]);
				alphacalc1=0;
				alphacalc2=0;
				if(j<M){				//We only do this for the variables that are arcs not the continuous artificial ones.
					for(i=0;i<num_feasact;i++){	
						for(k=countindfeas[i];k<=feascuts[feasact[i]].numcol;k++){
							if(feascuts[feasact[i]].ind[k]==j){ //if A[i][j] is different from 0								
								alphacalc1+=CGLPsol[pos_u[i]]*feascuts[feasact[i]].coeff[k];
								alphacalc2+=CGLPsol[pos_v[i]]*feascuts[feasact[i]].coeff[k];
								countindfeas[i]=k;
							}
							if(feascuts[feasact[i]].ind[k]>j) break;
						}
					}
				}
				for(i=0;i<num_optact;i++){	
					for(k=countindopt[i];k<=optcuts[optact[i]].numnz;k++){
						if(optcuts[optact[i]].ind[k]==j){ //if A[i][j] is different from 0
							alphacalc1+=CGLPsol[pos_u[num_feasact+i]]*optcuts[optact[i]].coeff[k];
							alphacalc2+=CGLPsol[pos_v[num_feasact+i]]*optcuts[optact[i]].coeff[k];
							countindopt[i]=k;
						}
						if(optcuts[optact[i]].ind[k]>j) break;
					}
				}
		/***********Now we will be adding the usual bounding constraints that are necessary**************************/
				if(j<M){		//This is if it falls within the arc variables
					//First it's the more than or equal to 0 of the arc variables
					alphacalc1+=-CGLPsol[pos_u[num_feasact+num_optact+j]]+CGLPsol[pos_u[num_feasact+num_optact+M+j]];
					alphacalc2+=-CGLPsol[pos_v[num_feasact+num_optact+j]]+CGLPsol[pos_v[num_feasact+num_optact+M+j]];
				}
				else{ //Now it's the more than or equal to 0 of the continuous variables
					alphacalc1+=-CGLPsol[pos_u[num_feasact+num_optact+M+j]];
					alphacalc2+=-CGLPsol[pos_v[num_feasact+num_optact+M+j]];
				}
		/**********************************************************************************************************/
				if(j==split_var){
					alphacalc1-=CGLPsol[pos_u0];
					alphacalc2+=CGLPsol[pos_v0];
				}
				if(j!=split_var){												//Cut strengthening
					mk=(alphacalc2-alphacalc1)/(CGLPsol[pos_u0]+CGLPsol[pos_v0]);
					alphacalc1=alphacalc1+CGLPsol[pos_u0]*ceil(mk);
					alphacalc2=alphacalc2-CGLPsol[pos_v0]*floor(mk);
					if(alphacalc1<alphacalc2 && ABS(alphacalc1)>0.00000001){	
						LifPcut[j]=alphacalc1*100;
						LPcuts[LifPcount].coeff[LPcuts[LifPcount].numnz]=alphacalc1;
						LPcuts[LifPcount].ind[LPcuts[LifPcount].numnz]=j;
						//optcuts[optcount].coeff[LPcuts[LifPcount].numnz]=alphacalc1;
						//optcuts[optcount].ind[LPcuts[LifPcount].numnz]=j;
						LPcuts[LifPcount].numnz++;
					}
					else if(ABS(alphacalc2)>0.00000001){
						LifPcut[j]=alphacalc2*100;
						LPcuts[LifPcount].coeff[LPcuts[LifPcount].numnz]=alphacalc2;
						LPcuts[LifPcount].ind[LPcuts[LifPcount].numnz]=j;
						//optcuts[optcount].coeff[LPcuts[LifPcount].numnz]=alphacalc2;
						//optcuts[optcount].ind[LPcuts[LifPcount].numnz]=j;
						LPcuts[LifPcount].numnz++;
					}
				}
				else{
					if(alphacalc1>alphacalc2 && ABS(alphacalc1)>0.00000001){	
						LifPcut[j]=alphacalc1*100;
						LPcuts[LifPcount].coeff[LPcuts[LifPcount].numnz]=alphacalc1;
						LPcuts[LifPcount].ind[LPcuts[LifPcount].numnz]=j;
						//optcuts[optcount].coeff[LPcuts[LifPcount].numnz]=alphacalc1;
						//optcuts[optcount].ind[LPcuts[LifPcount].numnz]=j;
						LPcuts[LifPcount].numnz++;
					}
					else if(ABS(alphacalc2)>0.00000001){ 
						LifPcut[j]=alphacalc2*100;
						LPcuts[LifPcount].coeff[LPcuts[LifPcount].numnz]=alphacalc2;
						LPcuts[LifPcount].ind[LPcuts[LifPcount].numnz]=j;
						//optcuts[optcount].coeff[LPcuts[LifPcount].numnz]=alphacalc2;     
						//optcuts[optcount].ind[LPcuts[LifPcount].numnz]=j;
						LPcuts[LifPcount].numnz++;
					}
				}
				
			}
			//optcuts[optcount].numnz=LPcuts[LifPcount].numnz;
			LPcuts[LifPcount].RHS=CGLPsol[pos_beta];
			//optcuts[optcount].RHS=LPcuts[LifPcount].RHS;
			//Creating the new arrays for optimality cuts
			//optcount++;
			//optcuts[optcount].coeff=create_double_vector(M+K+1);
			//optcuts[optcount].ind=create_int_vector(M+K+1);	
			LifPcut[M+K]=CGLPsol[pos_beta];			
			free(CGLPsol);
		}
		else value=0;
	}
	else{
		printf("Obtained a status of %d\n", status);
		value=0;
	}
	free(pos_u);
	free(pos_v);
	free(pos_alpha);
	free(alphacalc);
	free(countindfeas);
	free(countindopt);
	
	if ( lpaux != NULL ) {
			  status = CPXfreeprob (envaux, &lpaux);
			  if ( status ) {
				fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
			  }
	}
	if ( envaux != NULL ) {
			  status = CPXcloseCPLEX (&envaux);
			  if ( status ) {
				char  errmsg[1024];
				fprintf (stderr, "Could not close CPLEX environment.\n");
				CPXgeterrorstring (envaux, status, errmsg);
				fprintf (stderr, "%s", errmsg);
			  }
	}
	return value;
}