#include "def.h"

extern int        N, M, K ;
extern double      *d, *b;
extern double     MAX_DOUBLE;
//extern int        *open_facility;
//extern int        *customer_assign;
extern EDGE   *edges;
extern COMM   *com;



double solve_BendersMC(void)
{
	int i,j,k;
	int index,index1;  // auxiliar indices to fill in the constraint matrix
	
	//Variables to call cplex
	CPXLPptr  lp;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env;     // cplex environment.............................................
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
	double    value;   // objevtive value of solution ..................................
	double   *x;
	int		ArcsOpt;
	
	

	/*Declare the structures for the cut callbacks*/
	/************************************************/
	CUTINFO lazyconinfo;// These are the ones that are necessary to the model for the solution to actually be feasible. These will be verified at each integer node only.
	CUTINFO usercutinfo;     // These are the ones that are not necessary to the model for the solution to actually be feasible. These appear at every node in order to improve the lower bound
	/*user cuts info initialized*/
	/************************************************/

	average=create_double_vector(size_av); //added to calculate the average improvement
	for(i=0;i<size_av;i++) average[i]=0;
	/*****************************Going to initialize variables*************************/
    objsen = 1; //min
	lazycalled=0;
	Flag_Check=0;
	flag_keepcut=1;
	
	//Initialize CPLEX environment
	env = CPXopenCPLEX (&status); 
	if ( env == NULL ) {
		char  errmsg[1024];
		printf ("Could not open CPLEX. \n");
		CPXgeterrorstring (env, status, errmsg);
		printf ("%s", errmsg);
	}

    // Create the problem in CPLEX 
	strcpy(probname,"MP"); //Copy the name placed as the problem name
	lp = CPXcreateprob (env, &status, probname);  // Initialize the pointer that signals the problem 
	
	if ( env == NULL ) {
		char  errmsg[1024];
		printf ("Could not create LP. \n");
		CPXgeterrorstring (env, status, errmsg);
		printf ("%s", errmsg);
	}
	
	/************We entered into cplex the arrays with the design variable****************/
	status = CPXnewcols (env, lp, M, yobj, ylb, yub, yctype, NULL);
	if( status )  fprintf (stderr,"CPXnewcols failed.\n");
	/**********************************************************************************/
	/************We entered into cplex the arrays with the routing cost variables****************/
	status = CPXnewcols (env, lp, K, zobj, zlb, zub, zctype, NULL);

	/**********************************************************************************/	
	

	/**********************************We will add the active cuts we've obtained from LP*********************/
	/********************************First the feasibility cuts**************************************************/
	numrows = 1;
	d_vector(&rhs,numrows,"open_cplex:2");
	c_vector(&sense,numrows,"open_cplex:3");
	i_vector(&matbeg,numrows,"open_cplex:4");
	sense[0]='G';      
	matbeg[0] = 0;
	for(k=0;k<num_feasact;k++){	
		rhs[0]=1;
		if(k<0){
			status = CPXaddrows (env, lp, 0, 1, feascuts[feasact[k]].numcol, rhs, sense, matbeg, feascuts[feasact[k]].ind, feascuts[feasact[k]].coeff, NULL, NULL);
		}
		else{
			status=CPXaddusercuts (env,lp,1, feascuts[feasact[k]].numcol, rhs, sense, matbeg, feascuts[feasact[k]].ind, feascuts[feasact[k]].coeff,NULL);
			status=CPXaddlazyconstraints (env,lp,1, feascuts[feasact[k]].numcol, rhs, sense, matbeg, feascuts[feasact[k]].ind, feascuts[feasact[k]].coeff, NULL);
		}
		if( status ) fprintf (stderr,"CPXaddrows failed.\n");    
	}
	/********************************Next the optimality cuts that are active**************************************************/
	for(k=0;k<num_optact;k++){	
		rhs[0]=optcuts[optact[k]].RHS;
		if(k<0){
			status = CPXaddrows (env, lp, 0, 1, optcuts[optact[k]].numnz, rhs, sense, matbeg, optcuts[optact[k]].ind, optcuts[optact[k]].coeff, NULL, NULL);
		}
		else{
			status=CPXaddusercuts (env,lp,1, optcuts[optact[k]].numnz, rhs, sense, matbeg, optcuts[optact[k]].ind, optcuts[optact[k]].coeff,NULL);
			status=CPXaddlazyconstraints (env,lp,1, optcuts[optact[k]].numnz, rhs, sense, matbeg, optcuts[optact[k]].ind, optcuts[optact[k]].coeff, NULL);
		}
		if( status ) fprintf (stderr,"CPXaddrows failed.\n"); 
	}
	///***************************Calculating how many digits are in maxprior*****************/
	//digits=0;
	////printf("maxprior is %d ",maxprior);
	//while(maxprior >1)
	//{
	//	maxprior/=10;
	//	digits++;
	//}
	////printf("and has %d digits \n",digits);
	//for(j=0;j<M;j++){
	//	priority[j]=(int) (priority[j]/(powf(10,(digits-2))));//Giving priorities from 0-100 so that some coincide.
	//	//printf("variable %d has priority %d\n", j,priority[j]);
	//}

	///*******************************************************************************/
	end_time = clock();
	cputimeBD = (double)(end_time - start_time) / CLOCKS_PER_SEC;

	CPXchgobjsen(env, lp,CPX_MIN);
	CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use
	CPXsetintparam(env,CPX_PARAM_MIPDISPLAY,3); //different levels of output display
    CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,2);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfea
	//CPXsetintparam(env,CPX_PARAM_BRDIR,-1);
	CPXsetdblparam(env,CPX_PARAM_TILIM,Time-cputimeBD); // time limit
	//CPXsetintparam(env,CPX_PARAM_VARSEL, CPX_VARSEL_STRONG);
	/*Changing some of the precision*/
	CPXsetintparam(env,CPX_PARAM_NUMERICALEMPHASIS,1); //Numerical precision of the time
	//CPXsetintparam(env,CPX_PARAM_RINSHEUR, 5); //heuristic frequency of RINS
	//CPXsetintparam(env,CPX_PARAM_FPHEUR, 2); //heuristic frequency of Feasibility pump
	//CPXsetdblparam(env,CPX_PARAM_EPGAP, 0.0000000001); // e-optimal solution (%gap)
    //CPXsetdblparam(env,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)    	
	//CPXsetdblparam(env,CPX_PARAM_EPRHS, 0.0000001);// Feasibility tolerance
	//CPXsetdblparam(env,CPX_PARAM_EPINT, 0.0001); // integer precision


	//CPXsetintparam(env,CPX_PARAM_PREIND,0);
	CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON); //output display
	//CPXsetdblparam(env,CPX_PARAM_CUTSFACTOR, 1.0);
	//CPXsetintparam(env,CPX_PARAM_NODELIM,492);
    //CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.    
    //CPXsetdblparam(env,CPX_PARAM_TRELIM, 14000); // B&B memory limit    
	//CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
    //CPXsetintparam(env,CPX_PARAM_RINSHEUR, 5); //heuristic frequency and intensisty 	
	//CPXsetdblparam(env,CPX_PARAM_CUTUP,UpperBound+.01); // provide an initial upper bound
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(env,CPX_PARAM_PARALLELMODE, 1); 
	//CPXsetintparam(env,CPX_PARAM_PREIND,0);
	//CPXsetintparam(env,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_BESTBOUND);  // MIP emphasis: optimality, feasibility, moving best bound   
	
	/*Additional parameter settings for using the branch and cut*/
	/************************************************/
	status = CPXsetintparam (env, CPX_PARAM_PRELINEAR, 0);											/* Do not use presolve */
	//status = CPXsetintparam (env,CPXPARAM_Preprocessing_Reduce,1);									/*No dual reductions*/
	status = CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL); 					/* Turn on traditional search for use with control callbacks */
	status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);									 /* Let MIP callbacks work on the original model */   
   /************************************************/

	/* Code to use lazy constraints*/
	/************************************************/
	status = makelazyconstraint (env, lp, &lazyconinfo); 											//Linking the LP defined as the LP in the lazyconinfo structure and get sizes
	status = CPXsetlazyconstraintcallbackfunc (env, lazy_feas_callback1, &lazyconinfo); 

	/* Code to use usercut constraints*/
	/************************************************/
	status = makeusercuts (env, lp, &usercutinfo); 												//Linking the LP defined as the LP in the usercutinfo structure and get sizes
	status = CPXsetusercutcallbackfunc (env, user_cut_callback1 , &usercutinfo); 							//Now create the user cuts using our cutcallback which solves a separation problem
   /**********************************************************************************************************/
	/* Code to use Heuristic Callback*/
	/************************************************/
	status = CPXsetheuristiccallbackfunc (env, LagrangeHeur, NULL);
	/**********************************************************************************************************/
	/* Code to add branching priority*/
	/************************************************/
	//status = CPXcopyorder(env, lp, M, indbranch, priority, NULL);    
	/**********************************************************************************************************/

	heur_checkstart= clock();
	if(Flag_Print==1){printf("Started the clock\n"); /*getchar();*/}
	CPXmipopt(env,lp);  //solve the integer program
	status=CPXgetstat(env,lp);
	if(status==101|| status==102|| status==106 || status==105|| status==107|| status==108)
	{
		end_time = clock();
		CPXgetmipobjval(env,lp,&Upper_bound);		
		CPXgetbestobjval(env,lp,&best_lower_bound);
		nodecount=CPXgetnodecnt(env,lp);
		//Printing to the screen
		if(status==101|| status==102){
			printf("Optimal Solution found with value: %.2f\n",Upper_bound);
		}
		else{
		printf("Time limit reached best solution with value:%.2f\n",Upper_bound);
		}
		
	}
	else{
		end_time = clock();
		Upper_bound=-1;
		best_lower_bound=-1;
		nodecount=-1;
		printf("Unknown stopping criterion (%d) \n",status);		
	}
	free(average);
	if ( lp != NULL ) {
			  status = CPXfreeprob (env, &lp);
			  if ( status ) {
				fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
			  }
	}
	if ( env != NULL ) {
			  status = CPXcloseCPLEX (&env);
			  if ( status ) {
				char  errmsg[1024];
				fprintf (stderr, "Could not close CPLEX environment.\n");
				CPXgeterrorstring (env, status, errmsg);
				fprintf (stderr, "%s", errmsg);
			  }
	}
	return Upper_bound;
}

/*Defining make lazycuts*/
/************************************************/
static int makelazyconstraint (CPXENVptr  env, CPXLPptr   lp, CUTINFOptr lazyconinfo)
{
   int status = 0;
   cur_numcols = CPXgetnumcols (env, lp);
   lazyconinfo->lp = lp;
   lazyconinfo->numcols=cur_numcols;
   return (status);
} 

/*Defining make usercuts*/
/************************************************/
static int makeusercuts (CPXENVptr  env, CPXLPptr   lp, CUTINFOptr usercutinfo)
{
   int status = 0;
   cur_numcols = CPXgetnumcols (env, lp);
   usercutinfo->lp = lp;
   usercutinfo->numcols=cur_numcols;
   status=CPXgetbestobjval(env,lp,&best_lb);
   return (status);
} 

static int CPXPUBLIC lazy_feas_callback1 (CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p)
{
   int status = 0;
   /*Declaring the local cut structure*/
   /************************************************/
   CUTINFOptr cutinfo = (CUTINFOptr) cbhandle;
   int      numcols  = cutinfo->numcols;
   int cutnz;
   /************************************************/
   /*Variables used for creating the cut*/
   /************************************************/
   int      *cutind;
   double   *cutval;  
   int      *cutind2;
   double   *cutval2;
   double   cutvio;
   int      addcuts = 0;
   int		num_cuts=0;
   int		addcuts2=0;
   int      optvio=0;
   /************************************************/
   int      i, j, k,index=0;
   double      *temp_x;
   double	*objval;
   double   obval;
   int			depth;
   int			com_check; //Indicator as to what commodity will be brought to trial
	/*Obtaining the solution information necessary at the current node*/
	/************************************************/
   cutind=create_int_vector(M+K);
   cutval=create_double_vector(M+K);
   d_vector(&temp_x,cur_numcols,"open_cplex:0");	
   //One percent test
   /********************************************************/
   if(Flag_Print==1){printf("Entered lazy cutcallback\n"); }
   //printf("Entered lazy cutcallback\n"); getchar();
   flag_frac=0;
   last_feas=feascount;
   /*flag_MWeval=0;*/															//Create the variable that will be used to store the solution at the node
   objval= &obval;																							//Assign the memory slot to the pointer so that we may use it
   numcols = cur_numcols;	//Transfer the information on the number of variables which we obtained at the make lazy constraint callback
   lazycalled++; 
   *useraction_p = CPX_CALLBACK_DEFAULT;																	//At this point, don't add any cuts. We are only starting to solve the separation problem
   status = CPXgetcallbacknodex (env, cbdata, wherefrom, temp_x, 0, numcols-1);									//Obtain the solution at the node
   status = CPXgetcallbacknodeobjval (env, cbdata, wherefrom, objval);	//Obtain the objective value of the solution at the curent node
   status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
   globaldepth=depth;
   if(globaldepth>maxdepth) maxdepth=globaldepth;
   violation=0.0;
	for(j=0;j<M;j++){                  //fill up the values of ybar
				  if( temp_x[pos_y[edges[j].i][edges[j].j]] > 0.9){
					ybar[edges[j].i][edges[j].j] = 1;
				  }
				  else
					ybar[edges[j].i][edges[j].j] = 0;
	}
	//if(Flag_Print==1) printf("%.2lf routing cost at depth of %d ", obval, depth);
	if(1==1){
		start_feas=clock();
		for(i=0;i<num_origins;i++){
			addcuts=maxFlowAndCut(origins[i].s, origins[i].d, origins[i].dim, ybar, cutval,cutind);
			/******************************************************************************************/
			if(addcuts==1){
				if(Feas_Type==0 /*&& repeated_feas(cutval,cutind,feasnumnz)==0*/){
					status = CPXcutcallbackadd (env, cbdata, wherefrom,  feasnumnz, 1, 'G',  cutind, cutval, CPX_USECUT_FORCE);				
					//UpdFeas_pool(1,cutval, cutind, feasnumnz);
					cuts_SC[0]++;
				}
				com_check=origins[i].ind[com_check_ind][0];
				if(Feas_Type==2){
					StrongFeas_com(com_check);
					index = 0;//columns
					//printf("RHS=%.4lf\n", Righthand);
					for(j=0;j<M;j++){			 
						if(ABS(coeff[j])>0){
							cutind[index] = pos_y[edges[j].i][edges[j].j];
							cutval[index] = coeff[j];
							//if(coeff[j]>0)printf("%.4lf * y[%d][%d]\n", coeff[j],edges[j].i,edges[j].j);
							index++;	
						}
					}		
					status = CPXcutcallbackadd (env, cbdata, wherefrom,  index, Righthand, 'G',  cutind, cutval, CPX_USECUT_PURGE);
					if( status ){  fprintf (stderr,"CPXaddrows failed.\n"); getchar();}
					cuts_SC[0]=cuts_SC[0]+1; //Update feasibilty Cuts
					//UpdFeas_pool(1,cutval, cutind, index);
					//free(coeff);
					//getchar();
				}
				num_cuts=num_cuts+addcuts;	
			}
		}
		end_feas=clock();
		timeFeas=timeFeas+(double)(end_feas - start_feas) / CLOCKS_PER_SEC;
	}
	//if(num_cuts>0) printf("added %d cutsets at INTEGER solution\n", num_cuts);
	/* Tell CPLEX that feasibility cuts have been created */ 
	if ( num_cuts > 0 ) {
	  *useraction_p = CPX_CALLBACK_SET;
	  goto Terminate;
	 }
	/***************************************************************************/
	else{
		cutvio=solve_net_sub(temp_x);         //MIN_COST flow M&W
		for(k=optcount-vioptcuts;k<optcount;k++){					//Add the violated cuts
			status = CPXcutcallbackadd (env, cbdata, wherefrom,  optcuts[k].numnz, optcuts[k].RHS, 'G',  optcuts[k].ind, optcuts[k].coeff, CPX_USECUT_PURGE);
			cuts_SC[1]++;
		}
		if(vioptcuts==0){						//Try again except this time with the actual value
			solve_net_exact(temp_x);
			for(k=optcount-vioptcuts;k<optcount;k++){					//Add the violated cuts
				status = CPXcutcallbackadd (env, cbdata, wherefrom,  optcuts[k].numnz, optcuts[k].RHS, 'G',  optcuts[k].ind, optcuts[k].coeff, CPX_USECUT_PURGE);
				cuts_SC[1]++;
			}
		}
		if(vioptcuts>0){
			*useraction_p = CPX_CALLBACK_SET;
			//printf("added %d optimality cuts at integer solution\n", vioptcuts);
			goto Terminate;
		}
		goto Terminate;
	}

Terminate: 
	/******Updating the stabilizer only if there we found a new integer solution********/
	/*if(*useraction_p!=CPX_CALLBACK_SET){
		for(j=0;j<M;j++){
			stabilizer[j]=temp_x[j];
			core_point[edges[j].i][edges[j].j]=temp_x[j];
		}
	}
	flag_incumb=1;*/
	/****************************************/
	if(base_fact-decrement*(disc)==0)disc=0;
	free(temp_x);
	free(cutval);
	free(cutind);
return (status);

} 
/* End of lazycut callback */

/*Defining the usercuts that we are going to use for fractional solutions*/
/*********************************************/
static int CPXPUBLIC user_cut_callback1 (CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p)
{
   int status = 0;
   /*Declaring the local cut structure*/
   /************************************************/
   CUTINFOptr cutinfo = (CUTINFOptr) cbhandle;
   int      numcols  = cutinfo->numcols;
   /************************************************/
   /*Variables used for creating the cut*/
   /************************************************/
   int      *cutind;
   double   *cutval; 
   int      *cutind2;
   double   *cutval2;
   double   rhs;
   double   cutvio;
   int      addcuts = 0;
   int		num_cuts=0;
   int		cutnz;
   int		depth;
   int      optvio=0;
   /************************************************/
   int      i, j, k,index=0, SEQNUM;
   double      *temp_x;
   double	*objval;
   double   obval;
   double	change;
   int		Flag_feas=1;
   int		Found_Vio=0;
   int		com_check;
   double value;
   
   if(Flag_Print==1){ printf("Entered user cutcallback\n");}
	//printf("Entered user cutcallback\n"); getchar();/*Obtaining the solution information necessary at the current node*/
	/************************************************/
	d_vector(&temp_x,cur_numcols,"open_cplex:0");																//Create the variable that will be used to store the solution at the node
	objval= &obval;																							//Assign the memory slot to the pointer so that we may use it
	numcols = cur_numcols;																					//Transfer the information on the number of variables which we obtained at the make lazy constraint callback
	*useraction_p = CPX_CALLBACK_DEFAULT;																	//At this point, don't add any cuts. We are only starting to solve the separation problem
	status = CPXgetcallbacknodex (env, cbdata, wherefrom, temp_x, 0, numcols-1);									//Obtain the solution at the node
	status = CPXgetcallbacknodeobjval (env, cbdata, wherefrom, objval);	//Obtain the objective value of the solution at the curent node
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
	status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &SEQNUM);
	//printf("Best lb is %lf\n", best_lb);getchar();
	globaldepth=depth;
	if(globaldepth>maxdepth) maxdepth=globaldepth;
	CPXgetcallbacknodelb(env, cbdata, wherefrom,lbfix, 0, M+K-1);
	CPXgetcallbacknodeub(env, cbdata, wherefrom,ubfix, 0, M+K-1);	
	cutind=create_int_vector(M+K);
	cutval=create_double_vector(M+K);

	/*We begin to fill in the information needed to evaluate whether we will check for violated cuts or not*/
	/************************************************/
	//if(depth>=depthLP) max_cuts=LifPlimit+2;
			if(last_node!=SEQNUM){															//If we are at a new node, reset the counters of the nodes and the number of cuts added
				usercutsnum=0;
				//LifPcount=0;
				rep_node=0;
				flag_upd_av=1;
				Feas_unchanged=0;
				terminate=0;
				disc=0;
				for(i=0;i<size_av;i++){
					average[i]=0;                 
				}
			}

			if(depth==0){
				usercutsnum=0; 														//If we are at the rootnode, then act as if we have not added any cuts to avoid violating the max_cuts criteria
				LPval=obval;
			}
			/*Our tactical solution to deal with the numerical instability of the feasibility cuts*/
			/************************************************/
			if(ABS(obval-last_obj)>.01) Feas_unchanged=0;      //Has the solution changed
			else Feas_unchanged++;

			if(Feas_unchanged==3){
				Feas_Type=2;
				if(Flag_Print==1)printf("\n Now the Feas_type=%d\n", Feas_Type);
				//Feas_unchanged=0;
			}

			if(Feas_unchanged==4){
				Flag_feas=0;
			}
			if(Feas_unchanged==10){
				Flag_feas=0;
				terminate=1;
			}
			/************************************************/

			//Updating the last size_av solutions in the array to see the progress
			/************************************************/
			if(flag_upd_av==1){
				if(rep_node<size_av) average[rep_node]=obval;
				else {
					for(i=0;i<size_av-1;i++) average[i]=average[i+1];
					average[size_av-1]=obval;
				}
				flag_upd_av=0;
			}
			//Calculating the average improvement over the last size_av
			/************************************************/
			
				if(rep_node<size_av-1)	change=1;										//If we haven't obtained size_av different solutions at this node, then assure we don't violate increase restriction
				else change=(average[size_av-1]-average[0])/average[0];					//If we have obtained at least size_av different solutions at this node, then see how much has the solution increased compared to size_av solutions ago
				if(depth==0) change=0.05+1;
			//Updating the variables of last node and reps
			/************************************************/
			last_node=SEQNUM;	
			last_obj=obval;
			
	/************************************************/
	/******************************We now add the Lift and Project Cuts*******************************************************/
	//if(LifPcount>0 && depth==0 /*&& (Prev_incumbent-obval)/Prev_incumbent<.01*/ && added_LP==0){
	//	for(k=0;k<LifPcount;k++){
	//		status = CPXcutcallbackadd (env, cbdata, wherefrom,  LPcuts[k].numnz, LPcuts[k].RHS, 'G',  LPcuts[k].ind, LPcuts[k].coeff,CPX_USECUT_PURGE) ;
	//		added_LP++;
	//	}		
	//	printf("Added %d LP Cuts\n",LifPcount);
	//	//getchar();
	//	*useraction_p = CPX_CALLBACK_SET;
	//}

	//if(Flag_Print==1) printf("%.2lf routing cost at depth of %d with terminate=%d ", obval, depth, terminate);
	if( (depth<=depmult || (depth>depmult && SEQNUM%nodemult==0)) && usercutsnum<max_cuts && change>tol && terminate==0){
		flag_frac=1;
		last_feas=feascount;
		/**************************************/
		//Modifying the violation**
		//violation==(EPSILON)/(1/pow(10.0,2*(depth/depmult))); //Optimality Cut
		/*if(depth<=depmult)*/violation=0.00000001;
		/*else{ violation=LPval*.1;}*/
		cutsettol=0.9999; //Cutset inequality
		/**************************************/
		//printf("Opened the following arcs\n");
		for(j=0;j<M;j++){                  //fill up the values of ybar_frac				  
				ybar_frac[edges[j].i][edges[j].j] = temp_x[pos_y[edges[j].i][edges[j].j]];
				if(depth==0) LPsol[j]=temp_x[pos_y[edges[j].i][edges[j].j]];//Keeping the LP solution value
		}
		//printf("\n");
		if(1==1){
			start_feas=clock();
			for(i=0;i<num_origins;i++){
				for(j=0;j<origins[i].dim;j++){
					addcuts=maxFlowAndMinCut(origins[i].s, origins[i].d[j], ybar_frac, cutval,cutind);
					if(addcuts==1){
						com_check=origins[i].ind[j][0];
						if(Feas_Type==0 /*&& repeated_feas(cutval,cutind,feasnumnz)==0*/){
							status = CPXcutcallbackadd (env, cbdata, wherefrom,  feasnumnz, 1, 'G',  cutind, cutval, CPX_USECUT_FORCE);	
							//UpdFeas_pool(2,cutval, cutind, feasnumnz);
							cuts_SC[2]++;
						}
						if(Feas_Type==2){
							StrongFeas_com(com_check);
							cutnz= M;
							index = 0;//columns
							/*printf("RHS=%.4lf\n", Righthand);*/
							//contador=0;
							for(j=0;j<M;j++){			 
								if(ABS(coeff[j])>0){
									cutind[index] = pos_y[edges[j].i][edges[j].j];
									cutval[index] =coeff[j];
									index++;
								}
							}		
							status = CPXcutcallbackadd (env, cbdata, wherefrom,  index+1, Righthand, 'G',  cutind, cutval, CPX_USECUT_FORCE);
							if( status ){  fprintf (stderr,"CPXaddrows failed.\n");}
							//UpdFeas_pool(1,cutval, cutind, index);
							//free(coeff);
							cuts_SC[2]++;
						}
						num_cuts=num_cuts+addcuts;

					}
				}
			}
			end_feas=clock();
			timeFeas=timeFeas+(double)(end_feas - start_feas) / CLOCKS_PER_SEC;
		}
		//if(num_cuts>0) printf("added %d cutset cuts in fractional solution", num_cuts);
		/* Tell CPLEX that feasibility cuts have been created */ 
		if ( num_cuts > 0 && Flag_feas==1) {
			*useraction_p = CPX_CALLBACK_SET;
			//printf("Exited Callback\n");
			goto Terminate;
		}
		else{
			//printf("\nEntered to check for optimality\n"); getchar();
			rep_node++;
			flag_upd_av=1;
			/******Updating the stabilizer only if there has not been an incumbent yet********/
			if(depth>0 && flag_incumb==0){
				for(j=0;j<M;j++){
					stabilizer[j]=ini_core[j];
					core_point[edges[j].i][edges[j].j]=ini_core[j];
				}
				flag_incumb=1;
			}
			/****************************************/
			cutvio=solve_net_sub(temp_x);         //MIN_COST flow M&W
			if(vioptcuts>0){				
				for(k=optcount-vioptcuts;k<optcount;k++){					//Add the violated cuts
					status = CPXcutcallbackadd (env, cbdata, wherefrom,  optcuts[k].numnz, optcuts[k].RHS, 'G',  optcuts[k].ind, optcuts[k].coeff, CPX_USECUT_PURGE);
					//status = CPXcutcallbackaddlocal (env, cbdata, wherefrom,  optcuts[k].numnz, optcuts[k].RHS, 'G',  optcuts[k].ind, optcuts[k].coeff);
					cuts_SC[3]++;
				}
			}			
			if(depth==0 && (vioptcuts==0 ||Feas_Type==2)){						//Try again except this time with the actual value
				solve_net_exact(temp_x);
				for(k=optcount-vioptcuts;k<optcount;k++){					//Add the violated cuts
					status = CPXcutcallbackadd (env, cbdata, wherefrom,  optcuts[k].numnz, optcuts[k].RHS, 'G',  optcuts[k].ind, optcuts[k].coeff, CPX_USECUT_PURGE);
					//status = CPXcutcallbackaddlocal (env, cbdata, wherefrom,  optcuts[k].numnz, optcuts[k].RHS, 'G',  optcuts[k].ind, optcuts[k].coeff);
					cuts_SC[3]++;
				}
			}
			if(vioptcuts>0){
				*useraction_p = CPX_CALLBACK_SET;
				usercutsnum++;
				//printf("added %d optimality cuts at fractional solution\n", vioptcuts);
				goto Terminate;
			}
			if(LifPcount<LifPlimit && depth==0){
				split_var=frac_index(temp_x);
				printf("The most fractional index is %d with value %lf\n", split_var, temp_x[split_var]);
				//getchar();
				ConstrMatrix_LP(temp_x,1.0);
				if(LiftandP(split_var, temp_x)==1){
					if(depth==0)status = CPXcutcallbackadd (env, cbdata, wherefrom,  LPcuts[LifPcount].numnz, LPcuts[LifPcount].RHS, 'G',  LPcuts[LifPcount].ind, LPcuts[LifPcount].coeff, CPX_USECUT_PURGE);
					else status = CPXcutcallbackaddlocal (env, cbdata, wherefrom,  LPcuts[LifPcount].numnz, LPcuts[LifPcount].RHS, 'G',  LPcuts[LifPcount].ind, LPcuts[LifPcount].coeff);
					if( status )  fprintf (stderr,"CPXaddrows failed.\n");
					//printf("Added a Lift and project Cut\n");
					free(LifPcut);
					LifPcount++;
					TotalLPcuts++;
					*useraction_p = CPX_CALLBACK_SET;
					goto Terminate;
				}
				else{
					printf("No violated constraint foound\n");
					goto Terminate;
				}
			}
			goto Terminate;
		}
	}

Terminate:
	/*printf("Exited fractional separation\n");
	if(obval>=3088089.25)getchar();*/
	Feas_Type=0;
	if(*useraction_p != CPX_CALLBACK_SET){  /*printf("Exited without any violation \n");*/ /*getchar();*/ terminate=1;}
	if(Flag_Print==1){ printf("\n"); /*getchar();*/}
	free(temp_x);
	free(cutval);
	free(cutind);
return (status);
} 
/* End of usercut callback */



static void
free_and_null (char **ptr)
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
} /* END free_and_null */ 

static int CPXPUBLIC LagrangeHeur (CPXCENVptr env, void  *cbdata, int wherefrom, void *cbhandle, double  *objval_p, double *x, int *checkfeas_p,int *useraction_p)
{
   int status = 0; 
   int depth, SEQNUM;
   int       i,j, cols, addcuts=0,opened;
   double    roundobjval, lagrangeobjval,currincumbent;
   int       *feas = NULL;
   CPXCLPptr lp;
   double    *objcoefs = NULL;
   double    **cut_sets;
   int		limit=20,num_paths=3;
   double   *bestx;
   double	*temp_x;
   double  orig_ub;
   int flag_newinc=0;
   int divisor=8;
   int indexdiv;

   bestx=create_double_vector(M+K);
   temp_x=create_double_vector(M+K);
   cols = M+K;
   status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
   status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_SEQNUM, &SEQNUM);
   status = CPXgetcallbacknodex (env, cbdata, wherefrom, temp_x, 0, cols-1);
   globaldepth=depth;
   indexdiv=depth/divisor; //index in list
   *useraction_p = CPX_CALLBACK_DEFAULT;
	if(use_firstsolution==1 /*&& depth==4*/){
		use_firstsolution=0;  
		for(i=0;i<K;i++){
			x[pos_art[i]]=bestsol[M+i];
			//printf("Routing cost of commodity %d is %lf\n", i, Best_routes[i]);
		 }		   
		 for(j=0;j<M;j++){
			if(bestsol[j]>0.5){   
				x[pos_y[edges[j].i][edges[j].j]]=1;
				//printf("Value of integer (%d,%d) is %lf\n", edges[j].i, edges[j].j, Best_arcs[edges[j].i][edges[j].j]);

			}
			else x[pos_y[edges[j].i][edges[j].j]]=0;
		 }

         *objval_p =bestobj;   //Now telling it what the new objective value is
		 //printf("Writing the best found solution with value %lf\n",bestobj);
		 /*if(Flag_Print==1) printf("The lagrangean pre solution obtained is %lf\n", lagrangeobjval);*/
		 /* Have CPLEX check the solution for integer feasibility */
		 *checkfeas_p = 1;
		  /* Tell CPLEX that a solution is being returned */
		  *useraction_p = CPX_CALLBACK_SET;
	}
	/*********Here we are checking if there is any new incumbent******************/
	//CPXgetcallbackinfo (env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &bestobj);
	//if(bestobj<Prev_incumbent-0.0001){
	//	CPXgetcallbackincumbent (env, cbdata, wherefrom, bestsol, 0, cols-1);
	///******Updating the stabilizer only if there has not been an incumbent yet********/
	//	for(j=0;j<M;j++){
	//		if(bestsol[j]>0.5){
	//			stabilizer[j]=1.0;
	//			core_point[edges[j].i][edges[j].j]=1.0;
	//		}
	//		else{
	//			stabilizer[j]=.2;
	//			core_point[edges[j].i][edges[j].j]=0.2;
	//		}
	//	}
	//}
	//Prev_incumbent=bestobj;
	/****************************************/
	

	CPXgetcallbackinfo (env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &bestobj);
	Prev_incumbent=bestobj;
	if(SEQNUM%nodemult==0 && depth>lastdepth && lasttotinpath+N<totalinpaths && useheur==1 /*&&  heurchecktime>100*/){     /****Now we will apply our CG heuristic here*****/
		CPXgetcallbackinfo (env, cbdata, wherefrom, CPX_CALLBACK_INFO_MIP_FEAS, &status);
		CPXgetcallbackincumbent (env, cbdata, wherefrom, bestsol, 0, cols-1);
		CPXgetcallbacknodelb(env, cbdata, wherefrom,lbfix, 0, M+K-1);
		CPXgetcallbacknodeub(env, cbdata, wherefrom,ubfix, 0, M+K-1);
		
		lastdepth=depth;
		lasttotinpath=totalinpaths;
		heur_checkstart=clock();
		if(Prev_incumbent>master_p(ylb, indarcinpath)){ //Executing the heuristic again
			for(i=0;i<K;i++){
				x[pos_art[i]]=bestsol[M+i];
				//printf("Routing cost of commodity %d is %lf\n", i, Best_routes[i]);
			}		   
			for(j=0;j<M;j++){
				if(bestsol[j]>0.5){   
					x[pos_y[edges[j].i][edges[j].j]]=1;
					//printf("Value of integer (%d,%d) is %lf\n", edges[j].i, edges[j].j, Best_arcs[edges[j].i][edges[j].j]);

				}
				else x[pos_y[edges[j].i][edges[j].j]]=0;
			}

			*objval_p =bestobj;   //Now telling it what the new objective value is
			*checkfeas_p = 1;
			/* Tell CPLEX that a solution is being returned */
			*useraction_p = CPX_CALLBACK_SET;
		}
		
		currheurcheck=clock();
		heurchecktime += (double)(currheurcheck - heur_checkstart) / CLOCKS_PER_SEC;

		/******************Uncomment this to activate updating the corepoint when a new incumbent solution is found*******/
		//if(flag_newinc==1){
		//	/*printf("Updated corepoint\n");
		//	getchar();*/
		//	for(j=0;j<M;j++){
		//			core_point[edges[j].i][edges[j].j]=0.1*core_point[edges[j].i][edges[j].j]+0.9*bestx[j];
		//	}
		//}
		//
		/********************************************/
	}
	/***************************************************************************/
	

TERMINATE:
	free(bestx);
   return (0);
} /* END Lagrange_Heur */


double solve_BendersCS(int indic)
{
	int i,j,k;
	int index,index1;  // auxiliar indices to fill in the constraint matrix
	int nodecount;
	int globind;  //Complete counter of the first feasibility cuts
	//Variables to call cplex
	CPXLPptr  lp;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env;     // cplex environment.............................................
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
	double    value;   // objevtive value of solution ..................................
	double   *x;
	int		forlastrest; //for the last restriction

	/*Declare the structures for the cut callbacks*/
	/************************************************/
	CUTINFO lazyconinfo;// These are the ones that are necessary to the model for the solution to actually be feasible. These will be verified at each integer node only.
	CUTINFO usercutinfo;     // These are the ones that are not necessary to the model for the solution to actually be feasible. These appear at every node in order to improve the lower bound
	/*user cuts info initialized*/
	/************************************************/

	average=create_double_vector(size_av); //added to calculate the average improvement
	for(i=0;i<size_av;i++) average[i]=0;
	/*****************************Going to initialize variables*************************/
    objsen = 1; //min
	lazycalled=0;
	Flag_Check=0;
	flag_keepcut=1;
	
	//Initialize CPLEX environment
	env = CPXopenCPLEX (&status); 
	if ( env == NULL ) {
		char  errmsg[1024];
		printf ("Could not open CPLEX. \n");
		CPXgeterrorstring (env, status, errmsg);
		printf ("%s", errmsg);
	}

    // Create the problem in CPLEX 
	strcpy(probname,"MP"); //Copy the name placed as the problem name
	lp = CPXcreateprob (env, &status, probname);  // Initialize the pointer that signals the problem 
	
	if ( env == NULL ) {
		char  errmsg[1024];
		printf ("Could not create LP. \n");
		CPXgeterrorstring (env, status, errmsg);
		printf ("%s", errmsg);
	}
	
	d_vector(&ub_vec,M,"open_cplex:9");//upper bound
	for(j=0;j<M;j++){							//which are allowed
		ub_vec[j]=1*base_sol[j];
	}
	
	/************We entered into cplex the arrays with the design variable****************/
	if(indic==0){
		if(globalcounter>=1)status = CPXnewcols (env, lp, M, yobj, ylb, yub, yctype, NULL);
		else status = CPXnewcols (env, lp, M, yobj, ylb, yub, yctype, NULL);
	}
	else status = CPXnewcols (env, lp, M, yobj, ylb, yub, yctype, NULL);

	if( status )  fprintf (stderr,"CPXnewcols failed.\n");
	/**********************************************************************************/
	/************We entered into cplex the arrays with the routing cost variables****************/
	status = CPXnewcols (env, lp, K, zobj, zlb, zub, zctype, NULL);
	if( status )  getchar();
	/**********************************************************************************/	
	
	/**********************************Here we can choose to add these*********************/
	/**********************************************************************************/
	numrows = 1;
	numnz = M;
	globind=0;
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
			globind++;
			//if(globalcounter==0 && indic==0) UpdFeas_pool(1,matval, matind, index);
			status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
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
			globind++;
			//if(globalcounter==0 && indic==0)UpdFeas_pool(1,matval, matind, index);
			status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
			if( status ) fprintf (stderr,"CPXaddrows failed.\n");
		}
		free(matbeg);
		free(matind);
		free(matval);
		free(sense);
		free(rhs);
	

	/**********************************************************************************/

	/**********************************We will add the active cuts we've obtained from LP*********************/
	/********************************First the feasibility cuts**************************************************/
	numrows = 1;
	d_vector(&rhs,numrows,"open_cplex:2");
	c_vector(&sense,numrows,"open_cplex:3");
	i_vector(&matbeg,numrows,"open_cplex:4");
	sense[0]='G';      
	matbeg[0] = 0;
	for(k=0;k<num_feasact;k++){	
		rhs[0]=1;
		if(k<0){
			status = CPXaddrows (env, lp, 0, 1, feascuts[feasact[k]].numcol, rhs, sense, matbeg, feascuts[feasact[k]].ind, feascuts[feasact[k]].coeff, NULL, NULL);
		}
		else{
			status=CPXaddusercuts (env,lp,1, feascuts[feasact[k]].numcol, rhs, sense, matbeg, feascuts[feasact[k]].ind, feascuts[feasact[k]].coeff,NULL);
			status=CPXaddlazyconstraints (env,lp,1, feascuts[feasact[k]].numcol, rhs, sense, matbeg, feascuts[feasact[k]].ind, feascuts[feasact[k]].coeff, NULL);
		}
		if( status ) fprintf (stderr,"CPXaddrows failed.\n");    
	}
	/********************************Next the optimality cuts that are active**************************************************/
	for(k=0;k<num_optact;k++){	
		rhs[0]=optcuts[optact[k]].RHS;
		if(k<0){
			status = CPXaddrows (env, lp, 0, 1, optcuts[optact[k]].numnz, rhs, sense, matbeg, optcuts[optact[k]].ind, optcuts[optact[k]].coeff, NULL, NULL);
		}
		else{
			status=CPXaddusercuts (env,lp,1, optcuts[optact[k]].numnz, rhs, sense, matbeg, optcuts[optact[k]].ind, optcuts[optact[k]].coeff,NULL);
			status=CPXaddlazyconstraints (env,lp,1, optcuts[optact[k]].numnz, rhs, sense, matbeg, optcuts[optact[k]].ind, optcuts[optact[k]].coeff, NULL);
		}
		if( status ) fprintf (stderr,"CPXaddrows failed.\n"); 
	}

	printf("Added %d initial feasibility and %d optimality cuts\n", feascount, optcount);


	/*****************************New set of restrictions not exploring parts previously explored*****************************************************/
	d_vector(&matval,M,"open_cplex:6");
	d_vector(&rhs,1,"open_cplex:2");
	c_vector(&sense,1,"open_cplex:3");
	i_vector(&matbeg,1,"open_cplex:4");
	sense[0]='G';
	matbeg[0] = 0;
	for(j=0;j<M;j++) matval[j] = 1;
	for(i=0;i<globalcounter;i++){
		rhs[0]=(i+1)+1;
		if(indic==0){
			status =  CPXaddusercuts (env, lp, 1, base_soltrackind[i], rhs, sense, matbeg, base_soltrack[i], matval, NULL);
			status = CPXaddlazyconstraints (env, lp, 1, base_soltrackind[i], rhs, sense, matbeg, base_soltrack[i], matval, NULL);
		}
		else status = CPXaddrows(env, lp, 0, 1, base_soltrackind[i]+1, rhs, sense, matbeg, base_soltrack[i], matval, NULL, NULL);
	}
	free(matbeg);
	free(matval);
	free(sense);
	free(rhs);
	/**********************************************************************************/


	if(indic==0){				//Writing the restriction of variables not in the base solution solution to be 0
		numrows = 1;
		numnz = M;		
		d_vector(&rhs,numrows,"open_cplex:2");
		c_vector(&sense,numrows,"open_cplex:3");
		i_vector(&matbeg,numrows,"open_cplex:4");
		i_vector(&matind,numnz,"open_cplex:6");
		d_vector(&matval,numnz,"open_cplex:7");
		index = 0;//columns
		index1 = 0;//rows 	
		sense[index1]='L';      
		matbeg[index1] = index;
		if(Up_Extlim>0) rhs[index1++]=globalcounter+1+Up_Extlim-1;
		else rhs[index1++]=globalcounter+1;
		//printf("RHS in < Ext is %lf\n",rhs[0]);
		//getchar();
		for(j=0;j<M;j++){
			if(base_sol[j]==0){ //if node is origin
					 matind[index] = pos_y[edges[j].i][edges[j].j];
					 matval[index] = 1;			 
					 index++;
			}
		}
		status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
		if(Up_lim<size_sol-2){
			index = 0;//columns
			for(j=0;j<M;j++){
				if(base_sol[j]==1){ //if node is origin
						 matind[index] = pos_y[edges[j].i][edges[j].j];
						 matval[index] = 1;			 
						 index++;
				}
			}
			rhs[0]=Up_lim;
			status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
			//printf("RHS in < Int is %lf\n",rhs[0]);
			//getchar();
		}
		if(Low_lim>1){
			sense[0]='G'; 
			index = 0;//columns
			for(j=0;j<M;j++){
				if(base_sol[j]==1){ //if node is origin
						 matind[index] = pos_y[edges[j].i][edges[j].j];
						 matval[index] = 1;			 
						 index++;
				}
			}
			rhs[0]=Low_lim;
			status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
			//printf("RHS in > Int is %lf\n",rhs[0]);
			//getchar();
		}
		free(matbeg);
		free(matind);
		free(matval);
		free(sense);
		free(rhs);
	}
	else{					//Checking if there are any better solutions outside
		numrows = 1;
		numnz = M;		
		d_vector(&rhs,numrows,"open_cplex:2");
		c_vector(&sense,numrows,"open_cplex:3");
		i_vector(&matbeg,numrows,"open_cplex:4");
		i_vector(&matind,numnz,"open_cplex:6");
		d_vector(&matval,numnz,"open_cplex:7");
		index = 0;//columns
		index1 = 0;//rows 	
		sense[index1]='G';      
		matbeg[index1] = index;
		rhs[index1++]=(globalcounter+1)+1;
		for(j=0;j<M;j++){				
			if(base_sol[j]==0){ //if node is origin
					 matind[index] = pos_y[edges[j].i][edges[j].j];					 		 
					 matval[index] = 1;
					 index++;
					 
			}
		}
		status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
		free(matbeg);
		free(matind);
		free(matval);
		free(sense);
		free(rhs);
	}


	//printf("CutUP value is %lf\n",bestobj);
	CPXsetdblparam(env,CPX_PARAM_CUTUP,bestobj); // provide an initial upper bound
	CPXsetintparam(env,CPX_PARAM_VARSEL, CPX_VARSEL_PSEUDOREDUCED); //Selecting the variable to branch on
	CPXchgobjsen(env, lp,CPX_MIN);
	CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use
	CPXsetintparam(env,CPX_PARAM_MIPDISPLAY,3); //different levels of output display
	if(indic==1) CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,3);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfea
	if(indic==0) CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,2);//0:balanced; 1:feasibility; 2:optimality,3:bestbound, 4:hiddenfea
	CPXsetdblparam(env,CPX_PARAM_TILIM,Time-cputimeBD); // time limit

	/*Changing some of the precision*/
	//CPXsetintparam(env,CPX_PARAM_NUMERICALEMPHASIS,1); //Numerical precision of the time
	//CPXsetdblparam(env,CPX_PARAM_RINSHEUR,10.0); // e-optimal solution (%gap) 
	CPXsetdblparam(env,CPX_PARAM_EPGAP, 0.00000000000001); // e-optimal solution (%gap)
    //CPXsetdblparam(env,CPX_PARAM_EPAGAP, 0.0000000001); // e-optimal solution (absolute value)    	
	//CPXsetdblparam(env,CPX_PARAM_EPRHS, 0.0000001);// Feasibility tolerance
	//CPXsetdblparam(env,CPX_PARAM_EPINT, 0.0001); // integer precision


	//CPXsetintparam(env,CPX_PARAM_PREIND,0);
	CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON); //output display
	//CPXsetdblparam(env,CPX_PARAM_CUTSFACTOR, 1.0);
    //CPXsetintparam(env,CPX_PARAM_INTSOLLIM,1);    //stops after finding first integer sol.    
    //CPXsetdblparam(env,CPX_PARAM_TRELIM, 14000); // B&B memory limit    
	//CPXsetintparam(env,CPX_PARAM_REDUCE, 0);  // only needed when adding lazy constraints
    //CPXsetintparam(env,CPX_PARAM_HEURFREQ, -1); //heuristic frequency and intensisty 	
	
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_OPTIMALITY);  // MIP emphasis: optimality, feasibility, moving best bound
	//CPXsetintparam(env,CPX_PARAM_PARALLELMODE, 1); 
	//CPXsetintparam(env,CPX_PARAM_PREIND,0);
	//CPXsetintparam(env,CPX_PARAM_MIPORDIND,CPX_ON); // Turn on or off the use of priorities on bracnhing variables
	//CPXsetintparam(env,CPX_PARAM_MIPEMPHASIS,CPX_MIPEMPHASIS_BESTBOUND);  // MIP emphasis: optimality, feasibility, moving best bound   
	
	/*Additional parameter settings for using the branch and cut*/
	/************************************************/
	status = CPXsetintparam (env, CPX_PARAM_PRELINEAR, 0); 											 /* Do not use presolve */
	status = CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL); 					/* Turn on traditional search for use with control callbacks */
	status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);									 /* Let MIP callbacks work on the original model */   
   /************************************************/

	/* Code to use lazy constraints*/
	/************************************************/
	status = makelazyconstraint (env, lp, &lazyconinfo); 											//Linking the LP defined as the LP in the lazyconinfo structure and get sizes
	status = CPXsetlazyconstraintcallbackfunc (env, lazy_feas_callback1, &lazyconinfo); 
	/* Code to use usercut constraints*/
	/************************************************/
	status = makeusercuts (env, lp, &usercutinfo); 												//Linking the LP defined as the LP in the usercutinfo structure and get sizes
	status = CPXsetusercutcallbackfunc (env, user_cut_callback1 , &usercutinfo); 							//Now create the user cuts using our cutcallback which solves a separation problem
   /**********************************************************************************************************/
	/* Code to use Heuristic Callback*/
	/************************************************/
	/*if(globalcounter>0)*/status = CPXsetheuristiccallbackfunc (env, LagrangeHeur, NULL);
	/**********************************************************************************************************/
	printf("Finished declaring the problem\n");
	//getchar();
	if(Flag_Print==1){printf("Started the clock\n"); /*getchar();*/}
	CPXmipopt(env,lp);  //solve the integer program
	status=CPXgetstat(env,lp);
	/*************************************************************************************/
	if(status==101|| status==102){
		end_time = clock();
		printf("Solution found with ");
		CPXgetmipobjval(env,lp,&Upper_bound);
		if(indic==0) bestobj=Upper_bound;
		printf("Optimal Value: %.2f   \n",Upper_bound);
		printf("Added %d final feasibility and %d final optimality cuts\n", feascount, optcount);
		CPXgetbestobjval(env,lp,&best_lower_bound);
		nodecount=CPXgetnodecnt(env,lp);
		cputimeBD = (double)(end_time - start_time) / CLOCKS_PER_SEC;

		/*********************Printing*************/
		out = Open_File(Progress,"a+");
		if(globalcounter>=1)	fprintf(out,"; ; %d; %d; %d; %lf; %lf; %d; %d\n",feascount, optcount,nodecount,Upper_bound,cputimeBD, countimprovedheur, counttried);
		else fprintf(out,"%d; %d; %d; %lf; %lf;%d;%d\n",feascount, optcount,nodecount,Upper_bound,cputimeBD,countimprovedheur,counttried);
		fclose(out);
		/******************************************/
		if(indic==0){
			numcols = CPXgetnumcols (env, lp);
			d_vector(&x,numcols,"open_cplex:0");
			CPXgetmipx(env,lp,x,0, numcols-1);  // obtain the values of the decision variables
			//out = Open_File("CompNetworks.txt","a+");
			//fprintf(out,"%s; Optimal;",instance);
			for(j=0;j<M;j++){                  //fill up the values of ybar
			  if( x[pos_y[edges[j].i][edges[j].j]] > 0.9){
				//ybar[edges[j].i][edges[j].j] = 1;
				base_sol[j]=1;
				//core_point[edges[j].i][edges[j].j]=1;
				bestsol[j]=1;
				//fprintf(out,"%d; ", ybar[edges[j].i][edges[j].j]);
			  }
			  else{
				  //ybar[edges[j].i][edges[j].j] = 0;
				//core_point[edges[j].i][edges[j].j]=0;
				base_sol[j]=0;
				bestsol[j]=0;
				base_soltrack[globalcounter+1][base_soltrackind[globalcounter+1]]=j;  //Indices of variables not in solution
				base_soltrackind[globalcounter+1]++;
				//fprintf(out,"; ");
			  }
			}
			for(j=0;j<K;j++){
				base_solroute[j]=x[pos_art[j]];
				bestsol[M+j]=x[pos_art[j]];
			}
			//fprintf(out,"\n");
			//fclose(out);
			free(x);
		}
		/*********************************Printing the Solution******************************/
		if(Flag_PrintSol==1){
			for(i=0;i<M;i++){
						out = Open_File(Solution,"a+");					//Writing what instance we are solving
						fprintf(out,"%d; ", ybar[edges[i].i][edges[i].j]);
						fclose(out);
			}
			out = Open_File(Solution,"a+");					//Writing what instance we are solving
			fprintf(out,"\n");
			fclose(out);
		}
		
	}
	else if(status==103){
		end_time = clock();
		cputimeBD = (double)(end_time - start_time) / CLOCKS_PER_SEC;
		printf("Master Problem infeasible\n");
		out = Open_File(Progress,"a+");
		fprintf(out,"; ; %d; %d; %d; %lf; %lf; %d; %d\n",feascount, optcount,-1,Upper_bound,cputimeBD,countimprovedheur,counttried);
		fclose(out);
		//getchar();
	}
	else if(status==107|| status==108){
		printf("Time limit reached\n");
		end_time = clock();		
		CPXgetmipobjval(env,lp,&Upper_bound); 
		CPXgetbestobjval(env,lp,&best_lower_bound);
		nodecount=CPXgetnodecnt(env,lp);
		//printf("Going to print on %s\n",outfile);
		out = Open_File(Progress,"a+");
		fprintf(out,"%d; ",nodecount);
		fclose(out);
	}
	else{
		printf("Unknown stopping criterion (%d) \n",status);
	}
	

	free(average);
	if ( lp != NULL ) {
			  status = CPXfreeprob (env, &lp);
			  if ( status ) {
				fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
			  }
	}
	if ( env != NULL ) {
			  status = CPXcloseCPLEX (&env);
			  if ( status ) {
				char  errmsg[1024];
				fprintf (stderr, "Could not close CPLEX environment.\n");
				CPXgeterrorstring (env, status, errmsg);
				fprintf (stderr, "%s", errmsg);
			  }
	}
	return Upper_bound;
}