#include "def.h"

double solve_net_sub(double *temp_x){
	int		solstat, status = 0, inf_indicator; //last one is the infeasibility indicator
	int		i, k, j,index;
	double	obj_primal_sub, temp, value,x0,x0prime, lagrangperK;
	double	  sol_est; //Estimate of the routing cost.

	double *supply2;

	// Allocate memory for CPLEX PROBLEM DEFINITION
	CPXENVptr env2 = NULL;
	CPXNETptr net = NULL;
	
	// AND OUTPUT:
	double   *beta;    
	double   *alpha;
	double   acumub;
	double   *netX;
	double   *ubnew;
	double	 *ubnew2;
	double   *solnew;
	int		fakecount;
	///////////////////////////////////////////////////////////////////////////
	/****************************************************************************************/
	//Initialize the variables for everything
	vioptcuts=0;
	Lagrangeind=0;	
	if(globaldepth==0 && flag_frac==1){
		Update_core_point();
		alter=4;
	}
	if(globaldepth>0 && flag_frac==1){
		//Update_core_point();
		alter=4;
		//printf("Changed the corepoint and fixed to alter =%d\n", alter);
	}
	else{
		alter=4;
	}
	x0=0;	
	if(alter!=5){                                                        //If we are going to solve for usual M&W
		for(j=0;j<M;j++) {
			x0=x0+core_point[edges[j].i][edges[j].j];
			lagrangemult[j]=0;											//Initializing the variable we will use for the lagrangean heuristic
		}
	}
	//printf("x0=%lf\n", x0); getchar();
	//x0=MW_complete;
	sol_est=0;
	if(flag_frac==1){
		i=0;
	}
	/// OPEN CPLEX ////////////////////////////////////////////////////////////
	env2 = CPXopenCPLEX(&status);
	if ( env2 == NULL ) {
		char  errmsg[1024];
		printf ("Could not open CPLEX. \n");
		CPXgeterrorstring (env2, status, errmsg);
		printf ("%s", errmsg);
	}

	net = CPXNETcreateprob(env2, &status, "PrimalSubproblem");
	///////////////////////////////////////////////////////////////////////////
	for(j=0;j<M;j++){ 
		if(Flag_MaxCutPrint==1){
			if(temp_x[j]>0.5) ub[j]=x0+core_point[edges[j].i][edges[j].j];
			else ub[j]=core_point[edges[j].i][edges[j].j];
		}
		else{
			//Or else use the solution obtained from the MP
			if(flag_frac==0){
				if(ybar[edges[j].i][edges[j].j]>0.9) ub[j] = core_point[edges[j].i][edges[j].j]+x0;
				else ub[j] = core_point[edges[j].i][edges[j].j];
			}
			else{
				/*if(ybar_frac[edges[j].i][edges[j].j]>0.8) ub[j] = core_point[edges[j].i][edges[j].j]+x0;
				else if(ybar_frac[edges[j].i][edges[j].j]>0.5) ub[j] = inf_indicator*0.7 +(1-inf_indicator)*1.25*core_point[edges[j].i][edges[j].j]+x0*ybar_frac[edges[j].i][edges[j].j];
				else ub[j] = inf_indicator*0.5 +(1-inf_indicator)*2*core_point[edges[j].i][edges[j].j]+x0*ybar_frac[edges[j].i][edges[j].j];*/
				ub[j]=core_point[edges[j].i][edges[j].j]+x0*ybar_frac[edges[j].i][edges[j].j];
			}
		}
	}
	
	start_sep=clock();
	/// SOLVE K Subproblems ///////////////////////////////////////////////////
	k=0;
	inf_indicator=0;
	beta= create_double_vector(N);
	alpha=create_double_vector(M);
	netX=create_double_vector(M);
	supply2=create_double_vector(N);
	ubnew=create_double_vector(M);
	ubnew2=create_double_vector(M);
	solnew=create_double_vector(M);
	while(k<K) {
		for(j=0;j<M;j++){ 
			obj[j] =com[k].d*edges[j].c[k];
			if(obj[j]<0)getchar();
		}
		// Define the Supply:
		if(Flag_MaxCutPrint==1){
			for(i=0;i<N;i++){ supply[i]=0; supply2[i]=0;}
			supply[com[k].i] =  1.0+x0;
			supply2[com[k].i] =  1.0;
			supply[com[k].j] = -1.0-x0;
			supply2[com[k].j] =  -1.0;
		}
		else{
			for(i=0;i<N;i++){ supply[i]=0; supply2[i]=0;}
			supply[com[k].i] =  1.0+x0;
			supply2[com[k].i] =  1.0;
			supply[com[k].j] = -1.0-x0;
			supply2[com[k].j] =  -1.0;
		}

		// Delete existing network (if any). 
		if (CPXNETgetnumnodes(env2, net) > 0) {
			status = CPXNETdelnodes (env2, net, 0, CPXNETgetnumnodes (env2, net)-1);

			status= CPXNETdelarcs (env2, net, 0, CPXNETgetnumnodes (env2, net)-1);
		}
		// Set optimization sense:
		status = CPXNETchgobjsen(env2, net, CPX_MIN);
		// Add nodes to network along with their supply values:
		status = CPXNETaddnodes (env2, net, N, supply, NULL);
		// Add arcs to network along with their objective values and bounds:
		status = CPXNETaddarcs(env2, net, M, tail, head, lb, ub, obj, NULL);			
		// Optimize the problem and obtain solution:
		status = CPXNETprimopt(env2, net);
		// Get solution data:
		status = CPXNETsolution(env2, net, &solstat, &obj_primal_sub, netX, beta, NULL, alpha);


		if (solstat != CPX_STAT_OPTIMAL) { 
			//printf("Warning: infeasible! in the method of %d obtaining status of %d for commodity %d and fractional flag of %d\n", alter, solstat, k, flag_frac); /*getchar();*/
			inf_indicator++;
			Flag_Check++;
			if(inf_indicator>1){								//If we had already made adjustments and the subproblem is still infeasible, then go to the next commodity
				inf_indicator=1;
				k++;
			}
		} 
		else {												//We will start writing the cut
			/************Populating the paths structure for the CG heuristic****************/
			//First we are writing the new subnetwork we will solve over********/
			fakecount=0;
			//printf("Commodity (%d,%d) routed in M-W as:\n", com[k].i,com[k].j);
			for(i=0;i<M;i++){
				if(netX[i]>1.0001){
					ubnew[i]=1.0;
					ubnew2[i]=1.0;
					//printf("x(%d,%d)=%lf\n",tail[i],head[i],netX[i]);
					fakecount++;
				}
				else{
					ubnew2[i]=0.0;    //ubnew2 is the more restricted version where there are only some of them
					if(netX[i]>0.1){
						ubnew[i]=1.0;
					}
					else{
						ubnew[i]=0.0;
					}
				}
			}
			//printf("In total there were %d out of %d arcs greater than 0\n/**********/\n",fakecount, M);*/
			/*******First solve the less restrictive**************/
			status = CPXNETdelnodes (env2, net, 0, CPXNETgetnumnodes (env2, net)-1);
			status= CPXNETdelarcs (env2, net, 0, CPXNETgetnumnodes (env2, net)-1);
			// Set optimization sense:
			status = CPXNETchgobjsen(env2, net, CPX_MIN);
			// Add nodes to network along with their supply values:
			status = CPXNETaddnodes (env2, net, N, supply2, NULL);
			// Add arcs to network along with their objective values and bounds:
			status = CPXNETaddarcs(env2, net, M, tail, head, lb, ubnew, obj, NULL);			
			// Optimize the problem and obtain solution:
			status = CPXNETprimopt(env2, net);
			// Get solution data:
			status = CPXNETsolution(env2, net, &solstat, &obj_primal_sub, solnew, NULL, NULL, NULL);			
			if(solstat==CPX_STAT_OPTIMAL){
				/*fakecount=0;
				printf("Solution status is %d\nCommodity (%d,%d) routed in Dijkstra as:\n", solstat, com[k].i,com[k].j);
				for(i=0;i<M;i++){
					if(solnew[i]>0.0001){
						printf("x(%d,%d)=%lf\n",tail[i],head[i],solnew[i]);
						fakecount++;
					}
				}
				printf("In total there were %d out of %d arcs greater than 0\n\*****",fakecount, M);*/
				if(useheur==1)populate_pathperkSS(solnew, k);
			}					
			/****************************************************/
			
			/*******Now solve the more restrictive***************/
			status = CPXNETdelnodes (env2, net, 0, CPXNETgetnumnodes (env2, net)-1);
			status= CPXNETdelarcs (env2, net, 0, CPXNETgetnumnodes (env2, net)-1);
			// Add arcs to network along with their objective values and bounds:
			status = CPXNETaddnodes (env2, net, N, supply2, NULL);
			status = CPXNETaddarcs(env2, net, M, tail, head, lb, ubnew2, obj, NULL);			
			// Optimize the problem and obtain solution:
			status = CPXNETprimopt(env2, net);
			// Get solution data:
			status = CPXNETsolution(env2, net, &solstat, &obj_primal_sub, solnew, NULL, NULL, NULL);	
			if(solstat==CPX_STAT_OPTIMAL){
				/*fakecount=0;
				printf("Solution status is %d\nCommodity (%d,%d) routed in Dijkstra as:\n", solstat, com[k].i,com[k].j);
				for(i=0;i<M;i++){
					if(solnew[i]>0.0001){
						printf("x(%d,%d)=%lf\n",tail[i],head[i],solnew[i]);
						fakecount++;
					}
				}
				printf("In total there were %d out of %d arcs greater than 0\n",fakecount, M);*/
				if(useheur==1)populate_pathperkSS(solnew, k);
			}
			//getchar();
			/****************************************************/

			optcuts[optcount].RHS = -beta[com[k].j]+beta[com[k].i];
			for(i=0;i<N;i++){
				perturbnode[i][k]=beta[i];
				//printf("perturbnode[%d][%d]=%lf\n", i,k,perturbnode[i][k]);
			}
			value=optcuts[optcount].RHS;
			index=0;
			for (j=0; j<M; j++) {
				temp=MAX(0, -alpha[j]);	
				perturbarc[j][k]=temp;			//Possibly change to alpha directly
				//printf("perturbarc[%d][%d][%d]=%lf\n", edges[j].i,edges[j].j,k,perturbarc[j][k]);getchar();
				if(temp>0.1){
					optcuts[optcount].coeff[index]=temp;
					optcuts[optcount].ind[index]=pos_y[edges[j].i][edges[j].j];
					if(flag_frac==0)value=value-temp*ybar[edges[j].i][edges[j].j];
					else if(flag_frac==1)value=value-temp*ybar_frac[edges[j].i][edges[j].j];
					else if (Flag_MaxCutPrint==1)value=value-temp*temp_x[pos_y[edges[j].i][edges[j].j]];
					index++;
				}				
			}
			optcuts[optcount].coeff[index]=1;
			optcuts[optcount].ind[index]=pos_art[k];
			index++;
			optcuts[optcount].numnz=index;
			if(value-temp_x[pos_art[k]]>violation && flag_keepcut==1){						//Cut is violated so add to the cutpool
				if(optcount==optlim*M*K-1){
					//getchar();
					optlim++;
					optcuts= (CUTPOOL *) realloc(optcuts,optlim*M*K*sizeof *optcuts); 
				}
				vioptcuts++;
				optcount++;
				optcuts[optcount].coeff=create_double_vector(M+K+1);
				optcuts[optcount].ind=create_int_vector(M+K+1);				
			}
			k++;
			inf_indicator=0;
		}
		
	}
	end_sep=clock();
	timeSep=timeSep+(double)(end_sep - start_sep) / CLOCKS_PER_SEC;
	free(beta);
	free(alpha);
	free(netX);
	free(supply2);
	free(ubnew);
	free(ubnew2);
	free(solnew);
TERMINATE:
	status = CPXNETfreeprob(env2, &net);
	status = CPXcloseCPLEX(&env2);
	return sol_est;
}
double solve_net_exact(double *temp_x){
	int		solstat, status = 0, inf_indicator; //last one is the infeasibility indicator
	int		i, k, j,index;
	double	obj_primal_sub, temp, value,x0,x0prime, lagrangperK;
	double	  sol_est; //Estimate of the routing cost.

	// Allocate memory for CPLEX PROBLEM DEFINITION
	CPXENVptr env2 = NULL;
	CPXNETptr net = NULL;
	double	 *supply = create_double_vector(N);
	int		 *tail	 = create_int_vector(M);
	int	     *head	 = create_int_vector(M);
	double	 *obj	 = create_double_vector(M);
	double	 *ub	 = create_double_vector(M);
	double	 *lb	 = create_double_vector(M);
	// AND OUTPUT:
	double   *beta;    
	double   *alpha;
	double   acumub;
	///////////////////////////////////////////////////////////////////////////
	/****************************************************************************************/
	//Initialize the variables for everything
	if(Flag_MaxCutPrint==1){
		alter=1; //Are we preprocessing
	}
	vioptcuts=0;
	/// OPEN CPLEX ////////////////////////////////////////////////////////////
	env2 = CPXopenCPLEX(&status);
	if ( env2 == NULL ) {
		char  errmsg[1024];
		printf ("Could not open CPLEX. \n");
		CPXgeterrorstring (env2, status, errmsg);
		printf ("%s", errmsg);
	}

	net = CPXNETcreateprob(env2, &status, "PrimalSubproblem");
	///////////////////////////////////////////////////////////////////////////
	
	// Define network structure: //////////////////////////////////////////////
	for (j=0; j<M; j++) {
		tail[j] = edges[j].i;
		head[j] = edges[j].j;
		lb[j] = 0.0;
	}
	///////////////////////////////////////////////////////////////////////////

	start_sep=clock();
	/// SOLVE K Subproblems ///////////////////////////////////////////////////
	k=0;
	inf_indicator=0;
	beta= create_double_vector(N);
	alpha=create_double_vector(M);
	while(k<K) {
		for(j=0;j<M;j++){ 
			obj[j] =com[k].d*edges[j].c[k];
			if(obj[j]<0)getchar();		
			//Or else use the solution obtained from the MP
			if(flag_frac==0){
				if(ybar[edges[j].i][edges[j].j]>0.9) ub[j] = 1;
				else ub[j] = 0;
			}
			else{
				ub[j] = ybar_frac[edges[j].i][edges[j].j];
			}
		}
		// Define the Supply:
		for(i=0;i<N;i++) supply[i]=0;
		supply[com[k].i] =  1.0;
		supply[com[k].j] = -1.0;

		// Delete existing network (if any). 
		if (CPXNETgetnumnodes(env2, net) > 0) {
			status = CPXNETdelnodes (env2, net, 0, CPXNETgetnumnodes (env2, net)-1);
			status= CPXNETdelarcs (env2, net, 0, CPXNETgetnumnodes (env2, net)-1);
		}
		// Set optimization sense:
		status = CPXNETchgobjsen(env2, net, CPX_MIN);
		

		// Add nodes to network along with their supply values:
		status = CPXNETaddnodes (env2, net, N, supply, NULL);

		// Add arcs to network along with their objective values and bounds:
		status = CPXNETaddarcs(env2, net, M, tail, head, lb, ub, obj, NULL);
			
		// Optimize the problem and obtain solution:
		status = CPXNETprimopt(env2, net);

		// Get solution data:
		status = CPXNETsolution(env2, net, &solstat, &obj_primal_sub, NULL, beta, NULL, alpha);
		if (solstat != CPX_STAT_OPTIMAL) { 
			//printf("Warning: infeasible! in the method of %d obtaining status of %d for commodity %d and fractional flag of %d\n", alter, solstat, k, flag_frac); /*getchar();*/
			//getchar();
			inf_indicator++;
			Flag_Check++;
			if(inf_indicator>1){								//If we had already made adjustments and the subproblem is still infeasible, then go to the next commodity
				inf_indicator=1;
				k++;
			}
		} 
		else {												//We will start writing the cut
			optcuts[optcount].RHS = -beta[com[k].j]+beta[com[k].i];
			value=optcuts[optcount].RHS;
			for(i=0;i<N;i++){
				perturbnode[i][k]=beta[i];
				//printf("perturbnode[%d][%d]=%lf\n", i,k,perturbnode[i][k]);
			}			
			index=0;
			for (j=0; j<M; j++) {
				temp=MAX(0, -alpha[j]);	
				perturbarc[j][k]=temp;
				if(ABS(temp)>0.0001){
					optcuts[optcount].coeff[index]=temp;
					optcuts[optcount].ind[index]=pos_y[edges[j].i][edges[j].j];
					if(flag_frac==0)value=value-temp*ybar[edges[j].i][edges[j].j];
					else if(flag_frac==1)value=value-temp*ybar_frac[edges[j].i][edges[j].j];
					else if (Flag_MaxCutPrint==1)value=value-temp*temp_x[pos_y[edges[j].i][edges[j].j]];
					index++;
				}				
			}
			optcuts[optcount].coeff[index]=1;
			optcuts[optcount].ind[index]=pos_art[k];
			index++;
			optcuts[optcount].numnz=index;
			if(value-temp_x[pos_art[k]]>0.000001 || Flag_MaxCutPrint==1){						//Cut is violated so add to the cutpool
				if(optcount==optlim*M*K-1){
					//getchar();
					optlim++;
					optcuts= (CUTPOOL *) realloc(optcuts,optlim*M*K*sizeof *optcuts); 
				}
				vioptcuts++;
				optcount++;
				optcuts[optcount].coeff=create_double_vector(M+K+1);
				optcuts[optcount].ind=create_int_vector(M+K+1);
			}
			k++;
			inf_indicator=0;
		}
		
	}
	end_sep=clock();
	timeSep=timeSep+(double)(end_sep - start_sep) / CLOCKS_PER_SEC;
	free(beta);
	free(alpha);
TERMINATE:
	status = CPXNETfreeprob(env2, &net);
	status = CPXcloseCPLEX(&env2);
	free(supply);
	free(tail);
	free(head);
	free(obj);
	free(ub);
	free(lb);
	return value;
}

double Dijkstra ( int start, int *t, int **index, int *rep, int dim_dest, double *lambda){
	int current, u,v;
	int i,j,k,reps, count, lowindex;
	double newdist, temp,tot_dist;
	double		*distance;
	tot_dist=0;
	for(k=0;k<dim_dest;k++){		
		for(reps=0;reps<rep[k];reps++){
			current=start;
			count=1;
			list=create_int_vector(N);
			preced=create_int_vector(N);
			distance=create_double_vector(N);
			for(i=0;i<N;i++){
				distance[i]=MAX_double;
				list[i]=0;
			}
			preced[start]=0;
			distance[start]=0;
			list[start]=-1;
			while (count<N)
			{
				for(j=0;j<M;j++){//Check all the edges and update the value of predecesor and distance
					if(edges[j].i==current && (ybar[edges[j].i][edges[j].j]>0 || presolve_ind==1)){
						newdist=distance[current]+edges[j].c[index[k][reps]];
						if(newdist<distance[edges[j].j]){//add end node to candidate list
							distance[edges[j].j]=newdist;
							preced[edges[j].j]=current;
							list[edges[j].j]=1;
						}
					}
				}
				temp=MAX_double;
				for(i=0;i<N;i++){
					if(list[i]==1 && distance[i]<temp){
						temp=distance[i];//update smallest distance
						lowindex=i;//index of the following current
					}
				}
				list[lowindex]=-1;
				current=lowindex;
				count=count+1;
			}
			//Calculating the total routing distance				
			if(presolve_ind==0){
				tot_dist=tot_dist+com[index[k][reps]].d*distance[t[k]];
				lambda[index[k][reps]]=com[index[k][reps]].d*distance[t[k]];
			}
			else{
				shortest[index[k][reps]].costwdem=com[index[k][reps]].d*distance[t[k]];
				shortest[index[k][reps]].cost=distance[t[k]];
				lowest_routing=lowest_routing+shortest[index[k][reps]].costwdem;
			}
			v = t[k];
            u = preced[t[k]];
			//Here we are selecting which ones we had used
            while (v!=start) {
                    if(presolve_ind==1){
						shortest[index[k][reps]].path[shortest[index[k][reps]].len]=v;
						shortest[index[k][reps]].len++;
					}
					usedarc[u][v]++;
					v = u;
                    u = preced[v];					
            }
			if(presolve_ind==1){
				shortest[index[k][reps]].path[shortest[index[k][reps]].len]=start;
				shortest[index[k][reps]].len++;
			}
			/*printf("Routing commodity (%d,%d) which has distance %lf\n", com[index[k][reps]].i,com[index[k][reps]].j , distance[t[k]]);
			getchar();*/
			free(distance);
			free(preced);
			free(list);
		}
	}
	
	return tot_dist;
}



void StrongFeas_com(int commodity){
	int i,j,k;
	int index,index1;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound, best_lower_bound;
	int nodecount;
	//Variables to call cplex
	CPXLPptr  sublp;      // data strucutre to store a problem in cplex ...................
	CPXENVptr subenv;     // cplex environment.............................................
	int       numcols; // number of variables ..........................................
	int       numrows; // number of constraints.........................................
	int       numnz;   // number of non-zero elements in the matrix ....................
	int       objsen;  // optimization sense (min:1, max:-1 ) ..........................
	double    *obj;    // objective function coefficients ..............................
	double    *rhs;    // right hand side of constraints ................................
	double    *x;
	char      *sense;  // constraints sense (<=: 'L', =:'E', >=:'G') ...................
	int       *matbeg; // index of first non-zero element in each row...................
	int       *matind; // associated column of each non-zero element ...................
	double    *matval; // coefficient values for the non-zero elements of constraints....
	double    *lb;     // lower bounds of variables.....................................
	double    *ub;     // upper bounds of variables.....................................
	int       status;  // optimization status......................... .................
	char      probname[16]; // problem name for cplex .......................................
	char      *ctype;  // variable type ('C', 'I', 'B') only if integer.................
	double    value;   // objevtive value of solution ..................................
	int       num_lambda_var, num_lambda2_var, num_mu_var;
	int       *pos_lambda_plus;
	int       *pos_lambda_minus;
	int       opt_feas=0;
	int       stats;
	double	  right;
	/*char      its[14];
	char      strs[20];*/

	pos_lambda_plus = create_int_vector(N);//Assign them to the size of the Arc set
	pos_lambda_minus = create_int_vector(N);//Assign them to the size of the Arc set
	
	
	//Initialize CPLEX environment
	subenv = CPXopenCPLEX (&status); 
	if ( subenv == NULL ) {
		char  errmsg[1024];
		printf ("Could not open CPLEX. \n");
		CPXgeterrorstring (subenv, status, errmsg);
		printf ("%s", errmsg);
	}

    // Create the problem in CPLEX 
	strcpy(probname,"DSP"); //Copy the name placed as the problem name
	sublp = CPXcreateprob (subenv, &status, probname);  // Initialize the pointer that signals the problem 
	if ( subenv == NULL ) {
		char  errmsg[1024];
		printf ("Could not create LP. \n");
		CPXgeterrorstring (subenv, status, errmsg);
		printf ("%s", errmsg);
	}
	CPXchgobjsen(subenv, sublp, CPX_MAX);

	                                        //Define continuous positive lambda plus variables for each i in N and k in K 
     // index of columns
	index1 = 0;  
	numcols = N;
	d_vector(&obj,numcols,"open_cplex:1"); //objective function value
	d_vector(&lb,numcols,"open_cplex:8");//lower bound
	d_vector(&ub,numcols,"open_cplex:9");//upper bound
	c_vector(&ctype,numcols,"open_cplex:01");//type of the variable
	for(i=0;i<N;i++){
		   pos_lambda_plus[i]= index1;
		   if(i==com[commodity].j){
				obj[index1] =1;				
		   }
		   else {
			   if(i==com[commodity].i)  obj[index1]=-1;
			   else obj[index1]=0;
			}
			ctype[index1] ='C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
	}
	
	status = CPXnewcols (subenv, sublp, index1, obj, lb, NULL, NULL, NULL);
	if( status ) fprintf (stderr,"CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_lambda_var = index1;//number of Z variables

	                                        //Define continuous positive lambda minus variables for each i in N and k in K 
     // index of columns
	index1 = 0;  
	numcols = N;
	d_vector(&obj,numcols,"open_cplex:1"); //objective function value
	d_vector(&lb,numcols,"open_cplex:8");//lower bound
	d_vector(&ub,numcols,"open_cplex:9");//upper bound
	c_vector(&ctype,numcols,"open_cplex:01");//type of the variable
	for(i=0;i<N;i++){
		   pos_lambda_minus[i] =num_lambda_var+ index1;
		   if(i==com[commodity].j){
				obj[index1] =-1;				
		   }
		   else {
			   if(i==com[commodity].i)  obj[index1]=1;
			   else obj[index1]=0;
			}
			ctype[index1] ='C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
	}
	status = CPXnewcols (subenv, sublp, index1, obj, lb, NULL, NULL, NULL);
	if( status ) fprintf (stderr,"CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_lambda2_var =num_lambda_var+ index1;//number of Z variables
	
	                                        //Define Continuous mu_ijk variables
    index1 = 0;  // index of columns
	numcols = M;
	d_vector(&obj,numcols,"open_cplex:1");
	d_vector(&lb,numcols,"open_cplex:8");
	d_vector(&ub,numcols,"open_cplex:9");
	c_vector(&ctype,numcols,"open_cplex:01");
	for(j=0;j<M;j++){
			pos_mu[edges[j].i][edges[j].j] =  num_lambda2_var + index1; // position of the additional variables
			if(flag_frac==0)obj[index1] = -ybar[edges[j].i][edges[j].j];
			else obj[index1] = -ybar_frac[edges[j].i][edges[j].j];
			ctype[index1] = 'C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
	}
    status = CPXnewcols (subenv, sublp, index1, obj, NULL, NULL, NULL, NULL);
	if( status ) 
      fprintf (stderr,"CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_mu_var = index1;//number of x variables

	//Add dual restriction
	numrows = M;
	numnz = 5*M;
	d_vector(&rhs,numrows,"open_cplex:2");
	c_vector(&sense,numrows,"open_cplex:3");
	i_vector(&matbeg,numrows,"open_cplex:4");
	i_vector(&matind,numnz,"open_cplex:6");
	d_vector(&matval,numnz,"open_cplex:7");
    index = 0;//columns
    index1 = 0;//rows
	for(j=0;j<M;j++){
			 sense[index1]='L';
			 rhs[index1]=0;
			 matbeg[index1] = index;			 
			 matind[index] = pos_mu[edges[j].i][edges[j].j];  //dual variable associated to arc allocation
			 matval[index++] = -1;
			 matind[index]=pos_lambda_plus[edges[j].i];//origen del arco
			 matval[index++]=-1;
			 matind[index]=pos_lambda_plus[edges[j].j];//destino del arco
			 matval[index++]=1;
			 matind[index]=pos_lambda_minus[edges[j].i];//origen del arco
			 matval[index++]=1;
			 matind[index]=pos_lambda_minus[edges[j].j];//destino del arco
			 matval[index++]=-1;
			 index1++;
	}
	status = CPXaddrows (subenv, sublp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
    if( status ) 
      fprintf (stderr,"CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);

	//Adding the normalizing constraint
	numrows = 1;
	numnz = 2*N+M;
	d_vector(&rhs,numrows,"open_cplex:2");
	c_vector(&sense,numrows,"open_cplex:3");
	i_vector(&matbeg,numrows,"open_cplex:4");
	i_vector(&matind,numnz,"open_cplex:6");
	d_vector(&matval,numnz,"open_cplex:7");

    index = 0;//columns
    index1 = 0;//rows
    rhs[index1]=1;
	sense[index1]='L';
	matbeg[index1] = index;
	for(j=0;j<M;j++){		 
			 matind[index] = pos_mu[edges[j].i][edges[j].j];  //dual variable associated to arc allocation
			 matval[index++] = 1;
		}
	for(i=0;i<N;i++){
			matind[index]=pos_lambda_plus[i];//origen del arco
			matval[index++]=1;
			matind[index]=pos_lambda_minus[i];//origen del arco
			matval[index++]=1;
	}	 
	index1++;
	status = CPXaddrows (subenv, sublp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
    if( status ) 
      fprintf (stderr,"CPXaddrows failed.\n");
	free(matbeg);
	free(matind);
	free(matval);
	free(sense);
	free(rhs);
	//start_feas=clock();
	CPXlpopt(subenv,sublp);  //solve the linear program	
    stats=CPXgetstat(subenv,sublp);
	if(stats!=1){printf("status is %d\n",stats);getchar();}
	//else printf("Solved over the recesion cone\n");
	if(stats==1){
		CPXgetobjval(subenv,sublp,&value);
		numcols = CPXgetnumcols (subenv, sublp);
		d_vector(&x,numcols,"open_cplex:0");
		CPXgetx(subenv,sublp,x,0, numcols-1);  // obtain the values of the decision variables	
		for(j=0;j<M;j++) coeff[j]=0;                                      //initialize coeff and right hand side																				  
		//Calculating the Right hand side constant of the cut.	
		Righthand=-x[pos_lambda_plus[com[commodity].i]]+x[pos_lambda_minus[com[commodity].i]]+x[pos_lambda_plus[com[commodity].j]]-x[pos_lambda_minus[com[commodity].j]]; 
		for(j=0;j<M;j++){                                                                              //Calculating the coefficients of the variables of the cuts
				coeff[j]=coeff[j] + x[pos_mu[edges[j].i][edges[j].j]]; 
		}					
		free(x);
	}
	if ( sublp != NULL ) {
      status = CPXfreeprob (subenv, &sublp);
      if ( status ) {
        fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
      }
    }
    if ( subenv != NULL ) {
      status = CPXcloseCPLEX (&subenv);
      if ( status ) {
        char  errmsg[1024];
        fprintf (stderr, "Could not close CPLEX environment.\n");
        CPXgeterrorstring (subenv, status, errmsg);
        fprintf (stderr, "%s", errmsg);
      }
    }	
	free(pos_lambda_plus);
	free(pos_lambda_minus);
}

static void free_and_null (char **ptr)
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
} /* END free_and_null */ 
