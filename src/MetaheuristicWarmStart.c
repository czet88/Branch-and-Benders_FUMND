#include "def.h"

void populate_pathperkSS(double *x_sol, int com_ind){  //this function will receive as input the commodity index (com_ind) and the method SP was calculates (method) and populate the path obtained for it into paths[com_ind]
	int j,v,u,tempnumarcs, inddif,Flagadd;
	double temp_dist, temp_fixed;
	int *temp_increase;
	temp_increase=create_int_vector(M);	
	paths[totalinpaths].arcind=create_int_vector(M);
	paths[totalinpaths].distance=0;
	tempnumarcs=0;
	temp_dist=0;
	temp_fixed=0;
	for(j=0;j<M;j++){
		if(x_sol[j]>0.1){
			paths[totalinpaths].arcind[tempnumarcs]=j;
			path_indM[com_ind*M+j].Kpath_ind[path_indM[com_ind*M+j].dim]=num_paths[com_ind];
			temp_increase[j]=1;
			temp_dist+=edges[j].c[com_ind]; //Add the cost;
			temp_fixed+=edges[j].f; //Add to the fixed cost
			tempnumarcs++;
			indarcinpath[j]=1;
		}
	}
	paths[totalinpaths].distance=temp_dist;
	paths[totalinpaths].fixed=temp_fixed;
	paths[totalinpaths].PP=com[com_ind].d*paths[totalinpaths].distance; //Multiplied by the commodity amount
	paths[totalinpaths].numarcs=tempnumarcs;
	/******Now we evaluate if it's a different value******/
	inddif=0;
	for(v=0;v<num_paths[com_ind];v++){
		if(ABS(paths[totalinpaths].PP-paths[locationinpaths[com_ind][v]].PP)>1)inddif++;
		else if(ABS(tempnumarcs-paths[locationinpaths[com_ind][v]].numarcs)>0.5) inddif++;
		else if (ABS(temp_fixed-paths[locationinpaths[com_ind][v]].fixed)>1) inddif++;
		/*else{
			printf("New distance=%lf\nOld distance=%lf\n\nNewnumarcs=%d\nOldnumarcs=%d\nNewfixed=%lf\nOldnewfixed=%lf\n\n",paths[totalinpaths].PP,paths[locationinpaths[com_ind][v]].PP,tempnumarcs, paths[locationinpaths[com_ind][v]].numarcs,temp_fixed, paths[locationinpaths[com_ind][v]].fixed) ;
			getchar();
		}*/
	}
	if(inddif==num_paths[com_ind]){ 
		locationinpaths[com_ind][num_paths[com_ind]]=totalinpaths;
		num_paths[com_ind]++;
		totalinpaths++;
		Flagadd=1;
	}
	else{
		
		Flagadd=0;
		free(paths[totalinpaths].arcind);
	}
	for(j=0;j<M;j++){//Populating the path_indM array.
		if(Flagadd==1){
			path_indM[com_ind*M+j].dim=path_indM[com_ind*M+j].dim+temp_increase[j];
		}			
		temp_increase[j]=0;
	}
	/***************************************************/

	free(temp_increase);
}

double master_p(double *lbtree, double *ubtree)
{
	//Variables to call cplex
	int i,j,k,s;
	int index,index1,acum;  // auxiliar indices to fill in the constraint matrix
	double best_upper_bound, best_lower_bound, tempor;
	int nodecount, iterations, flag_addcol, flag_stop, currentcols;
	//Variables to call cplex
	CPXLPptr  lpclone = NULL;//Cloned problem for the LP relaxation
	CPXLPptr  lp;      // data strucutre to store a problem in cplex ...................
	CPXENVptr env;     // cplex environment.............................................
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
	double    value, valueMIP;   // objevtive value of solution ..................................
	double	  *xx; //here we'll store the solution of the CG.....................................
	double		lastobj=MAX_double;
	int		reptitions=0;
	 
	
	/*****************************************Here we will declare the variables****************************************/
	int			*pos_arc;  // positions of the design variables (double, even if the problem is integer) .....
	int			*pos_tree;  // positions of the routing variables
	int       num_arc_var, num_tree_var;
	int		currentcolsPO;

	currentcolsPO=0;
	for(s=0;s<K;s++){
		currentcolsPO+=num_paths[s];
	}

	pos_arc=create_int_vector(M);
	pos_tree=create_int_vector(currentcolsPO);


	//Initialize CPLEX environment
	env = CPXopenCPLEX (&status); 
	if ( env == NULL ) {
		char  errmsg[1024];
		printf ("Could not open CPLEX. \n");
		CPXgeterrorstring (env, status, errmsg);
		printf ("%s", errmsg);
	}

    // Create the problem in CPLEX 
	strcpy(probname,"CG4FCND"); //Copy the name placed as the problem name
	lp = CPXcreateprob (env, &status, probname);  // Initialize the pointer that signals the problem 
	if ( env == NULL ) {
		char  errmsg[1024];
		printf ("Could not create LP. \n");
		CPXgeterrorstring (env, status, errmsg);
		printf ("%s", errmsg);
	}
	
	
	///We begin declaring the arc variables of the design.
	index1 = 0;  // index of columns
	numcols = M;
	d_vector(&obj,numcols,"open_cplex:1"); //objective function value
	d_vector(&lb,numcols,"open_cplex:8");//lower bound
	d_vector(&ub,numcols,"open_cplex:9");//upper bound
	c_vector(&ctype,numcols,"open_cplex:3");
    for(i=0;i<M;i++){
	   pos_arc[i] = index1;
       obj[index1] = edges[i].f;
       lb[index1] = lbtree[index1];
       ub[index1] = ubtree[index1];
	   ctype[index1]='B';
       index1++;
	}
	status = CPXnewcols (env, lp, index1, obj, lb, ub, ctype, NULL);
	if( status ) 
      fprintf (stderr,"CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_arc_var = index1;//number of Z variables

	//Next are the variables indexing the paths.
	index1 = 0;  // index of columns
	numcols = currentcolsPO;
	d_vector(&obj,numcols,"open_cplex:1"); //objective function value
	d_vector(&lb,numcols,"open_cplex:8");//lower bound
	d_vector(&ub,numcols,"open_cplex:9");//upper bound
	c_vector(&ctype,numcols,"open_cplex:3");
    for(k=0;k<K;k++){
	   for(i=0;i<num_paths[k];i++){
		   pos_tree[index1] = num_arc_var+index1;
		   obj[index1] = paths[locationinpaths[k][i]].PP;
		   lb[index1] = 0;
		   ub[index1] = 1;
		   ctype[index1]='B';
		   comindofcolumn[index1]=k;
		   pathcomnumofcolumn[index1]=i;
		   index1++;
	   }
	}
	status = CPXnewcols (env, lp, index1, obj, lb, ub, ctype, NULL);
	if( status ) 
      fprintf (stderr,"CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_tree_var = num_arc_var+index1;//number of Z variables

	/**********************************************************************************************/
	//Now we begin writing the restrictions beginning with the assignment constraints must select one of them
	acum=M;
	for(k=0;k<K;k++){
		numrows = 1;
		numnz = num_paths[k];
		d_vector(&rhs,numrows,"open_cplex:2");
		c_vector(&sense,numrows,"open_cplex:3");
		i_vector(&matbeg,numrows,"open_cplex:4");
		i_vector(&matind,numnz,"open_cplex:6");
		d_vector(&matval,numnz,"open_cplex:7");
		index = 0;//columns
		index1 = 0;//rows    
		sense[index1]='E';
		matbeg[index1] = index;
		rhs[index1]=1;
		for(i=0;i<num_paths[k];i++){
			matind[index]=acum+index;
			matval[index++]=1;
		}
		index1++;
		acum=acum+num_paths[k];
		status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
		if( status )  fprintf (stderr,"CPXaddrows failed.\n");
		free(matbeg);
		free(matind);
		free(matval);
		free(sense);
		free(rhs);
	}

	
	//These are the constraints one for each path onto the design variables
	acum=M;
	for(k=0;k<K;k++){
		for(j=0;j<M;j++){
			//if(path_indM[k*M+j].dim>0){
					numrows = 1;
					numnz = path_indM[k*M+j].dim+1;
					d_vector(&rhs,numrows,"open_cplex:2");
					c_vector(&sense,numrows,"open_cplex:3");
					i_vector(&matbeg,numrows,"open_cplex:4");
					i_vector(&matind,numnz,"open_cplex:6");
					d_vector(&matval,numnz,"open_cplex:7");
					index = 0;//columns
					index1 = 0;//rows
					sense[index1]='L';
					matbeg[index1] = index;
					rhs[index1]=0;
					for(i=0;i<path_indM[k*M+j].dim;i++){					
						matval[index]=1;
						matind[index++]=acum+path_indM[k*M+j].Kpath_ind[i];
					}
					matval[index]=-1;
					matind[index++]=pos_arc[j];
					index1++;
					//status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
					status=CPXaddusercuts (env,lp,index1, index, rhs, sense, matbeg, matind, matval, NULL);
					status=CPXaddlazyconstraints (env,lp,index1, index, rhs, sense, matbeg, matind, matval, NULL);
					if( status )  fprintf (stderr,"CPXaddrows failed.\n");
					free(matbeg);
					free(matind);
					free(matval);
					free(sense);
					free(rhs);
			//}
		}
		acum=acum+num_paths[k];
	}

	//These are the constraints one for each path onto the design variables one by one (NOTE: WE HAD TRIED THIS BUT REALIZED THAT USING USER AND LAZY MAKES IT TAKE LONGER AND IT USING ADDROWS, THERE'S NOT MUCH OF AN ADVANTAGE OVER AGGREGATED FORMULATION
	//acum=M;
	//numrows = 1;
	//numnz = 2;
	//d_vector(&rhs,numrows,"open_cplex:2");
	//c_vector(&sense,numrows,"open_cplex:3");
	//i_vector(&matbeg,numrows,"open_cplex:4");
	//i_vector(&matind,numnz,"open_cplex:6");
	//d_vector(&matval,numnz,"open_cplex:7");
	//for(k=0;k<K;k++){
	//	for(j=0;j<M;j++){
	//		//if(path_indM[k*M+j].dim>0){
	//				
	//				for(i=0;i<path_indM[k*M+j].dim;i++){
	//					index = 0;//columns
	//					index1 = 0;//rows
	//					sense[index1]='L';
	//					matbeg[index1] = index;
	//					rhs[index1]=0;
	//					matval[index]=1;
	//					matind[index++]=acum+path_indM[k*M+j].Kpath_ind[i];
	//					matval[index]=-1;
	//					matind[index++]=pos_arc[j];
	//					index1++;
	//					status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
	//					//status=CPXaddusercuts (env,lp,index1, index, rhs, sense, matbeg, matind, matval, NULL);
	//					//status=CPXaddlazyconstraints (env,lp,index1, index, rhs, sense, matbeg, matind, matval, NULL);
	//					if( status )  fprintf (stderr,"CPXaddrows failed.\n");
	//				}						
	//		//}
	//	}
	//	acum=acum+num_paths[k];
	//}
	//free(matbeg);
	//free(matind);
	//free(matval);
	//free(sense);
	//free(rhs);

	/*****************************New set of restrictions not exploring parts previously explored for cut and solve*****************************************************/
	if(method==7){
		d_vector(&matval,M,"open_cplex:6");
		d_vector(&rhs,1,"open_cplex:2");
		c_vector(&sense,1,"open_cplex:3");
		i_vector(&matbeg,1,"open_cplex:4");
		sense[0]='G';
		matbeg[0] = 0;
		for(j=0;j<M;j++) matval[j] = 1;
		for(i=0;i<globalcounter;i++){
			rhs[0]=(i+1)+1;
			if(type_CS==0){
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
		if(type_CS==0){				//Writing the restriction of variables not in the base solution solution to be 0
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
	}

/*************************CPLEX parameters********************************************************/
	CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON); //output display
	CPXsetintparam(env,CPX_PARAM_MIPDISPLAY,3);
	CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use
	CPXsetdblparam(env,CPX_PARAM_CUTUP,bestobj-0.1);//nothing worse than bestobjective found so far
	CPXsetdblparam(env,CPX_PARAM_TILIM,30); // time limit
	//CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_OFF); //output display	
	//CPXsetintparam(env,CPX_PARAM_REDUCE, 0);
	//CPXsetintparam(env,CPX_PARAM_LPMETHOD, 4);	
	//status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	//CPXsetintparam (env, CPX_PARAM_PRELINEAR, 0);
	//CPXsetintparam(env,CPX_PARAM_PREIND,0);
/*************************Define Objective function, solve and record time to solve********************************************************/
	//printf("/*************************************************************************************/\n\n\n");getchar();
	met_start_time=clock();	
	CPXchgobjsen(env, lp, CPX_MIN);
	status= CPXmipopt(env,lp);
	end_time=clock();
	status=CPXgetstat(env,lp);
	printf("\n\nStatus is %d from cutup value %lf\n\n",status,bestobj);
	//printf("/*************************************************************************************/\n\n\n");getchar();
	counttried++;
	
	timeHeur+=(double)(end_time - met_start_time) / CLOCKS_PER_SEC;
	/*********************************************************************************************/
	
	if(status==101||status==102||status==107||status==109 || status==105){
		/*************************Record the information*********************************************/
		currentcols = CPXgetnumcols (env, lp);
		d_vector(&xx,currentcols,"open_cplex:7");
		CPXgetmipobjval(env,lp,&valueMIP);
		if(bestobj>valueMIP){
			countimprovedheur++;
			CPXgetmipx(env,lp,xx,0, currentcols-1);
			/*************Copying to the best solution************/
			for(j=0;j<M;j++){
				bestsol[j]=xx[j];  //copying the result of the design variables
			}
			acum=M;
			for(k=0;k<K;k++){		
				for(j=0;j<num_paths[k];j++){
					if(xx[acum+j]>0.5){
						bestsol[M+k]=paths[locationinpaths[k][j]].PP;
						break;
					}
				}
				acum=acum+num_paths[k];
			}
			bestobj=valueMIP;
			printf("We have an improved solution with objective value in Metaheur  %.2lf\n",valueMIP);/*getchar();*/
		}
		
		/*********************************************************************************************/
		free(xx);
	}
	
	
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
	free(pos_tree);
	free(pos_arc);
	return valueMIP;
}

