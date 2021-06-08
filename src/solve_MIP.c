#include "def.h"


extern int        N, M, K ;
extern double      *d, *b;
//extern int        *open_facility;
//extern int        *customer_assign;
extern EDGE   *edges;
extern COMM   *com;



double solve_MIP(void)
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
	double    *x;      //Variable for th esollution
	int       num_z_var, num_x_var;
	int		*ann_indices;
	CPXLONG		*ann_vals;
	int		val_count=0;
	int		ind_count=0;
	

	CUTINFO usercutinfo;     // These are the ones that are not necessary to the model for the solution to actually be feasible. These appear at every node in order to improve the lower bound

	pos_y=create_int_matrix(N,N);
	pos_x = (int ***) calloc(N, sizeof(int **));
	for(i=0;i<N;i++) pos_x[i] = create_int_matrix(N,K);//Assign them to the size of the Arc set
	/***For annotations for the Benders Deceomposition***/
	ann_indices=create_int_vector(M+M*K);
	ann_vals=(CPXLONG *) calloc (M+M*K, sizeof(CPXLONG));
	

	/*****************************Going to initialize the x and z variables at -1*************************/
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			pos_y[i][j]=-1;
			for(k=0;k<K;k++) pos_x[i][j][k]=-1;
		}
	}

    objsen = 1; //min

	
	//Initialize CPLEX environment
	env = CPXopenCPLEX (&status); 
	if ( env == NULL ) {
		char  errmsg[1024];
		printf ("Could not open CPLEX. \n");
		CPXgeterrorstring (env, status, errmsg);
		printf ("%s", errmsg);
	}

	

    // Create the problem in CPLEX 
	strcpy(probname,instance); //Copy the name placed as the problem name
	lp = CPXcreateprob (env, &status, probname);  // Initialize the pointer that signals the problem 
	if ( env == NULL ) {
		char  errmsg[1024];
		printf ("Could not create LP. \n");
		CPXgeterrorstring (env, status, errmsg);
		printf ("%s", errmsg);
	}

	if(Use_ben==1) status = CPXnewlongannotation (env, lp, "cpxBendersPartition", 0); //Will create a new Benders annotation 
	                                        //Define Binary y_ij variables for each ij in A 
    index1 = 0;  // index of columns
	numcols = M;
	d_vector(&obj,numcols,"open_cplex:1"); //objective function value
	d_vector(&lb,numcols,"open_cplex:8");//lower bound
	d_vector(&ub,numcols,"open_cplex:9");//upper bound
	c_vector(&ctype,numcols,"open_cplex:01");//type of the variable

    for(i=0;i<M;i++){
	   pos_y[edges[i].i][edges[i].j] = index1;
       obj[index1] = edges[i].f;
       ctype[index1] = 'B';
       lb[index1] = 0;
       ub[index1] = 1;
	   ann_indices[ind_count]=ind_count;
	   ann_vals[ind_count]=val_count;
	   ind_count++;
       index1++;
	}
	status = CPXnewcols (env, lp, index1, obj, lb, ub, ctype, NULL);
	if( status ) 
      fprintf (stderr,"CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_z_var = index1;//number of Z variables
	
	                                        //Define Continuous x_ijk variables (flows)
    index1 = 0;  // index of columns
	numcols = M*K;
	d_vector(&obj,numcols,"open_cplex:1");
	d_vector(&lb,numcols,"open_cplex:8");
	d_vector(&ub,numcols,"open_cplex:9");
	c_vector(&ctype,numcols,"open_cplex:01");
     for(k=0;k<K;k++){
		 val_count=val_count+1;		//For Benders Subproblem
		 for(j=0;j<M;j++){
			pos_x[edges[j].i][edges[j].j][k] = num_z_var + index1; // position of the additional variables
			obj[index1] = com[k].d*edges[j].c[k];
			ctype[index1] = 'C';
			lb[index1] = 0;
			ub[index1] = CPX_INFBOUND;
			index1++;
			ann_indices[ind_count]=ind_count;
			ann_vals[ind_count]=val_count;
			ind_count++;
		  }
	 }
    status = CPXnewcols (env, lp, index1, obj, lb, ub, ctype, NULL);
	if( status ) 
      fprintf (stderr,"CPXnewcols failed.\n");
	free(obj);
	free(lb);
	free(ub);
	free(ctype);
	num_x_var = index1;//number of x variables


	                                             //Add flow constraints
	for(k=0;k<K;k++){
		for(i=0;i<N;i++){
			numrows = 1;
			numnz = 2*N;
			d_vector(&rhs,numrows,"open_cplex:2");
			c_vector(&sense,numrows,"open_cplex:3");
			i_vector(&matbeg,numrows,"open_cplex:4");
			i_vector(&matind,numnz,"open_cplex:6");
			d_vector(&matval,numnz,"open_cplex:7");
			index = 0;//columns
			index1 = 0;//rows    
			sense[index1]='E';
			matbeg[index1] = index;
			if(com[k].i==i) rhs[index1]=-1;//if node is origin
			else{
					 if (com[k].j==i) rhs[index1]=1; //if node is destination
					 else
					{
						rhs[index1]=0;
					}
				}
			for(j=0;j<N;j++){             //Scanning all the arcs			 	  
					if(pos_x[j][i][k]!=-1){
						matind[index] = pos_x[j][i][k];// arc is an incoming arc of i
						matval[index++] = 1;
					}
					if(pos_x[i][j][k]!=-1){
						matind[index] = pos_x[i][j][k];// arc is an outgoing arc of i
						matval[index++] = -1;
					}
			}
			index1++;
			status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
			if( status )  fprintf (stderr,"CPXaddrows failed.\n");
			free(matbeg);
			free(matind);
			free(matval);
			free(sense);
			free(rhs);
		}
	}
	                                             //Add linking constraints  x_ik <= z_i  OK
	for(j=0;j<M;j++){
		for(k=0;k<K;k++){ 
			numrows = 1;
			numnz = 2;
			d_vector(&rhs,numrows,"open_cplex:2");
			c_vector(&sense,numrows,"open_cplex:3");
			i_vector(&matbeg,numrows,"open_cplex:4");
			i_vector(&matind,numnz,"open_cplex:6");
			d_vector(&matval,numnz,"open_cplex:7");
			index = 0;//columns
			index1 = 0;//rows
			sense[index1]='L';
			rhs[index1]= 0;
		    matbeg[index1] = index;
			matind[index] = pos_x[edges[j].i][edges[j].j][k];
			matval[index] = 1;
			index++;
			matind[index] = pos_y[edges[j].i][edges[j].j];
			matval[index] = -1;
			index++;
			index1++;
			status = CPXaddrows (env, lp, 0, index1, index, rhs, sense, matbeg, matind, matval, NULL, NULL);
			if( status ) fprintf (stderr,"CPXaddrows failed.\n");
			free(matbeg);
			free(matind);
			free(matval);
			free(sense);
			free(rhs);
		}
	}
        

    /* Set parameters */   
	/************************************************/
	//CPXsetdblparam(env,CPX_PARAM_CUTSFACTOR, 1);  //limit the number of cuts added by cplex 1.0002
	CPXsetdblparam(env,CPX_PARAM_TILIM,Time); // time limit
	/*For Benders*/
	/************************************************/
	if(Use_ben==1){
		CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use
		CPXsetintparam(env,CPX_PARAM_MIPDISPLAY,3); //different levels of output displ
		CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON); //output display
		status=CPXsetlongannotations (env, lp, 0, 1, M+M*K, ann_indices,ann_vals );
		status = CPXsetintparam (env, CPXPARAM_Benders_Strategy, 1);
		//status = CPXwriteannotations(env, lp, "myprob.ann");
		/*status = CPXwritebendersannotation (env, lp, "benders2.ann");*/
		/************************************************/
	}
	else{
		CPXsetintparam(env,CPX_PARAM_THREADS, 1); // Number of threads to use
		CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
		CPXsetintparam(env,CPX_PARAM_MIPDISPLAY,3); //different levels of output displ
		CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON); //output display
	}
   
    //strcat(probname,".lp");
	//CPXwriteprob(env,lp,probname,"lp");                          //write the model in .lp format if needed (to debug)
	/************************************************/
	start_time = clock(); 
	if(Use_ben==1) status = CPXbendersopt (env, lp);
	else status= CPXmipopt(env,lp);  //solve the integer program
	if ( status ) { 
			fprintf (stderr, "Failed to optimize MIP.\n");
			goto TERMINATE;
	}
    status=CPXgetstat(env,lp);
	printf("Status is %d\n", status);
	if(status==101|| status==102|| status==106 || status==105|| status==107|| status==108){
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
   TERMINATE:
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

    for(i=0;i<N;i++){
		for(j=0;j<N;j++) free(pos_x[i][j]);
		free(pos_y[i]);
		free(pos_x[i]);
	}
	free(pos_x);
	free(pos_y);
	return Upper_bound;
}


static void free_and_null (char **ptr)
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
} /* END free_and_null */ 
