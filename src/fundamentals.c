#include "def.h"
                                              //Definition of data structures

                                 // Main function
void main (int argc, char *argv[])
{

	char       path[50];
	char		LP[50];
	char       line[100];
	char		solver[100];
	int		num_inst, i, j, confirm, r,s, k, addcuts;
	FILE		*ini;
	double     extra;
	double		RelaxationVal;
	int		Bestsolind; //1 =currentbest, 2=routing solution, 3=modified dijkstra solution
	int		repcuts;

	MAX_double = 1000000000;
	/*************************************************************************/
	//List of Flags for our algorithm
	Flag_Print=0;
	Flag_PrintSol=0;
	Time=86400;
	LifPlimit=0;
	depthLP=10;
	LPval=-MAX_double;
	/***********************************************************************/

	if(argc < 3) {
		printf("Error: Input file not specified \n");
		exit(8);
	}
	ini=Open_File (argv[1],"r");
	t= time(NULL);
	tm = localtime(&t);

	fscanf(ini, "%d", &num_inst);
	//fscanf(ini, "%s", &outfile);
	//fscanf(ini, "%s", &Progress);
	//fscanf(ini, "%s", &Solution);

	sprintf(outfile,argv[2]);
	sprintf(Progress,"Pr");
	strcat(Progress,argv[1]);
	/******Fixing the parameters*****/
	type=1;
	alter=4;
	base_fact=0.5;
	decrement=0.1;
	depmult=5;
	max_cuts=1;
	size_av=100;
	useheur=0;
	LifPlimit=0;
	depthLP=0;
	/************************************/
	
	out = Open_File(outfile,"a+");
	fprintf(out,"\n%s",asctime(tm));
	fprintf(out,"input_name;instance;N;M;K;Algorithm;RootLP;TotalBenCuts;PrimalSol;GlobLB;NodeCount;CPUtime\n");
	fclose(out);
	for(i=0; i<num_inst;i++){
		OK=0;
		confirm=0;
		/***********************Now reading the information for the instance**************************/
		if(fscanf(ini,"%s",&instance)==1) confirm++;
		if(fscanf(ini,"%d",&method)==1) confirm++;
		if(method<1||method>7){
			OK=1;
		}
		terminate=0;// initializing the indicator for LP
		/**********************************************************************************************/
		//Reading the instance to be solved
		sprintf(path,"./NetDesign/");
		strcat(path,instance);
		Read_Instance(path);
		/**********************************************************************************************/
		//Begin providing output
		out = Open_File(outfile,"a+");					
		fprintf(out,"%s;%s;",argv[1], instance);
		fclose(out);
		Upper_bound=MAX_double;
		bestobj=MAX_double;

		if(OK==0){
			if(method<=2){
				sprintf(solver,"CPX");
				if(method==2){ 
					Use_ben=1;
					strcat(solver,"Ben");
				}
				printf("Solving %s MIP Model with %s\n",instance, solver);
				optMIP=solve_MIP();
				repcuts=-1;
				LPval=-1;
			}
			else if(method<=6){
				sprintf(solver,"Ben");
				switch(method){
				case 4:
					strcat(solver,"+H");
					useheur=1;
					break;
				case 5:
					strcat(solver,"+LP");
					LifPlimit=1;
					depthLP=5;
					break;				
				case 6:
					strcat(solver,"+H+LP");
					useheur=1;
					LifPlimit=1;
					depthLP=5;
				}
				//Now beginning the procedure
				flag_keepcut=1;
				Initialize_Variables();					
				Initialize_Corepoint();
				create_network();
				start_time=clock();
				
				printf("Solving %s MIP Model with %s\n",instance, solver);				
				LPval=solve_BendersPre(0); //Solve LP
				if(100*(bestobj-LPval)/ABS(LPval)>.01){					
					for(j=0;j<4;j++) cuts_SC[j]=0;
					optBD=solve_BendersMC();					
				}
				else{
					end_time=clock();
					Upper_bound=LPval;
					best_lower_bound=LPval;
					repcuts=cuts_SC[1]+cuts_SC[3]+cuts_SC[0]+cuts_SC[2];
					nodecount=-1;
				}
				repcuts=cuts_SC[1]+cuts_SC[3]+cuts_SC[0]+cuts_SC[2];
				Free_Variables();
			}
			else{//We are going to be using the Benders Cut and Solve
				sprintf(solver,"CS/LB");
				printf("Solving %s MIP Model with %s\n",instance, solver);
				Sort_Edge_Weights();
				Initialize_Variables();
				Initialize_Corepoint();
				create_network();
				start_time=clock();
				out = Open_File(Progress,"a+");					//Writing what instance we are solving
				fprintf(out,"%s ; %s ; ", argv[1], instance);
				fclose(out);
				flag_keepcut=1;
				//Solving the LP model
				LPval=solve_BendersPre(0); //Solve LP
				//Running heuristic
				newglobobj=0;
				for(k=0;k<K;k++){
					DijkstraSS(com[k].i,com[k].j,k,LPsol);
				}
				Bestsolind=1;
				if(bestobj>feas_bound2){
					Bestsolind=2;
					bestobj=feas_bound2;
				}
				if(bestobj>newglobobj){
					Bestsolind=3;
					bestobj=newglobobj;
				}
				/************************recording the base solution for our procedure**************/
				for(j=0;j<M;j++){
					switch(Bestsolind){
					case 1:
						base_sol[j]=bestsol[j];
						break;
					case 2:
						base_sol[j]=sol_route[j];    //making the base solution
						break;
					case 3:
						base_sol[j]=newinisol[j];
						break;
					default:
						break;
					}
					if(base_sol[j]<=0.00001){
						base_soltrack[globalcounter][base_soltrackind[globalcounter]]=j;   //record which ones should be 0
						base_soltrackind[globalcounter]++;
					}
				}
				for(k=0;k<K;k++){
					switch(Bestsolind){
					case 1:
						base_solroute[k]=bestsol[M+k];
						break;
					case 2:
						base_solroute[k]=shortest[k].costwdem;
						break;
					case 3:
						base_solroute[k]=newinisol[M+k];
						break;
					}
				}
				if(100*(bestobj-LPval)/ABS(LPval)>.01){					
					for(j=0;j<4;j++) cuts_SC[j]=0;
					type_CS=0;
					do{
						optBD_prev=optBD;
						optBD=solve_BendersCS(type_CS);
						globalcounter++;
						use_firstsolution=1; //reset so we use the solution from previous
						optcount=0;  //Add this because we're not using any of the old cuts
						feascount=0;
						counttried=0;
						countimprovedheur=0;
						lasttotinpath=0;
						cputimeBD = (double)(end_time - start_time) / CLOCKS_PER_SEC;
					}while(ABS(optBD_prev-optBD)>1);
					/**********Now solving for the rest of the feasible space***********/
					type_CS=1;
					optBD=solve_BendersCS(type_CS);
					/****************************************************************/
				}
				else{
					end_time=clock();
					out = Open_File(Progress,"a+");					//Writing what instance we are solving
					fprintf(out,"0;0;0; %lf; %lf; 0;0\n", bestobj,cputimeBD);
					fclose(out);
				}
				repcuts=cuts_SC[1]+cuts_SC[3]+cuts_SC[0]+cuts_SC[2];
				nodecount=-1;
				best_lower_bound=Upper_bound;
				Free_Variables();
			}
			//Preparing output
			cputimeBD = (double)(end_time - start_time) / CLOCKS_PER_SEC;
			out = Open_File(outfile,"a+");
			fprintf(out,"%d;%d;%d;%s;%lf;%d;%lf;%lf;%d;%lf\n", N,M,K,solver,LPval,repcuts, Upper_bound,best_lower_bound,nodecount,cputimeBD);
			fclose(out);
			/**********************************************************************************************/	 
			Free_Memory();
			//getchar();
		}
		else{
			ClassifyErrors();
		}
	}
}

void Initialize_Corepoint(void){
	int r,i,j,addcuts=0, num_cuts=0;
	int *cutind;
	double *cutval;
	double val;
	cutind=create_int_vector(M+1);
	cutval=create_double_vector(M+1);
	core_point=create_double_matrix(N,N);
	disc=0;	
	MW_complete=0;
	for(r=0;r<M;r++) {
		//stabilizer[r]=bestsol[r];
		//MW_complete=MW_complete+core_point[edges[r].i][edges[r].j];
		if(sol_route[r]>0.9){
			val=0.7;
			core_point[edges[r].i][edges[r].j]=val;
			stabilizer[r]=val;
			ybar_frac[edges[r].i][edges[r].j]=val;
			ini_core[r]=val;
			MW_complete=MW_complete+core_point[edges[r].i][edges[r].j];
		}
		else{
			val=0.2;
			core_point[edges[r].i][edges[r].j]=val;
			stabilizer[r]=val;
			ybar_frac[edges[r].i][edges[r].j]=val;
			ini_core[r]=val;
			MW_complete=MW_complete+core_point[edges[r].i][edges[r].j];			
		}
	}
	///*********Here we are just confirming that indeed the solution is correct We can remove this later on**********/
	cutsettol=0.99999;
	do{
		num_cuts=0;
		for(i=0;i<num_origins;i++){
			for(j=0;j<origins[i].dim;j++){
				addcuts=maxFlowAndMinCut(origins[i].s, origins[i].d[j], ybar_frac, cutval,cutind);
				if(addcuts==1){
					for(r=0;r<feasnumnz;r++){
						val=1/(float)feasnumnz+0.01;
						ybar_frac[edges[cutind[r]].i][edges[cutind[r]].j]=ybar_frac[edges[cutind[r]].i][edges[cutind[r]].j]+ val;
						core_point[edges[cutind[r]].i][edges[cutind[r]].j]=ybar_frac[edges[cutind[r]].i][edges[cutind[r]].j]+val;
						stabilizer[cutind[r]]+=val;
						ini_core[cutind[r]]+=val;
					}
					num_cuts=num_cuts+addcuts;
				}
			}
		}
	}while(num_cuts>0);
}



void Update_core_point(void){
	int j;
	double lambda2;
	if(flag_frac==1){
		if(alter==1) lambda2=0.0;
		else if(globaldepth==0){
			lambda2=.95;
		}
		else if(globaldepth>0){
			lambda2=.1;
		}
		for(j=0;j<M;j++) core_point[edges[j].i][edges[j].j]=lambda2*stabilizer[j]+(1-lambda2)*ybar_frac[edges[j].i][edges[j].j];
	}
	else{
		if(alter==1) lambda2=0.0;
		else lambda2=0.95;
		for(j=0;j<M;j++) core_point[edges[j].i][edges[j].j]=lambda2*stabilizer[j]+(1-lambda2)*ybar[edges[j].i][edges[j].j];
	}
}

void ClassifyErrors(void){
	switch(OK){
	case 1:
		out = Open_File(outfile,"a+");
		fprintf(out,"File does not exist\n  ");
		fclose(out);
		break;
	case 2:
		out = Open_File(outfile,"a+");
		fprintf(out,"Instance size does not coincide between input and instance file\n  ");
		fclose(out);
		break;
	case 3:
		out = Open_File(outfile,"a+");
		fprintf(out,"Missing data in instance file\n  ");
		fclose(out);
		Free_Memory();
		break;
	case 4:
		out = Open_File(outfile,"a+");
		fprintf(out,"Missing data in input file or bad definition of Tao 0r Gamma\n  ");
		fclose(out);
		break;
	default:
		out = Open_File(outfile,"a+");
		fprintf(out,"Unexplained error\n  ");
		fclose(out);
		break;
	}
}

void Initialize_Variables(void)
{
	int i,k,j,ind1,repxcomm,index1,numcols,mult;
	int *index;
	//fixed_costs=create_double_matrix(N,N);
	currentFlow = create_double_matrix(N,N);
	stabilizer=create_double_vector(M);
	LPsol=create_double_vector(M+K);
	branchedon=create_int_vector(M);
	ini_core=create_double_vector(M);
	fixed_0=create_int_vector(M);
	fixed_1=create_int_vector(M);
	cuts_SC=create_int_vector(4);
	lbfix=create_double_vector(M+K);
	ubfix=create_double_vector(M+K);
	/*****************What was previously in the fundamentals, setting up for the Benders********/
	nodemult=size_av;
	size_av= 2; 
	tol=0.05; 
	Flag_Pap_Mod=1; 
	Flag_Allcuts=0;
	decrement=0.1;
	feasact=create_int_vector(10000);
	optact=create_int_vector(10000);
	Feas_Type=0;
	Flag_MaxCutPrint=0;
	TotalLPcuts=0;
	flag_incumb=0;
	Prev_incumbent=10000000;
	globaldepth=-1;
	maxdepth=globaldepth;
	previndex=-1;
	added_LP=0;
	totalinpaths=0;
	timeHeur=0;
	/******************For the cutandsolve****************************/
	//if(method==3){
		base_sol=create_int_vector(M);
		base_solroute=create_double_vector(K);
		base_soltrack=create_int_matrix(20,M);
		base_soltrackind=create_int_vector(20);
		globalcounter=0;
		newinisol=create_double_vector(M+K);
	//}
	/*******************For our heuristic approach*****************/
	lastdepth=0;
	lasttotinpath=0;
	countimprovedheur=0;
	counttried=0;
	/**************************************************************/
	
	cuts_SC[0]=0;//counter for the feasibility cuts in  Integer Solutions
	cuts_SC[1]=0;//counter for the optimality cuts in Integer Solutions
	cuts_SC[2]=0;//counter of the feasibility  cuts in Fractional Solutions
	cuts_SC[3]=0;//counter of the optimality  cuts in Fractional Solutions
	feas_bound1=MAX_double;
	feas_bound2=MAX_double;
	timeSep=0;
	timeFeas=0;
	Lagrangeind=0;
	use_firstsolution=1;
	index=create_int_vector(N);
	lagrangemult=create_double_vector(M);
	origins = (Origin *) calloc(N, sizeof(Origin));
	/*****For initial solutions****/
	/***************************************************/
	sol_route=create_int_vector(M);
	sol_fixed=create_int_vector(M);
	truesol=create_int_vector(M);
	/***************************************************/
	/*****For network definitions later on****/
	ybar=create_int_matrix(N,N);
	ybar_frac=create_double_matrix(N,N);
	coeff=create_double_vector(M);
	usedarc=create_int_matrix(N,N);
	inisol=create_double_vector(M+K);
	routingcost=create_double_vector(K);
	bestsol=create_double_vector(M+K);
	perturbnode=create_double_matrix(N,K);																			//Matrix with the Lagrangean perturbations for the Slope scaling
	perturbarc=create_double_matrix(M,K);
	/******************For the variable elimination test*************/
	forelimintest=create_double_vector(M);
	indarcinpath=create_double_vector(M);
	indbranch=create_int_vector(M);
	priority=create_int_vector(M);

	/******************For the column generation heuristic*************/
	if(useheur==1){
		col_limits=create_int_vector(K);
		num_paths=create_int_vector(K);
		maxpathperK=M;
		locationinpaths=create_int_matrix(K,maxpathperK);
		comindofcolumn=create_int_vector(maxpathperK*K);//array of the commodity index in the column
		pathcomnumofcolumn=create_int_vector(maxpathperK*K);//array of the path number of the commodity index in the column
		path_indM=(Kpath_index *) calloc (K*M, sizeof(Kpath_index));
		if((paths=(path *) calloc (maxpathperK*K, sizeof(path *)) )== NULL) {
			printf ("\nError: Insuficient memory \n");
			exit (8);
		}
		//for(i=0;i<maxpathperK*K;i++){
		//	paths[i].arcind=create_int_vector(M); //initialized the memory for paths
		//}
		for(i=0;i<K;i++){
			col_limits[i]=maxpathperK;
			for(j=0;j<M;j++){
				path_indM[i*M+j].Kpath_ind=create_int_vector(maxpathperK);//For now we are only allowing M*20 paths
				path_indM[i*M+j].dim=0;
			}
		}
	}
	/***************************************************/
	/***************************************************/
	/*Initializing the Structure of the paths*/
	/***************************************************/
	shortest =(Routes *) calloc(K, sizeof(Routes)); //Structure containing the pool of optimality cuts we initialize it with MK
	for(j=0;j<K;j++){
		shortest[j].path=create_int_vector(N);
		shortest[j].len=0;
	}
	/***************************************************/
	/*for(j=0;j<M;j++){
		fixed_costs[edges[j].i][edges[j].j]=edges[j].f;
	}*/
	/****************************************************/
	/***************Initializing the cutpools***************/
	mult=M*K;
	optcuts =(CUTPOOL *) calloc(mult, sizeof(CUTPOOL)); //Structure containing the pool of optimality cuts we initialize it with MK
	feascuts =(FEASPOOL *) calloc(mult, sizeof(FEASPOOL)); //Structure containing the pool of optimality cuts we initialize it with MK
	LPcuts=(CUTPOOL *) calloc(mult, sizeof(CUTPOOL)); //structure containing the pool of lift and project cuts
	optcount=0;
	feascount=0;
	LifPcount=0;
	LifProunds=0;
	optlim=1;
	feaslim=1;
	optcuts[optcount].coeff=create_double_vector(M+K+1);
	optcuts[optcount].ind=create_int_vector(M+K+1);
	/*****Added to the code in Jan 2019(*******/
	splitvars=create_int_vector(M);
	indprior=create_int_vector(M);
	//Declaring and filling the data for the structure origins
	/***************************************************/
	for(i=0;i<N;i++){
		origins[i].d=create_int_vector(N);
		origins[i].ind=create_int_matrix(N,K);
		origins[i].rep=create_int_vector(N);
		origins[i].dim=0;
		index[i]=-1;
		for(j=0;j<N;j++) origins[i].rep[j]=0; //initialize repetitions at 0
	}
	num_origins=0;
	mindemand=MAX_double;
	maxdemand=0;
	aggdemand=0;
	flag_keepcut=0;
	for(k=0;k<K;k++){
		if(index[com[k].i]>-1) ind1=index[com[k].i];
		else{
			ind1=num_origins;
			index[com[k].i]=num_origins;
			origins[ind1].s=com[k].i;
			num_origins++;
		}
		repxcomm=-1;
		for(j=0;j<origins[ind1].dim;j++){			//Here we will determine whether it is a new destination or a previous one
			if(origins[ind1].d[j]==com[k].j){
				repxcomm=j;
				origins[ind1].ind[repxcomm][origins[ind1].rep[repxcomm]]=k;				
				break;
			}
		}
		if(repxcomm==-1){
			repxcomm=origins[ind1].dim;
			origins[ind1].dim++;
		}
		origins[ind1].d[repxcomm]=com[k].j;
		origins[ind1].ind[repxcomm][origins[ind1].rep[repxcomm]]=k;
		origins[ind1].rep[repxcomm]++;
		if(maxdemand<com[k].d)maxdemand=com[k].d;
		if(mindemand>com[k].d)mindemand=com[k].d;
		aggdemand=aggdemand+com[k].d;
	}
	avdemand=aggdemand/K;
	free(index);
	/*************Calculate the maximum number of destinations******/
	max_num_dest=0;
	min_num_dest=N;
	for(i=0;i<num_origins;i++){
		if(origins[i].dim>max_num_dest)	max_num_dest=origins[i].dim;
		if(origins[i].dim<min_num_dest)	min_num_dest=origins[i].dim;
	}
	/***************************************************/
	//Call on the routine to create the arrays to put into CPLEX for Bender's Model
	Model_Variables();
	/***************************************************/
}

void Model_Variables(void){
	int i,j,k,index1,numcols,index;
	/********************************************************************************************************************************/
	/******************Here we are Calculating the solution of routing on the shortest path of each commodity***********************/
	presolve_ind=1;
	lowest_routing=0;
	for(i=0;i<num_origins;i++){
			lowest_routing+=Dijkstra(origins[i].s, origins[i].d, origins[i].ind, origins[i].rep, origins[i].dim, NULL);
	}
	//printf("Pure Routing Cost is %lf\n",lowest_routing);
	presolve_ind=0;
	feas_bound2=lowest_routing;
	countsol_route=0;
	for(j=0;j<M;j++){
		if(usedarc[edges[j].i][edges[j].j]>=1){
			feas_bound2=feas_bound2+edges[j].f; 
			sol_route[j]=1;
			ybar_frac[edges[j].i][edges[j].j]=1;
			countsol_route++;
			inisol[j]=1;
		}
		else{
			sol_route[j]=0;
			ybar_frac[edges[j].i][edges[j].j]=0;
			inisol[j]=0;
		}
	}
	printf("Routing Cost solution is %lf\n",feas_bound2);
	/**************Initializing the global problem variables****************/ 
	pos_y=create_int_matrix(N,N);
	pos_art=create_int_vector(K);
	pos_mu = create_int_matrix(N,N);
	/****************************************************/

	/*********Define Binary y_ij variables for each ij in A*********/
    index1 = 0;  // index of columns
	numcols = M;
	d_vector(&yobj,numcols,"open_cplex:1"); //objective function value
	d_vector(&ylb,numcols,"open_cplex:8");//lower bound
	d_vector(&yub,numcols,"open_cplex:9");//upper bound
	c_vector(&yctype,numcols,"open_cplex:01");//type of the variable
    for(i=0;i<M;i++){
	   pos_y[edges[i].i][edges[i].j] = index1;
       yobj[index1] = edges[i].f;
       yctype[index1] = 'B';
       ylb[index1] = 0;
       yub[index1] = 1;
       index1++;
	   indbranch[i]=i;
	}
	/*************************************************************/
	/***************************Calculating how many digits are in maxprior*****************/
	digits=0;
	//printf("maxprior is %d ",maxprior);
	while(maxprior >1)
	{
		maxprior/=10;
		digits++;
	}
	//printf("and has %d digits \n",digits);

	/*******************************************************************************/
	/*********Define continuous artificial variables z***********/
	                                        //Define Continuous zs lower bounds
    index1 = 0;  // index of columns
	numcols = K;
	d_vector(&zobj,numcols,"open_cplex:1");
	d_vector(&zlb,numcols,"open_cplex:8");
	d_vector(&zub,numcols,"open_cplex:9");
	c_vector(&zctype,numcols,"open_cplex:01");
	lowestpath=MAX_double;
	for(k=0;k<K;k++){
    		pos_art[k] = M + index1; // position of the additional variables
			zobj[index1] = 1;
			zctype[index1] = 'C';
			zlb[index1] =shortest[k].costwdem;
			zub[index1] = CPX_INFBOUND;
			index1++;
			if(lowestpath>shortest[k].costwdem)
				lowestpath=shortest[k].costwdem;
	}		
	/*************************************************************/
}

void UpdFeas_pool(int cut_type, double *coeffic, int *position, int count)
{
	int m,k,j, stop=0;
	feascuts[feascount].numcol=count;
	feascuts[feascount].coeff=create_double_vector(count);
	feascuts[feascount].ind=create_int_vector(count);
	for(m=0;m<count;m++){
		feascuts[feascount].coeff[m]=coeffic[m];
		feascuts[feascount].ind[m]=position[m];
	}
	if(cut_type==1)feascuts[feascount].fraction=0;
	else feascuts[feascount].fraction=1;
	feascount++;
	if(feascount==feaslim*M*K){
		feaslim++;
		feascuts= (FEASPOOL *) realloc(feascuts,feaslim*M*K*sizeof *feascuts); 
	}
	
}

int create_network(void ){
	int j;
	supply = create_double_vector(N);
	tail	 = create_int_vector(M);
	head	 = create_int_vector(M);
	obj	 = create_double_vector(M);
	ub	 = create_double_vector(M);
	lb	 = create_double_vector(M);

	// Define network structure: //////////////////////////////////////////////
	for (j=0; j<M; j++) {
		tail[j] = edges[j].i;
		head[j] = edges[j].j;
		lb[j] = 0.0;
	}
	///////////////////////////////////////////////////////////////////////////
	return 0;

}



void Free_Variables(void)
{
	int i,k,m,ind1,j;
	/*****************free column generation stuff*************/
	if(useheur==1){
		for(i=0;i<K*M;i++){
			free(path_indM[i].Kpath_ind);
		}
		for(i=0;i<totalinpaths;i++){
			free(paths[i].arcind);
		}
		for(i=0;i<K;i++){
			free(locationinpaths[i]);
		}
		free(paths);
		free(col_limits);
		free(num_paths);
		free(comindofcolumn);
		free(pathcomnumofcolumn);
		free(locationinpaths);
		free(forelimintest);
		free(indarcinpath);
		free(path_indM);
	}
	/**********************************************************/

	
	for(i=0;i<num_origins;i++){
		free(origins[i].d);
	}
	free(origins);
	free(cuts_SC);
	free(lagrangemult);
	for(i=0;i<optcount;i++){
		free(optcuts[i].coeff);
		free(optcuts[i].ind);
	}
	for(i=0;i<LifPcount;i++){
		free(LPcuts[i].coeff);
		free(LPcuts[i].ind);
	}
	for(i=0;i<feascount;i++){
		free(feascuts[i].coeff);
		free(feascuts[i].ind);
	}
	for(i=0;i<N;i++){
		free(ybar[i]);
		free(ybar_frac[i]);
		free(usedarc[i]);
		free(pos_y[i]);
		free(perturbnode[i]);
		free(currentFlow[i]);
		free(pos_mu[i]);
	}
	for(j=0;j<M;j++){
		free(perturbarc[j]);
	}
	free(pos_mu);
	free(currentFlow);
	free(perturbnode);
	free(perturbarc);
	free(ybar);
	free(pos_y);
	free(ybar_frac);
	free(coeff);
	free(inisol);
	free(usedarc);
	free(routingcost);
	free(yobj);
    free(yctype);
    free(ylb);
	free(yub);
	free(zobj);
    free(zctype);
    free(zlb);
	free(zub);
	free(pos_art);
	free(bestsol);
	free(feasact);
	free(optact);
	free(stabilizer);
	free(fixed_0);
	free(fixed_1);
	free(LPsol);
	free(branchedon);
	free(lbfix);
	free(ubfix);
	free(supply);
	free(tail);
	free(head);
	free(obj);
	free(ub);
	free(lb);
	free(ini_core);
	free(sol_route);
	free(sol_fixed);
	free(indbranch);
	free(priority);
	/******************For the cutandsolve****************************/
	//if(method==3){
		free(base_sol);
		free(base_solroute);
		free(base_soltrack);
		free(base_soltrackind);
		free(newinisol);
	//}
/****Added to code in Janury 2019***************/
	free(splitvars);
	free(indprior);
	
}


void Read_Instance (const char *name )         // This function opens and external file and reads the input data from there
{
  int i,j,r,s;
  double art1,art2,art3, arccost;
  FILE *in;

  in = Open_File (name,"r");
  if(OK==0){
		  if(fscanf(in,"%d %d %d",&N, &M, &K) != 3){             // Read the first two integer numbers in the file (number of nodes and edges)
			fprintf(stderr,"ERROR: Cannot read number of nodes and edges \n");
			OK=2;
		  }

		  Initialize_Memory(); 
		  //create temporary matrices
		  /*******************************************************/
		  matrix_A=create_int_matrix(N,N);
		  matrix_K=create_int_matrix(N,N);
		  /*******************************************************/

		  for(r=0;r<M;r++){                                // Read the two endnodes, fixed cost and transportation cost of each edge
			  if(fscanf(in,"%d %d %lf %lf %lf %lf %lf",&edges[r].i,&edges[r].j, &arccost, &art1/*This would have been capacity*/, &edges[r].f, &art2, &art3) != 7){
				fprintf(stderr,"ERROR: Cannot read edges \n");
				OK=3;
			  }
			  edges[r].i=edges[r].i-1;
			  edges[r].j=edges[r].j-1;
			  edges[r].c=create_double_vector(K);
			  matrix_A[edges[r].i][edges[r].j]=r+1;//recording the indices of each in a matrix
			  for(s=0;s<K;s++){						//Here we are repeating the cost of the arc for all of the commodities
				  edges[r].c[s]=arccost;
			  }
			  
		  }

		  for(s=0;s<K;s++){
			  if(fscanf(in,"%d %d %lf",&com[s].i,&com[s].j,&com[s].d) != 3){            // Read the commodities to be routed 0(k),d(k),d
				fprintf(stderr,"ERROR: Cannot read distance matrix \n");
				OK=3;
			  }
			  com[s].i=com[s].i-1;
			  com[s].j=com[s].j-1;
			  matrix_K[com[s].i][com[s].j]=s+1;//recording the indices of each in a matrix
		  }
  
		  fclose(in);
  }
}



void Read_Instance_Canad (const char *name )         // This function opens and external file and reads the input data from there
{
  int i,j,r,s,temp1,temp2, total;
  double art1,art2,art3, arccost;
  FILE *in;

  in = Open_File (name,"r");
  if(OK==0){
		  if(fscanf(in,"%d %d %d",&N, &M, &K) != 3){             // Read the first two integer numbers in the file (number of nodes and edges)
			fprintf(stderr,"ERROR: Cannot read number of nodes and edges \n");
			OK=2;
		  }

		  Initialize_Memory();                       

		  for(r=0;r<M;r++){                                // Read the two endnodes, fixed cost and transportation cost of each edge
			  if(fscanf(in,"%d %d %lf %lf %d",&edges[r].i,&edges[r].j, &edges[r].f, &art2, &total) != 5){
				fprintf(stderr,"ERROR: Cannot read edges \n");
				OK=3;
			  }
			  edges[r].i=edges[r].i-1;
			  edges[r].j=edges[r].j-1;
			  if(edges[r].i < 0 ||edges[r].i>=N){ printf("Something Wrong! at r\n"); getchar();} 
			  if(total!=K){printf("Not a complete list of commodity prices at %d we got %d\n",r, total);getchar();}
			  for(s=0;s<total;s++){						//Here we are repeating the cost of the arc for all of the commodities
				  if(fscanf(in,"%lf %lf %lf",&art2,&art3, &art1) != 3){
					fprintf(stderr,"ERROR: Cannot read edges \n");
					OK=3;
				  }
				  edges[r].c[s]=art3;
			  }
		  }

		  for(s=0;s<K;s++){
			  fscanf(in,"%lf %d %lf",&art3,&temp1,&art1);
			  fscanf(in,"%lf %d %lf",&art3,&temp2,&art2);
			  if(art1>=0){
				  com[s].i=temp1-1;
				  com[s].j=temp2-1;				  
				  com[s].d=art1;
			  }
			  else{
				  com[s].i=temp2-1;
				  com[s].j=temp1-1;
				  com[s].d=art2;
			  }
		  }  
		  fclose(in);
  }
}

void Initialize_Memory (void)                 // This function initializes the memory required for the arrays (dynamic allocation of memory)
{
  int j;
  edges = (EDGE *) calloc(M, sizeof(EDGE));
  com=(COMM *) calloc(K, sizeof(COMM)); 
  for(j=0;j<M;j++){
	  edges[j].c=(double *) calloc (K, sizeof(double));
  }
}

void Free_Memory (void)                       // This function frees the memory that has been previously requested
{                                             // This function should always be called after you solved the problem and just 
 int i;	                                      // and jsut before terminating the program
 free(edges);
 free(com);
}

FILE *Open_File (const char *name, const char *mode)  // This function opens the file "name" under mode "mode" (reading, writing, etc)
{
 FILE *file;

 if((file=fopen(name,mode))==NULL) {
    printf("\nError: File cannot be opened \n");
    OK=1;
 }
 return file;
}

int **create_int_matrix (int rows, int Columns)
{
 int i;
 int **ptr;

 if((ptr=(int **) calloc (rows, sizeof(int *)) ) == NULL) {
   printf ("\nError: Insuficient memory \n");
   exit (8);
 }
 for(i=0;i<rows;i++)
   ptr[i] = create_int_vector(Columns);
 return ptr;
}

double **create_double_matrix (int rows, int Columns)
{
 int i;
 double **ptr;

 
 if((ptr=(double **) calloc (rows, sizeof(double *)))==NULL) {
    printf("\nError: Insuficient memory \n");
    exit(8);
  }
  for(i=0;i<rows;i++) {
	ptr[i]=create_double_vector(Columns);
  }
  return ptr;
}


int *create_int_vector (int dim)
{
 int *ptr;

 if((ptr=(int *) calloc (dim, sizeof(int)))==NULL) {
    printf("\nError: Insuficient memory \n");
    exit(8);
 }
 return ptr;
}

double *create_double_vector (int dim)
{
 double *ptr;

 if((ptr=(double *) calloc (dim, sizeof(double)))==NULL) {
    printf("\nError: Insuficient memory \n");
    exit(8);
 }
 return ptr;
}
// CPLEX functions to allocate memeory to arrays

void i_vector(int **vector,int n,char *s)
{
if((*vector=(int *)calloc(n,sizeof(int)))==NULL)
  //error(s);
  printf("Error: Insuficient memory \n");
return;
}

void d_vector(double **vector,int n,char *s)
{
if((*vector=(double *)calloc(n,sizeof(double)))==NULL)
 // error(s);
 printf("Error: Insuficient memory \n");
return;
}

void c_vector(char **vector,int n,char *s)
{
if((*vector=(char *)calloc(n,sizeof(char)))==NULL)
  //error(s);
  printf("Error: Insuficient memory \n");
return;
}

int maxFlowAndCut(int s, int *t, int dim_dest, int **Capacity, double *coeffi, int *ind){
    int		i,j;							// index variables for loops
    int		u,v;                            // index variables for nodes
	int		*visited;
	int		queueIndex, queueSize;          // just for use an array to store the connected part of s
	int		num_visit=0;
	int		*queue;
	int		retval=0;

	queue=create_int_vector(N);
	visited=create_int_vector(N);
	feasnumnz=0;
	// initialize everything else!
	for (i=0; i<N; i++) visited[i] = 0;     // nobody has been visited yet
    visited[s] = 1;                           //  ...nobody but 's'         
	queueSize = 0;                              // the queue contains 0 elements
	queue[queueSize++] = s;                     // now it contains s
    queueIndex = 0;                             // next element to leave the queue      
    // perform a BFS to find a augmenting path between 's' and 't'
    while (queueIndex < queueSize){
            u = queue[queueIndex++];                // get the next element
            // now look if you can extend the augmenting path from 'u' to any other 'v'
            for (v=0; v<N; v++) {
                if (visited[v]==0 && Capacity[u][v]==1) {
                    visited[v] = 1;               // if you can, visit the new node
                    queue[queueSize++] = v;           // put it in the queue           
                }
			}                              // avoid do extra work!           
	}  
    //We now have the connected component of s time to check if there exists a commodity that is not connected
	for(i=0;i<dim_dest;i++){
		if(visited[t[i]]==0){
				retval=1;
				com_check_ind=i;
				break;
		}
	}
	//if there was a destination not connected, then 
	for(j=0;j<M;j++){  
	  if(visited[edges[j].i]==1 && visited[edges[j].j]==0){
		coeffi[feasnumnz]=1;
		ind[feasnumnz]=pos_y[edges[j].i][edges[j].j];
		feasnumnz=feasnumnz+1;
	  }
	}	
	free(queue);
	free(visited);
	return retval;
}

////////////////////////////////////////////////////////////////////////////////
// Performs the Edmonds-Karp Algorithm and returns the maximum flow between
// 's' and 't'. The algorithm also remembers the connected component that
// contains 's' in the binary array 'visited'int 
int maxFlowAndMinCut (int s, int t, double **Capacity, double *coeffi, int *ind){
    int		i,j,retval,index;							// index variables for loops
    int		u,v;                            // index variables for nodes
	int		queueIndex, queueSize;          // just for use an array as a queue
    double	flow;							// the augmenting flow
    double	maxFlow = 0.0;                  // the total amount of flow
	double  *cutval;
	int		*cutind;
	int		*previous = create_int_vector(N);
	int		*visited = create_int_vector(N);
	int		*queue = create_int_vector(N);	// basic queue

	int DEBUG = s<0 ? YES : NO;
	s = ABS(s);

	feasnumnz=0;
    // initialize the 'currentFlow'    
    for (i=0; i<N; i++) { for (j=0; j<N; j++) { currentFlow[i][j] = 0.0; } }
    
    // perform as many flow augmentations as you can!
    flow = 1.0;                                     // this is just to start the loop
    while (flow > 0.0) {
        flow = 0.0;
        // initialize everything else!
		for (i=0; i<N; i++) { visited[i] = NO; }    // nobody has been visited yet
        visited[s] = YES;                           //  ...nobody but 's'
        queueSize = 0;                              // the queue contains 0 elements
		queue[queueSize++] = s;                     // now it contains s
        queueIndex = 0;                             // next element to leave the queue
      
        // perform a BFS to find a augmenting path between 's' and 't'
        while (queueIndex < queueSize){
            u = queue[queueIndex++];                // get the next element
            // now look if you can extend the augmenting path from 'u' to any other 'v'
            for (v=0; v<N; v++) {
                if (visited[v]==NO && Capacity[u][v]-currentFlow[u][v] > 0.0) {
                    visited[v] = YES;               // if you can, visit the new node
                    queue[queueSize++] = v;           // put it in the queue
                    previous[v] = u;                // remeber where you come from
                    if (v==t) break;                // and avoid do extra work! 
                }
			}

            // now, before continue, check if we have finished the BFS!
            if (visited[t] == YES) {
                // we are going to trace back our steps looking for how much flow we can send
                v = t;
                u = previous[t];
                flow = Capacity[u][v]-currentFlow[u][v];
                while (u!=s) {
                    v = u;
                    u = previous[v];
                    if (flow > Capacity[u][v]-currentFlow[u][v])
                        flow = Capacity[u][v]-currentFlow[u][v];
                }
                // now we can update the currentFlowMatrix with that amount of flow
                v = t;
                u = previous[t];
                currentFlow[u][v] += flow; // we add the flow in one direction
                currentFlow[v][u] -= flow; // and remove it from the other one!
                while (u!=s) {
                    v = u;
                    u = previous[v];
                    currentFlow[u][v] += flow; // we add the flow in one direction
                    currentFlow[v][u] -= flow; // and remove it from the other one!
                }        
                maxFlow += flow;                    // we have improved so we update
                break;                              // avoid do extra work!
            }
        }  
    }
	retval=0;

	if(Flag_MaxCutPrint==1){ //We have solved to obtain the LP solution now we will keep the cutset
		//contador=0;
		for(j=0;j<M;j++){  
			  if(visited[edges[j].i]==1 && visited[edges[j].j]==0){
				coeffi[feasnumnz]=1;
				ind[feasnumnz]=pos_y[edges[j].i][edges[j].j];
				feasnumnz++;
			  }
		}				
	}
	else{
		if (maxFlow<cutsettol/*-EPSILON2*/){
			//if(Flag_Print==1) printf("\n Maxflow is %.20lf and cutsettol is %.20lf \n", maxFlow, cutsettol);
			retval=1;
			for(j=0;j<M;j++){  
			  if(visited[edges[j].i]==1 && visited[edges[j].j]==0){
				coeffi[feasnumnz]=1;
				ind[feasnumnz]=pos_y[edges[j].i][edges[j].j];
				feasnumnz++;
			  }
			}
		}
	}
	/*out = Open_File("Cutset","w");
	for(i=0;i<feasnumnz;i++){
		fprintf(out,"(%d,%d); ",edges[ind[i]].i+1,edges[ind[i]].j+1 );
	}
	fclose(out);*/
	// clean:
	free(previous);
	free(queue);
	free(visited);
    return retval;
}

void Sort_Edge_Weights(void)                 // This funtion sorts the edges in nondecreasing order with respect to the weights
{
 qsort((COMM *) com, K, sizeof(com[0]), Compare_Weight);
}


int Compare_Weight(const void *a, const void *b)  // Auxiliar function used in qsort to explicitly state that the sorting is nondecreasing 
{                                                 // with respect to element "weight" of data structure EDGE
 if(((COMM *)a)->d>((COMM *)b)->d)
   return -1;
 if(((COMM *)a)->d<((COMM *)b)->d)
   return 1;
  return 0;
}




double DijkstraSS ( int start, int t, int comind, double * param){
	int current, u,v;
	int i,j,k,reps, count, lowindex;
	double newdist, temp,tot_dist;
	double		*distance;
	tot_dist=0;

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
			if(edges[j].i==current ){
				newdist=distance[current]+com[comind].d*/*(1-param[M+comind])**/edges[j].c[0]+(1-param[j])*edges[j].f;
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
	v = t;
	u = preced[t];
	newdist=0;
	//Here we are selecting which ones we had used
	while (v!=start) {
		if(newinisol[matrix_A[u][v]-1]<.99){
			newinisol[matrix_A[u][v]-1]=1;
			newglobobj=newglobobj+edges[matrix_A[u][v]-1].f;
		}
		newdist=newdist+com[comind].d*edges[matrix_A[u][v]-1].c[comind];
		v = u;
		u = preced[v];					
	}
	newglobobj+=newdist;
	newinisol[M+comind]=newdist;
	/*printf("Routing commodity (%d,%d) which has distance %lf\n", com[index[k][reps]].i,com[index[k][reps]].j , distance[t[k]]);
	getchar();*/
	free(distance);
	free(preced);
	free(list);

	
	return newdist;
}