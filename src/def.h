#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <string.h>
#define getrandom( min, max ) ((rand() % (int) (((max)+1)-(min)))+(min))
#define true 1
#define false 0
#define NOT_AVAIL -1
#define NONE -1
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) > (b)) ? (b) : (a))
#define ABS(x) (((x) > 0 ) ? (x) : -(x))	
#define ALL		2
#define YES     1	// True
#define NO      0	// False
#define NONE   -1	// Unknown
#define EPSILON	 0.000000000001    //Added to the fractional solution obtained at the node
#define EPSILON2 0.0000000001	////Added to check for feasibility cut violation


int create_network(void );
double solve_BendersPre(int);							//Declaration for solving the LP
FILE   *Open_File(const char *, const char *);
int ConstrMatrix_LP( double *x, double Tol_Act);
int frac_index(double *x);
int maxFlowAndCut(int , int *, int , int **, double *coeffi, int *ind);//For integer solutions. Delivers the cutset inequality that separates the current solution
int maxFlowAndMinCut(int , int , double **, double *, int *ind);//For fractional solutions. We will add one cut per common origin.
double DijkstraSS ( int start, int t, int comind, double * param);//Dijkstra to calculate the thing we want
double solve_net_exact(double *temp_x);					//Used to solve the exact Benders solution
void Model_Variables(void);								//Function to declare the variables of the Benders Model
void StrongFeas_com(int); //function that optimizes over the recession cone of a given commodity
void UpdFeas_pool(int cut_type, double *coeffic, int *position, int count); //function for updating the corresponding cutpool 1- Integer Optimality, 2- Integer Feasibility, 3-fractional Optimality, 4- Fractional Feasibility
double solve_net_sub(double *temp_x);
double solve_Dual(double *temp_x,int addedcuts);
double Dijkstra ( int , int *, int **, int *, int, double *);
void   Read_Instance(const char *);
void Initialize_Corepoint(void);
void Update_core_point(void);
void   Initialize_Memory (void) ;
void   Free_Memory(void);
void ClassifyErrors(void);	//Function that classifies the errors in input data that are detected but allows for the code to continue on to another instance or line
int    **create_int_matrix(int, int);
double **create_double_matrix (int, int);
int    *create_int_vector(int);
double *create_double_vector(int);
void i_vector(int **vector,int n,char *s);
void d_vector(double **vector,int n,char *s);
void c_vector(char **vector,int n,char *s);
double solve_MIP(void);
void Initialize_Variables(void);
void Free_Variables(void);
double solve_BendersMC(void);
int LiftandP(int split_var, double *soln);
double **currentFlow ;
double solve_BendersCS(int indic);

 /*****Columng generation heuristic********/
 void populate_pathperkSS(double *x_sol, int com_ind);
 double master_p(double *lbtree, double *ubtree);


 /*******Additional stuff for the commodity relaxations*******/
 void Sort_Edge_Weights(void);
 int Compare_Weight(const void *a, const void *b);

typedef struct Origin
{
	int s;       //The origin node shared by all these commodities
	int *d;		//An array containing all destinations for the commodities with this origin node
	int dim;	//The number of different destinations from this origin
	int **ind;	//An array containing the indices in the original commodity array of these commodities	
	int *rep;   //array with lengths of repetitions for each destination
} Origin;

typedef struct EDGE
{
int i;
int j;
double f;
double * c;
} EDGE;

typedef struct Routes //Structure that we will have for the shortest routes based on the whole network
{
int *path;				//The sequence of nodes on the path
int len;				//length of the path
double costwdem;			//cost of routing the demand
double cost;              //Costo of path
} Routes;

typedef struct COMM
{
int i;
int j;
double d;
} COMM;
                                              //Definition of global variables
typedef struct CUTPOOL {
  double *coeff;				//array of the coefficients
  int	 *ind;					//array of indices
  double  RHS;					//Right hand side
  int  numnz;				//non zero elements
  int	  fromfrac;				//1 if it's an optimality cut from fractional solution 0 from integer solution
} CUTPOOL;

typedef struct FEASPOOL {
  double *coeff;				//array of the coefficients
  int    *ind;					//array of index
  int	numcol;					//number of nonzeros
  int fraction;					//from fractional solution or not
} FEASPOOL;

struct cutinfo {
   CPXLPptr lp;
   int    numcols;
   double bestlb;
};

/******8New stuff for the CG heuristic*******/
                                              //Definition of global variables
typedef struct path {
 double value; //If >0 then it will stay
 double PP; //Value of the cost times demand
 double fixed; //value of the fixed cost
 double distance; //length of the path with respect to C
 int *arcind; //array of the indexes of arcs in the path
 int numarcs; //Number of arcs contained in the path
} path;

typedef struct Kpath_index
{
	int dim;	//number of arcs in the path of the commodity
	int *Kpath_ind;	//An array containing the indices in the original commodity array of these commodities	
} Kpath_index;
/*******************************************/

typedef struct cutinfo CUTINFO, *CUTINFOptr;  //Declaring the cutinfo and the cutinfoPtr 

static int makelazyconstraint (CPXENVptr env, CPXLPptr lp, CUTINFOptr cutinfo);
static int makeusercuts (CPXENVptr env, CPXLPptr lp, CUTINFOptr cutinfo);
static int CPXPUBLIC lazy_feas_callback1(CPXCENVptr env, void *cbdata, int wherefrom, void *info, int *useraction_p);
static int CPXPUBLIC user_cut_callback1(CPXCENVptr env, void *cbdata, int wherefrom, void *info, int *useraction_p);
static int CPXPUBLIC LagrangeHeur(CPXCENVptr env, void *cbdata, int wherefrom,  void *cbhandle, double *objval_p, double *x, int *checkfeas_p, int *useraction_p);

int			N, M,K;
int			type;
int			method;  //Variable we are going to use to distinguish what type of method to use: 1-CpleX; 2-Benders Single Cut; 3-Benders Multicut.
int			alter; //Variable that lets us know what type of alteration to the subproblem we want all guarantee feasibility: 1-Normal Iterative ; 2- M-W(LP); 3- Papadakos (LP); 4- M&W (MinCostFlow)  
int			disc;// Discount factor for corepoint update for Papadakos
double		decrement;//decrement amount for corepoint update for Papadakos
double		base_fact;// the base factor for corepoint update for Papadakos
int			Feas_Type;//Indicator: 0-Cutset inequalities (All violated); 1-Cutset Inequalities (Add only first found); 2- Use the strong feasibility cut separator
int			terminate; //Indicate when we have finished solving the LP.
double		**c;
double		**f;
EDGE		*edges;
COMM		*com;
COMM		*com2;
path		*paths; //Structure containing the pool of paths, each row is a different commodity and in each entry, is a path.
int			**matrix_A;//Matrix with indices of the arcs
int			**matrix_K;//Matrix with indices of the commodities
//Declaring all the cutpools
/*******************************************************************************************/
CUTPOOL     *optcuts; //Structure containing the pool of optimality cuts
CUTPOOL		*LPcuts;  //structure containing the pool of Lift and project cuts
FEASPOOL     *feascuts; //Structure containing the pool of optimality cuts
int			optlim;// Limit of the global optimality cut pool limit
int	        optcount;//Counting how many you integer optimality cuts have been added
int			feaslim;// Limit of the global optimality cut pool limit
int	        feascount;//Counting how many you integer optimality cuts have been added
int			feasnumnz;//Global counter of number of non zeros in feasibility cut
int			optnumnz;//Global counter of number of non zeros in feasibility cut
int			vioptcuts;//Number of violated optimality cuts found
/******************************************************************************************/
/**************Declaring the global arrays for the problem definition**********************/
int			**pos_y;  // positions of the design variables (double, even if the problem is integer) .....
int			***pos_x;  // positions of the routing variables
int			*pos_art;	//positions of the artificial variables of Benders Decomposition Multi cut 
double      *yobj;    // objective function coefficients of the design variables
double      *ylb;     //lower bound of the design variables
double      *yub;     //upper bound of the design variables
char		*yctype;  //Type of the design variables
int			**ybar;  //The y bar that we will share for integer solutions
double		**ybar_frac; //The y bar that we will share for fractional solutions
/******************************************************************************************/

/**************Declaring the global arrays used throughout the code**********************/
int			***pos_x;  // positions of the routing variables
int			*pos_art;	//positions of the artificial variables of Benders Decomposition Multi cut 
double      *yobj;    // objective function coefficients of the design variables
char		*yctype;  //Type of the design variables
double      *zobj;    // objective function coefficients of the design variables
double      *zlb;     //lower bound of the design variables
double      *zub;     //upper bound of the design variables
char		*zctype;  //Type of the design variables
/******************************************************************************************/
char		instance[20];           // External file containing all the input data of the problem
int			*list;
int			*preced;
int			*arco;
double		**lambda;
double		***mu;
int			cut_type;
double		*coeff;     //vector with coefficients for single cut
double      Righthand;  //Right hand side value for single cut
int			*cuts_SC;
FILE		*out;
time_t		t ;																														//Variables for taking the time and knowing when the code was exectued
struct tm	*tm;	
char		outfile[20];
char		Progress[20];
char		Solution[20];
clock_t		met_start_time;
clock_t		start_time, end_time;
clock_t		start_sep, end_sep, start_feas, end_feas;
double		timeMast, timeFeas, timeSep, timeHeur;
double		cputimeMIP,optMIP,optBD, optBDMC;
double		MAX_double;
int			OK;		//Declaring a variable to check if everything is being read correctly.
int			cur_numcols, num_origins;//number of origins and nuber of columns
int			flag_frac;//Flag that indicates whether I am at an integer or fractional separation problem
Origin		*origins;
double		**core_point;
int			lazycalled;                                              //Declaration of user-defined functions
/*Here are some of the variables needed to do the additional tricks of our branch and cut.*/
int		usercutsnum;																			//Number of user cuts added at the current node used for evaluating whether more should be added
int		last_node;																				//Variable stored the last node that was visited
int		rep_node;																				//The number of times this node has been already explored
double  last_obj;																				// The value of the last objective function
double  *average;																				//added to calculate the average improvement
/*Parameters set at their "almost fine tuned" values for our branch and cut*/
int		flag_upd_av;																			//Flag that lets us know that we need to update	
int		depmult;																					//Depth at witch we check to add user cuts
int		max_cuts;																					//Maximum number of cuts to be added at a node
int		size_av;																					//How many of the last solutions are to be evaluated at the root node to calculate the improvement 
double  tol;																					//Improvement Tolerance necessary to continue adding cuts at the root node 0.05%
int		Feas_unchanged;																			//Number of times the value at the node has gone unchanged.
int		com_check_ind;
/***************************************************************/
//Parameters for testing that we can play with
int		Flag_Print;																					//Flag if you want to print the progress at every separation problem
int		Flag_PrintSol;																			//Flag if you want to print the solution obtained at the end in a file.
int		Flag_Pap_Mod;																			//Flag that lets us know whether we are going to add only the cuts violated in the Papadakos
int		Flag_Allcuts;																			//Flag that defines whether in the Multicut version, we will add all of the cuts or only those that are violated.
int     Flag_Check;																				//Flag to tell us that there was an infeasibility in the net opt
int		Failed_Eval;																			//We use this to know whether we were successfully able to evaluate or not the fractional solution. If not, then we do not at the M-W restriction.
double	Time;																					//Value of the Time we give the algorithm
double	violation;
double  cutsettol;
double  bestobj;																				//best objective value from heuristic solution
double  *bestsol;																				//best solution found from heuristic
/***************************************************************/
//Parameters for our new heuristic
double *lagrangemult;																			//Array that will store the lagrangean multiplier to be used in the heuristic
int     Lagrangeind;																		    //Indicator that there exist lagrange multipliers
int		**usedarc;																					//Array where we will count the arcs that are being used.
int		last_feas;																					//Counter for the index of the last feasibility cut added																				
Routes	*shortest;																				//Shortest routes routing paths
Routes	*current;																				//Structure of the paths of current solution
double	Dom_coefficient;																			//Dominance coefficient to know which cost is more dominant over the other
double	*lb_art;																					//lowerbound on the artificial variable for commodity k-the shortest path for each commodity
int		presolve_ind;																				//Indicator=1 if in presolve meaning don't take into account ybar for Dijkstra 0 otherwise
double	lowest_routing;																			//Value of the lowest routing cost possible
int		degree_flag;																			//Flag to be used to know whether Dijkstra_FC is using cardinality or the actual fixed costs
double  min_conn;																				//Value of what is the minimum number of arcs to be open to ensure connectivity
double	**fixed_costs;																			//Matrix where fixed costs are stored as matrices
int		**Best_arcs;																				//Here we will store the best open arcs
double  *Best_routes;																			//Here we will store the routing cost of the current best arcs configuration
double  best_intval;																			//Best feasible integer objective value so far
double	*routingcost;                                                                         //Vector for the routing cost
int     use_firstsolution;                                                                     //Indicator to let you know that you will use the first solution found in the heuristic
double  *inisol;																				//Array for the initial solutions that we obtain
double  maxdemand;																				//Save the value of the maximum amount of demand to be sent
double  mindemand;																				//Save the value of the minimum amount of demand to be sent
double  avdemand;																				//arithmetic average of the demands
double	aggdemand;																				//total aggregate of the demands
double  **perturbnode;																			//Matrix with the Lagrangean perturbations for the Slope scaling
double	**perturbarc;																			//Matrix with the Lagrangean perturbations for the Slope scaling 
int	    flag_keepcut;																			//indicator that let's us know if we'll keep the cuts
/**************For Ivan's experiments**************************/
int  *sol_route;																				//open arcs of solution from efficient routing
int   *sol_fixed;																				//open arcs of solution from efficient fixed costs
int  *truesol;																					//open arcs of the real optimal solution
double **routheur;																				//routing values of the heuristic solution per commodity
double	feas_bound2;																		//Upper bound obtained basing ourselves on the routing costs
double   feas_bound1;																			//Upper bound obtained basing ourselves on the fixed costs
int countsol_route;																				//Number of open arcs from efficient routing																				
int countsol_fixed;																				//Number of open arcs from efficient fixed costs
int countsoltrue;																				//number of opern arcs at optimal solution
double fixedoptimalcost;																		//Part of the optimal solution given by the fixed cost
double routoptimalcost;                                                                         //Part of the optimal solution given by the routing cost
int routinopt;																					//Number of arcs ARE in routing solution that ARE in the optimal solution
int fixedinopt;																					//Number of arcs ARE in fixed cost solution that ARE in the optimal solution
int routoutopt;																					//Number of arcs ARE in routing solution that are NOT in the optimal solution
int fixedoutopt;																					//Number of arcs in fixed cost solution that are NOT in the optimal solution
int missrout;																					//Number of arcs NOT in routing solution that ARE in the optimal solution
int missfixed;																					//Number of arcs NOT in fixed cost solution that ARE not in the optimal solution																	
double MW_complete;
/**********************************************Parameters for Lift and Project*********************************************************************/
double **BigMatrix;																			   //Array that will store the active constraints
double *LifPcut;																			   //Array that will store the Lift and Project Cut
int    countBigMatrix;                                                                        //Number of active constraints in BigMatrix
int		Flag_MaxCutPrint;																	//Flag to tell that we will separate for cutsets despite having all of them placed already
double  *solin;																			//Aray where we'll keep the solution from the LP of the active constraints
int cutsetiter;																			//counter for cutset iteration
int optcutiter;																			//counter for optimality cut iteration
int completeactive;																		//1 if we will use the complete list of active constraints
int flag_entered;																		//Flag to let you know if you've entered into the routine of finding active constraints
int LifPcount;																			//Counter for the number of lift and project cuts done
double *alphacalc;																		//vector calculating the lift and project coefficients
double betacalc;																		//Right hand side of the lift and project cut
double origLP;																			//value of the original LP
double finallp;																			//value of the final LP	
int  bigMconstraints;																//Number of constraints at the moment
int  LandPdepth;																	//Depth up to which we will add Lift and Project cuts
int LifPlimit;																			//Maximum number of LifP permitted per node (-1 means add all cuts from fractional variables of solution)
int IndLifandP;																		//Global indicator stating when you haven't entered the Lift and Project procedure yet in a node
int num_feasact;																	//Active Feasibility cuts
int	num_optact;																	//Actuve Optimality Cuts
int *optact;																		//Array with the indices of the active optimality constraints
int *feasact;																		//Array with the indices of the active feasibility constraints
int split_var;																		//Variable on which we will split
int depthLP;																		//depth at which Lift and Project will be applied
int TotalLPcuts;																	//Total number of Lift and Project Cuts applied
//int size_act;																		//Size of the set of active optimality and feasiblity constraints
int Use_ben;				//Flag to know if we will use CPLEX Benders
double *stabilizer;

double Prev_incumbent;																//Value of the last incumbent solution
int       **pos_mu;

/**********************************************Parameters for Pre Processing*********************************************************************/
CPXLPptr  lpPre;      // data strucutre to store a problem in cplex ...................
CPXLPptr  lpclone;      // data strucutre to store a problem in cplex ...................
CPXENVptr envPre;     // cplex environment.............................................
int First_time;		  //Indicator for the presolve stating it's the first time enterin.
double *reducedcost;  //Value of the reduced cost so that we can perform variable elimination...
double true_opt;
double *ub_vec;   //array for the upper bound
double *ypre_obj;
double	RHSLPcutoff;

/*****************TO use within the branch and cut for the heuristic******************************/
int *fixed_0;
int *fixed_1;
int globaldepth;
double *LPsol, LPval;
int    *branchedon;  //Keep track of who's been branched on
double    *lbfix;  //lower bounds of variables at a given node
double	 *ubfix;// upper bounds of variables at a given node
int		previndex;
int maxdepth;

/********************Declaring the arrays of the network***********************************/
double	 *supply;
int		 *tail;
int	     *head;
double	 *obj;
double	 *ub;
double	 *lb;

int flag_incumb;
double *ini_core;
int		added_LP;
double best_lb;

/*************************For CG heuristic**********************/
Kpath_index	*path_indM;//array with the indexes of the paths to write the master problem constraint
int *col_limits;
int *num_paths;
int *comindofcolumn;//array of the commodity index in the column
int *pathcomnumofcolumn;//array of the path number of the commodity index in the column
int **locationinpaths;// This will contain the location of the different paths per commodity per num_path.
int totalinpaths; //counter of the number of paths in the paths structure.
int maxpathperK;
int lastdepth;//keep track of the last depth the heuristic entered in.
int lasttotinpath; //last totalin path
int countimprovedheur;//counter that says how many times an improved solution was found with our heuristic
int counttried;//count the number of times it was tried
int useheur;//use initial heuristic.
double  *forelimintest;//array to keep to check aong the way for the elimination test.
double    *indarcinpath;

/*******To test for Cut and solve******************************/
int   *base_sol;   //Base of the solution
double *base_solroute;
int		globalcounter;
int	   **base_soltrack; //tracking of the solutions
int    *base_soltrackind;
int		heurswitch;
int		count_rep; //Counter after the repetition of objective function value
double  LPbound;
int		rounds;//Number of increasings of RHS
int		Up_lim;
int		Low_lim;
int		Up_Extlim;
int		size_sol;
int		nodemult;
int       *pos_lambda_plus;
int       *pos_lambda_minus;
int       **pos_mu;
double		optBD_prev;
double      newobj;
int			indicateupdate;//indicate update in the cut-and-solve.
int		*indbranch; //indices for braching priority
int     *priority;//priority for branching calculated as the f/c
int     maxprior;//maximum value of priority
int     digits;
int    type_CS;


/******Additional stuff******/
 clock_t heur_checkstart;
 clock_t currheurcheck;
 double  heurchecktime;
 double lowestpath;

 /*************For Commodity relaxations***********/
 int  max_num_dest;				//Integer with the maximum number of different destinations for each commodity
 int  min_num_dest;				//Integer with the minimum number of different destinations for each commodity
 int  com_rel_count;			//Integer counter for the commodity relaxation
 double * com_z_obj;				//Array with the upperbounds considering the commodity relaxations.
 int   limit;
 /*************Additional stuff added on January 20 2019********/
 int LifProunds; //Number of rounds of lift-and-project cuts added
 int *splitvars; //order of split variables
 int *priority;//Integer of priorities.
 int *indprior;//index of priority.

 //For the one percent test
 double onePercentAway; //The value of the objective function at 1% away from optimum
 int  FlagOnepercent; //0 if not there yet 1 if there.

 /******For the new heuristic****/
 double *newinisol; //where we'll keep the new heuristic
 double newglobobj;//objective function value of the new heuristic solution

 /*******Global Reporting*******/
 double		Upper_bound;
 double		best_lower_bound;    //Value of the best lower bound so far
 int		nodecount;
 double		cputimeBD;


 