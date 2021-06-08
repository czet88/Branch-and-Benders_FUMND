LIB = ../../../../encs/pkg/cplex-12.7.1/root/cplex
CXX=gcc -O2
CPLXST=-DIL_STD  -L$(LIB)/lib/x86-64_linux/static_pic -lcplex -lm -lpthread 
LIB2=-I$(LIB)/include/
OBJ= fundamentals.o solve_MIP.o  solve_MasterMC.o solve_Subproblem.o MetaheuristicWarmStart.o Preprocess.o
NAME=Ben

$(NAME) : $(OBJ)
	$(CXX) -I$(LIB)/include/ilcplex $(OBJ) -o $(NAME) $(CPLXST)

fundamentals.o: fundamentals.c def.h
	$(CXX) -c fundamentals.c $< $(LIB2)
solve_MIP.o : solve_MIP.c def.h
	$(CXX) -c solve_MIP.c $< $(LIB2)
solve_MasterMC.o : solve_MasterMC.c def.h
	$(CXX) -c solve_MasterMC.c $< $(LIB2)
solve_Subproblem.o : solve_Subproblem.c def.h
	$(CXX) -c solve_Subproblem.c $< $(LIB2)
MetaheuristicWarmStart.o: MetaheuristicWarmStart.c def.h
	$(CXX) -c MetaheuristicWarmStart.c $< $(LIB2)
Preprocess.o: Preprocess.c def.h
	$(CXX) -c Preprocess.c $< $(LIB2)

clean :
	rm $(NAME) $(OBJ)

