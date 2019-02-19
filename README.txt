*Use of Benders executable ./Ben-  For execution in Linux
from the paper-"Exact algorithms based on Benders decomposition for the multicommodity uncapacitated fixed-charge network design"

Step 0
The Makefile must be updated to include the path to CPLEX in your local computer.

Step 1
Copy codes into your computer and type make to generate the Ben executable and other programming objects

Step 2
To execute, type "./Ben $NAME_OF_INPUT  $NAME_OF_OUTPUT"

/*********************************/
Structure of the input file

Number_of_instances
Name_of_instance   Solution_Method     

/*********************************/

Solution_Methods baremo
1-    Cplex's branch-and-cut on the disaggregated formulation
2-    Cplex's blackbox Benders implementation on the disaggregated formulation
3-    Branch-and-Benders cut algorithm without enhancements
4-    Branch-and-Benders cut algorithm with in-tree heuristic
5-    Branch-and-Benders cut algorithm with Lift-and-project cuts
6-    Branch-and-Benders cut algorithm with both in-tree heuristic and lift-and-project cuts
7-    Benders cut-and-solve. 