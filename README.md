# SURROGATE-ASF-code
     SURROGATEASF is a free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.
 
     SURROGATEASF is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
 
     You should have received a copy of the GNU General Public License
     along with SURROGATEASF.  If not, see <http://www.gnu.org/licenses/>.

  SURROGATE-ASF: An interactive surrogate-based method for solving 
  computationally expensive  multiobjective optimization problems of the form 
          min {f_1(x),…f_k(x)}    subject to   LB <= x <= UB
  where 
  f_i, i=1,..k (k=2,3), are computationally expensive objective function 
  whose closed-form formulas are not available (i.e., black-box functions). 
 
  LB and UB are the lower and upper bounds of the decision variables. 
  SURROAGTE-ASF aims at finding the most preferred solution for the 
  decision maker (user) within a  limited number of function evaluations. 
  SURROGATE-ASF consists of two phases: i.e., initialization and 
  decision making phases. In the initialization phase, SURROGATE-ASF employs 
  a surrogate-based single  objective optimization method to get some 
  information regarding the location of the solutions in the 
  decision/objective space. Then the original problem is decomposed into 
  a finite number of single-objective surrogate problems. In the decision 
  making phase, the interaction with the decision maker is conducted to find 
  the most preferred solution. The decision maker expresses his/her 
  preferences in the form of a reference point. 
   
  Input:
  The multiobjective optimization problem is defined in “P_objective.m”. 
  P_objective.m is a part of the implementation of RVEA method developed in 
 
  R. Cheng, Y. Jin, M. Olhofer and B. Sendhoff, A Reference Vector Guided 
  Evolutionary Algorithm for Many-objective Optimization, IEEE Transactions
  on Evolutionary Computation, 2016.
 
  In this implementation of SURROGATE-ASF, two surrogate-based methods 
  developed for computationally expensive single-objective optimization 
  problems are incorporated as listed below:
 
  MATSuMoTo developed in J. Muller and C. A. Shoemaker. Influence of ensemble 
  surrogate models and sampling strategy on the solution quality of algorithms 
  for computationally expensive black-box global optimization problems. 
  Journal of Global Optimization, 60(2):123-144, 2014.
 
  and
 
  StochasticRBF developed in R.G. Regis and C.A. Shoemaker, A Stochastic 
  Radial Basis Function Method for the Global Optimization of Expensive 
  Functions, INFORMS Journal on Computing, vol. 19, pp. 497-509, 2007 
 
  and 
 
  R.G. Regis and C.A. Shoemaker, Parallel Stochastic Global Optimization 
  Using Radial Basis Functions, INFORMS Journal on Computing, vol. 21, pp. 
  411-426, 2009.
 
  In order to solve single-objective surrogate problems, the DIRECT method 
  developed in 
 
  D. R. Jones, C. D. Perttunen, and B. E. Stuckman. Lipschitzian optimization 
  without the Lipschitz constant. Journal of Optimization Theory and Applications,
  79(1):157-181, 1993.” is incorporated. 

   Feel free to replace these methods with other methods. When using 
   SURROGATE-ASF, please cite the following paper:   
   
  Mohammad Tabatabaei, Markus Hartikainen, Karthik Sindhya, Jussi Hakanen, 
  Kaisa Miettinen, An Interactive Surrogate-based Method for Computationally 
  Expensive Multiobjective Optimization, submitted.
