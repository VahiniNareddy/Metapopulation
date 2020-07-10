//This code is written for metapopulation Model A with Ricker dynamics (growth parameter, r=2.2) in the sequence NDR. 

//inputs: {random_seed, lattice_size, noise:lambda, dispersal:epsilon}
//outputs: {initial_state configuration, 3 final states configurations (t, t+1, t+2), observables (lambda, order_parameter, Binder_cumulant)}
//functions: {neighbor table, run_simulation, generate_initial_state, generate_next_state, calculate order parameter, output observables, output state configurations}


#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

//variables
int warm_up; //Monte Carlo sweeps for the system to reach stationary state (warming up the system)
int max_mcs; // Total Monte Carlo sweeps
int nsize; // Lattice size
int xsize; // Length of the lattice
int ysize; // Breadth of the lattice
int seed; // Random seed for the run
gsl_rng * r;

//constants
double R=2.2; //Growth parameter of the Ricker map
double epsilon; // Dispersal
double lambda; // Noise 

//vectors
vector<double> basis_state_old; //Lattice configuration at time t-1
vector<double> basis_state_new; //Lattice configuration at time t
vector<double> basis_state_R; //Lattice configuration after applying Ricker process at time t
vector<double> basis_state_D; //Lattice configuration after applying Dispersal process at time t
vector<double> zeta; // Noise valuesat time t
vector< vector<int> > nn_table; //nearest neighbor table pointing to the four neighbors of each lattice site 
vector<double> lambda_values; //Lambda values for the observables
vector<double> orderparameter_values; //Order parameter values
vector<double> susceptibility_values; //Susceptibility values
vector<double> binder_values; //Binder cumulant values
 

//function prototypes
void gsl_rng_set (const gsl_rng * r, unsigned long int s);
void make_nntable();
void run_simulation(const double& lambda,const double& epsilon);
void generate_initial_state(void);  
void generate_next_state(const int& mcs,const double& lambda, const double& epsilon);
double get_orderparameter(const int& mcs);
void print_basisstates(const int& mcs, const double& lambda,const double& epsilon);
void print_observable_values();

int main(int argc, char** argv){

	gsl_rng_env_setup();
        r = gsl_rng_alloc(gsl_rng_default);

        seed=atoi(argv[1]);
	gsl_rng_set(r,seed);
	
	xsize = atoi(argv[2]); //Lattice size input
	ysize = xsize; //Length and breadth of the lattice is set to have same value
        nsize = xsize * ysize;
	
	nn_table.resize(nsize,vector<int>(4));
	basis_state_old.resize(nsize);
	basis_state_new.resize(nsize);
	basis_state_R.resize(nsize);
        basis_state_D.resize(nsize);
	zeta.resize(nsize);
	make_nntable();

	
	warm_up = 2000000;
	max_mcs = 3000005;
	
	
	lambda=atof(argv[3]);	//Noise input
	int no_of_values_lambda=1;
	double interval_lambda=0;

	epsilon = atof(argv[4]);  //Dispersal input
	int no_of_values_epsilon = 1;
	double interval_epsilon = 0;

	for(int i=0;i<no_of_values_lambda;++i){
	                                      for(int j=0;j<no_of_values_epsilon;++j){
	                                      					      run_simulation(lambda,epsilon);  //Monte Carlo simulation for given dispersal and noise
	   									      lambda+= interval_lambda;
	   									      epsilon+= interval_epsilon;
	   									      print_observable_values();  //Print observables at the end of the simulation
					      }
	}
	gsl_rng_free(r);	
}



//Function that runs metapopulation dynamics and calculates observables
void run_simulation(const double& lambda, const double& epsilon){
  //initialization
  double order_sum = 0;
  double order_sum_sqr = 0;
  double order_sum_four = 0;
  int loop_num = 0;
  int loop_num1 = 1;
  int loop_num2 = 2;
    
  
  //Creating initial state
  generate_initial_state();
  
  print_basisstates(0,lambda,epsilon);

  //Generating states with metapopualtion dynamics and calculating observables after the system reaches stationary state
  for (int mcs=1; mcs<max_mcs; ++mcs) {                                               
       				       generate_next_state(mcs,lambda,epsilon);
   
                                       if(mcs>=warm_up){
         						loop_num++;
	 						loop_num1++;
	 						loop_num2++;
         						double getorder = get_orderparameter(mcs);
         						order_sum+= getorder;
	          				        order_sum_sqr+=getorder*getorder;
	                                                order_sum_four+= getorder*getorder*getorder*getorder;
							// printing configurations at 0^6 runs after warm up
	 						if(loop_num%1000000==0){print_basisstates(mcs,lambda,epsilon);} 
         						if(loop_num1%1000000==0){print_basisstates(mcs,lambda,epsilon);}
	 						if(loop_num2%1000000==0){print_basisstates(mcs,lambda,epsilon);}
       				      }
  }
  double avg_order= order_sum/(loop_num);  //Averaged order parameter over time
  double binder = 1 - loop_num * order_sum_four/(3*order_sum_sqr*order_sum_sqr); //Averaged Binder cumulant
  double avg_suscep = order_sum_sqr*nsize/(loop_num)- nsize*avg_order*avg_order; //Averaged susceptibility

  //Saving observables in vectors
  lambda_values.push_back(lambda);
  orderparameter_values.push_back(avg_order);
  susceptibility_values.push_back(avg_suscep);
  binder_values.push_back(binder);
}



//initial state
void generate_initial_state(void){
	for (int j=0; j<nsize; ++j){basis_state_new[j]=2*gsl_rng_uniform(r);}
}



//Generating state at time t
void generate_next_state(const int& mcs, const double& lambda, const double& epsilon){
	double rdn=gsl_ran_gaussian(r, 1);

        for (int i=0; i<nsize; i++){basis_state_old[i] = basis_state_new[i];
 				    zeta[i]=gsl_ran_gaussian(r,1);
				    basis_state_R[i]=basis_state_old[i] *exp(R*(1-basis_state_old[i]));  //State after applying Ricker dynamics
        }
	
	for (int i=0; i<nsize; i++){
	    			    double nnfsum = 0;
            			    for(int j=0; j<4; j++){
	       			                           int k = nn_table[i][j];
	                                                   nnfsum +=  (epsilon/4)* basis_state_R[k];  //Summing over neighbors
	                            }
	    			    basis_state_D[i] = ((1-epsilon)* basis_state_R[i])+ nnfsum;    //State after Dispersal
	    			    basis_state_new[i]=basis_state_D[i]*exp(lambda*zeta[i]);    //New state at time t
        }
}

//Calculating order parameter 
double get_orderparameter(const int& mcs){
	double msum = 0;
	for(int i=0; i< nsize; ++i){msum+= basis_state_new[i]-basis_state_old[i];}
	return double(abs(msum/(2*nsize)));
}



//Printing state configuration at time mcd
void print_basisstates(const int& mcs, const double& lambda, const double& epsilon){
	ofstream myfile;
	myfile.open("state_size_"+to_string(nsize)+"_epsilon_"+to_string(epsilon)+"_seed_"+to_string(seed)+"_lambda_"+to_string(lambda)+"_mcs_"+to_string(mcs)+".txt");
	for (int i=0; i<nsize; i++){
                                    myfile << basis_state_new[i] << endl;
	}
	myfile.close();
}




//Printing observables
void print_observable_values(){
	ofstream observablefile;
	observablefile.open("observables_R_"+to_string(R)+"_epsilon_"+to_string(epsilon)+"_nsize_"+to_string(nsize)+"_seed_"+to_string(seed)+".txt");
	observablefile<< "size="<<nsize<<endl;
	for (int i=0; i<lambda_values.size();i++){
        					observablefile << lambda_values[i]<<"   "<<orderparameter_values[i]<<"   "<<binder_values[i]<<endl;
	}
	observablefile.close();
}


//Neighbor table to get the id of the four neighbors while the Dispersal process
//Two dimensional lattice is numbered as one dimensional chain

void make_nntable(){
	//neighbors:id {left:0 right:1 top:2 bottom:3}
	for (int i=0; i<nsize; i++){nn_table[i][0]=i-1;}
	for (int i=0; i<nsize; i++){nn_table[i][1]=i+1;}
	for (int i=0; i<nsize; i++){nn_table[i][2]=i+xsize;}
	for (int i=0; i<nsize; i++){nn_table[i][3]=i-xsize;}
	//periodic boundaries
	for (int i=0; i<xsize; i++){nn_table[i][3]=i-ysize+nsize;}
	for (int i=0; i<nsize; i += xsize){nn_table[i][0]=i+xsize-1;}
	for (int i=xsize-1; i<nsize; i += xsize){nn_table[i][1]=i+1-xsize;}
	for (int i=nsize-xsize; i<nsize; i++){nn_table[i][2]=i+ysize-nsize;}
}


