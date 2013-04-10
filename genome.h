#include <iostream>
#include <vector>
#include <utility>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_histogram.h>
#include <map>
#include <boost/lexical_cast.hpp>

#define VERBOSE 0
#define RNG gsl_rng_taus2
#define WRITEHIST 0

using namespace std;

struct fixation_event{
    int t_fix; //time of fixation
    int loci; //location 
    fixation_event(int t_fix_in, int loci_in);
};

struct block{
    int l; //Left index
    int r; //Right index
    int dat; //Interior of Block (integer)
    double block_fit;
    map<int,double> fit_loci; //<location, locus fitness>
    bool contains_locus(int locus);
    int get_width();
    bool operator==(const block &b2) const {return ((r == b2.r) && (l == b2.l) && (dat == b2.dat));}
    bool operator!=(const block &b2) const {return ((r != b2.r) && (l != b2.l) && (dat != b2.dat));}
    block();
    block(int l_in, int r_in, int dat_in);
    block(int l_in, int r_in, int dat_in, map<int,double> &block_fit_in);
};

class genome {
    protected:
        gsl_rng * rng;
        //int seed;
        //int get_random_seed();
        double genome_fit; //Genome Fitness
        //int L; //Number of Loci
        int find_block_index(int locus, int min, int max); 
    public:
        std::vector< block > g; //genome vector of isoancestral blocks
        int L;
        //vector < block > crossBlocks(block b1, block b2, int loc1, int loc2, int crossPoint);
        bool containsLocus(int ancestor, int locus);
        vector<double> cumulant_fit(); //returns cumulant fitness
        double sub_cumulant_fit(int start_locus, int l);
        double get_fit();
        //Operator Overload
        genome operator+(genome g2);
        friend ostream& operator<< (ostream &out,const genome &g);
        //Constructors & Destructors
        genome();
        genome(int c_in, int L_in, gsl_rng* rng_in);
        genome(int c_in, int L_in, gsl_rng* rng_in, map<int,double> &fit_in);
        genome(gsl_rng* rng_in, vector< block > g_in);
        ~genome();
};

class population {
    protected:
        gsl_rng * rng;
        int seed;
        double avgFit;
        vector <double> selective_weight;
        vector < vector<int> > genetic_weight;
        vector <vector<fixation_event> > fixations;
        int N; 
        int L;
        double r;
        int gen_pop;
        bool neutral;
        //Functions
        pair<int,int> selectParents(); //Randomly selects two parents from Ancestral Generation
        genome progeny(pair<int,int>); //From a given set of parents, calculates a recombinative offspring.
        int binarySearch_int(vector<double> &data, int min, int max, double &key);
        map<int,double> convert_vector_toMap(vector<double>);
        void update_selection_weights();
        void update_weights(); //Calculates the weight function for all loci for all ancestors.
        void updateBlockSizes(); //Get Histogram of Block sizes
        //int insert_fixation_loc(fixation_event &f,int min, int max); 
        int get_random_seed();
    public:
        gsl_histogram* blockHist;
        std::vector< genome > pop; 
        //Overloaded Constructors
        population(int N_in, int L_in, double r_in, int seed_in, gsl_rng* rng_in);
        population(int N_in, int L_in, double r_in, int seed_in, gsl_rng* rng_in, vector< vector< double > >& fit);
        //Destructors
        ~population();
        //Functions
        void evolve();
        void evolve(int gen); //Time evolves a population
        gsl_histogram* get_cumulant_fit_hist(int l, double sigma);
        double get_avg_fit();
        double get_variance();
        double get_allele_frequency(int ancestor, int locus);
        std::vector<int> get_weight(int ancestor); //Calculates the weight function for all loci for ancestral individual 'ancestor'.
        vector<vector<int> > get_all_weight();
        vector<pair<fixation_event, int> > get_flat_fixations();
        vector< vector<fixation_event> > get_fixations();
        void write_genetic_weight(const char* s);
        void write_fixations(const char* s);
        void writeBlockHist(const char* s);
};


