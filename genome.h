#include <iostream>
#include <vector>
#include <utility>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_histogram.h>
#include <map>

#define VERBOSE 0
#define RNG gsl_rng_taus2
#define WRITE 0

using namespace std;

struct block{
        int len; //Length 0<ln<L
        int dat; //Interior of Block (integer)
        block();
        block(int len_in, int dat_in);
        bool containsCross(int position, int crossPoint);
        bool operator==(const block &b2) const {return (len == b2.len);}
        bool operator!=(const block &b2) const {return (len != b2.len);}
        bool operator<(const block &b2) const {return (len < b2.len);}
        bool operator>(const block &b2) const {return (len > b2.len);}
};

class genome {
    protected:
        gsl_rng * rng;
        int seed;
        int get_random_seed();
        double fitness; //Genome Fitness
        int size;
    public:
        friend class population;
        std::vector< block > g; //genome vector of isoancestral blocks
        //Constructors & Destructors
        genome();
        genome(int c_in, int L_in);
        genome(int c_in, int L_in, double fit_in);
        genome(vector< block > g_in);
        ~genome();
        vector < block > crossBlocks(block b1, block b2, int loc1, int loc2, int crossPoint);
        bool containsLocus(int ancestor, int locus);
        double getFitness();
        //Operator Overload
        genome operator+(genome g2);
        friend ostream& operator<< (ostream &out,const genome &g);
};

class population {
    protected:
        gsl_rng * rng;
        int seed;
        int get_random_seed();
        //Functions
        pair<int,int> selectParents(); //Randomly selects two parents from Ancestral Generation
        genome progeny(pair<int,int>); //From a given set of parents, calculates a recombinative offspring.
        void writeBlockHist();
        double avgFit;
        map< pair<int,int>, double > fitnessL; //Static Fitness Landscape
        map< pair<int,int>, double>::iterator it; //Pointer to each hash entry.
        int size; //Size N
        int loci;
        gsl_histogram* blockHist;
        FILE *stream;
    public:
        std::vector< genome > pop; 
        //Overloaded Constructors
        population();
        population(int N);
        population(int N, int L);
        population(int N, int L, vector< vector< double > > fit);
        //Destructors
        ~population();
        //Functions
        void evolve();
        void evolve(int gen); //Time evolves a population
        double getAverageFit();
        std::vector<int> getWeight(int ancestor); //Calculates the weight function for all loci for ancestral individual 'ancestor'.
        std::vector< vector<int> > getAllWeights(); //Calculates the weight function for all loci for all ancestors.
        void updateBlockSizes(); //Get Histogram of Block sizes
};

