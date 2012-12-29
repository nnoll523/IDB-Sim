#include "genome.h"

// BLOCK FUNCTIONS
block::block() {
    len = 10;
    dat = 1;
}

block::block(int len_in, int dat_in) {
    len = len_in;
    dat = dat_in;
}

bool block::containsCross(int position, int crossPoint) {
    return (position <= crossPoint && crossPoint < len+position);
}

// GENOME FUNCTIONS

// Overloaded Constructors

genome::genome() {
    size = 10;
    fitness = 0;
    g.push_back(block(10,0));  //entire array of homogeneous genes
    seed = get_random_seed();
    rng = gsl_rng_alloc(RNG);
    gsl_rng_set(rng,seed);
}

genome::genome(int c_in, int L_in) {
    size = L_in;
    fitness = 0;
    g.push_back(block(L_in, c_in));  //entire array of homogeneous genes
    seed = get_random_seed();
    rng = gsl_rng_alloc(RNG);
    gsl_rng_set(rng,seed);
}

genome::genome(int c_in, int L_in, double fit_in) {
    size = L_in;
    fitness = fit_in;
    g.push_back(block(L_in, c_in));  //entire array of homogeneous genes
    seed = get_random_seed();
    rng = gsl_rng_alloc(RNG);
    gsl_rng_set(rng,seed);
}

genome::genome(vector < block > g_in) {
    g = g_in;
    size = 0;
    fitness = 0;
    for (int i=0; i<g.size(); i++) {
        size += g[i].len;
    }
    seed = get_random_seed();
    rng = gsl_rng_alloc(RNG);
    gsl_rng_set(rng,seed);

}
genome::~genome(){
    if (VERBOSE)
        cout << "Destructing genome...\n"; 
}

int genome::get_random_seed() {
	int seedtmp;
	ifstream urandom("/dev/urandom", ios::binary);
	if(urandom.bad()) {
		cerr<<"/dev/urandom gives bad stream, falling back to time + getpid + number_of_instances"<<endl;
		seedtmp = time(NULL) + getpid();
	} else {
		urandom.read(reinterpret_cast<char*>(&seedtmp),sizeof(seedtmp));
		urandom.close();
	}
	return seedtmp;
}

bool genome::containsLocus(int ancestor, int locus) {
    try {
        if (locus < 0 || locus >= size) {throw "Invalid Locus address";}
        int position = 0;
        for (int i = 0; i<g.size(); i++) {
            if (position > locus) {break;} //If we have traversed far enough and found no match, no sense in continuing walking.
            if (g[i].dat == ancestor && g[i].containsCross(position,locus)) {return true;}
            position += g[i].len;
        }
        return false;
        }
    catch (const char* Message) {
        cout << "Error:" << Message << "\n";
    }    
}

double genome::getFitness() {return fitness;}

genome genome::operator+(genome g2){
    try {
        if (g2.size != size) {
            throw "Genome sizes are not equal!";
        }
        else {
            vector<block> offspring; 
            gsl_rng_set(rng,seed);
            int crossPoint = (int)gsl_rng_uniform_int(rng,size-1); //Obtain the random crossover point-convention is that it occurs to right 
            //Decide which out of the two recombinative gametes will be returned
            int bin = (int)gsl_rng_uniform_int(rng,2);
            //Compute Crossover
            int loc1 = 0; int loc2 = 0; //Start at beginning.
            int i = 0; int j = 0; //Block Iterator
            if (bin == 0){
                while(!g[i].containsCross(loc1,crossPoint) || !g2.g[j].containsCross(loc2,crossPoint)) {
                    //Create offspring genome if i is still iterating.
                    if (i>=j && !g[i].containsCross(loc1,crossPoint)) { //If i stops iterating, it will immediately become less than j.
                        offspring.push_back(g[i]);
                    }
                    //Iterate through both genomes simultaneously, provided they don't individually contain the crossover point.
                    if (!g[i].containsCross(loc1,crossPoint)) {
                        loc1 += g[i].len;
                        i++; 
                    }
                    if (!g2.g[j].containsCross(loc2,crossPoint)) {
                        loc2 += g2.g[j].len;
                        j++;
                    }
                }
                vector<block> newBlocks = crossBlocks(g[i],g2.g[j], loc1, loc2, crossPoint); 
                for (int k = 0; k < newBlocks.size(); k++) {
                    offspring.push_back(newBlocks[k]);
                }   
                //Traverse the second genome until it ends. 
                for (int k = j+1; k < g2.g.size(); k++) {
                    offspring.push_back(g2.g[k]); //Transmit block from second parent.
                } 
                try {
                    if (genome(offspring).size != size) {
                        throw "Offspring genome is not equal to parental!";
                    }
                    return genome(offspring);
                }
                catch(const char* Message) {
                    cout << "Error:" << Message <<"\n";
                }
            }   
            else{
                //Same algorithm - flip roles.
                while(!g[j].containsCross(loc2,crossPoint) || !g2.g[i].containsCross(loc1,crossPoint)) {
                    //Create offspring genome if i is still iterating.
                    if (i>=j && !g2.g[i].containsCross(loc1,crossPoint)) { //If i stops iterating, it will immediately become less than j.
                        offspring.push_back(g2.g[i]);
                    }
                    //Iterate through both genomes simultaneously, provided they don't individually contain the crossover point.
                    if (!g2.g[i].containsCross(loc1,crossPoint)) {
                        loc1 += g2.g[i].len;
                        i++; 
                    }
                    if (!g[j].containsCross(loc2,crossPoint)) {
                        loc2 += g[j].len;
                        j++;
                    }
                }
                vector<block> newBlocks = crossBlocks(g2.g[i], g[j], loc1, loc2, crossPoint); 
                for (int k = 0; k < newBlocks.size(); k++) {
                    offspring.push_back(newBlocks[k]);
                } 
                //Traverse the second genome until it ends. push_back all blocks.
                for (int k = j+1; k < g.size(); k++) {
                    offspring.push_back(g[k]); //Transmit block for second parent.
                } 
                try {
                    if (genome(offspring).size != size) {
                        throw "Offspring genome is not equal to parental!";
                    }
                    return genome(offspring);
                }
                catch(const char* Message) {
                    cout << "Error:" << Message <<"\n";
                }
            }
        }
    }
    catch(const char* Message) {
        cout << "Error:" << Message << "\n";
    }
}

ostream& operator<< (ostream &out, const genome &gen) {
    for(int i=0; i<gen.g.size(); i++) {
        for(int j=0; j<gen.g[i].len; j++) {
            out << gen.g[i].dat;
        }
    }
    out << "\n";    
    return out;
}

vector< block > genome::crossBlocks(block b1, block b2, int loc1, int loc2, int crossPoint) {
    vector <block> v;
    if (b1.dat == b2.dat) {
        //Build single block. Corresponds to a merger event.
        v.push_back(block(loc2+b2.len-loc1, b1.dat));
    }
    else {
        //Build two blocks. crossPoint-loc1 >= 0. loc2+b2.len-1-crossPoint
        v.push_back(block((crossPoint+1-loc1),b1.dat));
        if (loc2+b2.len-1 > crossPoint) { //If crossover point is to the right of the entire second block - don't do anything.
            v.push_back(block((loc2+b2.len-1-crossPoint),b2.dat));
        }
    }
    return v;
}

//POPULATION METHODS

//Overloaded Constructors
population::population() {
    size = 1;
    loci = 10;
    avgFit = 0;
    stream = fopen("out.txt","w");
    pop.push_back(genome(1,10));
    seed = get_random_seed();
    rng = gsl_rng_alloc(RNG);
    gsl_rng_set(rng,seed);
    blockHist = gsl_histogram_alloc(loci);
    updateBlockSizes();
}
population::population(int N) {
    size = N;
    loci = 10;
    avgFit = 0;
    stream = fopen("out.txt","w");
    for(int i=0; i<size;i++) {
        pop.push_back(genome(i,10));
    }
    seed = get_random_seed();
    rng = gsl_rng_alloc(RNG);
    gsl_rng_set(rng,seed);
    blockHist = gsl_histogram_alloc(loci);
    updateBlockSizes();
}
population::population(int N, int L) { //Neutral Evolution!
    size = N;
    loci = L;
    avgFit = 0;
    stream = fopen("out.txt","w");
    for(int i=0;i<size;i++) {
        pop.push_back(genome(i,L));
    }
    seed = get_random_seed();
    rng = gsl_rng_alloc(RNG);
    gsl_rng_set(rng,seed);
    blockHist = gsl_histogram_alloc(loci);
    updateBlockSizes();
}
population::population(int N, int L, vector< vector<double> > fit) {
    size = N;
    loci = L;
    avgFit = 0;
    stream = fopen("out.txt","w");
    try {
        if (fit.size() != N || fit[0].size() != loci) {
            throw "Incorrect Fitness Landscape dimensions!";
        }
        vector< double > initFit (size); //Create a temp vector that will be discarded after class has been constructed that holds fitness.
        //Read in fitness landscape. Right now stored in a hash table - only beneficial for a sparse landscape.
        for (int i=0; i<fit.size(); i++){
            for (int j=0; j<fit[i].size(); j++) {
                if (fit[i][j] != 0) {
                    avgFit += fit[i][j];
                    initFit[i] += fit[i][j];
                    fitnessL.insert(pair< pair< int,int >,double >(pair<int,int>(i,j),fit[i][j])); //Create a map entry if non-zero.
                }
            }
        }
        avgFit/= size; //Total fitness divided by the number of individuals gives average fitness.
        //Initialize Population with corresponding fitness values.
        for(int i=0;i<size;i++) {
            pop.push_back(genome(i,L, initFit[i]));
        }
        seed = get_random_seed();
        rng = gsl_rng_alloc(RNG);
        gsl_rng_set(rng,seed);
        blockHist = gsl_histogram_alloc(loci);
        updateBlockSizes();
    }
    catch (const char* Message) {
        cout << "Error:" << Message << "\n";
    }
}
population::~population(){
    if (VERBOSE)
        cout << "Destructing Population...\n"; 
    gsl_histogram_free(blockHist); //Remove usage of histogram!
    fclose(stream);    
}

//Functions
void population::writeBlockHist(){
    gsl_histogram_fprintf(stream,blockHist,"%g","%g");
}

genome population::progeny(pair<int,int> parents) { //Creates an offspring genotype
    genome parent1 = pop[parents.first];
    genome parent2 = pop[parents.second];
    return parent1+parent2;
}
void population::updateBlockSizes(){
    gsl_histogram_set_ranges_uniform(blockHist,1,loci+1); //Questionable range choice. As it stands now, the bins are partitioned as follows {[1,2),[2,3)......[L,L+1)}. Thus if all have length L, then average will come out to be (L+.5). Will have to rescale.
    for (int i=0;i<size;i++){
       for (int j=0; j<pop[i].g.size(); j++) {    
            gsl_histogram_increment(blockHist,pop[i].g[j].len);
        }    
    }
}

vector<int> population::getWeight(int ancestor){
    try {
        if (ancestor >= size || ancestor < 0) {
            throw "Incorrect Ancestral Value!";
        }
        else {
            vector<int> weight (loci);
            for (int n=0; n<size; n++) {
                int loc = 0;
                for (int i=0; i<pop[n].g.size(); i++) {
                    if (pop[n].g[i].dat == ancestor) {
                        for (int j=loc; j<(loc+pop[n].g[i].len); j++) {
                            weight[j]++;
                        }
                    }
                    loc += pop[n].g[i].len - 1;
                }
            }
            return weight;
        }
    }
    catch (const char* Message) {
        cout << "Error:" << Message << "\n";
    }
}

vector< vector<int> > population::getAllWeights(){
    vector< vector<int> > temp (size);
    for (int i=0; i<size; i++){
        temp[i] = getWeight(i);
    }
    return temp;
}

void population::evolve(int gen) { //Wright-Fisher Model - Population size is held constant. Parental genotypes are chosen randomly. 
    for (int j=0; j<gen; j++){
        vector <genome> newPop(size);
        avgFit = 0;
        for(int i=0; i<size; i++) {
            newPop[i] = progeny(selectParents());
            for (map<pair<int,int>,double>::const_iterator it = fitnessL.begin(); it != fitnessL.end(); it++) {
                if (newPop[i].containsLocus(it->first.first, it->first.second)) {
                    newPop[i].fitness += it->second;
                    avgFit += it->second;
                }
            }
        }
        pop = newPop;
        avgFit/=size;
    }
    updateBlockSizes();
    if (WRITE) {writeBlockHist();}
}

void population::evolve() { //Wright-Fisher Model - Population size is held constant. Parental genotypes are chosen randomly. 
   evolve(1); 
}

double population::getAverageFit() {
    return avgFit;
}

pair<int,int> population::selectParents(){ //Parents are selected at random from population 
    pair<int,int> ancestors;
    if (fitnessL.empty()) { //Neutral Evolution
        ancestors.first = (int)gsl_rng_uniform_int(rng,size);
        ancestors.second = (int)gsl_rng_uniform_int(rng,size);
        //Recursive definition to ensure unique parents.
        if (ancestors.second == ancestors.first) {
            return selectParents();
        }
        else {
            return ancestors;
        }
    }
    else {
        //Need an algorithm that weights individuals based upon their fitness.
    }
}

int population::binarySearch_int(vector<double> &data, int min, int max, double key) { //Data is assumed to be sorted in ascending order!
    try {
        if (min > max) {throw "Encountered invalid search range";} //Something terrible has happened.
        int midpoint = min + ((max-min)/2); //Obtain the midpoint using integer arithmetic.
        if (midpoint > 0) {
            if (key > data[midpoint]) {binarySearch_int(data, midpoint+1, max, key);} //If the key is larger than the current bucket value, then the address must be to the right of our current position. 
            else if (key <= data[midpoint-1]) {binarySearch_int(data, min, midpoint-1, key);} //If the key is smaller than the bucket to our immediate left, then we know the address must be to the left of our current position.
            else {return midpoint;} //If other two cases fail, then the key must be located in our current bucket. 
        }
        else {return 0;} //If the midpoint is zero, search should terminate. 
    }
    catch (const char* Message) {
        cout << "Error:" << Message << "\n";
    }

}

int population::get_random_seed() {
	int seedtmp;
	ifstream urandom("/dev/urandom", ios::binary);
	if(urandom.bad()) {
		cerr<<"/dev/urandom gives bad stream, falling back to time + getpid + number_of_instances"<<endl;
		seedtmp = time(NULL) + getpid();
	} else {
		urandom.read(reinterpret_cast<char*>(&seedtmp),sizeof(seedtmp));
		urandom.close();
	}
	return seedtmp;
}
void dump_map(const std::map<pair<int,int>, double>& map) {
    for (std::map<pair<int,int>,double>::const_iterator it = map.begin(); it != map.end(); it++) {
        //cout << "Key: " << it->first << endl;
        cout << it-> second << endl;
    }
}
int main() {
    /*
    genome g0(0,10);
    genome g1(1,10);
    cout << g0;
    cout << g1;
    cout << g0.containsLocus(0,0) << "\n";
    cout << g1.containsLocus(0,0) << "\n";
    genome g2 = g0 + g1;
    cout << g2;
    cout << g2.containsLocus(0,6) << "\n";
    cout << g2.containsLocus(0,3) << "\n";
    
    genome g2 = g0 + g1;
    genome g3 = g2 + g0;
    genome g4 = g2 + g1;
    cout << g2;
    cout << g3; 
    cout << g4;
    
    cout << "(g1 + g2).g";
    for(int i=0; i<10;i++){
         cout << g2.g[i];
    }
    cout << '\n';
    cout << "(g3 + g2).g";
    for(int i=0; i<10;i++){
         cout << g4.g[i];
    }
    cout << '\n';
    
    population pop(10000,100000);
    cout << gsl_histogram_mean(pop.blockHist) << "\n";
    pop.evolve(2);
    vector<int> test = pop.getWeight(90);
    for (int i=0; i<test.size(); i++) {
        cout << test[i] << " ";
    }
    cout << "\n";
    cout << gsl_histogram_mean(pop.blockHist);
    cout << '\n';
    vector<int> test2 = pop.getWeight(50);
    for(int i=0; i<10; i++){
        cout << test2[i] << ",";
    }
    cout << '\n';
    
    for (int i=0;i<200;i++){
        pop.evolve(5);
    }
    */
    vector< vector<double> > test (10);
    for (int i=0; i<10; i++) {
        vector<double> temp (10);
        for (int j=0; j<10; j++) {
            temp[j] = i;
        }
        test[i] = temp;
    }
    population pop(10,10,test);
    /*
    pop.evolve(2);
    for (int i=0; i<10; i++) {
        cout << pop.pop[i] << " " << pop.pop[i].getFitness() << "\n";
    }
    cout << pop.getAverageFit() << "\n";
    */
    //dump_map(pop.fitnessL);
    //cout << "\n";
    return 0;
}

