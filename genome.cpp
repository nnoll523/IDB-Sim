#include "genome.h"

fixation_event::fixation_event(int t_fix_in, int loci_in) {
    t_fix = t_fix_in;
    loci = loci_in;
}

// BLOCK FUNCTIONS
block::block() {
    l = 0;
    r = 1;    
    dat = 1;
}

block::block(int l_in, int r_in, int dat_in) {
    try {
        if (l_in > r_in) {throw "Nonsensical position elements. L<R!";}
        l = l_in;
        r = r_in;
        dat = dat_in;
        block_fit = 0;
    }
    catch (const char* Message) {cout << "Error:" << Message << "\n";}
}

block::block(int l_in, int r_in, int dat_in, map<int,double> &block_fit_in) {
    try {
        if (l_in > r_in) {throw "Nonsensical position elements. L<R!";}
        else if (block_fit_in.size() > r_in - l_in + 1 || (block_fit_in.rbegin()->first > r_in) || (block_fit_in.begin()->first) < l_in) {throw "Incorrect fitness array dimensionality!";}
        else {
            l = l_in;
            r = r_in;
            dat = dat_in;
            block_fit = 0;
            if (!block_fit_in.empty()) {
                fit_loci.swap(block_fit_in); //Reads in fitness spatial landscape.
                for (map<int,double>::const_iterator it = fit_loci.begin(); it != fit_loci.end(); it++) {block_fit += it->second;}
            }
        /*for (int i=0; i<=(r-l); i++) {
            if (block_fit_in[i+l] != 0) {
                fit_loci.insert(pair<int,double> ((i+l),block_fit_in[i]));
                block_fit += block_fit_in[i];
            }
        }
        */
        }
    }
    catch (const char* Message) {cout << "Error:" << Message << "\n";}
}

bool block::contains_locus(int locus) {return (l <= locus && locus < r);}

int block::get_width() {return (r - l + 1);}

// GENOME FUNCTIONS

// Overloaded Constructors
genome::genome() {
    L = 10;
    genome_fit = 0;
    g.push_back(block(0,9,1));  //entire array of neutral homogeneous genes
    rng = gsl_rng_alloc(RNG);
    gsl_rng_set(rng,0);
}

genome::genome(int c_in, int L_in, gsl_rng* rng_in) {
    L = L_in;
    genome_fit = 0;
    g.push_back(block(0,L_in-1,c_in));  //entire array of neutral homogeneous genes
    //seed = seed_in ? seed_in : get_random_seed();
    rng = rng_in;
    //gsl_rng_set(rng,seed);
}

genome::genome(int c_in, int L_in, gsl_rng* rng_in, map<int,double> &fit_in) {
    L = L_in;
    g.push_back(block(0,L_in-1,c_in,fit_in));  //entire array of homogeneous genes
    genome_fit = g[0].block_fit;
    //seed = seed_in ? seed_in : get_random_seed();
    rng = rng_in;
    //gsl_rng_set(rng,seed);
}

genome::genome(gsl_rng* rng_in, vector < block > g_in) {
    g = g_in;
    L = 0;
    genome_fit = 0;
    for (int i=0; i<g.size(); i++) {
        try {
            if (i !=0 && g[i-1].r != (g[i].l-1)) {throw "Positions of blocks incorrect!";}
            L += g[i].get_width();
            genome_fit += g[i].block_fit;
        }
        catch (const char* Message) {cout<<"Error:"<<Message<<"\n";}
    }
    //seed = seed_in ? seed_in : get_random_seed();
    rng = rng_in;
    //gsl_rng_set(rng,seed);

}
genome::~genome(){
    if (VERBOSE)
        cout << "Destructing genome...\n"; 
}

int genome::find_block_index(int locus, int min, int max) {
    //Implements a binary search on the genome to find the address of the block that contains a given locus.
    try{
        if (locus < 0 || locus >= L) {throw "Invalid Locus address";}
        else if (min > max) {throw "Search for block failed!";}
        int midpoint = min + ((max-min)/2);
        if (locus < g[midpoint].l) {return find_block_index(locus,min,midpoint-1);}
        else if (locus > g[midpoint].r) {return find_block_index(locus,midpoint+1,max);}
        else {return midpoint;}
    }
    catch (const char* Message) {cout <<"Error:"<<Message<<"\n";}
}

bool genome::containsLocus(int ancestor, int locus) {
    try {
        if (locus < 0 || locus >= L) {throw "Invalid Locus address";}
        return (ancestor == g[find_block_index(locus,0,g.size())].dat);
    }
    catch (const char* Message) {
        cout << "Error:" << Message << "\n";
    }    
}

double genome::get_fit() {return genome_fit;}

vector<double> genome::cumulant_fit() {
    vector<double> cumulant_fit_genome (L);
    double total_fit = 0;
    for (int i=0; i<g.size(); i++) {
        for (map<int,double>::const_iterator it = g[i].fit_loci.begin(); it != g[i].fit_loci.lower_bound(g[i].r); it++) {
            total_fit += it->second;
            cumulant_fit_genome[it->first] = total_fit;
        }
    }
    //Fill in zeros.
    double last_non_zero = 0;
    for (int i=0; i<L; i++) {
        if (cumulant_fit_genome[i] == 0) {cumulant_fit_genome[i] = last_non_zero;}
        else {last_non_zero = cumulant_fit_genome[i];}
    }
    return cumulant_fit_genome;
}

double genome::sub_cumulant_fit(int start_locus, int l) {
    //Find block locations on the genome.
    int start_block_address = find_block_index(start_locus,0,L-1);
    int last_block_address = find_block_index(start_locus+l,0,L-1);
    double delta_f = 0;
    for (int i=start_block_address; i<last_block_address; i++) {delta_f += g[i].block_fit;} //Iterate over the intermediate blocks.
    //Iterate over the relevant subsection of the last block.
    for (map<int,double>::const_iterator it = g[last_block_address].fit_loci.begin(); it != g[last_block_address].fit_loci.lower_bound(l); it++)    {
        delta_f += it->second;   
    }
    return delta_f;
}

genome genome::operator+(genome g2){
    try {
        if (g2.L != L) {
            throw "Genome sizes are not equal!";
        }
        else {
            vector<block> offspring; 
            //gsl_rng_set(rng,seed);
            int crossPoint = (int)gsl_rng_uniform_int(rng,L-1); //Obtain the random crossover point. 
            //Compute Crossover
            int index1 = find_block_index(crossPoint,0,g.size());
            int index2 = g2.find_block_index(crossPoint,0,g2.g.size());
            for (int i=0; i<index1; i++) {
                offspring.push_back(g[i]);
            }
            if (g[index1].dat != g2.g[index2].dat){
                //If they are different -> create two blocks
                if (g[index1].fit_loci.empty() && g2.g[index2].fit_loci.empty()) { //If neutral, don't worry about updating fitness
                    if (crossPoint >= g[index1].l) {
                        offspring.push_back(block(g[index1].l,crossPoint,g[index1].dat));
                    }
                    if (crossPoint+1 <= g2.g[index2].r) {
                        offspring.push_back(block(crossPoint+1,g2.g[index2].r,g2.g[index2].dat));
                    }
                }
                else { //Partition up fitness hash tables amongst the composite blocks.
                    map<int,double> fit_loci_b1; map<int,double> fit_loci_b2;
                    if (!g[index1].fit_loci.empty()) {
                        fit_loci_b1.insert(g[index1].fit_loci.begin(),g[index1].fit_loci.lower_bound(crossPoint+1));
                    }
                    if (!g2.g[index2].fit_loci.empty()) {
                        fit_loci_b2.insert(g2.g[index2].fit_loci.lower_bound(crossPoint+1),g2.g[index2].fit_loci.end());
                    }
                    //Create two blocks - conditional on the fact that the crossover point didn't fall between blocks. 
                    //If so, then we don't need to create one or the other or both!
                    if (crossPoint >= g[index1].l) {
                        if (fit_loci_b1.empty()) {offspring.push_back(block(g[index1].l,crossPoint,g[index1].dat));}
                        else {offspring.push_back(block(g[index1].l,crossPoint,g[index1].dat, fit_loci_b1));}
                    }
                    if (crossPoint+1 <= g2.g[index2].r) {
                        if (fit_loci_b2.empty()) {offspring.push_back(block(crossPoint+1,g2.g[index2].r,g2.g[index2].dat));}
                        else {offspring.push_back(block(crossPoint+1,g2.g[index2].r,g2.g[index2].dat, fit_loci_b2));}
                    }
                }
            }
            else{
                //Else create one (merge event).
                block merge_block;
                if (g[index1].fit_loci.empty() && g2.g[index2].fit_loci.empty()) { //If neutral, don't worry about updating fitness
                    merge_block = block(g[index1].l,g2.g[index2].r,g[index1].dat);
                }
                else { //Stitch together the fitness hash tables for each block to create one composite hash table.
                    map<int,double> fit_loci;
                    if (!g[index1].fit_loci.empty()) {
                        fit_loci.insert(g[index1].fit_loci.begin(),g[index1].fit_loci.lower_bound(crossPoint+1));
                    }
                    if (!g2.g[index2].fit_loci.empty()) {
                        fit_loci.insert(g2.g[index2].fit_loci.lower_bound(crossPoint+1),g2.g[index2].fit_loci.end());
                    }
                    if (fit_loci.empty()) {merge_block = block(g[index1].l,g2.g[index2].r,g[index1].dat);}
                    else {merge_block = block(g[index1].l,g2.g[index2].r,g[index1].dat,fit_loci);}
                }
                offspring.push_back(merge_block);
            }
            for (int j=(index2+1); j<g2.g.size();j++) {
                offspring.push_back(g2.g[j]);
            }
            return genome(rng,offspring);
            /*
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
            */ 
        }
    }
    catch(const char* Message) {
        cout << "Error:" << Message << "\n";
    }
}

ostream& operator<< (ostream &out, const genome &gen) {
    for(int i=0; i<gen.g.size(); i++) {
        for(int j=0; j<(gen.g[i].r-gen.g[i].l+1); j++) {
            out << gen.g[i].dat;
        }
    }
    out << "\n";    
    return out;
}

/*
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
*/

//POPULATION METHODS

//Overloaded Constructors
population::population(int N_in, int L_in, double r_in, int seed_in, gsl_rng* rng_in) { //Neutral Evolution!
    try {
        if (r_in < 0 || r_in > 1 || N_in < 0 || L_in <0) {throw "Bad inputs. Check arguments.";}
        else {
            N = N_in;
            L = L_in;
            r = r_in;
            gen_pop = 0;
            neutral = true;
            avgFit = 0;
            selective_weight = vector<double>(N); //All zeros
            genetic_weight = vector< vector<int> >(N,vector<int>(L));
            fixations = vector< vector<fixation_event> > (N);
            seed = seed_in ? seed_in:get_random_seed();
            rng = rng_in;
            gsl_rng_set(rng,seed);
            blockHist = gsl_histogram_alloc(L);
            gsl_histogram_set_ranges_uniform(blockHist,1,L+1);
            for(int i=0;i<N;i++) {
                pop.push_back(genome(i,L,rng));
            }
            updateBlockSizes();
        }
    }
    catch(const char* Message) {cout << "Error: " << Message << "\n";}
}

population::population(int N_in, int L_in, double r_in, int seed_in, gsl_rng* rng_in, vector< vector<double> >& fit) {
    try {
        if (fit.size() != N_in || fit[0].size() != L_in) {throw "Incorrect Fitness Landscape dimensions!";}
        else if (r_in < 0 || r_in > 1 || N_in < 0 || L_in <0) {throw "Bad inputs. Check arguments.";}
        else {
            N = N_in;
            L = L_in;
            r = r_in;
            gen_pop = 0;
            neutral = false;
            avgFit = 0;
            selective_weight = vector<double>(N);
            genetic_weight = vector< vector<int> >(N,vector<int>(L));
            fixations = vector< vector<fixation_event> > (N);
            seed = seed_in ? seed_in : get_random_seed();
            rng = rng_in;
            gsl_rng_set(rng,seed);
            blockHist = gsl_histogram_alloc(L);
            gsl_histogram_set_ranges_uniform(blockHist,1,L+1); //Questionable range choice. As it stands now, the bins are partitioned as follows {[1,2),[2,3)......[L,L+1)}. Thus if all have length L, then average will come out to be (L+.5). Will have to rescale.
            map<int,double> genome_fit;
            for(int i=0;i<N;i++) {
                genome_fit = convert_vector_toMap(fit[i]);
                if (!genome_fit.empty()) {pop.push_back(genome(i, L, rng, genome_fit));}
                else {pop.push_back(genome(i,L,rng));}
                avgFit += pop[i].get_fit();
            }
            avgFit/=N;
            update_selection_weights();
            updateBlockSizes();
        }
    }
    catch (const char* Message) {
        cout << "Error:" << Message << "\n";
    }
}

population::~population(){
    if (VERBOSE)
        cout << "Destructing Population...\n"; 
    gsl_histogram_free(blockHist); //Remove usage of histogram!
}

//Functions
map<int,double> population::convert_vector_toMap(vector<double> data) {
    map<int,double> out;
    for (int i=0; i<data.size(); i++) { 
        if (data[i] != 0) {out.insert(pair<int,double> (i,data[i]));}
    }
    return out;
}

void population::update_selection_weights() { 
    double sum_weight = 0;
    for(int i=0;i<N;i++) {
        sum_weight += exp(pop[i].get_fit()-avgFit); //Increment our temp variable by the weight of the ith ind.
        selective_weight[i] = sum_weight; // The ith bin should include the sum from 0 to i.
    }
}

void population::update_weights(){
    vector< vector<int> > all_weight (N, vector<int>(L)); //Initialize the necessary matrix.
    for (int i=0; i<N; i++){ //For all members of current population
        for (int j=0; j<pop[i].g.size(); j++) { //Iterate over the ith genome
            for (int k=pop[i].g[j].l; k<=pop[i].g[j].r; k++) { //Increment the necessary tiles of the matrix.
                all_weight[pop[i].g[j].dat][k]++;
                if (all_weight[pop[i].g[j].dat][k] == N && genetic_weight[pop[i].g[j].dat][k] != N) { //Keep track of fixations.
                    fixations[pop[i].g[j].dat].push_back(fixation_event(gen_pop,k));
                } 
            }
        }
    }
    genetic_weight.swap(all_weight); //swap contents with the temporary vector. 
}

void population::writeBlockHist(const char* s){
    FILE *stream = fopen(s,"w");
    gsl_histogram_fprintf(stream,blockHist,"%g","%g");
    fclose(stream);
}

void population::write_genetic_weight(const char* s) {
    ofstream out_weight(s, ios::out);
    for (vector<vector<int> >::const_iterator it_row = genetic_weight.begin(); it_row != genetic_weight.end(); it_row++) {
        vector<int> row = *it_row;
        for (vector<int>::const_iterator it_col = row.begin(); it_col != row.end(); it_col++) {
            out_weight << *it_col << " ";  
        }
        out_weight << "\n";
    }
    out_weight.close();
}
void population::write_fixations(const char* s) {
    ofstream out_fix(s, ios::out);
    for (int i=0; i<fixations.size(); i++) {
        if (!fixations[i].empty()) {
            for (int j=0; j<fixations[i].size(); j++) {
                out_fix << fixations[i][j].t_fix << " " << fixations[i][j].loci << " " << i;
                out_fix << "\n";
            }
        }
    }
    out_fix.close();
}

genome population::progeny(pair<int,int> parents) { //Creates an offspring genotype
    genome parent1 = pop[parents.first];
    genome parent2 = pop[parents.second];
    return parent1+parent2;
}

void population::updateBlockSizes(){
    gsl_histogram_reset(blockHist);
    for (int i=0;i<N;i++){
       for (int j=0; j<pop[i].g.size(); j++) {    
            gsl_histogram_increment(blockHist,pop[i].g[j].get_width());
        }    
    }
}

vector<int> population::get_weight(int ancestor){
    try {
        if (ancestor >= N || ancestor < 0) {throw "Incorrect Ancestral Value!";}
        else {return genetic_weight[ancestor];}
    }
    catch (const char* Message) {
        cout << "Error:" << Message << "\n";
    }
}

vector<vector<int> > population::get_all_weight() {
    return genetic_weight;
}

vector< vector <fixation_event> > population::get_fixations() {
    return fixations;
}

vector< pair<fixation_event, int> > population::get_flat_fixations() {
    vector< pair<fixation_event, int> > flat_fixations;
    for (int i = 0; i < fixations.size(); i++) {
        if (!fixations[i].empty()) {
            for (int j =0; j <fixations[i].size(); j++) {
                flat_fixations.push_back(pair<fixation_event,int>(fixations[i][j],i));
            }
        }
    }
    return flat_fixations;
}

void population::evolve(int gen) { //Wright-Fisher Model - Population size is held constant. Parental genotypes are chosen randomly. 
    for (int j=0; j<gen; j++){
        vector<genome> newPop (N);
        avgFit = 0;
        for(int i=0; i<N; i++) {
            if (r == 1) {   
                newPop[i] = progeny(selectParents());
                avgFit += newPop[i].get_fit();
            }
            else if ((r != 0) && gsl_rng_uniform(rng) <= r) {
                newPop[i] = progeny(selectParents());
                avgFit += newPop[i].get_fit();
            }
            else {
                if (neutral) {newPop[i] = pop[(int)gsl_rng_uniform_int(rng,N)];}
                else {
                    double fit_rand = selective_weight[N-1]*gsl_rng_uniform(rng);
                    newPop[i] = pop[binarySearch_int(selective_weight, 0, N-1, fit_rand)];
                    avgFit += newPop[i].get_fit();
                }            
            }
        }
        avgFit/= N;
        pop = newPop;
        gen_pop++;
        update_weights();
        //updateBlockSizes();
        if (!neutral) {update_selection_weights();}
    }
    update_weights();
    updateBlockSizes();
    //if (WRITEHIST) {writeBlockHist();}
}

void population::evolve() { //Wright-Fisher Model - Population size is held constant. Parental genotypes are chosen randomly. 
   evolve(1); 
}

double population::get_avg_fit() {
    return avgFit;
}

double population::get_variance() {
    double var = 0;
    for (int i=0; i<N; i++) {
        var += pow(pop[i].get_fit() - avgFit,2);
    }
    return var/=(N-1);
}

double population::get_allele_frequency(int ancestor, int locus) {
    double count = genetic_weight[ancestor][locus];
    return count/N;
}

pair<int,int> population::selectParents(){ //Parents are selected at random from population 
    pair<int,int> ancestors;
    if (neutral) { //Neutral Evolution
        ancestors.first = (int)gsl_rng_uniform_int(rng,N);
        ancestors.second = (int)gsl_rng_uniform_int(rng,N);
    }
    else { //Selection!
        double fit_rand1 = selective_weight[N-1]*gsl_rng_uniform(rng);
        double fit_rand2 = selective_weight[N-1]*gsl_rng_uniform(rng);
        ancestors.first = binarySearch_int(selective_weight, 0, N-1, fit_rand1);
        ancestors.second = binarySearch_int(selective_weight, 0, N-1, fit_rand2);
    }
    //Recursive definition to ensure unique parents.
    if (ancestors.second == ancestors.first) {
        return selectParents();
    }
    else {
        return ancestors;
    }
}

int population::binarySearch_int(vector<double> &data, int min, int max, double &key) { //Data is assumed to be sorted in ascending order!
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

gsl_histogram* population::get_cumulant_fit_hist(int l, double sigma) {
    //Initialize histogram.
    gsl_histogram* fit_hist = gsl_histogram_alloc(1000); //Think more about bin size.
    gsl_histogram_set_ranges_uniform(fit_hist, -10*l*sigma, 10*l*sigma);
    for (int n=0; n<N; n++) { //Iterate over the entire population.
        for (int i=0; i<L-l; i+l/2) { //Center Locus @ i+l/2
            gsl_histogram_increment(fit_hist,pop[n].sub_cumulant_fit(i,l)); //Histogram delta_f values. 
        }   
    }
    return fit_hist;
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

void dump_map(const std::map<int, double>& map) {
    for (std::map<int,double>::const_iterator it = map.begin(); it != map.end(); it++) {
        cout << (it-> first) << "," << (it-> second) << endl;
    }
}

int get_random_seed() {
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

vector< vector<double> > get_rand_fitScape(int N, int L, double mu, double s, gsl_rng* rng) {
    //Decide how many non-dark matter alleles are sorting through the population. 
    unsigned int polymx = gsl_ran_poisson(rng,N*L*mu); 
    if (polymx > N*L) {polymx =N*L;} 
    cout << polymx << "\n";
    //Decide placement of founding polymorphisms.
    pair<int,int>* polymx_loc = new pair<int,int>[polymx];
    pair<int,int>* genome_loc = new pair<int,int>[N*L];
    for (int i=0; i<N; i++) {
        for (int j=0; j<L; j++) {
            genome_loc[i*L+j]= pair<int,int> (i,j);
        }
    } 
    gsl_ran_choose(rng, polymx_loc, polymx, genome_loc, N*L, sizeof(pair<int,int>));
    delete [] genome_loc;
    //Decide on selective strength
    vector<double> polymx_s (polymx);
    for (int i=0; i<polymx; i++) {polymx_s[i] = gsl_ran_gaussian(rng,s);}
    //Build fitness landscape
    vector< vector<double> > fitScape_pop (N);
    for (int i=0; i<N; i++) {fitScape_pop[i] = vector<double> (L);}
    for (int i=0; i<polymx; i++) {fitScape_pop[polymx_loc[i].first][polymx_loc[i].second] = polymx_s[i];}
    delete [] polymx_loc;
    return fitScape_pop;
}

vector< vector<double> > seed_fitScape(vector<vector<double> > fitScape, int num_seeds, double s, int width) {
    for (int i=0; i<num_seeds; i++) {
        for (int j=0; j<width; j++) {
            fitScape[i][j] = s/width;
        }
    }
    return fitScape;
}

void find_fixation_prob(int N, int L, double r, double s, gsl_rng* rng, int runs, vector<vector <double> > fit_scape, ofstream& out) {
    fit_scape = seed_fitScape(fit_scape, 5, s, 1); //seed 5 SNP's
    for (int i=0; i<5; i++) { //Output fitness background.
            for (int j=0; j<L; j++) {
                out << fit_scape[i][j] << " ";
            }
            out << "\n";
    }
    for (int n=0; n<runs; n++) {
        population pop(N, L, r, 0, rng, fit_scape);
        pop.evolve(N); //Evolve on neutral time scale.
        for (int m=0; m<5; m++) {
            out << pop.get_allele_frequency(m, 0) << " ";
        }
        out << "\n";
    }
}

void update_fit_cumulant_histogram(gsl_histogram* h, int l, vector<pair<fixation_event,int> > fixations, vector<vector<double> > fitness_landscape) {
    for (vector<pair<fixation_event, int> >::const_iterator it = fixations.begin(); it != fixations.end(); it++) {
        //Must consider boundary conditions.
        int L = fitness_landscape[0].size();
        int min; int max;
        if (it->first.loci < (l/2)) {
            min = 0;
            max = l-1;
        }
        else if (it->first.loci > (L-1)-l/2) {
            max = L-1;
            min = L-1-l+1; 
        }
        else {
            min = it->first.loci - l/2 + (1-l%2); //Last term used to correct for parity.
            max = it->first.loci + l/2;
        }
        double delta_fit = 0;
        for (int i=min; i<=max; i++) { //Sum over relevant loci.
            delta_fit += fitness_landscape[it->second][i];
        }
        gsl_histogram_increment(h,delta_fit); //Increment relevant bin.
    }
}

void update_fit_cumulant_histogram(gsl_histogram* h, int l, vector<vector<double> >& fitness_landscape) {
    //If no fixation vector is passed, assume you want unconditioned histogram.
    double delta_fit = 0; double last_delta_fit = 0;
    int L = fitness_landscape[0].size();
    for (int n=0; n<fitness_landscape.size(); n++) {
        for (int i=0; i<(L-l); i++) { //Iterate over entire genome.
            if (i == 0) {
                delta_fit = 0;
                for (int j=0; j<l; j++) {
                    delta_fit += fitness_landscape[n][j];
                }
                last_delta_fit = delta_fit;
            }
            else {
                delta_fit = delta_fit - fitness_landscape[n][i-1] + fitness_landscape[n][i+l];
                last_delta_fit = delta_fit;
            }
            gsl_histogram_increment(h,delta_fit);
        }
    }
}

void update_fit_cumulant_histogram_2(gsl_histogram* h, int l, vector<vector<double> >& fitness_landscape) {
    int L = fitness_landscape[0].size();
    for (int n=0; n<fitness_landscape.size(); n++) {
        for (int i=0; i<(L-l); i+=((l/2)+(l%2))) { //Iterate over entire genome.
            double delta_fit = 0;
            for (int j=i; j<=i+l; j++) {
                delta_fit += fitness_landscape[n][j];
            }
            gsl_histogram_increment(h,delta_fit);
        }
    }
}
            //double delta_fit = 0;
            /*
            int min; int max;
            if (i < (l/2)) {
                min = 0;
                max = l-1;
            }
            else if (i > (L-1)-l/2) {
                max = L-1;
                min = L-l; 
            }
            else {
                min = i - l/2 + (1-l%2); //Last term used to correct for parity.
                max = i + l/2;
            }
            */
            //for (int j=min; j<=max; j++) {
            //    delta_fit += fitness_landscape[n][j];
            //}
vector< pair<int,double> > calc_correlation(vector<vector <fixation_event> > fix, vector<vector <double> > fit_land) {
    vector < pair<int,double> > pt_correlation;
    for (int i=0; i<fix.size(); i++) {
        if (!fix[i].empty()) {
            for (int j=0; j<fix[i].size(); j++) {
                for (int k=j; k<fix[i].size(); k++) {
                    int min = (fix[i][j].loci>fix[i][k].loci)?fix[i][k].loci:fix[i][j].loci;
                    int max = (fix[i][j].loci<fix[i][k].loci)?fix[i][k].loci:fix[i][j].loci;
                    double delta_f = 0;
                    for (int l=min; l<=max; l++) {delta_f += fit_land[i][l];}
                    pt_correlation.push_back(pair<int,double>(max-min+1,delta_f));
                }
            }
        }
    }
    return pt_correlation;
}

void write_histogram(gsl_histogram* h, const char* s){
    FILE *stream = fopen(s,"w");
    gsl_histogram_fprintf(stream,h,"%g","%g");
    fclose(stream);
}

int main() {
    /*
    for (int i=0; i<g4.g.size(); i++) {
        dump_map(g4.g[i].fit_loci);
    }
    cout << g4.get_fit() << "\n";
    
    genome g6 = g5 + g3;
    cout << g6 << "Size:" << g6.g.size() << "\n";
    for (int i=0; i<g6.g.size(); i++) {
        dump_map(g6.g[i].fit_loci);
    }
     
    int x = g0.find_block_index(20,0,g0.g.size());
    cout << g0.g[x].l <<"," << g0.g[x].r<< "\n";
    
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
    */ 
    /*
    double sigma = .001;
    vector< vector<double> > fit_landscape = get_rand_fitScape(100,1000,1,sigma, rng_in);
        //Calculte length dependent histogram
        for (int l=1; l<=101; l+=1) {
            gsl_histogram* hist = gsl_histogram_alloc(5000);
            gsl_histogram_set_ranges_uniform(hist, -l*1.5*sigma, l*1.5*sigma);
            update_fit_cumulant_histogram(hist, l, fit_landscape);
            s << "hist_full(" << l << "-l).txt";
            write_histogram(hist, s.str().c_str());
            s.str(std::string());

            gsl_histogram* hist2 = gsl_histogram_alloc(5000);
            gsl_histogram_set_ranges_uniform(hist2, -l*1.5*sigma, l*1.5*sigma);
            update_fit_cumulant_histogram_2(hist2, l, fit_landscape);
            s << "hist_block(" << l << "-l).txt";
            write_histogram(hist2, s.str().c_str());
            s.str(std::string());

            gsl_histogram_free(hist);
            gsl_histogram_free(hist2);
        }
        */
    gsl_rng* rng_in = gsl_rng_alloc (RNG);
    std::stringstream s;
    //double r = 1;
    //int max_window = 500;
    for (double r = 0; r <= 0; r += .01) {
    for (double sigma = .003; sigma <= .003; sigma += .004) {
        vector< vector<double> > fit_landscape = get_rand_fitScape(1000,10000,.05,sigma, rng_in);
        /*
        gsl_histogram** hist_l = new gsl_histogram*[max_window];   
        for (int l=1; l<=max_window; l+=1) {
            hist_l[l-1] = gsl_histogram_alloc(10000);
            gsl_histogram_set_ranges_uniform(hist_l[l-1],-1.5*sigma*l,1.5*sigma*l);
        }
        */
        /*
        s.str(std::string());
        s << "coalescent_data(1000-N,10000-L," << std::setprecision(6) << sigma << "-z," << r << "-r).txt";
        ofstream out_l(s.str().c_str(), ios::out);
        */
        s << "selection_test(1000-N,10000-L," << std::setprecision(6) << sigma << "-z," << r << "-r).txt";
        ofstream out_s(s.str().c_str(), ios::out);
        for (int iter = 0; iter <= 15; iter += 1) {
            //cout << iter << "\n";
            //out_s << iter << "\n";
            //Initialize Population
            population pop(1000, 10000, r, 0, rng_in, fit_landscape);
            int T_max = 0;
            if (sigma == .003) {T_max = 300;}
            else if (sigma == .005) {T_max = 250;}
            else {T_max = 100;}
            for (int t=0; t < T_max; t += 1) {
                out_s << t << " " << pop.get_avg_fit() << " " << pop.get_variance() << "\n";
                pop.evolve(); 
            }
            /*
            if (sigma == .01) {
                pop.evolve(650);
            }
            else {
                pop.evolve(450);
            }
            
            //Calculate length distribution of fixed blocks. 
            for (int i = 0; i < pop.get_fixations().size(); i++) {
                if (!pop.get_fixations()[i].empty()) {
                    vector<fixation_event> ancestor = pop.get_fixations()[i];
                    int width = 1; 
                    int start_locus = ancestor[0].loci; int last_locus = start_locus;
                    double fix_time_avg = ancestor[0].t_fix;
                    double local_fit_avg = fit_landscape[i][start_locus];
                    for (int j = 1; j < ancestor.size(); j++) {
                        if (ancestor[j].loci == last_locus + 1) {
                            //Encountered piece of continuous segment. Increment everything.
                            width++;
                            last_locus++;
                            fix_time_avg += ancestor[j].t_fix;
                            local_fit_avg += fit_landscape[i][start_locus];
                        }
                        else {
                            //End of contiguous segment found.  
                            local_fit_avg /= width; fix_time_avg /= width;
                            double local_fit_std = 0;
                            double fix_time_std = 0;
                            if (width > 1) {
                                //Calculate local fitness variance
                                for (int locus = start_locus; locus <= last_locus; locus++) {
                                    local_fit_std += pow(fit_landscape[i][locus] - local_fit_avg,2);
                                }
                                local_fit_std /= (width-1);
                                local_fit_std = sqrt(local_fit_std);
                                //Calculate fixation time variance
                                for (int k = j-width; k < j; k++) {
                                    fix_time_std += pow(ancestor[k].t_fix - fix_time_avg,2);
                                }
                                fix_time_std /= (width-1);
                                fix_time_std = sqrt(fix_time_std);
                            }
                            //Print data to file.
                            out_l << width << " " << local_fit_avg << " " << local_fit_std << " " << fix_time_avg << " " << fix_time_std << "\n";
                            //Reset vanguard.
                            width = 1;
                            start_locus = ancestor[j].loci; last_locus = start_locus;
                            fix_time_avg = ancestor[j].t_fix;
                            local_fit_avg = fit_landscape[i][start_locus];
                        }
                    }
                }
            }
            */
            //Evolve the population to obtain fixation events.
            /*
            for (int i = 1; i <= 50; i++){
                pop.evolve(8);
                out_l << gsl_histogram_mean(pop.blockHist) << " " << gsl_histogram_sigma(pop.blockHist) << "\n";
            } 
            //s << "fix_event(1000-N,10000-L," << std::setprecision(6) << sigma << "-z-" << iter << "-iter).txt";
            //pop.write_fixations(s.str().c_str());
            //s.str(std::string());
            
            s << "block_hist(100-N,5000-L," << std::setprecision(6) << sigma << "-z-" << iter << "-iter).txt";
            pop.writeBlockHist(s.str().c_str());
            */
           
            //Calculate two point differences.
            /*
            s.str(std::string()); 
            vector< pair<int,double> > correlation =  calc_correlation(pop.get_fixations(), fit_landscape);               
            s << "fix_pairs(100-N,5000-L," << std::setprecision(6) << sigma << "-z-" << iter << "-iter).txt";
            ofstream out_fit(s.str().c_str(), ios::out);
            for (vector< pair<int,double> >::const_iterator it = correlation.begin(); it != correlation.end(); it++) {
                out_fit << it->first << " " << it->second;
                out_fit << "\n";
            }
            out_fit.close();
            s.str(std::string());
            
            //Calculte length dependent histogram
            for (int l=1; l<=max_window; l+=1) {
                gsl_histogram* hist_fix_temp = gsl_histogram_alloc(10000);
                gsl_histogram_set_ranges_uniform(hist_fix_temp, -1.5*sigma*l, 1.5*sigma*l);
                update_fit_cumulant_histogram(hist_fix_temp, l, pop.get_flat_fixations(), fit_landscape);
                gsl_histogram_add(hist_l[l-1],hist_fix_temp);
                gsl_histogram_free(hist_fix_temp);
            }
            
            index = 0;
            s << "dist_diff(100-N,5000-L," << std::setprecision(6) << sigma << "-z-" << iter << "-iter).txt";
            ofstream out_diff(s.str().c_str(), ios::out);
            for (vector <double>::const_iterator it = diff_dist.begin(); it != diff_dist.end(); it++) {
                out_diff << *it << " " << 2*index+1 << " ";
                out_diff << "\n";
                index++;
            }
            out_diff.close();
            s.str(std::string());
            */
        }
        out_s.close();
        s.str(std::string());
    }
    }
        /*
        s << "dist_diff(500-N,10000-L," << std::setprecision(6) << sigma << "-z).txt";
        ofstream out_diff(s.str().c_str(), ios::out);
        for (int l=1; l<=max_window; l++) {
            gsl_histogram* hist_fit = gsl_histogram_alloc(10000);
            gsl_histogram_set_ranges_uniform(hist_fit, -1.5*sigma*l, 1.5*sigma*l);
            update_fit_cumulant_histogram(hist_fit, l, fit_landscape);
            gsl_histogram_scale(hist_fit,1/gsl_histogram_sum(hist_fit));

            gsl_histogram_scale(hist_l[l-1],1/gsl_histogram_sum(hist_l[l-1]));
            vector <double> fit_cumulant (10000);
            vector <double> fix_cumulant (10000);
            for (int i=0; i<10000; i++) {
                if (i==0) {
                    fit_cumulant[i] = 0;
                    fix_cumulant[i] = 0;
                }
                else {
                    fit_cumulant[i] = fit_cumulant[i-1] + gsl_histogram_get(hist_fit,i);
                    fix_cumulant[i] = fix_cumulant[i-1] + gsl_histogram_get(hist_l[l-1],i);
                }
            }
            double diff = 0;
            for (int i=0; i<10000; i++) {
                diff += fix_cumulant[i] - fit_cumulant[i];
            }
            diff /= 10000;
            out_diff << l << " " << diff << "\n";
        }
        out_diff.close();
        s.str(std::string());
        delete [] hist_l;
    }
    */
        //Calculate the Histogram of Various Lengths on the fixed landscape. Write to text files.
        /*    
        for (int l=1; l<=201; l+=2) {
            gsl_histogram* hist_cumulant_fit_fixed = gsl_histogram_alloc(5000);
            gsl_histogram_set_ranges_uniform(hist_cumulant_fit_fixed, -l*.6, l*.6);
            update_fit_cumulant_histogram(hist_cumulant_fit_fixed, l, pop.get_fixations(), fit_landscape);
            s << "hist_fix(N1000,L20000," << l << "-l," << std::setprecision(6)<< r <<"-r).txt";
            write_histogram(hist_cumulant_fit_fixed, s.str().c_str());
            s.str(std::string());
            gsl_histogram_free(hist_cumulant_fit_fixed);
        }
        
        }
        */
    //FILE *stream = fopen("out.txt","w");
    //FILE *stream2 = fopen("out2.txt", "w");
    //gsl_histogram_fprintf(stream,hist_cumulant_fit,"%g","%g");
    //gsl_histogram_fprintf(stream2,hist_cumulant_fit_fixed,"%g","%g");
    //fclose(stream2);
    //fclose(stream);
    /*
    for (double r=0; r<=.5; r+=.05) {
        s << "fix_prob(1000N256Ld1s" << std::setprecision(6) << r << "r).txt";
        ofstream output(s.str().c_str(), ios::out);
        for (int n=0; n<20; n++) {
            vector< vector<double> > fit_landscape = get_rand_fitScape(1000,256,.01,.1, rng_in); //Initialize fitscape.
            find_fixation_prob(1000, 256, r, .1, rng_in, 100, fit_landscape, output); 
        }
        s.str(std::string());
        output.close();
    }
    */
    //vector< vector<double> > fit_landscape = get_rand_fitScape(1000,256,.01,.1, rng_in);
    //std::stringstream s; 
    /*
    for (int t=0; t<50; t++) {
        for (double r=0; r<=.5; r+=.005) {
            population pop(1000,256,r, get_random_seed(), rng_in, fit_landscape);
            pop.evolve(1000);
            s << "fix(" << std::setprecision(6) << r <<  ")L256-" << t << ".txt";
            pop.write_fixations(s.str().c_str());
            s.str(std::string());
        }
    }
    */
    /*
    population pop(1000,256,0,get_random_seed(), rng_in, fit_landscape);
    pop.evolve(10);
    pop.write_fixations("test.txt");
    population pop2(1000,256,0,get_random_seed(), rng_in, fit_landscape);
    pop.evolve(10);
    pop.write_fixations("test2.txt");
    */
    //for (int i=0; i<25; i++) {
        //cout << "Average Fitness: " << pop.get_avg_fit() << ", Variance: " << pop.get_variance() << "\n";;
        //pop.evolve(40);
    //}
    //pop.write_genetic_weight("out_weight.txt");
    //pop.write_fixations("fix.txt");
    //cout << gsl_histogram_mean(pop.blockHist) << "\n";
    /* 
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
    
    vector< vector<double> > test (1000);
    for (int i=0; i<1000; i++) {
        vector<double> temp (1000);
        for (int j=0; j<1000; j++) {
            temp[j] = i * .000001;
        }
        test[i] = temp;
    }
    
    population pop(1000,1000);
    pop.evolve(1000);
    
    for (int t=0; t<10; t++) {
        cout << pop.getAverageFit() << "\n";
        cout << pop.getVariance() << "\n";
        pop.evolve();
    }
    cout << pop.getAverageFit() << "\n";
    cout << pop.getVariance() << "\n";
    */
    return 0;
}

