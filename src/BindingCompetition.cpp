//
// Created by Smith, Maureen on 04.06.18.
//

#include "BindingCompetition.hpp"
#include "../extern/CppNumericalSolvers/cppoptlib/solver/lbfgsbsolver.h"
#include <chrono>
#include <random>

std::valarray<kd_type> UnboundProtein::getSpeciesKds(const species::species_map& spec) {
    std::valarray<kd_type> kds(spec.size());
    unsigned i=0;
    for(auto& s:spec) {
    //for(int i=0; i < spec.size(); ++i) {
        kds[i] = (s.second).getKd();
        ++i;
    }
    return kds;
}

std::valarray<frequency_type> UnboundProtein::getSpeciesFrequencies(const species::species_map& spec){
    std::valarray<frequency_type> freq(spec.size());
    unsigned i=0;
    for(auto& s:spec) {
    //for(int i=0; i < spec.size(); ++i) {
        freq[i] = s.second.getFreq();
        ++i;
    }
    return freq;
}

//TODO: Umbau nach counts
std::valarray<count_type> UnboundProtein::getSpeciesCounts(const species::species_map &spec) {
    std::valarray<count_type> count(spec.size());
    unsigned i=0;
    for(auto& s:spec) {
        //for(int i=0; i < spec.size(); ++i) {
        count[i] = s.second.getCount();
        ++i;
    }
    return count;
}

//TODO Beschreibe/Umbenennen: drawing number of sequences of a binomial distribution
count_type drawNumberOfSequences(const unsigned int N, const double p) {
    const auto seed = static_cast<count_type>(std::chrono::system_clock::now().time_since_epoch().count());
    std::default_random_engine generator(seed);
    std::binomial_distribution<int> bino(N, p);
    return bino(generator);
}


//TODO: Umbau nach counts
//double UnboundProtein::solve(std::valarray<double>& S_bound, std::valarray<double>& S_unbound) {
double UnboundProtein::solve(std::valarray<count_type>& S_bound, std::valarray<count_type>& S_unbound) {

    // choose a starting point (the amount of unbound protein, in the beginning = total amount of protein)
    //TODO: Umbau nach counts
    //TVector B(1); B << double(B_TOT);
    TVector B(1); B << B_TOT;
    // set boundaries for the amount of free protein (either all free or all bound)
    TVector lo(1); lo << 0;
    //TODO: Umbau nach counts
    //TVector up(1); up  << double(B_TOT);
    TVector up(1); up << B_TOT;
    this->setLowerBound(lo);
    this->setUpperBound(up);

    // choose a solver
    cppoptlib::LbfgsbSolver<UnboundProtein> solver;
    // and minimize the function
    solver.minimize(*this, B);

    //TODO: Umbau nach counts
    //S_bound = frequencies/(1.0+(kds/B[0]));
    //S_unbound = frequencies - S_bound;

    //TODO: ich will ja nur die anzahl von der species samples, nicht alles f_s*f_s_b
    //std::valarray<double> f_bound = frequencies/(1.0+(kds/B[0]));
    std::valarray<frequency_type> f_bound = 1.0/(1.0+(kds/B[0]));
    //auto f_unbound = frequencies - f_bound;

    //TODO: Sampling statt runden, vielleicht beeser in BindinCompetition
    //first sample the number of unbound sequences of the sequence variant
    auto& constants = constants::Constants::get_instance();

    // sample for each simulated frequency the number of actual bound counts from the total count of that species
//    std::transform(std::begin(f_bound), std::end(f_bound), std::begin(counts), std::begin(S_bound), [m = constants.M](const auto p, const auto s) {
//        //auto blub =  std::min(drawNumberOfSequences(m, p),s);
//        auto blub =  std::floor(m*p)
//                + std::min(drawNumberOfSequences(m, p),s);
//        std::cout << "bound p" << p << " Stot "  << " " << s << " " << blub << " " << (s-blub) << std::endl;
//        return blub;
//    });
    //TODO catch Error: if S_bound is not empty
    std::transform(std::begin(f_bound), std::end(f_bound), std::begin(counts), std::begin(S_bound), [](const auto p, const auto s) {
        auto blub =  drawNumberOfSequences(s, p);
        //auto blub = floor(s*p);
        //sample die rundung anhand des rests
        //blub += drawNumberOfSequences(1, (s*p)-floor(s*p));
        //TODO output entfernen
        std::cout << "bound p" << p << " Stot "  << " " << s << " " << blub << " " << (s-blub) << " " << round(s*p) - blub<< std::endl;
        return blub;
    });

    //The remaining counts of per species are the unbound
    S_unbound = counts - S_bound;


    //std::cout << S_bound.size() << std::endl;
    //for(int i = 0; i < kds.size(); ++i ){
        //S_bound[i] = round((B[0]*counts[i])/(kds[i]+double(B[0])));
        //S_unbound[i] = counts[i] - S_bound[i];
        //std::cout << "freq " << frequencies[i] << " ";
        //std::cout << i << " " << S_bound[i] << " " << S_unbound[i] << std::endl;

    //}
//    std::cout << "B + bound" << B[0] << " " <<  S_bound.sum() << " " <<  B[0] + S_bound.sum() << " " << B_TOT << std::endl;

    return B[0] ;
}

const count_type UnboundProtein::getB_tot() {
    return B_TOT;
}

const std::valarray<kd_type> &UnboundProtein::getKds() const {
    return kds;
}

const std::valarray<frequency_type> &UnboundProtein::getFrequencies() const {
    return frequencies;
}


