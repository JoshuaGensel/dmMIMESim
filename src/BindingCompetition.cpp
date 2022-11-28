//
// Created by Smith, Maureen on 04.06.18.
//

#include "BindingCompetition.hpp"
#include "Generator.hpp"
#include "cppoptlib/solver/lbfgsbsolver.h"

std::valarray<kd_type> UnboundProtein::getSpeciesKds(const species::species_map& spec)
{
    std::valarray<kd_type> kds(spec.size());
    unsigned i = 0;
    for (auto& s : spec)
    {
        kds[i] = (s.second).getKd();
        ++i;
    }
    return kds;
}

std::valarray<frequency_type> UnboundProtein::getSpeciesFrequencies(const species::species_map& spec)
{
    std::valarray<frequency_type> freq(spec.size());
    unsigned i = 0;
    for (auto& s : spec)
    {
        freq[i] = s.second.getFreq();
        ++i;
    }
    return freq;
}

std::valarray<count_type> UnboundProtein::getSpeciesCounts(const species::species_map& spec)
{
    std::valarray<count_type> count(spec.size());
    unsigned i = 0;
    for (auto& s : spec)
    {
        count[i] = s.second.getCount();
        ++i;
    }
    return count;
}

/*
 * Draw number of sequence from Binomial distribution with
 *  @param N the total number of a particular sequence
 *  @param p the probability being bound
 */
count_type drawBinomialNumber(const unsigned int N, const double p)
{
    std::default_random_engine& generator = Generator::get_instance()->engine;
    std::binomial_distribution<int> bino(N, p);
    return bino(generator);
}

double UnboundProtein::solve(std::valarray<count_type>& S_bound, std::valarray<count_type>& S_unbound)
{

    // choose a starting point (the amount of unbound protein, in the beginning = total amount of protein)
    TVector B(1);
    B << B_TOT;
    // set boundaries for the amount of free protein (either all free or all bound)
    TVector lo(1);
    lo << 0;
    TVector up(1);
    up << B_TOT;
    this->setLowerBound(lo);
    this->setUpperBound(up);

    // choose a solver
    cppoptlib::LbfgsbSolver<UnboundProtein> solver;
    // and minimize the function
    solver.minimize(*this, B);

    std::valarray<frequency_type> f_bound = 1.0 / (1.0 + (kds / B[0]));

    // sample for each simulated frequency the number of actual bound counts from the total count of that species
    // TODO catch Error: if S_bound is not empty
    std::transform(std::begin(f_bound), std::end(f_bound), std::begin(counts), std::begin(S_bound),
                   [](const auto p, const auto s) {
                       auto nSeq = drawBinomialNumber(s, p);
                       return nSeq;
                   });

    // The remaining counts of per species are the unbound
    S_unbound = counts - S_bound;

    return B[0];
}

const count_type UnboundProtein::getB_tot()
{
    return B_TOT;
}

const std::valarray<kd_type>& UnboundProtein::getKds() const
{
    return kds;
}

const std::valarray<frequency_type>& UnboundProtein::getFrequencies() const
{
    return frequencies;
}
