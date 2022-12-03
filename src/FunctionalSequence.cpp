//
// Created by Smith, Maureen on 31.05.18.
//

#include "Constants.hpp"
#include "FunctionalSequence.hpp"
#include "Generator.hpp"
#include "Mutation.hpp"
#include "Species.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>

// Initialize static member instance
FunctionalSequence* FunctionalSequence::instance = NULL;

FunctionalSequence* FunctionalSequence::get_instance()
{
    if (instance == NULL)
    {
        std::cerr << "FunctionalSequence not correctly instantiated (get before create or after release)." << std::endl;
        throw std::runtime_error("FunctionalSequence not correctly instantiated (get before create or after release).");
    }
    return instance;
}
// FunctionalSequence(FunctionalSequence const &);

FunctionalSequence* FunctionalSequence::create_instance(const constants::Constants& params)
{
    if (instance != NULL)
    {
        std::cerr << "FunctionalSequence was already instantiated (create after create)." << std::endl;
        throw std::runtime_error("FunctionalSequence was already instantiated (create after create).");
    }
    instance = new FunctionalSequence(params);
    return instance;
}

FunctionalSequence* FunctionalSequence::create_instance(const std::string& inputPath)
{
    if (instance != NULL)
    {
        std::cerr << "FunctionalSequence was already instantiated (create after create)." << std::endl;
        throw std::runtime_error("FunctionalSequence was already instantiated (create after create).");
    }
    instance = new FunctionalSequence(inputPath);
    return instance;
}

void FunctionalSequence::release_instance()
{
    if (instance == NULL)
    {
        std::cerr << "Deletion of not correctly instantiated FunctionalSequence (release after release)." << std::endl;
        throw std::runtime_error("Deletion of not correctly instantiated FunctionalSequence (release after release).");
    }
    delete instance;
    instance = NULL;
}

std::vector<double> FunctionalSequence::drawKdValues()
{
    std::vector<double> kds(this->params.SVal);
    std::default_random_engine& generator = Generator::get_instance()->engine;
    std::bernoulli_distribution bd(this->params.P_EFFECT);
    std::lognormal_distribution<double> lnd(0, 1);
    for (int i = 0; i < kds.size(); ++i)
    {
        // first sample the positions which have an effect (binomial/bernoulli distributed), then the strength of
        // the effect(log normal distributed)
        kds[i] = bd(generator) ? lnd(generator) : this->params.KD_WT;
    }
    return kds;
}

std::vector<double> FunctionalSequence::readKdValues(const std::string& inputPath)
{
    fs::path kdFile(fs::canonical(inputPath) / "single_kds.txt");

    // TODO: make more efficient by using the param file in inputPath to know
    //  the size of the vector beforehand
    std::vector<double> kds;

    std::ifstream infile(kdFile);
    std::string line;
    while (std::getline(infile, line))
    {
        // each line in this file contains one item (1 column file)
        kds.push_back(std::stod(line));
    }
    return kds;
}

bool FunctionalSequence::checkValidReadValues()
{

    if (this->kds.size() != this->params.SVal)
    {
        std::cerr << "Provided single_kds.txt does not match experiment params (sequence length L)." << std::endl;
        return false;
    }
    if (this->epistasis.size() != this->params.PWVal)
    {
        std::cerr << "Provided pairwise_epistasis.txt does not match exp params (sequence length L)." << std::endl;
        return false;
    }

    return true;
}

std::vector<double> FunctionalSequence::drawEpistasis()
{
    if (this->params.EPIMUTEXCL)
        return drawEpistasis_restricted();
    else
        return drawEpistasis_unrestricted();
}

// TODO: remove
std::vector<double> FunctionalSequence::drawEpistasis_unrestricted()
{
    std::vector<double> epistasis(this->params.PWVal);
    std::default_random_engine& generator = Generator::get_instance()->engine;
    std::bernoulli_distribution bd(this->params.P_EPISTASIS);
    std::lognormal_distribution<double> lnd(0, 1);
    for (int i = 0; i < epistasis.size(); ++i)
    {
        // first sample if position pair has epistatic effect (bernoulli) and then the value of the epistasis (log
        // normal distributed)
        epistasis[i] = bd(generator) ? lnd(generator) : this->params.NO_EPISTASIS;
    }
    return epistasis;
}

// mututally exclusive epistasis of positions that individually have an effect
// each position only interacts with one other position
std::vector<double> FunctionalSequence::drawEpistasis_restricted()
{
    std::vector<double> epistasis(this->params.PWVal);
    std::fill(epistasis.begin(), epistasis.end(), this->params.NO_EPISTASIS);

    std::default_random_engine& generator = Generator::get_instance()->engine;
    std::bernoulli_distribution bd(this->params.P_EPISTASIS);
    std::lognormal_distribution<double> lnd(0, 1);
    std::uniform_int_distribution<unsigned int> symb_unif(0, (this->params.Q - 1) - 1); // [low, high]

    std::vector<int> sym_pos_effect;   // symbol at position have effect
    std::vector<int> selected_effects; // single position has effect

    int nSymb = this->params.Q - 1;
    for (int i = 0; i < this->params.L; ++i)
    {
        // select only one mutation per position
        // TODO: allow multiple? If so, add check in while loop for pos1=pos2
        for (int j = i * nSymb; j < i * nSymb + nSymb; j++)
        {
            if (this->kds[j] != this->params.KD_WT)
            {
                sym_pos_effect.push_back(j);
            }
        }
        std::shuffle(sym_pos_effect.begin(), sym_pos_effect.end(), generator);
        selected_effects.push_back(sym_pos_effect[0]);
        sym_pos_effect.clear();
    }
    std::shuffle(selected_effects.begin(), selected_effects.end(), generator);

    int counter = 0;
    while (counter < selected_effects.size() - 1)
    {
        if (bd(generator))
        {
            int epi_index = FunctionalSequence::getPairIndex(getMutationFromVectorIndex(selected_effects[counter]),
                                                             getMutationFromVectorIndex(selected_effects[counter + 1]));

            epistasis[epi_index] = lnd(generator);
            counter += 2; // next pair
        }
        ++counter;
    }

    return epistasis;
}

std::vector<double> FunctionalSequence::readEpistasis(const std::string& inputPath)
{
    fs::path epiFile(fs::canonical(inputPath) / "pairwise_epistasis.txt");

    // TODO: make more efficient by using the param file in inputPath to know
    //  the size of the vector beforehand
    std::vector<double> epistasis;

    std::ifstream infile(epiFile);
    std::string line;
    while (std::getline(infile, line))
    {
        // each line in this file contains one item (1 column file)
        epistasis.push_back(std::stod(line));
    }
    return epistasis;
}

const double FunctionalSequence::getKd(const Mutation& p) const
{
    return kds.at(getVectorIndex(p));
}

unsigned int FunctionalSequence::getVectorIndex(const Mutation& m) const
{
    // all positions x symbols before (+ 1 for actual position -1 as index starts at 0) + the symbol of the actual
    // position
    return ((m.getPosition() - 1) * (this->params.Q - 1) + m.getSymbol());
}

Mutation FunctionalSequence::getMutationFromVectorIndex(const int index) const
{
    int pos = std::floor(index / (this->params.Q - 1)) + 1; // +1 because pos is 1-indexed
    int sym = index % (this->params.Q - 1);

    return Mutation(pos, sym);
}

std::vector<int> FunctionalSequence::constructPairIndex() const
{
    /*
    0: (1,2)
    1: (1,3)
    2: (1,4)
    3: (1,5)
    4: (2,3)
    5: (2,4)
    6: (2,5)
    7: (3,4)
    8: (3,5)
    9: (4,5)

    (1,X)	1	->	4
    (2,X)	2	->	3
    (3,X)	3	->	2
    (4,X)	4	->	1

    n_pairs = 1 2 3 4
    n_pairs = 4 3 2 1

    tmp = 4 7 9 10
    preceedingPosPairs = 0 4 7 9 10

    (3,5) --> preceedingPosPairs[3-1] + (5-3-1) = preceedingPosPairs[2] + 1 = 7 + 1 = 8 (0-based index)
    */

    // pre-calculate number of pairs preceeding a pair(pos_i, ...) (triangle numbers)
    std::vector<int> n_pairs(this->params.L - 1); // (1,2), (1,3), (1,4), ..., (L-1,L)
    std::iota(n_pairs.begin(), n_pairs.end(), 1); // 1, 2, 3, ..., L-1
    std::reverse(n_pairs.begin(), n_pairs.end()); // L-1, L-2, ..., 1
    int tmp[this->params.L - 1];
    std::partial_sum(n_pairs.begin(), n_pairs.end(), tmp); // L-1, (L-1)+(L-2), ..., L*(L-1)/2
    std::vector<int> preceedingPosPairs(tmp, tmp + (params.L - 1));
    preceedingPosPairs.insert(preceedingPosPairs.begin(), 0);

    return preceedingPosPairs;
}

int FunctionalSequence::getPairIndex(const int p1, const int p2, const int sym1, const int sym2) const
{
    if (p1 >= p2)
    {
        std::cerr << "p1=" << p1 << ", p2=" << p2 << std::endl;
        throw std::invalid_argument("p1 should be smaller than p2");
    }

    int nPreceedingPos1Pairs = this->pairIndex[p1 - 1]; // getPosition is 1-based
    int nPreceedingPos2Pairs = p2 - p1 - 1;
    int index = nPreceedingPos1Pairs + nPreceedingPos2Pairs;
    index = index * std::pow(this->params.Q - 1, 2); // for number of mutation pairs in alphabet
    index += (this->params.Q - 1) * sym1 + sym2;

    return index;
}

int FunctionalSequence::getPairIndex(const Mutation& a, const Mutation& b) const
{
    Mutation m1 = std::min(a, b);
    Mutation m2 = std::max(a, b);

    return FunctionalSequence::getPairIndex(m1.getPosition(), m2.getPosition(), m1.getSymbol(), m2.getSymbol());
}

void FunctionalSequence::writeKdsToFile(const std::string& filename)
{
    std::ofstream outfile(filename);

    if (outfile.good())
    {
        std::for_each(kds.cbegin(), kds.cend(), [&outfile](const auto& entry) { outfile << entry << '\n'; });
    }
}
void FunctionalSequence::writeEpistasisToFile(const std::string& filename)
{
    std::ofstream outfile(filename);

    if (outfile.good())
    {
        std::for_each(epistasis.cbegin(), epistasis.cend(),
                      [&outfile](const auto& entry) { outfile << entry << '\n'; });
    }
}

const double& FunctionalSequence::getEpistasis(const Mutation& a, const Mutation& b) const
{
    return epistasis.at(FunctionalSequence::getPairIndex(a, b));
}

const std::vector<double>& FunctionalSequence::getKd() const
{
    return kds;
}

const std::vector<double>& FunctionalSequence::getEpistasis() const
{
    return epistasis;
}