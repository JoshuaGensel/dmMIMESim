//
// Created by Smith, Maureen on 31.05.18.
//

#include "Constants.hpp"
#include "FunctionalSequence.hpp"
#include "Generator.hpp"
#include "Species.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>

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

std::vector<double> FunctionalSequence::drawEpistasis()
{
    std::vector<double> epistasis(this->params.PWVal);
    std::default_random_engine& generator = Generator::get_instance()->engine;
    std::bernoulli_distribution bd(this->params.P_EPISTASIS);
    std::lognormal_distribution<double> lnd(0, 1);
    // TODO: Kann epistasis auch vorhanden sein, wenn eine einzelmutation keinen effect hat? (hatte es in R so
    // realisiert dass beide einene Effekt haben müssen.
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
std::vector<double> FunctionalSequence::drawEpistasis(const std::vector<double> kds)
{
    std::vector<double> epistasis(this->params.PWVal);
    std::fill(epistasis.begin(), epistasis.end(), this->params.NO_EPISTASIS);

    std::default_random_engine& generator = Generator::get_instance()->engine;
    std::bernoulli_distribution bd(this->params.P_EPISTASIS);
    std::lognormal_distribution<double> lnd(0, 1);

    // pre-calculate number of pairs preceeding pairs(pos_i, ...) (triangle numbers)
    std::vector<int> n_pairs(this->params.L - 1); // (1,2), (1,3), (1,4), ..., (L-1,L)
    std::iota(n_pairs.begin(), n_pairs.end(), 1); // 1, 2, 3, ..., L-1
    std::reverse(n_pairs.begin(), n_pairs.end()); // L-1, L-2, ..., 1
    int index_offset[this->params.L - 1];
    std::partial_sum(n_pairs.begin(), n_pairs.end(), index_offset); // L-1, (L-1)+(L-2), ..., L*(L-1)/2

    std::vector<int> has_effect; // single position has effect

    for (int i = 0; i < kds.size(); ++i)
    {
        if (kds[i] != this->params.KD_WT)
        {
            has_effect.push_back(i);
        }
    }
    std::shuffle(has_effect.begin(), has_effect.end(), generator);

    int counter = 0;
    while (counter < has_effect.size() - 1)
    {
        if (bd(generator))
        {
            int pos1 = std::min(has_effect[counter], has_effect[counter + 1]);
            int pos2 = std::max(has_effect[counter], has_effect[counter + 1]);

            int pair_index = pos2 + (pos1 > 0 ? index_offset[pos1 - 1] : 0);
            epistasis[pair_index] = lnd(generator);
            ++counter;
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

// TODO test
unsigned int FunctionalSequence::getVectorIndex(const Mutation& m) const
{
    // all positions x symbols before (+ 1 for actual position -1 as index starts at 0) + the symbol of the actual
    // position
    return ((m.getPosition() - 1) * (this->params.Q - 1) + m.getSymbol());
}

// TODO test
unsigned long long FunctionalSequence::getMatrixVectorIndex(const Mutation& a, const Mutation& b) const
{
    // TODO weg damit nach test?
    // i must always be smaller than j
    //    if(pos1 > pos2) {
    //        auto k = pos2;
    //        auto l = mut2;
    //        pos2=pos1;
    //        pos1=k;
    //        mut2=mut1;
    //        mut1=l;
    //    }
    // same as in R mut -1 because index starts at 0
    // for the index calculation, the positions have to be in ascending order
    // unsigned int res = (c.L*(c.L-1)/2) - ((c.L-i+1)*((c.L-i+1)-1)/2) + j - i - 1;
    // TODO oha... determine the id of the sequence with pairwise mutations and substract the ID range for the
    // sequences with no or 1 mutation
    unsigned long long res;
    if (a.getPosition() < b.getPosition())
        // TODO diese Funktionen woanders hin als in Species? Eigentlich hat Species ja nichts mit Functional
        // Sequence zu tun
        res = species::mutPosToSpecIdx({a, b}, this->params) - this->params.NMUT_RANGE.at(1) - 1;
    else if (a.getPosition() > b.getPosition())
        res = species::mutPosToSpecIdx({b, a}, this->params) - this->params.NMUT_RANGE.at(1) - 1;
    else
        // TODO throw exception for i == j?
        std::cerr << "Two mutations at the same position are not possible.";

    return res;
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
    return epistasis.at(FunctionalSequence::getMatrixVectorIndex(a, b));
}

const std::vector<double>& FunctionalSequence::getKd() const
{
    return kds;
}

const std::vector<double>& FunctionalSequence::getEpistasis() const
{
    return epistasis;
}