//
// Created by Smith, Maureen on 31.05.18.
//

#include "Constants.hpp"
#include "FunctionalSequence.hpp"
#include "Species.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>

std::vector<double> FunctionalSequence::drawKdValues()
{
    auto& constants = constants::Constants::get_instance();
    std::vector<double> kds(constants.SVal);
    // const auto seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
    const auto seed = static_cast<unsigned int>(0);
    std::default_random_engine generator(seed);
    std::bernoulli_distribution bd(constants.P_EFFECT);
    std::lognormal_distribution<double> lnd(0, 1);
    for (int i = 0; i < kds.size(); ++i)
    {
        // first sample the positions which have an effect (binomial/bernoulli distributed), then the strength of the
        // effect(log normal distributed)
        kds[i] = bd(generator) ? lnd(generator) : constants.KD_WT;
    }
    return kds;
}

std::vector<double> FunctionalSequence::drawEpistasis()
{
    auto& constants = constants::Constants::get_instance();
    std::vector<double> epistasis(constants.PWVal);
    // const auto seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
    const auto seed = static_cast<unsigned int>(0);
    std::default_random_engine generator(seed);
    std::bernoulli_distribution bd(constants.P_EPISTASIS);
    std::lognormal_distribution<double> lnd(0, 1);
    // TODO: Kann epistasis auch vorhanden sein, wenn eine einzelmutation keinen effect hat? (hatte es in R so
    // realisiert dass beide einene Effekt haben müssen.
    for (int i = 0; i < epistasis.size(); ++i)
    {
        // first sample if position pair has epistatic effect (bernoulli) and then the value of the epistasis (log
        // normal distributed)
        epistasis[i] = bd(generator) ? lnd(generator) : constants.NO_EPISTASIS;
    }
    return epistasis;
}

FunctionalSequence& FunctionalSequence::get_instance()
{
    static FunctionalSequence instance;
    return instance;
}

const double FunctionalSequence::getKd(const Mutation& p)
{
    return kds.at(getVectorIndex(p));
}

// TODO test
unsigned int FunctionalSequence::getVectorIndex(const Mutation& m)
{
    auto& c = constants::Constants::get_instance();
    // all positions x symbols before (+ 1 for actual position -1 as index starts at 0) + the symbol of the actual
    // position
    return ((m.getPosition() - 1) * (c.Q - 1) + m.getSymbol());
}

// TODO test
unsigned int FunctionalSequence::getMatrixVectorIndex(const Mutation& a, const Mutation& b)
{
    auto& c = constants::Constants::get_instance();

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
    // TODO oha... determine the id of the sequence with pairwise mutations and substract the ID range for the sequences
    // with no or 1 mutation
    unsigned int res;
    if (a.getPosition() < b.getPosition())
        // TODO diese Funktionen woanders hin als in Species? Eigentlich hat Species ja nichts mit Functional Sequence
        // zu tun
        res = species::mutPosToSpecIdx({a, b}) - c.NMUT_RANGE.at(1) - 1;
    else if (a.getPosition() > b.getPosition())
        res = species::mutPosToSpecIdx({b, a}) - c.NMUT_RANGE.at(1) - 1;
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

const double& FunctionalSequence::getEpistasis(const Mutation& a, const Mutation& b)
{
    return epistasis.at(FunctionalSequence::getMatrixVectorIndex(a, b));
}

const std::vector<double>& FunctionalSequence::getKd()
{
    return kds;
}

const std::vector<double>& FunctionalSequence::getEpistasis()
{
    return epistasis;
}
