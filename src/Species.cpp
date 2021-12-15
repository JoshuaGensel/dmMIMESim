//
//  Species.cpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 03.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#include "Constants.hpp"
#include "Count.hpp"
#include "FunctionalSequence.hpp"
#include "Generator.hpp"
#include "Species.hpp"
#include "Utils.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <valarray>

namespace species
{
    Species::Species(const unsigned long long id, const constants::Constants& param)
        : specId(id), params(param), numMut(getNumberOfMutationsById()), mutatedPositions(specIdxToMutPos()), count(0),
          mutCountBound(0), mutCountUnbound(0), errorCountBound(0.0), errorCountUnbound(0.0)
    {
    }

    unsigned int Species::getNumberOfMutationsById()
    {
        return species::getNumberOfMutationsById(this->specId, this->params);
    }

    mutVector Species::specIdxToMutPos()
    {
        return species::specIdxToMutPos(this->specId, this->params);
    }

    // TODO weg
    //    ref::ref_map Species::createRead() {
    //        constants::Constants& constants = constants::Constants::get_instance();
    //        ref::ref_map read;
    //        read.reserve(constants.L);
    //        for(unsigned i=1; i <= constants.L; ++i) {
    //            //TODO Workaround to set wildtype to A and mutations to C
    //            read.add({i, nucleotid::nucleobase{1}});
    //        }
    //
    //        for(auto& mutPos:mutatedPositions) {
    //            auto pos = std::find_if(read.begin(), read.end(), [&mutPos](const auto& val)
    //            {
    //                return val.first == mutPos.getPosition();
    //            });
    //            if(pos != read.end())
    //                read.remove(pos);
    //            //TODO workaround: add the mutation as nucleobase, here simply  + 1 (->wt=A=1), create a wt to mut
    //            interpretation, or make it more general also for AA...and +1 since muts start at 0
    //            read.add({mutPos.getPosition(), nucleotid::nucleobase{int(mutPos.getSymbol()+1+1)}});
    //        }
    //        return read;
    //    }

    const unsigned long long Species::getSpecId() const
    {
        return specId;
    }

    const constants::Constants& Species::getParams() const
    {
        return params;
    }

    unsigned int Species::getCount() const
    {
        return count;
    }

    const unsigned int Species::getNumMut() const
    {
        return numMut;
    }

    const mutVector& Species::getMutatedPositions() const
    {
        return mutatedPositions;
    }

    double Species::getFreq() const
    {
        return this->count / double(this->params.M);
    }

    double Species::getKd() const
    {
        return kd;
    }

    void Species::setCount(unsigned int count)
    {
        Species::count = count;
    }

    void Species::incrementCount()
    {
        ++this->count;
    }

    void Species::decrementCount()
    {
        --this->count;
    }
    //    const ref::ref_map &Species::getRead() const {
    //        return read;
    //    }

    std::valarray<char> Species::getSequence(std::set<Mutation> errors)
    {
        std::map<const unsigned int, char> mutMap = {
            {0, 'C'},
            {1, 'G'},
            {2, 'T'},
        };

        std::valarray<char> sequence('A', this->params.L);
        for (Mutation mut : this->mutatedPositions)
        {
            sequence[mut.getPosition() - 1] = mutMap[mut.getSymbol()];
        }

        if (errors.size() > 0)
        {
            // if a real mutation has error, the according symbol need to be updated. In case it turns into wild
            // type delete it
            for (const auto& mut : this->mutatedPositions)
            {
                auto it = errors.find(mut);
                if (it != errors.end())
                {
                    auto mutSymbol = mut.getSymbol();
                    auto errorSymbol = it->getSymbol();
                    // in case the real mutations symbol is drawn, read it as wild type and correct the count (erase
                    // from rest of errors)
                    if (errorSymbol == mutSymbol)
                    {
                        sequence[mut.getPosition() - 1] = 'A';
                    }
                    else
                    {
                        // otherwise replace with new mutation
                        sequence[mut.getPosition() - 1] = mutMap[errorSymbol];
                    }
                    errors.erase(it);
                }
            }

            // add remaining errors
            for (const auto& err : errors)
            {
                sequence[err.getPosition() - 1] = mutMap[err.getSymbol()];
            }
        }
        return sequence;
    }

    unsigned int Species::getMutCountBound() const
    {
        return mutCountBound;
    }

    void Species::setMutCountBound(unsigned int mutCountBound)
    {
        Species::mutCountBound = mutCountBound;
    }

    unsigned int Species::getMutCountUnbound() const
    {
        return mutCountUnbound;
    }

    void Species::setMutCountUnbound(unsigned int mutCountUnbound)
    {
        Species::mutCountUnbound = mutCountUnbound;
    }

    int Species::getErrorCountBound() const
    {
        return errorCountBound;
    }

    void Species::setErrorCountBound(int errorCountBound)
    {
        Species::errorCountBound = errorCountBound;
    }

    void Species::addErrorCountBound(int errorCountBound)
    {
        Species::errorCountBound += errorCountBound;
    }

    int Species::getErrorCountUnbound() const
    {
        return errorCountUnbound;
    }

    void Species::setErrorCountUnbound(int errorCountUnbound)
    {
        Species::errorCountUnbound = errorCountUnbound;
    }

    void Species::addErrorCountUnbound(int errorCountUnbound)
    {
        Species::errorCountUnbound += errorCountUnbound;
    }

    double Species::getTotalFractionBound()
    {
        return mutCountBound / double(this->params.M);
    }

    double Species::getTotalFractionUnbound()
    {
        return mutCountUnbound / double(this->params.M);
    }

    double Species::getFractionBound()
    {
        return getTotalFractionBound() / (getTotalFractionBound() + getTotalFractionUnbound());
    }

    double Species::getFractionUnbound()
    {
        return getTotalFractionUnbound() / (getTotalFractionBound() + getTotalFractionUnbound());
    }

    void Species::computeSpeciesKd()
    {
        FunctionalSequence* effects = FunctionalSequence::get_instance();
        Species::kd = 1.0;
        // TODO nachfragen: wie setzt sich der Gesamteffekt zusammen? prod(Kd_i)*prod(e_ij) ?  oder prod(e_ij^2)?
        // additive effect of epistasis (since we have the exponential of the epistasis here, it is multiplicative
        for (auto mutPos1_it = begin(Species::mutatedPositions); mutPos1_it != end(Species::mutatedPositions);
             ++mutPos1_it)
        {
            Species::kd *= effects->getKd(*mutPos1_it);
            // std::cout << "mut pos1 " << *mutPos1_it << std::endl;
            for (auto mutPos2_it = mutPos1_it + 1;
                 mutPos1_it != Species::mutatedPositions.end() && mutPos2_it != Species::mutatedPositions.end();
                 ++mutPos2_it)
            {
                // std::cout << "mut pos2 " << *mutPos2_it << std::endl;
                // compute the kd for a species by multiplying all single kds of the mutations and add the pairwise
                // epistasis factor totalEpistasisPerPos[*mutPos1_it-1] *= effects->getEpistasis(*mutPos1_it-1,
                // *mutPos2_it-1); totalEpistasisPerPos[*mutPos2_it-1] *= effects->getEpistasis(*mutPos1_it-1,
                // *mutPos2_it-1);
                Species::kd *= effects->getEpistasis(*mutPos1_it, *mutPos2_it);
            }
        }
    }

    species_map drawSpeciesIds(const constants::Constants& params)
    {
        std::default_random_engine& generator = Generator::get_instance()->engine;

        // contains the map with all sequence species
        species::species_map species_map;

        // Break down the drawing of all possible (allowed) species ids into 2 smaller ones:
        // first draw a the number of mutations from 0 to MAX_MUT, with the given probabilities...
        std::discrete_distribution<> d(begin(params.P_NMUT), end(params.P_NMUT));
        // then draw uniformly the id from the id range for this particular number of mutations
        std::vector<std::uniform_int_distribution<>> unif(params.MAX_MUT + 1);

        // create distributions for all numbers of mutations beforehand
        for (int numMut = 1; numMut <= params.MAX_MUT + 1; ++numMut)
        {
            unif[numMut - 1] =
                std::uniform_int_distribution<>(params.NMUT_RANGE[numMut - 1] + 1, params.NMUT_RANGE[numMut]);
        }
        // count the given species
        for (int n = 0; n < params.M; ++n)
        {
            // draw number of mutations
            const int numMut = d(generator);
            // if no mutations, the id is always 1
            int id = 1;
            if (numMut > 0)
            {
                id = unif[numMut - 1](generator);
            }
            // create new object if not yet present (return value gives iterator and flag if insertion happened)
            // the id is the key for the map, and also the parameter for the constructor for the species class
            auto currentEntry = species_map.try_emplace(id, id, params);
            if (currentEntry.second)
                currentEntry.first->second.computeSpeciesKd();
            currentEntry.first->second.incrementCount();
            //            if(species_map.find(id) == species_map.end()) {
            //                //the id is the key for the map, and also the parameter for the constructor for the
            //                species class species_map.emplace(id, id); species_map.at(id).computeSpeciesKd();
            //
            //            }
            //            species_map.at(id).incrementCount();
        }
        return species_map;
    }

    species_map readFromFile(const std::string& inputPath, utils::SampleID id, const constants::Constants& params)
    {
        fs::path speciesFile(fs::canonical(inputPath) / "species" / (utils::SampleIDStr(id) + ".txt"));

        species::species_map species_vec;

        std::ifstream infile(speciesFile);
        std::string line;

        std::getline(infile, line); // header row
        while (std::getline(infile, line))
        {
            // per row: SpeciesID \t #instances
            std::vector<std::string> row_values;
            utils::split_string(line, '\t', row_values);

            unsigned long long id = std::stoull(row_values.at(0));
            int count = std::stoi(row_values.at(1));

            auto currentObj = species_vec.emplace(id, Species{id, params});
            while (count > 0)
            {
                currentObj.first->second.incrementCount();
                --count;
            }
            currentObj.first->second.computeSpeciesKd();
        }
        return species_vec;
    }

    unsigned int speciesMapSum(species_map species_vec)
    {
        unsigned int count = 0;
        for (auto it = species_vec.begin(); it != species_vec.end(); ++it)
        {
            count = count + it->second.getCount();
        }
        return count;
    }

    void writeSpeciesToFile(const std::string& out_file, species_map& spec_map, std::valarray<unsigned int>& S_pool,
                            const std::string& header)
    {
        std::ofstream outfile(out_file);

        if (outfile.good())
        {
            outfile << header;

            int specIdx = 0;
            for (auto it = spec_map.begin(); it != spec_map.end(); ++it)
            {

                // only include non-zero occurences
                if (S_pool[specIdx] > 0)
                {
                    // print SpeciesID
                    outfile << it->first;
                    // print number of sequences in pool
                    outfile << '\t' << S_pool[specIdx] << '\n';
                }
                ++specIdx;
            }
        }
    }

    void writeSequencesToFile(const std::string& out_file, species_map& spec_map,
                              std::vector<std::set<Mutation>> errors, std::valarray<unsigned int> S_pool)
    {
        std::ofstream outfile(out_file);
        if (outfile.good())
        {
            int specIdx = 0;
            int errIdx = 0;
            for (auto it = spec_map.begin(); it != spec_map.end(); ++it)
            {
                // only include non-zero occurences
                while (S_pool[specIdx] > 0)
                {
                    // print sequence string
                    std::string sequence;
                    for (char c : it->second.getSequence(errors[errIdx]))
                    {
                        sequence.push_back(c);
                    }
                    outfile << sequence << '\n';
                    --S_pool[specIdx];
                    ++errIdx;
                }
                ++specIdx;
            }
        }
    }

    // TODO testen
    std::set<Mutation> drawError_2(const constants::Constants& params)
    {
        std::default_random_engine& generator = Generator::get_instance()->engine;
        std::binomial_distribution<int> bino(params.L, params.P_ERR);
        // random generator for the position with an error
        std::uniform_int_distribution<> unif_err(1, params.L);
        // random generator for the mutated symbol (if sequence symbol is wt, the errror is one of the mutations;
        //  if the sequence symbol is mutated, the according symbol is the wild type or one of the other mutations
        std::uniform_int_distribution<> unif_sym(0, params.Q - 1 - 1);

        // draw number of errors
        int numErrors = bino(generator);

        // TODO: mal weg lassen, da ich ja eh keine species berechne
        // numErrors = std::min<int>(numErrors, params.MAX_MUT*2 - mutations.size());

        // containing numErrors Mutations with unique position
        std::set<Mutation> uniquePositions;
        if (numErrors > 0)
        {
            while (uniquePositions.size() < numErrors)
            {
                // draw position with error
                auto position = unif_err(generator);
                //...and the symbol
                auto symbol = unif_sym(generator);
                // first time we sample this position
                uniquePositions.emplace(position, symbol);
            }
        }
        return uniquePositions;
    }

    // TODO: rekursiver Aufruf? Aber dafür müsste jedesmal für irgendein L' (Restlänger nach aktueller Position) die ID
    // ranges berechnet werden, oder mache ich das eh?
    mutVector specIdxToMutPos(const unsigned long long specId, const constants::Constants& params)
    {
        auto numMut = getNumberOfMutationsById(specId, params);
        // collect the mutated position with the respective mutation symbol
        mutVector mutPos;
        mutPos.reserve(numMut);

        // check if Id is valid
        if (specId <= params.NMUT_RANGE.back() && numMut > 0)
        {
            unsigned int Lact = params.L;
            unsigned int numMutAct = numMut;
            auto mSymbols = params.Q - 1;

            // get the id within the range of number of mutations (substract the ids for the sequences with less
            // mutations)
            unsigned int idAct = specId - params.NMUT_RANGE[numMut - 1];

            // determine each mutations position seen from the mutations position before...
            for (unsigned int m = 0; m < numMut; ++m)
            {
                // for each possible positions within the length the actual mutation covers a range of ids depending on
                // the residual mutations to follow
                std::vector<unsigned long> cumSumRange(Lact - (numMutAct - 1));
                // bool indexFound = false;
                unsigned int i = 0;
                // initialise first value of the vector for cummulative sum (do it so complicated to not compute the
                // whole range if not necessary)
                cumSumRange[0] = utils::nChoosek(Lact - 1, numMutAct - 1) * std::pow(mSymbols, numMut);
                // find the id within the ranges and get the index (=position)
                while (idAct > cumSumRange[i] && i < Lact - (numMutAct - 1))
                {
                    ++i;
                    // for(; i<Lact-(numMutAct-1) && !indexFound; ++i) {
                    // TODO test
                    cumSumRange[i] = utils::nChoosek(Lact - i - 1, numMutAct - 1) * std::pow(mSymbols, numMut);
                    if (i > 0)
                    {
                        cumSumRange[i] += cumSumRange[i - 1];
                    }
                }

                // TODO weg
                // mutPos[m] = i + 1;
                //  the symbol of a posisition  is given for (q-1) ^ numMut-1 times (e.g. with three mutations and 2
                //  symbols 2 ^2 times) : AAA, AAB, ABA, ABB
                int symbolCombiPerPos = std::pow(mSymbols, numMutAct - 1);
                // TODO getestet? wenn ja, comments weg
                // find symbol
                unsigned mut = (int)std::floor((idAct - 1) / symbolCombiPerPos) % mSymbols;

                unsigned int pos = i + 1; // TODO Utils hat noch das- *cumSumRange.begin()
                unsigned int prePos = mutPos.begin() == mutPos.end() ? 0 : mutPos.rbegin()->getPosition();

                // arguments: the two pair_constructor parameter pos and mut, adding the last cummulative position
                //  (= position seend frim the beginning of the sequence)
                mutPos.emplace_back(pos + prePos, mut);

                // indexFound = true;

                // the residual length after the actual mutation
                Lact = Lact - pos;
                // the redsiudal number of mutations after the actual mutation
                --numMutAct;
                // the id within the residual length (-1 because the indices start at 0 and another -1 because we
                // substract the ids of the preceeding mutations range
                idAct = idAct - (pos == 1 ? 0 : cumSumRange[i - 1]); // hier stimmt was nicht
            }
        }
        //.... and get the correct positions within the sequence with cumsum
        // TODO in die schleife rein, wie in Utils..
        // std::partial_sum(mutPos.begin(), mutPos.end(), mutPos.begin());

        return (mutPos);
    }

    // TODO testen! (vorallem das mit partial sum)
    unsigned long long mutPosToSpecIdx(const mutVector& mutPos, const constants::Constants& params)
    {
        unsigned numMut = mutPos.size();
        // id for 0 mutations is 1
        unsigned long long specId = 1;
        if (numMut > 0)
        {
            // add the ids for the sequences with less mutations
            specId = params.NMUT_RANGE[numMut - 1];
            mutVector mutPos_new = mutPos;
            // get the indices for the individual length segments for each position
            std::partial_sum(mutPos.begin(), mutPos.end(), mutPos_new.begin(),
                             [](const Mutation& x, const Mutation& y) {
                                 return Mutation(y.getPosition() - x.getPosition(), y.getSymbol());
                             });
            // TODO weg get the indices for the individual length segments for each position
            //             for(unsigned i = 1; i<mutPos_new.size(); ++i) {
            //                 mutPos_new[i].getPosition() -= mutPos[i-1].getPosition();
            //             }
            unsigned Lact = params.L;
            unsigned numMutAct = numMut;

            // Notiz an mich selbst: enforcing const elements in range iteration (C++17)
            for (auto const& mutation : std::as_const(mutPos_new))
            {
                if (numMutAct == 1)
                {
                    // TODO überall checken: the possible mutations are represented from 0 upwards
                    specId += mutation.getPosition() + mutation.getSymbol();
                }
                else
                {
                    // if the position of the actual mutations is != 1, the id is depending on the next position
                    if (mutation.getPosition() != 1)
                    {
                        for (int i = 1; i < mutation.getPosition(); ++i)
                        {
                            // specId += utils::nChoosek(Lact-i, numMutAct-1);
                            specId += utils::nChoosek(Lact - i, numMutAct - 1) * pow(params.Q - 1, numMutAct);
                        }
                    }
                    Lact -= mutation.getPosition();
                    --numMutAct;
                }
            }
        }

        return specId;
    }

    unsigned getNumberOfMutationsById(const unsigned long long specId, const constants::Constants& params)
    {
        // gives the index where the content is still lower than the given id
        auto low_it = std::lower_bound(std::begin(params.NMUT_RANGE), std::end(params.NMUT_RANGE), specId);
        return (low_it - std::begin(params.NMUT_RANGE));
    }

    void countErrors(const unsigned int S, const mutVector& mutatedPositions, const constants::Constants& params,
                     count::counter_1& counter_1d, count::counter_2& counter_2d)
    {
        // sample error for all sequences of the bound / unbound fraction
        for (int b = 0; b < S; ++b)
        {
            auto uniquePositions = drawError_2(params);
            if (uniquePositions.size() > 0)
            {
                // if a real mutation has error, the according symbol need to be updated. In case it turns into wild
                // type delete it
                for (const auto& mut : mutatedPositions)
                {
                    auto it = uniquePositions.find(mut);
                    if (it != uniquePositions.end())
                    {
                        //+1 because symbols start from 1, +1 because wt is 1
                        auto mutSymbol = mut.getSymbol() + 2;
                        auto errorSymbol = it->getSymbol() + 2;
                        counter_1d.decrement(mut.getPosition(), mutSymbol);
                        counter_2d.count(mut.getPosition(), mutSymbol, -1);
                        // in case the real mutations symbol is drawn, read it as wild type and correct the count (erase
                        // from rest of errors)
                        if (errorSymbol == mutSymbol)
                        {
                            counter_1d.increment(mut.getPosition(), 1);
                            counter_2d.count(mut.getPosition(), 1, 1);
                        }
                        else
                        {
                            // otherwise replace with new mutation
                            counter_1d.increment(mut.getPosition(), errorSymbol);
                            counter_2d.count(mut.getPosition(), errorSymbol, 1);
                        }
                        uniquePositions.erase(it);
                    }
                }

                // add errors to count
                for (const auto& err : uniquePositions)
                {
                    auto errorSymbol = err.getSymbol() + 2;

                    counter_1d.decrement(err.getPosition(), 1);
                    counter_1d.increment(err.getPosition(), errorSymbol);

                    counter_2d.count(err.getPosition(), 1, -1);
                    counter_2d.count(err.getPosition(), errorSymbol, 1);
                }
            }
        }
    }

    // TODO neu: count directlly all mutations and adding errors
    count::counter_collection countMutationsWithErrors(const std::valarray<unsigned int>& SBound,
                                                       const std::valarray<unsigned int>& SUnbound,
                                                       const species_map& spec_map, const constants::Constants& params)
    {
        // to get the correct counts from the valarrays increment the index

        // containg counter for the bound and unbound fractions
        auto counters = count::counter_collection(params.L, params.Q);

        auto boundSum = SBound.sum();
        auto unboundSum = SUnbound.sum();

        // Because the majority is wildtype, count all als wildtype (symbol = 1)....
        for (unsigned i = 1; i <= params.L; ++i)
        {
            counters.counter_bound_1d.count(i, 1, boundSum);
            counters.counter_unbound_1d.count(i, 1, unboundSum);
            for (unsigned j = i + 1; j <= params.L && i < params.L; ++j)
            {
                counters.counter_bound_2d.count(i, j, 1, 1, boundSum);
                counters.counter_unbound_2d.count(i, j, 1, 1, unboundSum);
            }
        }

        int specIdx = 0;
        // for (auto it = spec_map.begin(); it != spec_map.end(); ) {
        for (auto it = spec_map.begin(); it != spec_map.end(); ++it)
        {

            // add the counts for the mutated positions (and substract it from the initial wild type value
            for (auto mutIt = it->second.getMutatedPositions().begin(); mutIt != it->second.getMutatedPositions().end();
                 ++mutIt)
            {
                counters.counter_bound_1d.count((*mutIt).getPosition(), (*mutIt).getSymbol() + 2, SBound[specIdx]);
                counters.counter_bound_1d.count((*mutIt).getPosition(), 1, -SBound[specIdx]);

                counters.counter_unbound_1d.count((*mutIt).getPosition(), (*mutIt).getSymbol() + 2, SUnbound[specIdx]);
                counters.counter_unbound_1d.count((*mutIt).getPosition(), 1, -SUnbound[specIdx]);

                //.... substract the ones where position 1 is mutated, and count as wt mut
                counters.counter_bound_2d.count((*mutIt).getPosition(), (*mutIt).getSymbol() + 2, SBound[specIdx]);
                counters.counter_bound_2d.count((*mutIt).getPosition(), 1, -SBound[specIdx]);

                counters.counter_unbound_2d.count((*mutIt).getPosition(), (*mutIt).getSymbol() + 2, SUnbound[specIdx]);
                counters.counter_unbound_2d.count((*mutIt).getPosition(), 1, -SUnbound[specIdx]);

                // change the counts for the double mutants
                for (auto mut2It = mutIt + 1; mut2It != it->second.getMutatedPositions().end(); ++mut2It)
                {
                    counters.counter_bound_2d.count((*mutIt).getPosition(), (*mut2It).getPosition(),
                                                    (*mutIt).getSymbol() + 2, (*mut2It).getSymbol() + 2,
                                                    SBound[specIdx]);
                    counters.counter_bound_2d.count((*mutIt).getPosition(), (*mut2It).getPosition(),
                                                    (*mutIt).getSymbol() + 2, 1, -SBound[specIdx]);

                    counters.counter_unbound_2d.count((*mutIt).getPosition(), (*mut2It).getPosition(),
                                                      (*mutIt).getSymbol() + 2, (*mut2It).getSymbol() + 2,
                                                      SUnbound[specIdx]);
                    counters.counter_unbound_2d.count((*mutIt).getPosition(), (*mut2It).getPosition(),
                                                      (*mutIt).getSymbol() + 2, 1, -SUnbound[specIdx]);
                }
            }

            // sample error for bound sequences
            countErrors(SBound[specIdx], it->second.getMutatedPositions(), it->second.getParams(),
                        counters.counter_bound_1d, counters.counter_bound_2d);
            countErrors(SUnbound[specIdx], it->second.getMutatedPositions(), it->second.getParams(),
                        counters.counter_unbound_1d, counters.counter_unbound_2d);

            ++specIdx;
            // spec_map.erase(it++);
        }
        return counters;
    }

    species_map combineSpecies(species_map firstSpec, species_map secondSpec, const constants::Constants& firstParams,
                               const constants::Constants& secondParams)
    {
        assert((firstParams.L == secondParams.L) && "Both species sets should have the same length!");
        assert((firstParams.Q == secondParams.Q) && "Both species sets should have the same alphabet!");

        constants::Constants* combined_params = new constants::Constants(firstParams, secondParams); // 2MAX_MUT

        assert((speciesMapSum(firstSpec) == speciesMapSum(secondSpec)) &&
               "Both species sets should have the same number of sequences!");
        int total = speciesMapSum(firstSpec);

        // contains the map with all sequence species
        species::species_map species_map;

        int done_counter = 0;
        // loop through both maps simultaneously
        for (auto it1 = firstSpec.begin(), end1 = firstSpec.end(), it2 = secondSpec.begin(), end2 = secondSpec.end();
             it1 != end1 || it2 != end2;)
        {
            assert((it1->second.getCount() != 0) && (it2->second.getCount() != 0) &&
                   "This should not happen if both species sets should have the same number of sequences");
            assert((done_counter < total) && "This should not happen if iterator loop is well-defined");

            mutVector combined = combineMutations(it1->second.getMutatedPositions(), it2->second.getMutatedPositions());
            unsigned int id = mutPosToSpecIdx(combined, *combined_params);

            // create new object if not yet present (return value gives iterator and flag if insertion happened)
            // the id is the key for the map, and also the parameter for the constructor for the species class
            auto currentEntry = species_map.try_emplace(id, id, *combined_params);
            if (currentEntry.second)
                currentEntry.first->second.computeSpeciesKd();
            currentEntry.first->second.incrementCount();

            it1->second.decrementCount();
            it2->second.decrementCount();
            ++done_counter;

            if ((it1->second.getCount() == 0) && (it1 != end1))
            {
                ++it1;
            }

            if ((it2->second.getCount() == 0) && (it2 != end2))
            {
                ++it2;
            }
        }
        return species_map;
    }

    mutVector combineMutations(mutVector firstMut, mutVector secondMut)
    {
        mutVector comb_mut(secondMut);
        // newest mutations appear before old ones in vector
        comb_mut.insert(comb_mut.end(), firstMut.begin(), firstMut.end());

        // A duplicate mutation means the position is mutated back to the wildtype
        // and hence, that this mutation should be removed from the mutVector.

        // Create a map to store the frequency of each mutations in vector
        std::map<Mutation, int> countMutMap;
        for (auto it = comb_mut.begin(); it != comb_mut.end();)
        {
            auto result = countMutMap.insert(std::pair<Mutation, int>(*it, 1));
            if (result.second == false) // element already exists
            {
                it = comb_mut.erase(it);
            }
            else
            {
                ++it;
            }
        }

        // When two different mutations exist for the same position, keep the newest one.

        // Create a map to store the frequency of each position in vector
        std::map<unsigned int, int> countPosMap;
        for (auto it = comb_mut.begin(); it != comb_mut.end();)
        {
            auto result = countPosMap.insert(std::pair<unsigned int, int>(it->getPosition(), 1));
            if (result.second == false) // element already exists
            {
                it = comb_mut.erase(it);
            }
            else
            {
                ++it;
            }
        }

        std::sort(comb_mut.begin(), comb_mut.end());
        return comb_mut;
    }
}
