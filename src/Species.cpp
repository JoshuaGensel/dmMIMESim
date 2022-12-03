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
#include <unordered_map>
#include <valarray>

namespace species
{
    Species::Species(const utils::ids ids, const constants::Constants& param)
        : specIds{ids}, params{param}, numMuts{getNumberOfMutationsByIds()}, mutatedPositions{specIdsToMutPos()},
          count{0}, mutCountBound{0}, mutCountUnbound{0}, errorCountBound{0}, errorCountUnbound{0}
    {
        if (specIds.size() != params.NCHUNKS)
        {
            std::cerr << "specIds vector does not match params object." << std::endl;
            throw std::runtime_error("Mismatch between specIds and params.");
        }
    }

    std::vector<unsigned int> Species::getNumberOfMutationsByIds()
    {
        return species::getNumberOfMutationsByIds(this->specIds, this->params);
    }

    mutVector Species::specIdsToMutPos()
    {
        return species::specIdsToMutPos(this->specIds, this->params);
    }

    const utils::ids Species::getSpecIds() const
    {
        return specIds;
    }

    const constants::Constants& Species::getParams() const
    {
        return params;
    }

    unsigned int Species::getCount() const
    {
        return count;
    }

    const std::vector<unsigned int> Species::getNumMuts() const
    {
        return numMuts;
    }

    const unsigned int Species::getNumMut_total() const
    {
        return std::accumulate(numMuts.begin(), numMuts.end(), 0);
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

        // additive effect of epistasis (since we have the exponential of the epistasis here, it is multiplicative
        for (auto mutPos1_it = begin(Species::mutatedPositions); mutPos1_it != end(Species::mutatedPositions);
             ++mutPos1_it)
        {
            Species::kd *= effects->getKd(*mutPos1_it);

            for (auto mutPos2_it = mutPos1_it + 1; mutPos2_it != end(Species::mutatedPositions); ++mutPos2_it)
            {
                // compute the kd for a species by multiplying all single kds of the mutations and add the pairwise
                // epistasis factor totalEpistasisPerPos[*mutPos1_it-1] *= effects->getEpistasis(*mutPos1_it-1,
                // *mutPos2_it-1); totalEpistasisPerPos[*mutPos2_it-1] *= effects->getEpistasis(*mutPos1_it-1,
                // *mutPos2_it-1);
                Species::kd *= effects->getEpistasis(*mutPos1_it, *mutPos2_it);
            }
        }
    }

    // TODO: parallelize
    species_map drawSpeciesIds(const constants::Constants& params)
    {
        std::default_random_engine& generator = Generator::get_instance()->engine;

        // contains the map with all sequence species
        species::species_map species_map;

        // Break down the drawing of all possible (allowed) species ids into 3 smaller ones:
        // first draw a the number of mutations from 0 to MAX_MUT, with the given probabilities...
        std::discrete_distribution<> d(begin(params.P_NMUT), end(params.P_NMUT));
        // then draw uniformly the chunk indices where the mutations reside in
        std::uniform_int_distribution<unsigned int> chunk_unif(0, params.NCHUNKS - 1);
        // then draw uniformly, per chunk, the id from the id range for the particular number of mutations
        std::vector<std::uniform_int_distribution<utils::id>> id_unif(params.MAX_MUT + 1);

        id_unif[0] = std::uniform_int_distribution<utils::id>(1, 1);
        // create distributions for all numbers of mutations beforehand
        for (int numMut = 1; numMut <= params.MAX_MUT; ++numMut)
        {
            id_unif[numMut] =
                std::uniform_int_distribution<utils::id>(params.NMUT_RANGE[numMut - 1] + 1, params.NMUT_RANGE[numMut]);
        }
        // count the given species
        for (int n = 0; n < params.M; ++n)
        {
            // draw number of mutations
            const int numMut = d(generator);
            // if no mutations, the id is always 1
            utils::ids ids(params.NCHUNKS, 1);

            if (numMut > 0)
            {
                // draw the chunk indices that contain mutations
                std::unordered_map<int, int> numMutPerChunkId;
                std::unordered_map<int, int>::iterator it;

                for (int i = 0; i < numMut; ++i)
                {
                    // check if key `c` exists in the map or not
                    unsigned int chunk_index = chunk_unif(generator);

                    it = numMutPerChunkId.find(chunk_index);
                    // key already present on the map
                    if (it != numMutPerChunkId.end())
                    {
                        it->second++; // increment map's value for key `chunk_index`
                    }
                    // key not found
                    else
                    {
                        numMutPerChunkId.insert(std::make_pair(chunk_index, 1));
                    }
                }

                // draw the chunk ids
                for (it = numMutPerChunkId.begin(); it != numMutPerChunkId.end(); it++)
                {
                    ids[it->first] = id_unif[it->second](generator);
                }
            }
            // create new object if not yet present (return value gives iterator and flag if insertion happened)
            // the id is the key for the map, and also the parameter for the constructor for the species class
            auto currentEntry = species_map.emplace(ids, species::Species(ids, params));
            if (currentEntry.second)
            {
                currentEntry.first->second.computeSpeciesKd();
            }
            currentEntry.first->second.incrementCount();
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
            // per row: SpeciesID_chunk0 \t SpeciesID_chunk1 \t ... \t SpeciesID_chunkN-1 \t #instances
            std::vector<std::string> row_values;
            utils::split_string(line, '\t', row_values);

            utils::ids ids(params.NCHUNKS, 1);
            for (int i = 0; i < params.NCHUNKS; ++i)
            {
                ids[i] = std::stoull(row_values.at(i));
            }
            int count = std::stoi(row_values.at(params.NCHUNKS)); // second-to-last-entry in row

            auto currentObj = species_vec.emplace(ids, Species{ids, params});
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

    void writeSpeciesToFile(const std::string& out_file, unsigned int& n_chunks, species_map& spec_map,
                            std::valarray<unsigned int>& S_pool)
    {

        std::string header = "";
        for (int i = 0; i < n_chunks; ++i)
        {
            header += "speciesID_" + std::to_string(i) + "\t";
        }
        header += "count\tKd\n";
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
                    // print SpeciesID chunks
                    for (utils::id id : it->first)
                        outfile << std::to_string(id) << '\t';
                    // print number of sequences in pool
                    outfile << '\t' << S_pool[specIdx];
                    // print species Kd
                    outfile << '\t' << it->second.getKd() << '\n';
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

    std::vector<std::set<Mutation>> drawErrors(const constants::Constants& params)
    {
        std::vector<std::set<Mutation>> errors;
        errors.reserve(params.M);

        for (int i = 0; i < params.M; ++i)
        {
            errors.push_back(drawError_2(params));
        }
        return errors;
    }

    mutVector specIdToMutPos(const utils::id specId, const constants::Constants& params, unsigned int offset)
    {
        unsigned int numMut = getNumberOfMutationsById(specId, params);
        mutVector mutPos;
        mutPos.reserve(numMut);

        // check if Id is valid
        if (specId <= params.NMUT_RANGE.back() && numMut > 0)
        {
            unsigned int Lact = params.chunkL;
            unsigned int numMutAct = numMut;
            auto mSymbols = params.Q - 1;

            // get the id within the range of number of mutations (substract the ids for the sequences with less
            // mutations)
            unsigned int idAct = specId - params.NMUT_RANGE[numMut - 1];

            // determine each mutations position seen from the mutations position before...
            for (unsigned int m = 0; m < numMut; ++m)
            {
                // for each possible positions within the length the actual mutation covers a range of ids
                // depending on the residual mutations to follow
                std::vector<utils::id> cumSumRange(Lact - (numMutAct - 1));
                unsigned int i = 0;
                // initialise first value of the vector for cummulative sum (do it so complicated to not compute
                // the whole range if not necessary)
                cumSumRange[0] = std::roundl(utils::nChoosek(Lact - 1, numMutAct - 1) * std::pow(mSymbols, numMut));
                // find the id within the ranges and get the index (=position)
                while (idAct > cumSumRange[i] && i < Lact - (numMutAct - 1))
                {
                    ++i;
                    cumSumRange[i] =
                        std::roundl(utils::nChoosek(Lact - i - 1, numMutAct - 1) * std::pow(mSymbols, numMut));
                    if (i > 0)
                    {
                        cumSumRange[i] += cumSumRange[i - 1];
                    }
                }

                //  the symbol of a posisition  is given for (q-1) ^ numMut-1 times (e.g. with three mutations
                //  and 2 symbols 2 ^2 times) : AAA, AAB, ABA, ABB
                int symbolCombiPerPos = std::pow(mSymbols, numMutAct - 1);
                // find symbol
                unsigned mut = (int)std::floor((idAct - 1) / symbolCombiPerPos) % mSymbols;

                unsigned int pos = i + 1;
                unsigned int prePos = mutPos.begin() == mutPos.end() ? offset : mutPos.rbegin()->getPosition();

                // arguments: the two pair_constructor parameter pos and mut, adding the last cummulative
                // position
                //  (= position seend frim the beginning of the sequence)
                mutPos.emplace_back(pos + prePos, mut);

                // the residual length after the actual mutation
                Lact = Lact - pos;
                // the redsiudal number of mutations after the actual mutation
                --numMutAct;
                // the id within the residual length (-1 because the indices start at 0 and another -1 because
                // we substract the ids of the preceeding mutations range
                idAct = idAct - (pos == 1 ? 0 : cumSumRange[i - 1]); // hier stimmt was nicht
            }
        }

        return mutPos;
    }

    mutVector specIdsToMutPos(const utils::ids specIds, const constants::Constants& params)
    {
        // collect the mutated position with the respective mutation symbol
        mutVector mutPos;
        for (int i = 0; i < params.NCHUNKS; i++)
        {
            int offset = i * params.chunkL;
            mutVector mutPos_chunk = specIdToMutPos(specIds[i], params, offset);
            mutPos.reserve(mutPos.size() + mutPos_chunk.size());
            mutPos.insert(mutPos.end(), mutPos_chunk.begin(), mutPos_chunk.end());
        }

        return mutPos;
    }

    utils::id mutPosToSpecId(const mutVector& mutPos, const constants::Constants& params)
    {
        unsigned numMut = mutPos.size();
        // id for 0 mutations is 1
        utils::id specId = 1;
        if (numMut > 0)
        {
            // add the ids for the sequences with less mutations
            specId = params.NMUT_RANGE[numMut - 1];
            mutVector mutPos_new = mutPos;
            // get the indices for the individual length segments for each position
            if (numMut > 1)
            {
                std::adjacent_difference(mutPos.begin(), mutPos.end(), mutPos_new.begin(),
                                         [](const Mutation& x, const Mutation& y) {
                                             return Mutation(x.getPosition() - y.getPosition(), x.getSymbol());
                                         });
            }
            unsigned Lact = params.chunkL;
            unsigned numMutAct = numMut;

            for (auto const& mutation : std::as_const(mutPos_new))
            {
                if (numMutAct == 1)
                {
                    for (int i = 1; i < mutation.getPosition(); ++i)
                    {
                        specId += pow(params.Q - 1, numMut);
                    }
                    specId += mutation.getSymbol() + 1;
                }
                else
                {
                    // if the position of the actual mutations is != 1, the id is depending on the next position
                    for (int i = 1; i < mutation.getPosition(); ++i)
                    {
                        if (numMut > numMutAct)
                        {
                            specId += pow(params.Q - 1, numMut - numMutAct) * (params.Q - 1) *
                                      std::roundl(utils::nChoosek(Lact - i, numMutAct - 1)) *
                                      pow(params.Q - 1, numMutAct - 1);
                        }
                        else
                            specId += (params.Q - 1) * std::roundl(utils::nChoosek(Lact - i, numMutAct - 1)) *
                                      pow(params.Q - 1, numMutAct - 1);
                    }
                    // if position ==1, and for symbol remainder if position>1
                    specId += (mutation.getSymbol()) * pow(params.Q - 1, numMutAct - 1);
                    Lact -= mutation.getPosition();
                    --numMutAct;
                }
            }
        }

        return specId;
    }

    utils::ids mutPosToSpecIds(const mutVector& mutPos, const constants::Constants& params)
    {
        // collect the mutated position with the respective mutation symbol
        utils::ids specIds(params.NCHUNKS, 1);
        mutVector mutPos_new = mutPos;

        for (int i = 0; i < params.NCHUNKS; ++i)
        {
            mutVector mutPos_chunk;

            int posMax = (i + 1) * params.chunkL;
            while ((mutPos_new.size() > 0) && (mutPos_new[0].getPosition() < posMax))
            {
                mutPos_chunk.emplace_back(mutPos_new[0]);
                mutPos_new.erase(mutPos_new.begin());
            }
            specIds[i] = mutPosToSpecId(mutPos_chunk, params);
        }
        return specIds;
    }

    unsigned getNumberOfMutationsById(const utils::id specId, const constants::Constants& params)
    {
        // gives the index where the content is still lower than the given id
        auto low_it = std::lower_bound(std::begin(params.NMUT_RANGE), std::end(params.NMUT_RANGE), specId);

        // will incorrectly return numMut = 2*MAX_MUT+1 for all specId's > params.NMUT_RANGE[-1]
        if (low_it == std::end(params.NMUT_RANGE))
        {
            std::cerr << "specId " << specId << "is out of range ";
            for (auto i : params.NMUT_RANGE)
                std::cerr << i << ' ';
            std::cerr << std::endl;
            throw std::out_of_range("specId is out of range.");
        }
        return (low_it - std::begin(params.NMUT_RANGE));
    }

    std::vector<unsigned int> getNumberOfMutationsByIds(const utils::ids specIds, const constants::Constants& params)
    {
        std::vector<unsigned int> numMuts(params.NCHUNKS, 0);

        for (int i = 0; i < params.NCHUNKS; i++)
        {
            numMuts[i] = getNumberOfMutationsById(specIds[i], params);
        }
        return numMuts;
    }

    void countErrors(std::vector<std::set<Mutation>> errors, species_map& species_vec,
                     std::valarray<unsigned int>& S_bound, std::valarray<unsigned int>& S_unbound,
                     count::counter_collection& counter)
    {
        int errIdx = 0;
        int specIdx = 0;
        for (auto it = species_vec.begin(); it != species_vec.end(); ++it)
        {
            // species occurs in bound pool
            if (S_bound[specIdx] > 0)
            {
                mutVector mutatedPositions = it->second.getMutatedPositions();

                // loop over #occurences in bound pool
                for (int i = 0; i < S_bound[specIdx]; ++i)
                {
                    auto uniquePositions = errors[errIdx];
                    countErrors_oneSeq(uniquePositions, mutatedPositions, counter.counter_bound_1d,
                                       counter.counter_bound_2d);
                    ++errIdx;
                }
            }

            ++specIdx;
        }

        // reset specIdx and loop through species_vec one more time for unbound pool
        specIdx = 0;
        for (auto it = species_vec.begin(); it != species_vec.end(); ++it)
        {
            // species occurs in undbound pool
            if (S_unbound[specIdx] > 0)
            {
                mutVector mutatedPositions = it->second.getMutatedPositions();

                // loop over #occurences in unbound pool
                for (int i = 0; i != S_unbound[specIdx]; ++i)
                {
                    std::set<Mutation> uniquePositions = errors[errIdx];
                    countErrors_oneSeq(uniquePositions, mutatedPositions, counter.counter_unbound_1d,
                                       counter.counter_unbound_2d);
                    ++errIdx;
                }
            }
            ++specIdx;
        }
    }

    void countErrors_oneSeq(std::set<Mutation> uniquePositions, const mutVector& mutatedPositions,
                            count::counter_1& counter_1d, count::counter_2& counter_2d)
    {
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

    count::counter_collection countMutations(const std::valarray<unsigned int>& SBound,
                                             const std::valarray<unsigned int>& SUnbound, const species_map& spec_map,
                                             const constants::Constants& params)
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
            ++specIdx;
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
            utils::ids ids = mutPosToSpecIds(combined, *combined_params);

            // create new object if not yet present (return value gives iterator and flag if insertion happened)
            // the id is the key for the map, and also the parameter for the constructor for the species class
            auto currentEntry = species_map.try_emplace(ids, ids, *combined_params);
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
