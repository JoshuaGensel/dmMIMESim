//
//  Species.hpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 03.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#ifndef Species_hpp
#define Species_hpp

#include "Constants.hpp"
#include "Count.hpp"
#include "Mutation.hpp"

#include <array>
#include <map>
#include <random>
#include <set>
#include <stdio.h>
#include <unordered_map>
#include <vector>

namespace species
{
    // TODO weg nach test
    // typedef std::array<unsigned int, Constants::MAX_MUT> mutArr;
    // typedef std::vector<unsigned int> mutVector;
    // TODO ist ja eigentlich generell, nicht nur Species. Irgedwo anders hinschieben?
    // TODO stattdessen einfach Mutation?
    // typedef std::pair<unsigned, unsigned> posMutPair;
    // contain each mutated position with the respective mutation
    typedef std::vector<Mutation> mutVector;
    // typedef std::map<int, std::array<int, 2>> idCountMap;
    // TODO umändern in unordered_map und vergleichen
    typedef std::unordered_map<int, int> idCountMap;

    class Species
    {

      private:
        const unsigned long long specId;
        const constants::Constants& params;
        unsigned int count;
        const unsigned int numMut;
        // mutatated positions need to be in ascending order
        // TODO where to test this, exception?
        const mutVector mutatedPositions;
        // double freq;
        // KD of the the given sequence, adding all single Kds of the mutations and the epistatic effects for pairs
        double kd = 1.0;
        // TODO for later count purposes directly create a read. nee weg damit, nimmt zu viel platz ein
        // const ref::ref_map read;

        /**** after ODE: save bound and unbound fraction information ****/

        unsigned mutCountBound;
        unsigned mutCountUnbound;

        int errorCountBound;
        int errorCountUnbound;

        unsigned int getNumberOfMutationsById();

        mutVector specIdxToMutPos();

      public:
        Species(const unsigned long long id, const constants::Constants& param);

        const unsigned long long getSpecId() const;

        const constants::Constants& getParams() const;

        unsigned int getCount() const;

        const unsigned int getNumMut() const;

        const mutVector& getMutatedPositions() const;

        double getFreq() const;

        double getKd() const;

        void setCount(unsigned int count);

        void incrementCount();

        void decrementCount();

        void computeSpeciesKd();

        // generate nt sequence based on mutated positions and optional errors
        std::valarray<char> getSequence(std::set<Mutation> errors = {});

        // TODO weg
        // const ref::ref_map &getRead() const;
        // ref::ref_map createRead();

        unsigned int getMutCountBound() const;

        void setMutCountBound(unsigned int mutCountBound);

        unsigned int getMutCountUnbound() const;

        void setMutCountUnbound(unsigned int mutCountUnbound);

        int getErrorCountBound() const;

        void setErrorCountBound(int errorCountBound);

        void addErrorCountBound(int errorCountBound);

        int getErrorCountUnbound() const;

        void setErrorCountUnbound(int errorCountUnbound);

        void addErrorCountUnbound(int errorCountUnbound);

        // Fraction of all M sequences
        double getTotalFractionBound();

        double getTotalFractionUnbound();

        // fraction of the one particular species
        double getFractionBound();

        double getFractionUnbound();
    };

    using species_map = std::unordered_map<unsigned long long, Species>;

    species_map drawSpeciesIds(const constants::Constants& params);

    mutVector specIdxToMutPos(const unsigned long long specId, const constants::Constants& params);

    unsigned long long mutPosToSpecIdx(const mutVector& mutPos, const constants::Constants& params);

    unsigned getNumberOfMutationsById(const unsigned long long specId, const constants::Constants& params);

    species_map readFromFile(const std::string& inputPath, utils::SampleID id, const constants::Constants& params);

    unsigned int speciesMapSum(species_map species_vec);

    /**
     *
     * @param out_file
     * @param species_map the drawn Species IDs
     * @param S_pool counts per drawn Species IDs present in pool of interest
     * @param header optional argument. If it is not given, it is assumed that we have species and #occurences
     */
    void writeSpeciesToFile(const std::string& out_file, species_map& spec_map, std::valarray<unsigned int>& S_pool,
                            const std::string& header = "speciesID\tcount\n");

    void writeSequencesToFile(const std::string& out_file, species_map& spec_map,
                              std::vector<std::set<Mutation>> errors, std::valarray<unsigned int> S_pool);

    std::set<Mutation> drawError_2(const constants::Constants& params);

    // draws params.M error Mutation sets
    std::vector<std::set<Mutation>> drawErrors(const constants::Constants& params);

    // first S_bound.size() entries in errors are applied to the bound counters, last S_unbound.size() to the unbound
    // counters
    void countErrors(std::vector<std::set<Mutation>> errors, species_map& species_vec,
                     std::valarray<unsigned int>& S_bound, std::valarray<unsigned int>& S_unbound,
                     count::counter_collection& counter);

    void countErrors_oneSeq(std::set<Mutation> uniquePositions, const mutVector& mutatedPositions,
                            count::counter_1& counter_1d, count::counter_2& counter_2d);

    count::counter_collection countMutations(const std::valarray<unsigned int>& SBound,
                                             const std::valarray<unsigned int>& SUnbound, const species_map& spec_map,
                                             const constants::Constants& params);

    species_map combineSpecies(species_map firstSpec, species_map secondSpec, const constants::Constants& firstParams,
                               const constants::Constants& secondParams);

    mutVector combineMutations(mutVector firstMut, mutVector secondMut);
}
#endif /* Species_hpp */
