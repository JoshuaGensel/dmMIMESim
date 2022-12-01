//
//  test_species.cpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 24.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#include "Constants.hpp"
#include "FunctionalSequence.hpp"
#include "Generator.hpp"
#include "Species.hpp"

#include <gtest/gtest.h>
#include <iostream>
#include <valarray>

class SpeciesTest : public testing::Test
{
  protected:
    // Per-test-suite set-up.
    // Called before the first test in this test suite.
    void SetUp()
    {
        unsigned int length = 50;
        unsigned int q = 4;
        unsigned int m = 1000;
        double p_mut = 0.1;
        double p_error = 0.01;
        double p_effect = 0.5;
        double p_epistasis = 0.75;
        unsigned int seed = 0;
        double B_tot = 2.0;
        unsigned max_mut = 4;
        unsigned epi_restrict = 2;
        fs::path outputDir = std::filesystem::temp_directory_path();

        params = new constants::Constants(length, q, m, p_mut, p_error, p_effect, p_epistasis, seed, B_tot, max_mut,
                                          epi_restrict, outputDir);
        Generator::create_instance(params->SEED);
        effects = FunctionalSequence::create_instance(*params);
    }

    // Per-test-suite tear-down.
    // Called after the last test in this test suite.
    void TearDown()
    {
        Generator::release_instance();
        FunctionalSequence::release_instance();
    }

    const constants::Constants* params;
    FunctionalSequence* effects;
};

TEST_F(SpeciesTest, specIdToMutPos)
{

    unsigned int length = 25;
    unsigned int q = 4;
    unsigned int m = 1000;
    double p_mut = 0.2;
    double p_error = 0.01;
    double p_effect = 0.5;
    double p_epistasis = 0.75;
    unsigned int seed = 0;
    double B_tot = 2.0;
    unsigned epi_restrict = 2;
    fs::path outputDir = std::filesystem::temp_directory_path();

    constants::Constants* small_params = new constants::Constants(length, q, m, p_mut, p_error, p_effect, p_epistasis,
                                                                  seed, B_tot, epi_restrict, outputDir);

    auto chunkL = small_params->chunkL;
    auto numSymbols = q - 1;

    // Test id 1 = no  mutation
    utils::id id = 1;
    auto mutations = species::specIdToMutPos(id, *small_params);
    EXPECT_EQ(mutations.size(), 0);
    EXPECT_TRUE(mutations.empty());

    // test Id range for one mutation
    for (auto pos = 1; pos <= chunkL; ++pos)
    {
        for (auto symb = 0; symb < numSymbols; ++symb)
        {
            ++id;
            mutations = species::specIdToMutPos(id, *small_params);
            EXPECT_EQ(mutations.size(), 1);
            EXPECT_EQ(mutations.back().getPosition(), pos);
            EXPECT_EQ(mutations.back().getSymbol(), symb);
        }
    }

    // test Id range for two mutations
    for (auto pos1 = 1; pos1 <= chunkL - 1; ++pos1)
    {
        for (auto pos2 = pos1 + 1; pos2 <= chunkL; ++pos2)
        {
            for (auto symb1 = 0; symb1 < numSymbols; ++symb1)
            {
                for (auto symb2 = 0; symb2 < numSymbols; ++symb2)
                {
                    ++id;
                    mutations = species::specIdToMutPos(id, *small_params);
                    EXPECT_EQ(mutations.size(), 2);
                    EXPECT_EQ(mutations.at(0).getPosition(), pos1);
                    EXPECT_EQ(mutations.at(1).getPosition(), pos2);
                    EXPECT_EQ(mutations.at(0).getSymbol(), symb1);
                    EXPECT_EQ(mutations.at(1).getSymbol(), symb2);
                }
            }
        }
    }

    // test Id range for three mutations
    for (auto pos1 = 1; pos1 <= chunkL - 2; ++pos1)
    {
        for (auto pos2 = pos1 + 1; pos2 <= chunkL - 1; ++pos2)
        {
            for (auto pos3 = pos2 + 1; pos3 <= chunkL; ++pos3)
            {
                for (auto symb1 = 0; symb1 < numSymbols; ++symb1)
                {
                    for (auto symb2 = 0; symb2 < numSymbols; ++symb2)
                    {
                        for (auto symb3 = 0; symb3 < numSymbols; ++symb3)
                        {
                            ++id;
                            mutations = species::specIdToMutPos(id, *small_params);
                            EXPECT_EQ(mutations.size(), 3);
                            EXPECT_EQ(mutations.at(0).getPosition(), pos1);
                            EXPECT_EQ(mutations.at(1).getPosition(), pos2);
                            EXPECT_EQ(mutations.at(2).getPosition(), pos3);
                            EXPECT_EQ(mutations.at(0).getSymbol(), symb1);
                            EXPECT_EQ(mutations.at(1).getSymbol(), symb2);
                            EXPECT_EQ(mutations.at(2).getSymbol(), symb3);
                        }
                    }
                }
            }
        }
    }

    // test Id range for four mutations
    for (auto pos1 = 1; pos1 <= length - 3; ++pos1)
    {
        for (auto pos2 = pos1 + 1; pos2 <= length - 2; ++pos2)
        {
            for (auto pos3 = pos2 + 1; pos3 <= length - 1; ++pos3)
            {
                for (auto pos4 = pos3 + 1; pos4 <= length; ++pos4)
                {
                    for (auto symb1 = 0; symb1 < numSymbols; ++symb1)
                    {
                        for (auto symb2 = 0; symb2 < numSymbols; ++symb2)
                        {
                            for (auto symb3 = 0; symb3 < numSymbols; ++symb3)
                            {
                                for (auto symb4 = 0; symb4 < numSymbols; ++symb4)
                                {
                                    ++id;
                                    mutations = species::specIdToMutPos(id, *small_params);
                                    EXPECT_EQ(mutations.size(), 4);
                                    EXPECT_EQ(mutations.at(0).getPosition(), pos1);
                                    EXPECT_EQ(mutations.at(1).getPosition(), pos2);
                                    EXPECT_EQ(mutations.at(2).getPosition(), pos3);
                                    EXPECT_EQ(mutations.at(3).getPosition(), pos4);
                                    EXPECT_EQ(mutations.at(0).getSymbol(), symb1);
                                    EXPECT_EQ(mutations.at(1).getSymbol(), symb2);
                                    EXPECT_EQ(mutations.at(2).getSymbol(), symb3);
                                    EXPECT_EQ(mutations.at(3).getSymbol(), symb4);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

TEST_F(SpeciesTest, specIds)
{
    utils::ids ids{1, 200};
    species::Species spec = species::Species(ids, *params);
    EXPECT_EQ(spec.getSpecIds().size(), 2);

    utils::ids invalid(3, 1);

    EXPECT_THROW(
        {
            try
            {
                species::Species(invalid, *params);
            }
            catch (const std::runtime_error& e)
            {
                // and this tests that it has the correct message
                EXPECT_STREQ("Mismatch between specIds and params.", e.what());
                throw;
            }
        },
        std::runtime_error);
}

TEST_F(SpeciesTest, getNumberOfMutationsById)
{
    // upper bound ID per numMutations
    // MAXMUT is 4, there should be room for IDs up to 8 mutations

    EXPECT_EQ(species::getNumberOfMutationsById(1, *params), 0);
    EXPECT_EQ(
        species::getNumberOfMutationsById(utils::nChoosek(params->chunkL, 1) * std::pow(params->Q - 1, 1) + 2, *params),
        1);
    EXPECT_EQ(
        species::getNumberOfMutationsById(utils::nChoosek(params->chunkL, 2) * std::pow(params->Q - 1, 2) + 2, *params),
        2);
    EXPECT_EQ(
        species::getNumberOfMutationsById(utils::nChoosek(params->chunkL, 3) * std::pow(params->Q - 1, 3) + 2, *params),
        3);
    EXPECT_EQ(
        species::getNumberOfMutationsById(utils::nChoosek(params->chunkL, 4) * std::pow(params->Q - 1, 4) + 2, *params),
        4);
    EXPECT_EQ(
        species::getNumberOfMutationsById(utils::nChoosek(params->chunkL, 5) * std::pow(params->Q - 1, 5) + 2, *params),
        5);
    EXPECT_EQ(
        species::getNumberOfMutationsById(utils::nChoosek(params->chunkL, 6) * std::pow(params->Q - 1, 6) + 2, *params),
        6);
    EXPECT_EQ(
        species::getNumberOfMutationsById(utils::nChoosek(params->chunkL, 7) * std::pow(params->Q - 1, 7) + 2, *params),
        7);
    EXPECT_EQ(
        species::getNumberOfMutationsById(utils::nChoosek(params->chunkL, 8) * std::pow(params->Q - 1, 8) + 2, *params),
        8);

    try
    {
        // utils::ids invalid_ids({1, llround(utils::nChoosek(params->chunkL, 9) * (params->Q - 1) + 2)});
        auto invalid =
            species::Species({1, llround(utils::nChoosek(params->chunkL, 9) * std::pow(params->Q - 1, 9) + 2)}, *params)
                .getNumMuts();
        FAIL() << "Expected specId is out of range.";
    }
    catch (std::out_of_range const& err)
    {
        EXPECT_EQ(err.what(), std::string("specId is out of range."));
    }
    catch (...)
    {
        FAIL() << "Expected specId is out of range.";
    }
}

TEST_F(SpeciesTest, DrawSpeciesId)
{

    species::species_map speciesIDs = species::drawSpeciesIds(*params);
    std::valarray<int> numMutCounts(0, params->MAX_MUT + 1);
    for (auto it = speciesIDs.begin(); it != speciesIDs.end(); ++it)
    {
        // Test if the id range is valid
        for (auto chunkId : it->first)
        {
            EXPECT_GT(chunkId, 0);
            EXPECT_LE(chunkId, params->NMUT_RANGE[params->MAX_MUT]);
        }
        numMutCounts[it->second.getNumMut_total()] += it->second.getCount();
    }
    // Test that all number of mutations were drawn
    for (auto& count : numMutCounts)
    {
        EXPECT_GT(count, 0);
    }
    EXPECT_EQ(params->M, numMutCounts.sum());
}

TEST_F(SpeciesTest, DrawSpeciesId_big)
{

    // this is a big param setting that should still work
    // when increasng L or Q, one should reduce M!

    unsigned int length = 500;
    unsigned int q = 4;
    unsigned int m = 1000000;
    double p_mut = 0.01;
    double p_error = 0.001;
    double p_effect = 0.5;
    double p_epistasis = 0.75;
    unsigned int seed = 0;
    double B_tot = 2.0;
    unsigned epi_restrict = 0;
    fs::path outputDir = std::filesystem::temp_directory_path();

    constants::Constants* big_params = new constants::Constants(length, q, m, p_mut, p_error, p_effect, p_epistasis,
                                                                seed, B_tot, epi_restrict, outputDir);

    FunctionalSequence::release_instance();
    auto effects = FunctionalSequence::create_instance(*big_params);

    species::species_map speciesIDs = species::drawSpeciesIds(*big_params);
    std::valarray<int> numMutCounts(0, big_params->MAX_MUT + 1);

    for (auto it = speciesIDs.begin(); it != speciesIDs.end(); ++it)
    {
        // Test if the id range is valid
        for (auto chunkId : it->first)
        {

            EXPECT_GT(chunkId, 0);
            EXPECT_LE(chunkId, big_params->NMUT_RANGE[big_params->MAX_MUT]);
        }
        numMutCounts[it->second.getNumMut_total()] += it->second.getCount();
    }

    // Test that all number of mutations were drawn
    for (auto& count : numMutCounts)
    {
        EXPECT_GT(count, 0);
    }

    EXPECT_EQ(big_params->M, numMutCounts.sum());
}

TEST_F(SpeciesTest, mutPosToSpecId)
{
    for (utils::id id = 1; id <= params->NMUT_RANGE[5]; ++id)
    {
        auto mutPos = species::specIdToMutPos(id, *params);

        utils::id computedId = species::mutPosToSpecId(mutPos, *params);

        EXPECT_EQ(id, computedId);
    }
}

TEST_F(SpeciesTest, DrawErrors)
{
    auto errors_lib1 = species::drawErrors(*params);
    auto errors_lib2 = species::drawErrors(*params);

    EXPECT_EQ(errors_lib1.size(), errors_lib2.size());

    // errors should not be the same
    EXPECT_FALSE(std::equal(errors_lib1.begin(), errors_lib1.end(), errors_lib2.begin()));
}

TEST_F(SpeciesTest, DrawError_2)
{
    double nSeqs = 10e6;
    std::vector<int> numErrors(nSeqs, 0);
    std::vector<int> posCount(params->L, 0);

    for (int i = 0; i < nSeqs; ++i)
    {
        std::set<Mutation> errors = species::drawError_2(*params);
        numErrors[i] = errors.size();

        for (auto err : errors)
        {
            posCount[err.getPosition() - 1] += 1;
        }
    }

    // average number of errors should reflect the error probability
    EXPECT_NEAR(params->P_ERR,
                static_cast<double>(std::accumulate(numErrors.begin(), numErrors.end(), 0)) / nSeqs / params->L,
                2 * params->P_ERR);

    // errors should be approximately uniform distributed over the positions
    int tolerance = std::round(params->P_ERR * nSeqs);
    for (int i = 0; i < params->L - 1; ++i)
    {
        for (int j = i + 1; j < params->L; ++j)
        {
            EXPECT_NEAR(posCount[i], posCount[j], 2 * tolerance);
        }
    }
}

// specIdsToMutPos

// getNumberOfMutationsByIds

// getSequence

// computeSpeciesKd

// drawSpeciesIds

// readFromFile

// # --- species

// speciesMapSum

// writeSpeciesToFile

// writeSequencesToFile

// drawError_2

// drawError

// specIdsToMutPos

// mutPosToSpecIds

// getNumberOfMutationsById

// getNumberOfMutationsByIds

// countErrors

// countErrors_oneSeq

// countMutations

// combineSpecies

// combineMutations
