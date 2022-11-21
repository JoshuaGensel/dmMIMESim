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
        params = new constants::Constants(50, 2, 12000, 0.01, 0.001, 0.5, 0.75, 0, 2.0, false,
                                          std::filesystem::temp_directory_path());
        Generator::create_instance(params->SEED);
        effects = FunctionalSequence::create_instance(*params);
    }

    // Per-test-suite tear-down.
    // Called after the last test in this test suite.
    void TearDown()
    {
        FunctionalSequence::release_instance();
        Generator::release_instance();
    }

    const constants::Constants* params;
    FunctionalSequence* effects;
};

TEST_F(SpeciesTest, SpecIdxToMutPos)
{
    constants::Constants* small_params = new constants::Constants(10, 2, 12000, 0.01, 0.001, 0.5, 0.75, 0, 2.0, false,
                                                                  std::filesystem::temp_directory_path());

    auto q = small_params->Q;
    auto length = small_params->L;
    auto numSymbols = q - 1;

    // Test id 1 = no  mutation
    unsigned int id = 1;
    auto positions = species::specIdxToMutPos(id, *small_params);
    EXPECT_EQ(positions.size(), 0);
    EXPECT_TRUE(positions.empty());

    // test Id range for one mutation
    for (auto pos = 1; pos <= length; ++pos)
    {
        for (auto symb = 1; symb < q; ++symb)
        {
            ++id;
            positions = species::specIdxToMutPos(id, *small_params);
            EXPECT_EQ(positions.size(), 1);
            EXPECT_EQ(positions.back().getPosition(), pos);
        }
    }

    // test Id range for two mutations
    for (auto pos1 = 1; pos1 <= length - 1; ++pos1)
    {
        for (auto pos2 = pos1 + 1; pos2 <= length; ++pos2)
        {
            for (auto symb1 = 1; symb1 < q; ++symb1)
            {
                // for (auto symb2 = 1; symb2 < q; ++symb2) {
                ++id;
                positions = species::specIdxToMutPos(id, *small_params);
                EXPECT_EQ(positions.size(), 2);
                EXPECT_EQ(positions.at(0).getPosition(), pos1);
                EXPECT_EQ(positions.at(1).getPosition(), pos2);
                //}
            }
        }
    }

    // test Id range for three mutations
    for (auto pos1 = 1; pos1 <= length - 1; ++pos1)
    {
        for (auto pos2 = pos1 + 1; pos2 <= length; ++pos2)
        {
            for (auto pos3 = pos2 + 1; pos3 <= length; ++pos3)
            {
                for (auto symb1 = 1; symb1 < q; ++symb1)
                {
                    // for (auto symb2 = 1; symb2 < q; ++symb2) {
                    ++id;
                    positions = species::specIdxToMutPos(id, *small_params);
                    EXPECT_EQ(positions.size(), 3);
                    EXPECT_EQ(positions.at(0).getPosition(), pos1);
                    EXPECT_EQ(positions.at(1).getPosition(), pos2);
                    EXPECT_EQ(positions.at(2).getPosition(), pos3);
                    //}
                }
            }
        }
    }
    // test Id range for three mutations
    for (auto pos1 = 1; pos1 <= length - 1; ++pos1)
    {
        for (auto pos2 = pos1 + 1; pos2 <= length; ++pos2)
        {
            for (auto pos3 = pos2 + 1; pos3 <= length; ++pos3)
            {
                for (auto pos4 = pos3 + 1; pos4 <= length; ++pos4)
                {
                    for (auto symb1 = 1; symb1 < q; ++symb1)
                    {
                        ++id;
                        positions = species::specIdxToMutPos(id, *small_params);
                        EXPECT_EQ(positions.size(), 4);
                        EXPECT_EQ(positions.at(0).getPosition(), pos1);
                        EXPECT_EQ(positions.at(1).getPosition(), pos2);
                        EXPECT_EQ(positions.at(2).getPosition(), pos3);
                        EXPECT_EQ(positions.at(3).getPosition(), pos4);
                    }
                }
            }
        }
    }
}

TEST_F(SpeciesTest, GetNumberOfMutationsById)
{
    // upper bound ID per numMutations
    // MAXMUT is 4, there should be room for IDs up to 8 mutations
    EXPECT_EQ(species::Species(1, *params).getNumMut(), 0);
    EXPECT_EQ(species::Species(51, *params).getNumMut(), 1);
    EXPECT_EQ(species::Species(1276, *params).getNumMut(), 2);
    EXPECT_EQ(species::Species(20876, *params).getNumMut(), 3);

    EXPECT_EQ(species::Species(251176, *params).getNumMut(), 4);
    EXPECT_EQ(species::Species(2369936, *params).getNumMut(), 5);
    EXPECT_EQ(species::Species(18260636, *params).getNumMut(), 6);

    EXPECT_EQ(species::Species(118145036, *params).getNumMut(), 7);
    EXPECT_EQ(species::Species(655023686, *params).getNumMut(), 8);

    try
    {
        species::Species(3160457386, *params).getNumMut();
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
        EXPECT_GT(it->first, 0);
        EXPECT_LE(it->first, params->NMUT_RANGE[params->MAX_MUT]);
        numMutCounts[it->second.getNumMut()] += it->second.getCount();
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
    constants::Constants* big_params = new constants::Constants(100, 2, 1200000, 0.01, 0.001, 0.5, 0.3, 0, 2.0, false,
                                                                std::filesystem::temp_directory_path());

    FunctionalSequence::release_instance();
    auto effects = FunctionalSequence::create_instance(*big_params);

    species::species_map speciesIDs = species::drawSpeciesIds(*big_params);
    std::valarray<int> numMutCounts(0, big_params->MAX_MUT + 1);

    for (auto it = speciesIDs.begin(); it != speciesIDs.end(); ++it)
    {
        // Test if the id range is valid
        EXPECT_GT(it->first, 0);
        EXPECT_LE(it->first, big_params->NMUT_RANGE[big_params->MAX_MUT]);
        numMutCounts[it->second.getNumMut()] += it->second.getCount();
    }

    // Test that all number of mutations were drawn
    for (auto& count : numMutCounts)
    {
        EXPECT_GT(count, 0);
    }

    EXPECT_EQ(big_params->M, numMutCounts.sum());
}

TEST_F(SpeciesTest, MutPosToSpecIdx)
{
    int max = params->NMUT_RANGE[params->MAX_MUT];
    for (int i = 1; i <= params->NMUT_RANGE[params->MAX_MUT]; ++i)
    {
        unsigned id = i;
        auto mutPos = species::specIdxToMutPos(id, *params);
        unsigned computedId = species::mutPosToSpecIdx(mutPos, *params);
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
