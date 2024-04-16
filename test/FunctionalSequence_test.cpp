#include "Constants.hpp"
#include "FunctionalSequence.hpp"
#include "Generator.hpp"
#include "Species.hpp"

#include <gtest/gtest.h>
#include <gtest/gtest_prod.h>

class FunctionalSequenceTest : public ::testing::Test
{
  protected:
    void SetUp()
    {
        unsigned int length = 25;
        unsigned int q = 4;
        unsigned int m = 100;
        double p_mut = 0.01;
        double p_error = 0.01;
        double p_effect = 0.5;
        double p_epistasis = 0.5;
        unsigned int seed = 0;
        double B_tot = 2.0;
        // unsigned max_mut = -1;
        unsigned epi_restrict = 2;
        fs::path outputDir = std::filesystem::temp_directory_path();

        params = new constants::Constants(length, q, m, p_mut, p_error, p_effect, p_epistasis, seed, B_tot,
                                          epi_restrict, outputDir);
        Generator::create_instance(params->SEED);
    }
    void TearDown()
    {
        Generator::release_instance();
    }

    constants::Constants* params;
};

TEST_F(FunctionalSequenceTest, NoCreation)
{
    EXPECT_THROW(
        {
            try
            {
                FunctionalSequence::get_instance();
            }
            catch (const std::runtime_error& e)
            {
                EXPECT_STREQ("FunctionalSequence not correctly instantiated (get before create or after release).",
                             e.what());
                throw;
            }
        },
        std::runtime_error);
}

TEST_F(FunctionalSequenceTest, CreationAndRelease)
{

    FunctionalSequence::create_instance(*params);
    FunctionalSequence::release_instance();
}

TEST_F(FunctionalSequenceTest, AllEpistasis)
{
    auto all_epi_params =
        new constants::Constants(params->L, params->Q, params->M, params->P_MUT, params->P_ERR, params->P_EFFECT, 1.0,
                                 params->SEED, params->BTOT, params->EPI_RESTRICT, params->OUTPUT_DIR);
    FunctionalSequence::create_instance(*all_epi_params);
    FunctionalSequence::release_instance();
}

TEST_F(FunctionalSequenceTest, NoEpistasis)
{
    auto no_epi_params =
        new constants::Constants(params->L, params->Q, params->M, params->P_MUT, params->P_ERR, params->P_EFFECT, 0.0,
                                 params->SEED, params->BTOT, params->EPI_RESTRICT, params->OUTPUT_DIR);
    FunctionalSequence::create_instance(*no_epi_params);
    FunctionalSequence::release_instance();
}

class FunctionalSequenceTestPrivate : public ::testing::Test
{
  public:
    FunctionalSequence* effects;
    fs::path inputPath;

  protected:
    void SetUp()
    {
        inputPath = "./test/data/FunctionalSequence_test";
        const constants::Constants cons = constants::readParameters(inputPath);
        Generator::create_instance(0);
        effects = FunctionalSequence::create_instance(cons);
    }
    void TearDown()
    {
        Generator::release_instance();
        FunctionalSequence::release_instance();
    }
};

TEST_F(FunctionalSequenceTestPrivate, ReadKdValues)
{
    std::vector<double> kds = effects->readKdValues(inputPath);

    for (float i = 0; i < kds.size(); ++i)
        EXPECT_EQ(i + 1.1, kds[i]);
}

TEST_F(FunctionalSequenceTestPrivate, ReadEpistasis)
{
    std::vector<double> epi = effects->readEpistasis(inputPath);

    for (float i = 0; i < epi.size(); ++i)
        EXPECT_EQ(i + 1.1, epi[i]);
}

TEST_F(FunctionalSequenceTestPrivate, GetMutationFromVectorIndex)
{
    auto mut = effects->getMutationFromVectorIndex(49);
    EXPECT_EQ(mut.getPosition(), 25);
    EXPECT_EQ(mut.getSymbol(), 1);
}

TEST_F(FunctionalSequenceTestPrivate, PairIndex)
{
    EXPECT_EQ(effects->pairIndex.back(), (effects->params.L * (effects->params.L - 1)) / 2);
    EXPECT_EQ(effects->pairIndex.front(), 0);
    EXPECT_EQ(effects->pairIndex[2], (0 + (effects->params.L - 1) + (effects->params.L - 2)));
}

TEST_F(FunctionalSequenceTestPrivate, GetPairIndex)
{
    EXPECT_EQ(effects->getPairIndex(1, 2, 0, 0), 0);
    EXPECT_EQ(effects->getPairIndex(1, 2, 0, 1), 1);
    EXPECT_EQ(effects->getPairIndex(2, 3, 0, 0), (effects->params.L - 1) * std::pow(effects->params.Q - 1, 2));
    EXPECT_EQ(effects->getPairIndex(2, 3, 1, 1),
              (effects->params.L - 1) * std::pow(effects->params.Q - 1, 2) + (effects->params.Q - 1) * 1 + 1);
    EXPECT_EQ(effects->getPairIndex((effects->params.L - 1), effects->params.L, (effects->params.Q - 1) - 1,
                                    (effects->params.Q - 1) - 1),
              effects->params.PWVal - 1);
}

TEST_F(FunctionalSequenceTestPrivate, DrawEpistasis_restricted)
{
    int counter = 0;
    for (auto epiVal : effects->drawEpistasis_restricted())
    {
        if (epiVal != effects->params.NO_EPISTASIS)
            ++counter;
    }

    EXPECT_EQ(effects->params.P_EPISTASIS, 1.0);
    EXPECT_EQ(counter, std::floor(effects->params.L / 2));
}

TEST_F(FunctionalSequenceTestPrivate, DrawEpistasis_semi_restricted)
{
    int counter = 0;
    for (auto epiVal : effects->drawEpistasis_semi_restricted())
    {
        if (epiVal != effects->params.NO_EPISTASIS)
            ++counter;
    }

    EXPECT_EQ(effects->params.P_EPISTASIS, 1.0);
    EXPECT_EQ(counter, effects->params.SVal / 2);
}

TEST_F(FunctionalSequenceTestPrivate, DrawEpistasis_unrestricted)
{
    int counter = 0;
    for (auto epiVal : effects->drawEpistasis_unrestricted())
    {
        if (epiVal != effects->params.NO_EPISTASIS)
            ++counter;
    }

    EXPECT_EQ(effects->params.P_EPISTASIS, 1.0);
    EXPECT_EQ(counter, effects->params.PWVal);
}
