#include "Constants.hpp"

#include <gtest/gtest.h>
#include <numeric>

class ConstantsTest : public testing::Test
{
  protected:
    void SetUp()
    {
        const unsigned int length = 50;
        const unsigned int q = 2;
        const unsigned int m = 12000;
        const unsigned int s = 0;
        const double p_mut = 0.01;
        const double p_error = p_mut / 10;
        const double p_effect = 0.5;
        const double p_epistasis = 0.75;
        const double B_tot = 2.0;
        const bool epi_mut_excl = false;

        params = new constants::Constants(length, q, m, p_mut, p_error, p_effect, p_epistasis, s, B_tot, epi_mut_excl,
                                          std::filesystem::path());
    }

    // Some expensive resource shared by all tests.
    constants::Constants* params;
};

TEST_F(ConstantsTest, ConstructFromValues)
{
    const unsigned int length = 12;
    const unsigned int q = 4;
    const unsigned int m = 100;
    const unsigned int s = 356;
    const double p_mut = 0.1;
    const double p_error = p_mut / 10;
    const double p_effect = 0.5;
    const double p_epistasis = 0.3;
    const double B_tot = 2.0;
    const bool epi_mut_excl = true;
    constants::Constants cons_val = constants::Constants(length, q, m, p_mut, p_error, p_effect, p_epistasis, s, B_tot,
                                                         epi_mut_excl, std::filesystem::path());
    EXPECT_EQ(cons_val.L, length);
    EXPECT_EQ(cons_val.Q, q);
    EXPECT_EQ(cons_val.M, m);
    EXPECT_EQ(cons_val.SEED, s);
    EXPECT_EQ(cons_val.P_MUT, p_mut);
    EXPECT_EQ(cons_val.P_EFFECT, p_effect);
    EXPECT_EQ(cons_val.P_EPISTASIS, p_epistasis);
    EXPECT_TRUE(cons_val.OUTPUT_DIR.empty());
    EXPECT_TRUE(cons_val.EPIMUTEXCL);
}

TEST_F(ConstantsTest, ConstructFromPath)
{
    constants::Constants cons_path = constants::readParameters(std::filesystem::path("./testParameters"));
    EXPECT_EQ(cons_path.L, 15);
    EXPECT_EQ(cons_path.Q, 4);
}

TEST_F(ConstantsTest, ComputeMaxMut)
{
    EXPECT_EQ(params->MAX_MUT, 4);

    constants::Constants* big_params = new constants::Constants(100, 2, 1200000, 0.01, 0.001, 0.5, 0.75, 0, 2.0, false,
                                                                std::filesystem::temp_directory_path());
    EXPECT_EQ(big_params->MAX_MUT, 8);

    constants::Constants* bigger_params = new constants::Constants(500, 2, 120000, 0.001, 0.0001, 0.5, 0.75, 0, 2.0,
                                                                   false, std::filesystem::temp_directory_path());
    EXPECT_EQ(bigger_params->MAX_MUT, 5);
}

TEST_F(ConstantsTest, SetNMutRange)
{
    const std::vector<utils::id> nMutRange = params->NMUT_RANGE;
    EXPECT_EQ(nMutRange.size(), 2 * 4 + 1);

    const std::vector<utils::id> nMutRange_true{1, 51, 1276, 20876, 251176, 2369936, 18260636, 118145036, 655023686};
    EXPECT_EQ(nMutRange, nMutRange_true);
}

TEST_F(ConstantsTest, SetP_NMut)
{
    std::vector<long double> pNMut = params->P_NMUT;
    EXPECT_EQ(std::accumulate(pNMut.begin(), pNMut.end(), double(0)), 1);
    EXPECT_GT(pNMut[0], pNMut[1]);
    EXPECT_GT(pNMut[1], pNMut[2]);
    EXPECT_GT(pNMut[2], pNMut[3]);
    EXPECT_GT(pNMut[3], pNMut[4]);
}
