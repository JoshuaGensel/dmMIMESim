#include "Constants.hpp"

#include <gtest/gtest.h>
#include <numeric>

class ConstantsTest : public testing::Test
{
  protected:
    void SetUp()
    {
        length = 100;
        q = 4;
        m = 100;
        s = 0;
        p_mut = 0.01;
        p_error = p_mut / 10;
        p_effect = 0.5;
        p_epistasis = 0.75;
        B_tot = 2.0;
        max_mut = 4;
        epi_restrict = 2;

        params = new constants::Constants(length, q, m, p_mut, p_error, p_effect, p_epistasis, s, B_tot, max_mut,
                                          epi_restrict, std::filesystem::path());
    }

    // Some expensive resource shared by all tests.
    unsigned int length;
    unsigned int q;
    unsigned int m;
    unsigned int s;
    double p_mut;
    double p_error;
    double p_effect;
    double p_epistasis;
    double B_tot;
    unsigned max_mut;
    unsigned epi_restrict;

    constants::Constants* params;
};

TEST_F(ConstantsTest, ConstructFromValues)
{
    constants::Constants cons_val = constants::Constants(length, q, m, p_mut, p_error, p_effect, p_epistasis, s, B_tot,
                                                         max_mut, epi_restrict, std::filesystem::path());
    EXPECT_EQ(cons_val.L, length);
    EXPECT_EQ(cons_val.chunkL, 25);
    EXPECT_EQ(cons_val.NCHUNKS, length / 25);
    EXPECT_EQ(cons_val.Q, q);
    EXPECT_EQ(cons_val.M, m);
    EXPECT_EQ(cons_val.SEED, s);
    EXPECT_EQ(cons_val.P_MUT, p_mut);
    EXPECT_EQ(cons_val.P_EFFECT, p_effect);
    EXPECT_EQ(cons_val.P_EPISTASIS, p_epistasis);
    EXPECT_EQ(cons_val.MAX_MUT, 4);
    EXPECT_EQ(cons_val.EPI_RESTRICT, 2);
    EXPECT_TRUE(cons_val.OUTPUT_DIR.empty());
}

TEST_F(ConstantsTest, NoMUTMAX)
{
    m = 100;
    p_mut = 0.1;
    length = 100;
    // should make for max_mut = 14;
    constants::Constants cons_val = constants::Constants(length, q, m, p_mut, p_error, p_effect, p_epistasis, s, B_tot,
                                                         epi_restrict, std::filesystem::path());

    EXPECT_EQ(cons_val.MAX_MUT, 14);
}

TEST_F(ConstantsTest, ConstructFromPath)
{
    constants::Constants cons_path = constants::readParameters(std::filesystem::path("./test/data/Constants_test"));
    EXPECT_EQ(cons_path.L, 25);
    EXPECT_EQ(cons_path.Q, 4);
}

TEST_F(ConstantsTest, InvalidValues)
{
    const unsigned int invalid_length = 30; // not a multiplative of 25
    unsigned invalid_epi_restrict = 3;

    EXPECT_THROW(
        {
            try
            {
                constants::Constants(invalid_length, q, m, p_mut, p_error, p_effect, p_epistasis, s, B_tot, max_mut,
                                     epi_restrict, std::filesystem::path());
            }
            catch (const std::runtime_error& e)
            {
                EXPECT_STREQ("L should be multiplative of chunkL=25.", e.what());
                throw;
            }
        },
        std::runtime_error);

    EXPECT_THROW(
        {
            try
            {
                constants::Constants(length, q, m, p_mut, p_error, p_effect, p_epistasis, s, B_tot, max_mut,
                                     invalid_epi_restrict, std::filesystem::path());
            }
            catch (const std::runtime_error& e)
            {
                EXPECT_STREQ("epi_restrict should be 0, 1 or 2.", e.what());
                throw;
            }
        },
        std::runtime_error);
}

TEST_F(ConstantsTest, NegativeMAXMUTVal)
{
    const int neg_max_mut = -1;

    constants::Constants no_mutmax_params = constants::Constants(length, q, m, p_mut, p_error, p_effect, p_epistasis, s,
                                                                 B_tot, epi_restrict, std::filesystem::path());

    constants::Constants neg_mutmax_params =
        constants::Constants(length, q, m, p_mut, p_error, p_effect, p_epistasis, s, B_tot, neg_max_mut, epi_restrict,
                             std::filesystem::path());

    EXPECT_EQ(no_mutmax_params.MAX_MUT, neg_mutmax_params.MAX_MUT);
}

TEST_F(ConstantsTest, SetP_NMut)
{
    std::vector<long double> pNMut = params->P_NMUT;
    EXPECT_EQ(std::accumulate(pNMut.begin(), pNMut.end(), double(0)), 1);
}
