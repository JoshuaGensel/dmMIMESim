//
// Created by Smith, Maureen on 31.05.18.
//

#include "Constants.hpp"
#include "FunctionalSequence.hpp"
#include "catch.hpp"

#include <iostream>

TEST_CASE("Testing FunctionalSequence Class")
{
    constants::Constants& cons = constants::Constants::get_instance();

    SECTION("test function getMatrixIndex which converts 2 indicices into the corresponding index  of a vectoral "
            "symmetric matrix")
    {

        FunctionalSequence& effect = FunctionalSequence::get_instance();
        unsigned int i;
        unsigned int j;
        unsigned int idx = 0;

        auto epistasis = effect.getEpistasis();

        for (i = 1; i < cons.L; ++i)
        {
            for (j = i + 1; j <= cons.L; ++j)
            {
                REQUIRE(effect.getEpistasis(i, j) == epistasis[idx]);
                ++idx;
            }
        }
    }

    SECTION("Test if its singleton")
    {
        FunctionalSequence& effect = FunctionalSequence::get_instance();
        // Kds are picked randomly. If it would not be a singleton, the second drawn effects would be different
        FunctionalSequence& effect2 = FunctionalSequence::get_instance();
        for (int i = 1; i <= cons.L; ++i)
        {
            REQUIRE(effect.getKd(i) == effect2.getKd(i));
        }
    }

    SECTION("Test total species Kd computation")
    {
        // TODO Test schreiben
    }
}