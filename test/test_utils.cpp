//
// Created by Smith, Maureen on 07.06.18.
//

#include "Constants.hpp"
#include "Species.hpp"
#include "Utils.hpp"
#include "catch.hpp"

#include <iostream>
#include <valarray>

TEST_CASE("Testing Utils", "[UtilsTest]")
{

    // TODO test: no output path given
    // TODO test: output path given, no parameter file given -> take
    // TODO test: output path give
    SECTION("Test function read parameter: no outputpath given")
    {
        REQUIRE_THROWS_AS(utils::readParameters(fs::path()), std::invalid_argument);
    }

    SECTION("Test function read parameter: outputpath not given")
    {
        fs::path outputPath = "./results";

        REQUIRE_NOTHROW(utils::readParameters(outputPath));
        REQUIRE(fs::exists(outputPath));
    fs:
        remove(outputPath);
    }

    // TODO: WEG DAMIT wie kann ich für verschiedene Test cases verschieden set ups erstellen so dass sie innerhalb des
    // cases für alle sections gilt, aber nicht für alle cases?
    constants::Constants& cons = constants::Constants::get_instance();

    SECTION("Test  function n choose k (binomial coefficient)") {}
}