//
// Created by Smith, Maureen on 06.06.18.
//

#include "BindingCompetition.hpp"
#include "Constants.hpp"
#include "FunctionalSequence.hpp"
#include "Species.hpp"
#include "catch.hpp"

#include <iostream>

TEST_CASE("Testing Binding Competition")
{
    // TODO: WEG DAMIT wie kann ich für verschiedene Test cases verschieden set ups erstellen so dass sie innerhalb des
    // cases für alle sections gilt, aber nicht für alle cases?
    //    const unsigned int length = 10;
    //    const unsigned int q = 2;
    //    const double p_mut = 0.1;
    //    const double p_error = p_mut/10;
    //    const double p_effect = 0.5;
    //    const double p_epistasis = 0.3;

    constants::Constants& cons = constants::Constants::get_instance();

    std::cout << "****** Sample mutational effects *******" << std::endl;
    // Create Ground Truth: Effects of each mutated position and epistatic effects and sequencing noise
    FunctionalSequence& effects = FunctionalSequence::get_instance();
    std::vector<Mutation> mutationsPerPos_vec;
    mutationsPerPos_vec.reserve(cons.L);
    for (int i = 1; i <= cons.L; ++i)
    {
        for (int j = 1; j < cons.Q; ++j)
        {
            // create instance for each
            mutationsPerPos_vec.emplace_back(i, j, effects.getKd(i));
        }
    }

    std::cout << "****** Create species *******" << std::endl;
    // Create M species
    auto specId_map = species::drawSpeciesIds();
    species::species_map species_vec;
    for (auto it = specId_map.begin(); it != specId_map.end(); ++it)
    {
        auto currentObj = species_vec.emplace(it->first, it->first);
        // first is a pointer to just constructed pair, second is the species object to call the methods
        currentObj.first->second.setCount(it->second);
        currentObj.first->second.computeSpeciesKd();
    }
    UnboundProtein f(species_vec);

    SECTION("Test ODE")
    {
        std::valarray<count_type> S_bound_tot(species_vec.size());
        ;
        std::valarray<count_type> S_unbound_tot(species_vec.size());
        ;
        double B = f.solve(S_bound_tot, S_unbound_tot);
        // check if the frequencies are summing up to 1
        REQUIRE(Approx(f.getFrequencies().sum()) == 1.0);
        // Test if the given equations for the ODE are fullfilled:
        // for (int i = 0; i < species_vec.size(); ++i) {
        int i = 0;
        for (auto& spec : species_vec)
        {
            REQUIRE(S_bound_tot[i] + S_unbound_tot[i] == Approx(spec.second.getCount()));
            // TODO because S_bound is sampled from the binomial distribution, the values have to be roughly similar, b
            // but cannot be tested here
            // REQUIRE(S_bound_tot[i] == Approx(S_unbound_tot[i]*B/spec.second.getKd()));

            REQUIRE((S_bound_tot[i] + S_unbound_tot[i]) / (double)cons.M == Approx(f.getFrequencies()[i]));
            // TODO s.o
            //            REQUIRE(S_bound_tot[i] == Approx(S_unbound_tot[i]*B/f.getKds()[i]));

            ++i;
        }
        REQUIRE(Approx(S_bound_tot.sum() / (double)cons.M + B).epsilon(0.001) == f.getB_tot());
    }
}