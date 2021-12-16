//
//  main.cpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 03.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#include "BindingCompetition.hpp"
#include "Constants.hpp"
#include "Count.hpp"
#include "FunctionalSequence.hpp"
#include "Generator.hpp"
#include "Species.hpp"
#include "Utils.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>

namespace fs = std::filesystem;

species::species_map firstRoundMutagenesis(fs::path outputPath)
{
    const constants::Constants& params = constants::readParameters(outputPath);

    std::cout << "****** Create new species *******" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    // Create M species, the map contains the counts for each sampled sequence id
    species::species_map species_vec = species::drawSpeciesIds(params);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    return species_vec;
}

std::tuple<species::species_map, species::species_map, constants::Constants> secondRoundMutagenesis(fs::path inputPath,
                                                                                                    fs::path outputPath)
{
    // TODO: fix printing in function readParameters
    const constants::Constants& paramFirst = constants::readParameters(inputPath);
    const constants::Constants& paramSecond = constants::readParameters(outputPath);

    std::cout << "****** Load old species *******" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    species::species_map specOldBound = species::readFromFile(inputPath, utils::SampleID::mut_bound, paramSecond);
    species::species_map specOldUnbound = species::readFromFile(inputPath, utils::SampleID::mut_unbound, paramSecond);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    unsigned int nBound = species::speciesMapSum(specOldBound);
    unsigned int nUnbound = species::speciesMapSum(specOldUnbound);

    // change M parameter, all other constants stay the same
    constants::Constants paramSecondBound =
        constants::Constants(paramSecond.L, paramSecond.Q, nBound, paramSecond.P_MUT, paramSecond.P_ERR,
                             paramSecond.P_EFFECT, paramSecond.P_EPISTASIS, paramSecond.SEED, paramSecond.OUTPUT_DIR);
    constants::Constants paramSecondUnbound =
        constants::Constants(paramSecond.L, paramSecond.Q, nUnbound, paramSecond.P_MUT, paramSecond.P_ERR,
                             paramSecond.P_EFFECT, paramSecond.P_EPISTASIS, paramSecond.SEED, paramSecond.OUTPUT_DIR);

    std::cout << "****** Create new species *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // Create M species, the map contains the counts for each sampled sequence id
    auto fut_specNewBound = std::async(std::launch::async, species::drawSpeciesIds, paramSecondBound);
    auto fut_specNewUnbound = std::async(std::launch::async, species::drawSpeciesIds, paramSecondUnbound);
    species::species_map specNewBound = fut_specNewBound.get();
    species::species_map specNewUnbound = fut_specNewUnbound.get();
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Combine old and new species *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    auto fut_specCombBound = std::async(std::launch::async, species::combineSpecies, specOldBound, specNewBound,
                                        paramFirst, paramSecondBound);
    auto fut_specCombUnbound = std::async(std::launch::async, species::combineSpecies, specOldUnbound, specNewUnbound,
                                          paramFirst, paramSecondBound);
    auto result = std::make_tuple(fut_specCombBound.get(), fut_specCombUnbound.get(), paramSecond);
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    return result;
}

std::tuple<count::counter_collection, std::valarray<unsigned int>, std::valarray<unsigned int>>
runSelection(species::species_map species_vec, const constants::Constants& params)
{
    std::cout << "****** Solve ODE to infer bound and unbound fraction *******" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    // TODO: Umbau nach counts
    // std::valarray<double> f_bound_tot;
    // std::valarray<double> f_unbound_tot;
    std::valarray<unsigned int> S_bound(species_vec.size());
    std::valarray<unsigned int> S_unbound(species_vec.size());

    // std::valarray<int> f_bound_tot(species_vec.size());
    // std::valarray<int> f_unbound_tot(species_vec.size());
    //  set up the ODE (binding competition) and solve it to get the bound and unbound fractions (from the total
    //  amount M) in equilibrium
    UnboundProtein f(species_vec);
    f.solve(S_bound, S_unbound);

    // stimmt ja so nicht mehr, da unrdered map
    //    std::cout << "wt bound unbound freq. " << S_bound[0] << " " << S_unbound[0] << std::endl;
    //    std::cout << "mut bound unbound freq. " << S_bound[1] << " " << S_unbound[1] << std::endl;
    //    std::cout << "mut bound unbound freq. " << S_bound[2] << " " << S_unbound[2] << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Count *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // Carefull: The map is extended by species that occur only because of sequencing error, hence the length of
    // S_bound and S_unbound dont fit any more with the length of the map
    // TODO: Umbau nach counts
    // TODO weg
    // species::addCountsWithError(S_bound, S_unbound, species_vec);

    auto counters = species::countMutations(S_bound, S_unbound, species_vec, params);

    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    return std::make_tuple(counters, S_bound, S_unbound);
}

void writeToFile(fs::path outputPath, count::counter_collection& counter, species::species_map& species_vec,
                 std::vector<std::set<Mutation>>& errors, std::valarray<unsigned int>& S_bound,
                 std::valarray<unsigned int>& S_unbound, utils::SampleID boundID, utils::SampleID unboundID)
{
    counter.counter_bound_1d.write_to_file(outputPath / "1d" / (utils::SampleIDStr(boundID) + ".txt"));
    counter.counter_unbound_1d.write_to_file(outputPath / "1d" / (utils::SampleIDStr(unboundID) + ".txt"));
    counter.counter_bound_2d.write_to_file(outputPath / "2d" / (utils::SampleIDStr(boundID) + ".txt"));
    counter.counter_unbound_2d.write_to_file(outputPath / "2d" / (utils::SampleIDStr(unboundID) + ".txt"));

    species::writeSpeciesToFile(outputPath / "species" / (utils::SampleIDStr(boundID) + ".txt"), species_vec, S_bound);
    species::writeSpeciesToFile(outputPath / "species" / (utils::SampleIDStr(unboundID) + ".txt"), species_vec,
                                S_unbound);
}

int main(int argc, const char* argv[])
{

    std::cout << "****** Set up constants *******" << std::endl;
    // measuring realtime duration (see std::clock for cpu time)
    auto start = std::chrono::high_resolution_clock::now();

    fs::path outputPath("../results");
    fs::path inputPath("");
    bool use_inputPath = false;

    if (argc > 1)
    {
        outputPath = argv[1];
    }

    // TODO: Cleaner solution for this?
    // TODO: Document that input path can only be used when output path is specified as well.
    if (argc > 2)
    {
        inputPath = argv[2];
        use_inputPath = true;
    }

    if (!fs::exists(outputPath))
    {
        fs::create_directory(outputPath);
        std::cout << "Create output directory " << fs::canonical(outputPath) << std::endl;
    }
    else
    {
        std::cout << "Using output directory " << fs::canonical(outputPath) << std::endl;
    }

    if (use_inputPath)
    {
        if (!fs::exists(inputPath))
        {
            std::cerr << "Input directory does not exist:\n" << fs::canonical(inputPath) << std::endl;
            return 1;
        }
        else
        {
            std::cout << "Using input directory " << fs::canonical(inputPath) << std::endl;
        }

        if (!fs::exists(inputPath / "single_kds.txt"))
        {
            std::cerr << "Input file does not exist:\n" << fs::canonical(inputPath) / "single_kds.txt" << std::endl;
            return 2;
        }
        else if (!fs::exists(inputPath / "pairwise_epistasis.txt"))
        {
            std::cerr << "Input file does not exist:\n"
                      << fs::canonical(inputPath) / "pairwise_epistasis.txt" << std::endl;
            return 3;
        }
        else
        {
            std::cout << "All required files are present in input directory. " << std::endl;
        }
    }

    // get the constants of the simulation from the parameters file
    const constants::Constants& cons = constants::readParameters(outputPath);
    constants::writeParameters(cons);
    std::cout << "MaxMut " << cons.MAX_MUT << std::endl;

    // create an instance of the random number generator
    Generator::create_instance(cons.SEED);

    // create subdirectories for the single and double mutant counts
    fs::create_directory(outputPath / "2d");
    fs::create_directory(outputPath / "1d");
    fs::create_directory(outputPath / "species");

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Sample mutational effects *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();

    // Create or load Ground Truth: Effects of each mutated position and epistatic effects and sequencing noise
    auto* effects =
        (use_inputPath) ? FunctionalSequence::create_instance(inputPath) : FunctionalSequence::create_instance(cons);

    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    // Draw errors for both libraries
    auto errors_lib1 = species::drawErrors(cons);
    auto errors_lib2 = species::drawErrors(cons);

    // first round of mutagenesis
    if (!use_inputPath)
    {
        species::species_map species_vec = firstRoundMutagenesis(outputPath);
        auto results = runSelection(species_vec, cons);
        count::counter_collection counters = std::get<0>(results);
        std::valarray<unsigned int> S_bound = std::get<1>(results);
        std::valarray<unsigned int> S_unbound = std::get<2>(results);

        std::cout << "****** Create unmutated wild type library  *******" << std::endl;
        start = std::chrono::high_resolution_clock::now();
        // The "control expereriment / wild type library" contains only wildtype sequences
        species::species_map wtSpecies_vec;
        auto currentObj = wtSpecies_vec.emplace(std::make_pair<int, species::Species>(1, {1, cons}));
        currentObj.first->second.setCount(cons.M);
        currentObj.first->second.computeSpeciesKd();
        end = std::chrono::high_resolution_clock::now();
        diff = end - start;
        std::cout << "Duration: " << diff.count() << " s\n";

        auto wtResults = runSelection(wtSpecies_vec, cons);
        count::counter_collection counters_wt = std::get<0>(wtResults);
        std::valarray<unsigned int> wtS_bound = std::get<1>(wtResults);
        std::valarray<unsigned int> wtS_unbound = std::get<2>(wtResults);

        std::cout << "****** Add errors to counts *******" << std::endl;
        start = std::chrono::high_resolution_clock::now();
        countErrors(errors_lib1, species_vec, S_bound, S_unbound, counters);
        countErrors(errors_lib2, wtSpecies_vec, wtS_bound, wtS_unbound, counters_wt);
        end = std::chrono::high_resolution_clock::now();
        diff = end - start;
        std::cout << "Duration: " << diff.count() << " s\n";

        std::cout << "****** Write output to file *******" << std::endl;
        start = std::chrono::high_resolution_clock::now();
        writeToFile(outputPath, counters, species_vec, errors_lib1, S_bound, S_unbound, utils::SampleID::mut_bound,
                    utils::SampleID::mut_unbound);
        writeToFile(outputPath, counters, wtSpecies_vec, errors_lib2, wtS_bound, wtS_unbound, utils::SampleID::wt_bound,
                    utils::SampleID::wt_unbound);
        end = std::chrono::high_resolution_clock::now();
        diff = end - start;
        std::cout << "Duration: " << diff.count() << " s\n";
    }
    // second round of mutagenesis
    else
    {
        auto results = secondRoundMutagenesis(inputPath, outputPath);
        species::species_map species_vec_b = std::get<0>(results);
        species::species_map species_vec_u = std::get<1>(results);
        const constants::Constants& cons_in = std::get<2>(results);

        // TODO: multithreading for results_bound, results_unbound; first change printing behavior
        auto results_bound = runSelection(species_vec_b, cons_in);
        count::counter_collection counters_bound = std::get<0>(results_bound);
        std::valarray<unsigned int> S_bound_bound = std::get<1>(results_bound);
        std::valarray<unsigned int> S_bound_unbound = std::get<2>(results_bound);
        auto results_unbound = runSelection(species_vec_u, cons_in);
        count::counter_collection counters_unbound = std::get<0>(results_unbound);
        std::valarray<unsigned int> S_unbound_bound = std::get<1>(results_unbound);
        std::valarray<unsigned int> S_unbound_unbound = std::get<2>(results_unbound);

        std::cout << "****** Add errors to counts *******" << std::endl;
        start = std::chrono::high_resolution_clock::now();
        countErrors(errors_lib1, species_vec_b, S_bound_bound, S_bound_unbound, counters_bound);
        countErrors(errors_lib2, species_vec_u, S_unbound_bound, S_unbound_unbound, counters_unbound);
        end = std::chrono::high_resolution_clock::now();
        diff = end - start;
        std::cout << "Duration: " << diff.count() << " s\n";

        std::cout << "****** Write output to file *******" << std::endl;
        start = std::chrono::high_resolution_clock::now();
        writeToFile(outputPath, counters_bound, species_vec_b, errors_lib1, S_bound_bound, S_bound_unbound,
                    utils::SampleID::mut_bound_bound, utils::SampleID::mut_bound_unbound);
        writeToFile(outputPath, counters_unbound, species_vec_u, errors_lib2, S_unbound_bound, S_unbound_unbound,
                    utils::SampleID::mut_unbound_bound, utils::SampleID::mut_unbound_unbound);
        end = std::chrono::high_resolution_clock::now();
        diff = end - start;
        std::cout << "Duration: " << diff.count() << " s\n";
    }

    std::cout << "****** Write true values to files *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::cout << "Write epistasis" << std::endl;
    effects->writeEpistasisToFile(outputPath / "pairwise_epistasis.txt");
    std::cout << "Write KD" << std::endl;
    effects->writeKdsToFile(outputPath / "single_kds.txt");

    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";
    start = std::chrono::high_resolution_clock::now();
    // write pairwise effects
    std::cout << "Write pairwise effects" << std::endl;
    std::ofstream outfile(outputPath / "pairwise_kds.txt");
    if (outfile.good())
    {
        for (unsigned pos1 = 1; pos1 < cons.L; ++pos1)
        {
            for (unsigned pos2 = pos1 + 1; pos2 <= cons.L; ++pos2)
            {
                for (unsigned sym1 = 0; sym1 < cons.Q - 1; ++sym1)
                {
                    Mutation mut1{pos1, sym1};
                    for (unsigned sym2 = 0; sym2 < cons.Q - 1; ++sym2)
                    {
                        Mutation mut2{pos2, sym2};
                        auto doubleKd = effects->getKd(mut1) * effects->getKd(mut2) * effects->getEpistasis(mut1, mut2);
                        outfile << doubleKd << '\n';
                    }
                }
            }
        }
        outfile.close();
    }

    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Write parameter into File *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();

    constants::writeParameters(cons.OUTPUT_DIR, cons);

    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Done *******" << std::endl;

    return 0;
}
