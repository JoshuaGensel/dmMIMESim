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
#include "Species.hpp"
#include "Utils.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace fs = std::filesystem;

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

    utils::readParameters(outputPath);

    // get the newly created instance of the constants
    auto& cons = constants::Constants::get_instance();

    std::cout << "MaxMut " << cons.MAX_MUT << std::endl;

    // The 4 output files are saved with their ids where wild_type_bound = firstId, wild_type_unbound = firstId+1,
    // mut_bound = firstId+2, mut_unbound = firstId+3
    int firstId = 1;
    std::string wt_bound_id(std::to_string(firstId));
    std::string wt_unbound_id(std::to_string(firstId + 1));
    std::string mut_bound_id(std::to_string(firstId + 2));
    std::string mut_unbound_id(std::to_string(firstId + 3));

    // TODO in einen Test mit einbringen
    //     std::cout << "Maximal number of mutations: " << cons.MAX_MUT << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Sample mutational effects *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();

    // Create Ground Truth: Effects of each mutated position and epistatic effects and sequencing noise
    FunctionalSequence& effects =
        use_inputPath ? FunctionalSequence::get_instance(inputPath) : FunctionalSequence::get_instance();

    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Create species *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // Create M species, the map contains the counts for each sampled sequence id
    species::species_map species_vec = species::drawSpeciesIds();

    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    // TODO: diese Art der Abfrage in die Tests packen
    //     std::cout << "wt species count + freq. " << species_vec.at(1).getCount() << " " <<
    //     species_vec.at(1).getFreq() << std::endl; std::cout << "mut species count + freq. " <<
    //     species_vec.at(2).getCount() << " " << species_vec.at(2).getFreq() << std::endl; std::cout << "mut species
    //     count + freq. " << species_vec.at(3).getCount() << " " << species_vec.at(3).getFreq() << std::endl;
    //     //std::cout << "mut bound unbound freq. " << species_vec.at(20877).getCount() << " " <<
    //     species_vec.at(20877).getFreq() << std::endl;

    std::cout << "****** Create unmutated wild type library  *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // auto wtSpecId_map = species::drawWildtypeErrors();
    // std::vector<species::Species> wtSpecies_vec;
    // The "control expereriment / wild type library" contains only wildtype sequences
    species::species_map wtSpecies_vec;
    auto currentObj = wtSpecies_vec.emplace(1, 1);
    currentObj.first->second.setCount(cons.M);
    currentObj.first->second.computeSpeciesKd();

    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Solve ODE to infer bound and unbound fraction of species *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // TODO: Umbau nach counts
    // std::valarray<double> f_bound_tot;
    // std::valarray<double> f_unbound_tot;
    std::valarray<unsigned int> S_bound(species_vec.size());
    std::valarray<unsigned int> S_unbound(species_vec.size());

    // std::valarray<int> f_bound_tot(species_vec.size());
    // std::valarray<int> f_unbound_tot(species_vec.size());
    //  set up the ODE (binding competition) and solve it to get the bound and unbound fractions (from the total amount
    //  M) in equilibrium
    UnboundProtein f(species_vec);
    f.solve(S_bound, S_unbound);

    // stimmt ja so nicht mehr, da unrdered map
    //    std::cout << "wt bound unbound freq. " << S_bound[0] << " " << S_unbound[0] << std::endl;
    //    std::cout << "mut bound unbound freq. " << S_bound[1] << " " << S_unbound[1] << std::endl;
    //    std::cout << "mut bound unbound freq. " << S_bound[2] << " " << S_unbound[2] << std::endl;

    //*************solve for wt
    // TODO: Umbau nach counts
    // std::valarray<double> f_bound_tot_wt;
    // std::valarray<double> f_unbound_tot_wt;
    std::valarray<unsigned int> S_bound_wt(wtSpecies_vec.size());
    std::valarray<unsigned int> S_unbound_wt(wtSpecies_vec.size());
    // set up the ODE (binding competition) and solve it to get the bound and unbound fractions (from the total amount
    // M) in equilibrium
    UnboundProtein f_wt(wtSpecies_vec);
    f_wt.solve(S_bound_wt, S_unbound_wt);

    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Add noise and Count *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // Carefull: The map is extended by species that occur only because of sequencing error, hence the length of S_bound
    // and S_unbound dont fit any more with the length of the map
    // TODO: Umbau nach counts
    // TODO weg
    // species::addCountsWithError(S_bound, S_unbound, species_vec);

    // species::addCountsWithError(S_bound_wt, S_unbound_wt, wtSpecies_vec);

    std::cout << "Mutation Library" << std::endl;
    ;
    auto counters = species::countMutationsWithErrors(S_bound, S_unbound, species_vec);

    std::cout << "Wild type Library" << std::endl;
    ;
    auto counters_wt = species::countMutationsWithErrors(S_bound_wt, S_unbound_wt, wtSpecies_vec);

    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Write output to file *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();

    if (!fs::exists(outputPath))
        ;
    fs::create_directory(outputPath);

    // create subdirectories for the single and double mutant counts
    fs::create_directory(outputPath / "2d");
    fs::create_directory(outputPath / "1d");
    // create for each barcode a textfile and write the counts into the files
    counters.counter_bound_1d.write_to_file(outputPath / "1d" / (mut_bound_id + ".txt"));
    counters.counter_unbound_1d.write_to_file(outputPath / "1d" / (mut_unbound_id + ".txt"));
    counters.counter_bound_2d.write_to_file(outputPath / "2d" / (mut_bound_id + ".txt"));
    counters.counter_unbound_2d.write_to_file(outputPath / "2d" / (mut_unbound_id + ".txt"));

    counters_wt.counter_bound_1d.write_to_file(outputPath / "1d" / (wt_bound_id + ".txt"));
    counters_wt.counter_unbound_1d.write_to_file(outputPath / "1d" / (wt_unbound_id + ".txt"));
    counters_wt.counter_bound_2d.write_to_file(outputPath / "2d" / (wt_bound_id + ".txt"));
    counters_wt.counter_unbound_2d.write_to_file(outputPath / "2d" / (wt_unbound_id + ".txt"));

    fs::create_directory(outputPath / "sequences");
    species::write_to_file(outputPath / "sequences" / (mut_bound_id + ".txt"), species_vec, S_bound);
    species::write_to_file(outputPath / "sequences" / (mut_unbound_id + ".txt"), species_vec, S_unbound);

    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Write true values to files *******" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::cout << "Write epistasis" << std::endl;
    effects.writeEpistasisToFile(outputPath / "pairwise_epistasis.txt");
    std::cout << "Write KD" << std::endl;
    effects.writeKdsToFile(outputPath / "single_kds.txt");

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
                        auto doubleKd = effects.getKd(mut1) * effects.getKd(mut2) * effects.getEpistasis(mut1, mut2);
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

    utils::writeParameters(cons.OUTPUT_DIR);

    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "Duration: " << diff.count() << " s\n";

    std::cout << "****** Done *******" << std::endl;

    return 0;
}
