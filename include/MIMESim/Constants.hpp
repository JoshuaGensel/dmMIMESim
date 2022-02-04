//
//  Constants.hpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 04.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#ifndef Constants_hpp
#define Constants_hpp

#include "Utils.hpp"

#include <array>
#include <filesystem>
#include <math.h>
#include <stdio.h>

namespace fs = std::filesystem;

namespace constants
{
    struct Constants
    {
        // static = lifetime during execution (run time);
        // constexpr value available at compile time -> constant,  (not necessarily run time)

        /***** Constants regarding output *******/
        // directory where the results and parameter file are stored
        const fs::path OUTPUT_DIR = "./results";
        // TODO lösche static const std::string OUTPUT_DIR;
        // parameter file name never changes
        static constexpr auto PARAMETER_FILE = "parameters.txt";

        /***** Constants regarding kds sampling *******/
        // absolute wildtype Kd
        const double KD_WT = 1.0;
        // probability for each position to have an effect when mutated
        // static constexpr double P_EFFECT = 0.5;
        const double P_EFFECT = 0.5;
        // no multiplicative effect, when there is no epistasis
        const double NO_EPISTASIS = 1.0;
        // probability for each position pair to have epistatic effects
        // static constexpr double P_EPISTASIS = 0.75;
        const double P_EPISTASIS = 0.75;

        /**** Constants regarding species sampling *******/
        // maximal number of mutations in one sequence
        // static constexpr unsigned int MAX_MUT = 3;
        const unsigned int MAX_MUT;
        // number of total sequences
        const unsigned int M = 12 * pow(10, 6);
        // sequence length
        const unsigned int L = 50;
        // number of single values with al symbols
        const unsigned int SVal;
        // number of pairwise values with all combinations of symbols
        const unsigned int PWVal;
        // symbols per position
        const unsigned int Q = 2;
        // probability for a mutation
        const double P_MUT = 0.01;
        // probability for sequencing error
        const double P_ERR = 0.001;
        // id range for 0..MAX_MUT mutations for a sequence length L: (0 -> 1, 1 -> 2..L+1, etc)
        // const std::array<unsigned int, MAX_MUT+1> NMUT_RANGE;
        const std::vector<unsigned long long> NMUT_RANGE;
        // const std::array<double, MAX_MUT+1> P_NMUT;
        const std::vector<long double> P_NMUT;

        // seed for random number generator
        const unsigned int SEED;

        // total amount protein, relative to M
        const double BTOT = 2.0;

        std::vector<unsigned long long> setNMutRange(const unsigned int maxMut, const unsigned int L,
                                                     const unsigned int q);
        const std::vector<long double> setP_NMut(unsigned int max_mut, unsigned int l, double p_mut);
        unsigned int computeMaxMut(unsigned int m, unsigned int l, double p_mut);

        Constants(unsigned int length, unsigned int q, unsigned int m, double p_mut, double p_error, double p_effect,
                  double p_epistasis, unsigned int seed, double B_tot, fs::path outputDir)
            : L{length}, M{m}, SVal{length * (q - 1)}, PWVal{uint((length * (length - 1) / 2) * std::pow(q - 1, 2))},
              Q{q}, P_MUT{p_mut}, P_ERR{p_error}, P_EFFECT{p_effect}, P_EPISTASIS{p_epistasis}, MAX_MUT{computeMaxMut(
                                                                                                    m, length, p_mut)},
              NMUT_RANGE{setNMutRange(computeMaxMut(m, length, p_mut), length, q)},
              P_NMUT{setP_NMut(computeMaxMut(m, length, p_mut), length, p_mut)}, SEED{seed}, BTOT{B_tot},
              OUTPUT_DIR{outputDir} {};


        // Constructor for combining species sets; most params are the same as in params2, except MAX_MUT related params
        Constants(Constants const& params1, Constants const& params2)
            : L{params2.L}, Q{params2.Q}, M{params2.M}, MAX_MUT{params1.MAX_MUT + params2.MAX_MUT},
              NMUT_RANGE{setNMutRange(std::max(params1.MAX_MUT, params2.MAX_MUT) +
                                          std::round(params2.P_ERR / params2.P_MUT * params2.MAX_MUT),
                                      params2.L, params2.Q)},
              SVal{params2.SVal}, PWVal{params2.PWVal}, P_MUT{params2.P_MUT}, P_ERR{params2.P_ERR},
              P_EFFECT{params2.P_EFFECT}, P_EPISTASIS{params2.P_EPISTASIS}, P_NMUT{params2.P_NMUT}, SEED{params2.SEED},
              BTOT{params2.BTOT}, OUTPUT_DIR{params2.OUTPUT_DIR} {};
    };

    /**
     * Read in parameters from a given parameter file in the given result directory. If there is no file, use
     * default parameters
     * TODO: entweder result ordner angeben als Muss, wo ggf die parameter liste drin ist
     */
    const Constants& readParameters(const fs::path& outputPath);

    /**
     * TODO noch abfragen, dass nur bei "" in den cout geschrieben werden soll?
     * Write parameters into a parameter file in the given output path.
     * If the path does not exist, the parameters are printed into cout
     */
    void writeParameters(const fs::path& outputPath, const Constants& params);

    /**
     * Write parameters into cout
     */
    void writeParameters(const Constants& params);
}
#endif /* Constants_hpp */
