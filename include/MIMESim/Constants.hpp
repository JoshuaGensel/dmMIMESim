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
        const std::vector<unsigned int> NMUT_RANGE;
        // const std::array<double, MAX_MUT+1> P_NMUT;
        const std::vector<double> P_NMUT;

        // std::array<unsigned int, MAX_MUT+1> setNMutRange();
        // std::array<double, MAX_MUT+1> setP_NMut();
        std::vector<unsigned int> setNMutRange(const unsigned int maxMut, const unsigned L);
        // TODO: austauschen
        std::vector<unsigned int> setNMutRange(const unsigned int maxMut, const unsigned int L, const unsigned int q);
        std::vector<double> setP_NMut(const unsigned int MaxMut, const unsigned L, const double pMut);
        unsigned int computeMaxMut(const unsigned int L, const double pMut);

        Constants(unsigned int length, unsigned int q, double p_mut, double p_error, double p_effect,
                  double p_epistasis, fs::path outputDir)
            : L(length), SVal(length * (q - 1)), PWVal((length * (length - 1) / 2) * std::pow(q - 1, 2)), Q(q),
              P_MUT(p_mut), P_ERR(p_error), P_EFFECT(p_effect), P_EPISTASIS(p_epistasis),
              MAX_MUT(computeMaxMut(length, p_mut / (q - 1))),
              NMUT_RANGE(setNMutRange(computeMaxMut(length, p_mut / (q - 1)), L, q)),
              P_NMUT(setP_NMut(computeMaxMut(length, p_mut / (q - 1)), length, p_mut)), OUTPUT_DIR(outputDir){};

        // Constructor for combining species sets; most params are the same as in params2, except MAX_MUT related params
        Constants(Constants const &params1, Constants const &params2)
            : L(params2.L), Q(params2.Q), M(params2.M), MAX_MUT(params1.MAX_MUT + params2.MAX_MUT),
              NMUT_RANGE(setNMutRange(params1.MAX_MUT + params2.MAX_MUT, params2.L, params2.Q)), SVal(params2.SVal),
              PWVal(params2.PWVal), P_MUT(params2.P_MUT), P_ERR(params2.P_ERR), P_EFFECT(params2.P_EFFECT),
              P_EPISTASIS(params2.P_EPISTASIS), P_NMUT(params2.P_NMUT), OUTPUT_DIR(params2.OUTPUT_DIR){};

        // Constants(unsigned int length, unsigned int q, double p_mut) : L(length),PWVal(L*(L-1)/2), Q(q),
        // P_MUT(p_mut), NMUT_RANGE(setNMutRange()), P_NMUT(setP_NMut()) {};

      private:
        // TODO: workaround, 3 mal die gleiche routine (computeMaxMut) aufrufen, weil es erst am  ende alles gespeichert
        // wird, anders lösen?
        // Lösung siehe Species: indexToMutPos, erst aufruf ohne parameter um dann den member parameter mitugeben.
        // TODO: computeMaxMut mit p_mut/q-1, da ja p_mut die wkeit ist dass eine position überhaupt mutiert, aber nicht
        // fpr jede einzelne muation
    };

    /**
     * Read in parameters from a given parameter file in the given result directory. If there is no file, use default
     * parameters
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
