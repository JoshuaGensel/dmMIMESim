//
//  Constants.cpp
//  DCABenchmark
//
//  Created by Smith, Maureen on 04.05.18.
//  Copyright © 2018 Smith, Maureen. All rights reserved.
//

#include "Constants.hpp"

#include <exception>
#include <fstream>
#include <iostream>
#include <math.h>

namespace constants
{

    /*
     * determination of the maximal number of mutations (Expected number of occurrence >5), must be minimum 3
     */
    // TODO: Is it correct to enforce the minimum of 3?
    unsigned int Constants::computeMaxMut(unsigned int m, unsigned int l, double p_mut)
    {
        unsigned int maxMut = 2;
        do
        {
            ++maxMut;
        } while ((static_cast<long double>(m) * utils::nChoosek(l, maxMut + 1) * pow(p_mut, maxMut + 1) *
                      pow(1 - p_mut, l - maxMut - 1) >
                  5) &&
                 (maxMut < l));

        return maxMut;
    }

    std::vector<unsigned long long> Constants::setNMutRange(const unsigned int maxMut, const unsigned int L,
                                                            const unsigned int q)
    {
        // compute the number of possible sequence for 0..MAX_MUT mutations, the cumulative sum gives the id range for
        // each number of mutations larger range for adding errors (i.e. more mutated positions)
        std::vector<unsigned long long> nMutRange(maxMut * 2 + 1);
        nMutRange[0] = 1;
        for (unsigned int i = 1; i <= maxMut * 2; ++i)
        {
            // cummulative sum of L choose i position combinations and (q-1)^i symbol combinations
            nMutRange[i] = nMutRange[i - 1] + llround(utils::nChoosek(L, i) * pow(q - 1, i));

            if (nMutRange[i] < nMutRange[i - 1])
            {
                throw std::overflow_error("Overflow error in setNMutRange.");
            }
        }
        return (nMutRange);
    }

    /*
     * computing the probability for each number of mutations
     */
    const std::vector<long double> Constants::setP_NMut(unsigned int max_mut, unsigned int l, double p_mut)
    {
        std::vector<long double> p_nmut(max_mut + 1);
        long double p_sum = 0;
        for (unsigned i = 0; i < max_mut; ++i)
        {
            p_nmut[i] = utils::nChoosek(l, i) * pow(p_mut, i) * pow(1 - p_mut, l - i);
            p_sum += p_nmut[i];
        }
        // the probabilities for n > n_max are added to n_max
        p_nmut[max_mut] = 1 - p_sum;
        return p_nmut;
    };

    const Constants& readParameters(const fs::path& outputPath)
    {

        if (!outputPath.empty())
        {
            // if output directory does not exist yet, create it
            if (!fs::exists(outputPath))
            {
                fs::create_directory(outputPath);
                std::cout << "Create output directory " << fs::canonical(outputPath) << std::endl;
            }

            fs::path paraFile(fs::canonical(outputPath));

            // the path can either be a directory, where the parameter file has the standard name, or it the path
            // contains a given parameter file
            if (fs::is_directory(outputPath))
            {
                paraFile = (fs::canonical(outputPath) / Constants::PARAMETER_FILE);
            }

            // dafault parameters, in case no file is given, or paramater are not set in the file
            unsigned int L = 50;
            // symbols per position
            unsigned int q = 2;
            // mutation probability
            double p_mut = 0.1;
            double p_error = p_mut / 10.0;
            double p_effect = 0.5;
            double p_epistasis = 0.3;
            unsigned int M = 12 * pow(10, 6);
            unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());

            int max_mut = -1;
            double B_tot = 2.0;
            bool epiMutExcl = false;

            // if the output directory contains the parameter file, read it
            if (fs::exists(paraFile) && fs::is_regular_file(paraFile))
            {
                // Paths are implicitly convertible to and from std::basic_strings, so we can call
                std::ifstream infile(paraFile);
                if (infile.good())
                {
                    std::cout << "Reading file " << paraFile << std::endl;
                    std::string line;
                    std::string param;
                    std::string val;

                    try
                    {

                        while (std::getline(infile, line))
                        {

                            // stream through each line to read the parameter name and its value, seperated by tabular
                            std::istringstream lineSS(line);

                            std::getline(lineSS, param, '\t');
                            std::getline(lineSS, val, '\t');
                            // TODO make variable
                            // if(param == "kd_wt")
                            //     kd_wt = std::stoi(val);
                            if (param == "M")
                                M = std::stoi(val);
                            if (param == "L")
                                L = std::stoi(val);
                            if (param == "q")
                                q = std::stoi(val);
                            if (param == "p_mut")
                                p_mut = std::stod(val);
                            if (param == "p_error")
                                p_error = std::stod(val);
                            if (param == "p_effect")
                                p_effect = std::stod(val);
                            if (param == "p_epistasis")
                                p_epistasis = std::stod(val);
                            if (param == "seed")
                                seed = std::stoi(val);
                            if (param == "max_mut")
                                max_mut = std::stoi(val);
                            if (param == "B_tot")
                                B_tot = std::stod(val);
                            if (param == "epistasis_mut_excl")
                            {
                                if ((val == "true") || (val == "True"))
                                    epiMutExcl = true;
                            }
                        }
                        infile.close();
                        std::cout << " ... successful." << std::endl;
                    }
                    catch (const std::invalid_argument& ia)
                    {
                        std::cerr << "Error in parameter file. Invalid argument for parameter " << param << " ("
                                  << ia.what() << ")" << '\n';
                        std::cout << "Using default parameters." << std::endl;
                    }
                }
            }
            else
            {
                //...if no config file is given, create the constants with the default parameter values
                std::cout << "No parameter file given. Using default parameters." << std::endl;
            }

            // Create constants which are used through out this test set
            if (max_mut == -1)
            {
                // no fixed max_mut
                Constants* cons =
                    new Constants(L, q, M, p_mut, p_error, p_effect, p_epistasis, seed, B_tot, epiMutExcl, outputPath);
                return *cons;
            }
            else
            {
                // with fixed max_mut
                Constants* cons = new Constants(L, q, M, p_mut, p_error, p_effect, p_epistasis, seed, B_tot, max_mut,
                                                epiMutExcl, outputPath);
                return *cons;
            }
        }
        else
        {
            std::cerr << "Output path is a mandatory parameter" << std::endl;
            throw std::invalid_argument("Output path is a mandatory parameter");
        }
    }

    void writeParameters(const fs::path& outputPath, const Constants& params)
    {
        // initialize the stream with the
        std::ostream* paraStream = &std::cout;
        std::ofstream ofs;

        if (fs::exists(outputPath) && fs::is_directory(outputPath))
        {
            auto paraFile(outputPath / params.PARAMETER_FILE);
            std::cout << "Parameter File: " << fs::canonical(outputPath) << std::endl;
            ofs.open(paraFile);
            paraStream = &ofs;
        }

        if ((*paraStream).good())
        {
            std::string epistasis_mut_excl = params.EPIMUTEXCL ? "true" : "false";
            (*paraStream) << "seed\t" << params.SEED << '\n';
            (*paraStream) << "### paratemers regarding kd sampling ###\n";
            (*paraStream) << "kd_wt\t" << params.KD_WT << '\n';
            (*paraStream) << "p_effect\t" << params.P_EFFECT << '\n';
            (*paraStream) << "p_epistasis\t" << params.P_EPISTASIS << '\n';
            (*paraStream) << "epistasis_mut_excl\t" << epistasis_mut_excl << '\n';
            (*paraStream) << "### paramters regarding sequence sampling ###\n";
            (*paraStream) << "L\t" << params.L << '\n';
            (*paraStream) << "q\t" << params.Q << '\n';
            (*paraStream) << "M\t" << params.M << '\n';
            (*paraStream) << "p_mut\t" << params.P_MUT << '\n';
            (*paraStream) << "p_error\t" << params.P_ERR << '\n';
            (*paraStream) << "max_mut\t" << params.MAX_MUT << '\n';
            (*paraStream) << "### paramters regarding binding competition ###\n";
            (*paraStream) << "B_tot\t" << params.BTOT << '\n';
        }
    }

    void writeParameters(const Constants& params)
    {
        writeParameters("", params);
    }
}
