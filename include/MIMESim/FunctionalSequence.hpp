//
// Created by Smith, Maureen on 31.05.18.
//

#ifndef DCABENCHMARK_FUNCTIONALSEQUENCE_HPP
#define DCABENCHMARK_FUNCTIONALSEQUENCE_HPP

#include "Mutation.hpp"

#include <iostream>
#include <random>
#include <vector>

class FunctionalSequence
{

  private:
    // original Kd for each site
    const std::vector<double> kds;
    // original epistasis values for each position pair
    const std::vector<double> epistasis;

    // singleton class has a private constructor which is called by a public method
    FunctionalSequence() : kds(drawKdValues()), epistasis(drawEpistasis()){};

    // singleton class has a private constructor which is called by a public method
    FunctionalSequence(const std::string &inputPath) : kds(readKdValues(inputPath)), epistasis(readEpistasis(inputPath)){};

    /**
     * Draws the Kds value according to log normal distribution with a probability p_kd for each L postions and q-1
     * possible mutations
     * @return the vector with all Kd values
     */
    std::vector<double> drawKdValues();

    /**
     * Reads saved Kds values from disk from the inputPath
     * @return the vector with all Kd values
     */
    std::vector<double> readKdValues(const std::string &inputPath);

    /**
     * Draws the epistasis value for each PWVal (= each pair of positions and each combination of possible mutation
     * symbols) according to log normal distribution with a probability p_epi
     * @return the vector with all Epistasis values
     */
    std::vector<double> drawEpistasis();

    /**
     * Reads saved Epistasis values from disk from the inputPath
     * @return the vector with all Epistasis values
     */
    std::vector<double> readEpistasis(const std::string &inputPath);

    /*
     * computes the index within the vector representation of a symmetric matrix
     * @param a Mutation 1
     * @param b Mutation 2
     * @return the Index of the Vector of the the pairwise values in the matrix (upper or lower triangle)
     */
    unsigned int getMatrixVectorIndex(const Mutation& a, const Mutation& b);

    /**
     * computes index within vector of position and symbol
     * @param m Mutation = pair of position and symbol
     * @return the index within the vector
     */
    unsigned int getVectorIndex(const Mutation& m);

  public:
    /**
     * @return Vector of Kd values for or all possible mutations in ascending order.
     */
    const std::vector<double>& getKd();

    /**
     * @param p Mutated position and symbol
     * @return Kd value for the given mutation.
     */
    const double getKd(const Mutation& p);

    /**
     * @return Vector of epistasis values for or all pairwise mutations in ascending order with position i<j.
     */
    const std::vector<double>& getEpistasis();

    //
    /**
     * The position of mutations a has to be smaller than the position of mutations b.
     * @param a Mutation 1
     * @param b Mutation 2
     * @return Epistasis value for the given pair of mutations
     */
    const double& getEpistasis(const Mutation& a, const Mutation& b);

    /**
     *
     * @return singleton instance of the true KD and Epistasis information
     */
    static FunctionalSequence& get_instance();

    /**
     *
     * @return singleton instance of the true KD and Epistasis information as provided in inputPath
     */
    static FunctionalSequence &get_instance(const std::string &inputPath);

    /**
     * Write all Kd values in consecutive order into the given outputfile, delimited by a new line
     * @param filename Outputfile
     */
    void writeKdsToFile(const std::string& filename);

    /**
     * Write all Epistasis values in consecutive order into the given outputfile, delimited by a new line
     * @param filename Outputfile
     */
    void writeEpistasisToFile(const std::string& filename);

    FunctionalSequence(FunctionalSequence const&) = delete;

    FunctionalSequence(FunctionalSequence&&) = delete;
};

#endif // DCABENCHMARK_FUNCTIONALSEQUENCE_HPP
