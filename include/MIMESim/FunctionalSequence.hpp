//
// Created by Smith, Maureen on 31.05.18.
//

#ifndef DCABENCHMARK_FUNCTIONALSEQUENCE_HPP
#define DCABENCHMARK_FUNCTIONALSEQUENCE_HPP

#include "Constants.hpp"
#include "Mutation.hpp"

#include <gtest/gtest_prod.h>
#include <iostream>
#include <random>
#include <vector>

class FunctionalSequence
{

  private:
    FRIEND_TEST(FunctionalSequenceTestPrivate, ReadKdValues);
    FRIEND_TEST(FunctionalSequenceTestPrivate, ReadEpistasis);
    FRIEND_TEST(FunctionalSequenceTestPrivate, GetMutationFromVectorIndex);
    FRIEND_TEST(FunctionalSequenceTestPrivate, PairIndex);
    FRIEND_TEST(FunctionalSequenceTestPrivate, DrawEpistasis_restricted);
    FRIEND_TEST(FunctionalSequenceTestPrivate, DrawEpistasis_semi_restricted);
    FRIEND_TEST(FunctionalSequenceTestPrivate, DrawEpistasis_unrestricted);

    // constants used for the simulation
    const constants::Constants params;
    // original Kd for each site
    const std::vector<double> kds;
    // index used in getPairIndex to index epistasis values
    const std::vector<int> pairIndex;
    // original epistasis values for each position pair
    const std::vector<double> epistasis;
    // The one, single instance
    static FunctionalSequence* instance;

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
    std::vector<double> readKdValues(const std::string& inputPath);

    /**
     * Draws the epistasis value for each PWVal (= each pair of positions and each combination of possible mutation
     * symbols) according to log normal distribution with a probability p_epi
     * @return the vector with all Epistasis values
     */
    std::vector<double> drawEpistasis();

    std::vector<double> drawEpistasis_restricted();
    std::vector<double> drawEpistasis_unrestricted();

    /**
     * Reads saved Epistasis values from disk from the inputPath
     * @return the vector with all Epistasis values
     */
    std::vector<double> readEpistasis(const std::string& inputPath);

    std::vector<int> constructPairIndex() const;

    int getPairIndex(const Mutation& a, const Mutation& b) const;

    int getPairIndex(const int p1, const int p2, const int sym1, const int sym2) const;

    /**
     * computes index within vector of position and symbol
     * @param m Mutation = pair of position and symbol
     * @return the index within the vector
     */
    unsigned int getVectorIndex(const Mutation& m) const;

    Mutation getMutationFromVectorIndex(const int index) const;

    // Constructors are private, to prevent client code from creating new
    // instances of this class. The only instance shall be retrieved through the
    // create_instance(...) function.
    FunctionalSequence(const constants::Constants& params)
        : params{params}, kds{drawKdValues()}, epistasis{drawEpistasis()}, pairIndex{constructPairIndex()} {};
    FunctionalSequence(const std::string& inputPath)
        : params{constants::readParameters(inputPath)}, kds{readKdValues(inputPath)},
          epistasis{readEpistasis(inputPath)}, pairIndex{constructPairIndex()} {};

    FunctionalSequence(const FunctionalSequence&) = delete;
    FunctionalSequence(FunctionalSequence&&) = delete;
    void operator=(const FunctionalSequence&) = delete;

  public:
    static FunctionalSequence* create_instance(const constants::Constants& params);
    static FunctionalSequence* create_instance(const std::string& inputPath);
    static FunctionalSequence* get_instance();
    static void release_instance();
    bool checkValidReadValues();

    /**
     * @return Vector of Kd values for or all possible mutations in ascending order.
     */
    const std::vector<double>& getKd() const;

    /**
     * @param p Mutated position and symbol
     * @return Kd value for the given mutation.
     */
    const double getKd(const Mutation& p) const;

    /**
     * @return Vector of epistasis values for or all pairwise mutations in ascending order with position i<j.
     */
    const std::vector<double>& getEpistasis() const;

    //
    /**
     * The position of mutations a has to be smaller than the position of mutations b.
     * @param a Mutation 1
     * @param b Mutation 2
     * @return Epistasis value for the given pair of mutations
     */
    const double& getEpistasis(const Mutation& a, const Mutation& b) const;

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
};
#endif // DCABENCHMARK_FUNCTIONALSEQUENCE_HPP
