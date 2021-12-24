//
// Created by Smith, Maureen on 31.05.18.
//

#ifndef Mutation_hpp
#define Mutation_hpp

class Mutation
{
  private:
    // position from 1 .. L
    unsigned int position; // const! not declared here, to be able to put it into a set
    // TODO ändern? von 1 .. Q-1?
    // mutation symbol from 0 .. Q-1 -1
    unsigned int symbol; // const! not declared here, to be able to put it into a set

  public:
    Mutation(const unsigned int p, const unsigned int s);

    const unsigned int getPosition() const;

    const unsigned int getSymbol() const;

    // to be able to build a set of Mutations (where the position is unique), the following operators have to be
    // overloaded.
    bool operator<(const Mutation& mut) const;

    bool operator==(const Mutation& mut) const;
};

#endif // DCABENCHMARK_MUTATION_H
