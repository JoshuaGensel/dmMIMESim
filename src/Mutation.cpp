//
// Created by Smith, Maureen on 31.05.18.
//

#include "Mutation.hpp"

Mutation::Mutation(const unsigned int p, const unsigned int s) : position(p), symbol(s) {}

const unsigned int Mutation::getPosition() const
{
    return position;
}

const unsigned int Mutation::getSymbol() const
{
    return symbol;
}

// TODO: fix symbol comparison
bool Mutation::operator<(const Mutation& mut) const
{
    return (this->getPosition() < mut.getPosition()); // && (this->getSymbol() < mut.getSymbol());
}

bool Mutation::operator==(const Mutation& mut) const
{
    return (this->getPosition() == mut.getPosition()) && (this->getSymbol() == mut.getSymbol());
}