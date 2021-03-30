/** \file include/pomerol/StatesClassification.h
** \brief Declaration of BlockNumber and StatesClassification classes
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/


#ifndef __INCLUDE_STATESCLASSIFICATION_H
#define __INCLUDE_STATESCLASSIFICATION_H

#include "Misc.h"
#include "Index.h"
#include "IndexClassification.h"
#include "Operator.h"
#include "Symmetrizer.h"

namespace Pomerol{

/** A small wrapper around int to hold a number of block.
 *  A Block is a sub-matrix of Hamiltonian which is separated from the others ( i.e. the Hamiltonian is block-diagonal). */
struct BlockNumber;

/** InnerQuantumState labels the states inside of the block of Fock States. Has no physical meaning. */
typedef unsigned long InnerQuantumState;

/** This class handles all information about Fock states.
 *  It makes a classification of Fock states into blocks.
 */
template<bool Complex = false>
class StatesClassification : public ComputableObject {
    /** Total number of states = 2^(IndexInfo.size()) */
    unsigned long StateSize;
    /** Total number of modes of the system. Equal to IndexInfo.size() */
    ParticleIndex IndexSize;

    /** A map between all BlockNumbers and QuantumNumbers */
    std::map<BlockNumber,QuantumNumbers<Complex>> BlockToQuantum;
    /** A map between all QuantumNumbers and BlockNumbers */
    std::map<QuantumNumbers<Complex>,BlockNumber> QuantumToBlock;
    /** A storage for all FockStates, each subvector correspond to the states which belong to a block with a given BlockNumber. */
    std::vector<std::vector<FockState> > StatesContainer;
    /** Index all states to belong to a block. */
    std::vector<BlockNumber> StateBlockIndex;

    /** A reference to an IndexClassification object */
    const IndexClassification<Complex> &IndexInfo;
    /** A reference to a Symmetrizer object. This will be used for classification of the states. */
    const Symmetrizer<Complex> &Symm;
public:
    /** Constructor
     * \param[in] IndexInfo A reference to an IndexClassification object
     */
    StatesClassification(const IndexClassification<Complex>& IndexInfo, const Symmetrizer<Complex> &Symm);

    /** Perform a classification of all FockStates */
    void compute();

   /** get total number of Quantum States ( 2^IndexInfo.size() ) */
    const unsigned long getNumberOfStates() const;

    /** get a vector of all FockStates with a given set of QuantumNumbers
     * \param[in] in A set of quantum numbers to get a vector of FockStates
     */
    const std::vector<FockState>& getFockStates( QuantumNumbers<Complex> in ) const;
    const std::vector<FockState>& getFockStates( BlockNumber in ) const;
    const size_t getBlockSize( BlockNumber in ) const;

    /** get a FockState, corresponding to an internal InnerQuantumState
     * \param[in] QuantumNumbers of block in which the InnerQuantumState is located
     * \param[in] m InnerQuantumState for which the correspondence is required
     */
    const FockState getFockState( BlockNumber in, InnerQuantumState m) const;
    const FockState getFockState( QuantumNumbers<Complex> in, InnerQuantumState m) const;
    /** get InnerQuantumState of a given FockState. Since FockState is associated with
     * the Block number no explicit BlockNumber or QuantumNumbers is required
     * \param[in] state FockState for which the correspondence is required
     */
    const InnerQuantumState getInnerState( FockState state) const;
    const InnerQuantumState getInnerState( QuantumState state) const;

    /** Returns a number of Block which corresponds to given Quantum Numbers
     * \param[in] in A set of QuantumNumbers to find corresponding BlockNumber
     */
    BlockNumber getBlockNumber(QuantumNumbers<Complex> in) const;

    /** Returns QuantumNumbers for a given BlockNumber
     * \param[in] in A BlockNumber to find a set of corresponding QuantumNumbers
     */
    QuantumNumbers<Complex> getQuantumNumbers(BlockNumber in) const;
    /** Returns total amount of non-vanishing blocks */
    BlockNumber NumberOfBlocks() const;

    /** Returns QuantumNumbers of a given FockState
     * \param[in] in A FockState for which the QuantumNumbers are requested
     */
    QuantumNumbers<Complex> getQuantumNumbers(FockState in) const;
    QuantumNumbers<Complex> getQuantumNumbers(QuantumState in) const;
    /** Returns BlockNumber of a given FockState
     * \param[in] in A FockState for which the BlockNumber is requested
     */
    BlockNumber getBlockNumber(FockState in) const;
    BlockNumber getBlockNumber(QuantumState in) const;

    /** Checks that a block with a given QuantumNumbers does not vanish
     * \param[in] in A set of QuantumNumbers to check
     */
    //bool checkQuantumNumbers(QuantumNumbers in) const;


    /** Exception - wrong state. */
    class exWrongState : public std::exception { virtual const char* what() const throw(); };
};


/** This class represents a number of a current block of FockStates which have
 * the same quantum numbers. If such block can't exist ( can't be created by anyone else )
 * it's value is assigned to an ERROR_BLOCK_NUMBER.
 * The classification of blocks is now done by StatesClassification class
 */
struct BlockNumber {
    /** The number of current block */
    int number;
    /** Empty constructor */
    BlockNumber(){};
    /** Copy constructor from int number
     * \param[in] number A number to copy from
     */
    BlockNumber(int number):number(number){};
    /** Type conversion to integer */
    operator int() const {return number;}
    /** post-increment operator */
    BlockNumber& operator ++(int unused){number++; return *this;}
    /** Returns true if such block exists */
    bool isCorrect(){return number >= 0;}
    /** Operator < */
    bool operator<(const BlockNumber& rhs) const ;
    /** Operator > */
    bool operator==(const BlockNumber& rhs) const;
};

/** All blocks with this number are treated as nonexistent */
const BlockNumber ERROR_BLOCK_NUMBER = -1;

extern template class StatesClassification<false>;
extern template class StatesClassification<true>;

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_STATESCLASSIFICATION_H
