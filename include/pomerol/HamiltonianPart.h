/** \file include/pomerol/HamiltonianPart.h
** \brief Declaration of the HamiltonianPart class.
**
** \author Andrey Antipov (Andrey.E.Antipov@gmail.com)
** \author Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/

#ifndef __INCLUDE_HAMILTONIANPART_H
#define __INCLUDE_HAMILTONIANPART_H

#include "Misc.h"
#include "IndexClassification.h"
#include "StatesClassification.h"

#ifdef ENABLE_SAVE_PLAINTEXT
#include <boost/filesystem/path.hpp>
#endif

namespace Pomerol{

template<bool Complex> class Hamiltonian;

/** HamiltonianPart is a class, which stores and diagonalizes the block of the Hamiltonian, which corresponds to a set of given quantum numbers. */
template<bool Complex = false>
class HamiltonianPart : public ComputableObject {

public:

  using MelemT = MelemType<Complex>;
  using MatrixT = MatrixType<Complex>;

private:

    /** A reference to the IndexClassification object. */
    const IndexClassification<Complex> &IndexInfo;
    /** A reference to the IndexHamiltonian object. */
    const IndexHamiltonian<Complex> &F;
    /** A reference to the StateClassification object. */
    const StatesClassification<Complex> &S;

    /** The number of Block. Defined in StatesClassification. */
    BlockNumber Block;
    /** QuantumNumbers of the block. Consructed in Symmetrizer and defined in StatesClassification. */
    QuantumNumbers<Complex> QN;

    /** A matrix filled with matrix elements of HamiltonianPart in the space of FockState's.
     *  After diagonalization it stores the eigenfunctions of the problem in a rows of H. */
    MatrixT H;
    /** A vector of eigenvalues of the HamiltonianPart. */
    RealVectorType Eigenvalues;

    friend class Hamiltonian<Complex>;

public:

    /** Constructor.
     * \param[in] IndexInfo IndexClassification object. Provides information about the indices in the problem.
     * \param[in] F IndexHamiltonian object. Provides all Terms required to fill the Hamiltonian part.
     * \param[in] S StatesClassification object. Provides information about Fock States of the problem.
     * \param[in] Block The BlockNumber of current part. It is a genuine id of the part. */
    HamiltonianPart(const IndexClassification<Complex> &IndexInfo,
                    const IndexHamiltonian<Complex> &F,
                    const StatesClassification<Complex> &S,
                    const BlockNumber& Block);

    /** Fill in the H matrix. */
    void prepare(void);
    /** Diagonalize the H matrix and get EigenValues. */
    void compute(void);

    bool reduce(RealType ActualCutoff); // Useless now

    /** Return the total dimensionality of the H matrix. This corresponds to the one in StatesClassfication. */
    InnerQuantumState getSize(void) const;

    /** Get the matrix element of the HamiltonianPart by the number of states inside the part. */
    MelemT getMatrixElement(InnerQuantumState m, InnerQuantumState n) const; //return H(m,n)
    /** Get the matrix element of the Hamiltonian within two given FockStates. */
    MelemT getMatrixElement(FockState m, FockState n) const; //return H(m,n)

    /** Get the eigenvalue of the H matrix.
     * \param[in] Number of eigenvalue. */
    RealType getEigenValue(InnerQuantumState state) const;

    /** Returns calculated eigenvalues. */
    const RealVectorType& getEigenValues() const;

    /** Return the hamiltonian part matrix. */
    const MatrixT& getMatrix() const;

    /** Return the lowest Eigenvalue of the current part. */
    RealType getMinimumEigenvalue() const;
    /** Return the eigenstate of the H matrix.
     * \param[in] Number of eigenvalue. */
    VectorType<Complex> getEigenState(InnerQuantumState state) const;

    /** Return the QuantumNumbers associated with the Hamiltonian part. */
    QuantumNumbers<Complex> getQuantumNumbers() const;
    /** Return the BlockNumber associated with the Hamiltonian part. */
    BlockNumber getBlockNumber() const;

    /** Print the part of hamiltonian to screen. */
    void print_to_screen() const;

    /** Save the data to the file.
     * \param[in] path Path to the file.
     */
    #ifdef ENABLE_SAVE_PLAINTEXT
    bool savetxt(const boost::filesystem::path &path1);
    #endif
};

extern template class HamiltonianPart<false>;
extern template class HamiltonianPart<true>;

} // end of namespace Pomerol
#endif // endif :: #ifndef __INCLUDE_HAMILTONIANPART_H
