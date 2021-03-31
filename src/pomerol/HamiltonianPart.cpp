#include"pomerol/HamiltonianPart.h"
#include"pomerol/StatesClassification.h"
#include<sstream>
#include<Eigen/Eigenvalues>

#ifdef ENABLE_SAVE_PLAINTEXT
#include<boost/filesystem.hpp>
#include<boost/filesystem/fstream.hpp>
#endif

// class HamiltonianPart

namespace Pomerol{

HamiltonianPart::HamiltonianPart(const IndexClassification& IndexInfo, const IndexHamiltonian &F, const StatesClassification &S, const BlockNumber& Block):
    ComputableObject(),
    IndexInfo(IndexInfo),
    F(F), S(S),
    Block(Block), QN(S.getQuantumNumbers(Block))
{
}

void HamiltonianPart::prepare()
{
    size_t BlockSize = S.getBlockSize(Block);

    H.resize(BlockSize,BlockSize);
    H.setZero();
    std::map<FockState,ComplexType>::const_iterator melem_it;

    for(InnerQuantumState right_st=0; right_st<BlockSize; right_st++)
    {
        FockState ket = S.getFockState(Block,right_st);
        std::map<FockState,ComplexType> mapStates = F.actRight(ket);
        for (melem_it=mapStates.begin(); melem_it!=mapStates.end(); melem_it++) {
            FockState bra = melem_it -> first;
            ComplexType melem = melem_it -> second;
            //DEBUG("<" << bra << "|" << melem << "|" << F << "|" << ket << ">");
            InnerQuantumState left_st = S.getInnerState(bra);
//            if (left_st > right_st) { ERROR("!"); exit(1); };
            H(left_st,right_st) = melem;
        }
    }

//    H.triangularView<Eigen::Lower>() = H.triangularView<Eigen::Upper>().transpose();
//    assert(MatrixType(H.triangularView<Eigen::Lower>()) == MatrixType(H.triangularView<Eigen::Upper>().transpose()));
    assert((H.adjoint() - H).array().abs().maxCoeff() < 100*std::numeric_limits<RealType>::epsilon());
    Status = Prepared;
}

void HamiltonianPart::compute()		//method of diagonalization classificated part of Hamiltonian
{
    if (Status >= Computed) return;
    if (H.rows() == 1) {
        assert (std::abs(H(0,0) - std::real(H(0,0))) < std::numeric_limits<RealType>::epsilon());
        Eigenvalues.resize(1);
        Eigenvalues << std::real(H(0,0));
        H(0,0) = 1;
        }
    else {
	    Eigen::SelfAdjointEigenSolver<MatrixType> Solver(H,Eigen::ComputeEigenvectors);
	    H = Solver.eigenvectors();
	    Eigenvalues = Solver.eigenvalues();	// eigenvectors are ready
    }
    Status = Computed;
}


ComplexType HamiltonianPart::getMatrixElement(InnerQuantumState m, InnerQuantumState n) const	//return  H(m,n)
{
    return H(m,n);
}

RealType HamiltonianPart::getEigenValue(InnerQuantumState state) const // return Eigenvalues(state)
{
    if ( Status < Computed ) throw (exStatusMismatch());
    return Eigenvalues(state);
}

const RealVectorType& HamiltonianPart::getEigenValues() const
{
    if ( Status < Computed ) throw (exStatusMismatch());
    return Eigenvalues;
}

InnerQuantumState HamiltonianPart::getSize(void) const
{
    return S.getBlockSize(Block);
}

BlockNumber HamiltonianPart::getBlockNumber(void) const
{
    return S.getBlockNumber(QN);
}

QuantumNumbers HamiltonianPart::getQuantumNumbers(void) const
{
    return QN;
}



void HamiltonianPart::print_to_screen() const
{
    INFO(H << std::endl);
}

const MatrixType& HamiltonianPart::getMatrix() const
{
    return H;
}

VectorType HamiltonianPart::getEigenState(InnerQuantumState state) const
{
    if ( Status < Computed ) throw (exStatusMismatch());
    return H.col(state);
}

RealType HamiltonianPart::getMinimumEigenvalue() const
{
    if ( Status < Computed ) throw (exStatusMismatch());
    return Eigenvalues.minCoeff();
}

bool HamiltonianPart::reduce(RealType ActualCutoff)
{
    if ( Status < Computed ) throw (exStatusMismatch());
    InnerQuantumState counter=0;
    for (counter=0; (counter< (unsigned int)Eigenvalues.size() && Eigenvalues[counter]<=ActualCutoff); ++counter){};
    std::cout << "Left " << counter << " eigenvalues : " << std::endl;
    if (counter)
	{std::cout << Eigenvalues.head(counter) << std::endl << "_________" << std::endl;
	Eigenvalues = Eigenvalues.head(counter);
	H = H.topLeftCorner(counter,counter);
	return true;
    }
    else return false;
}

#ifdef ENABLE_SAVE_PLAINTEXT
bool HamiltonianPart::savetxt(const boost::filesystem::path &path1)
{
    boost::filesystem::create_directory(path1);
    boost::filesystem::fstream out;
    if (Status >= Computed) {
        out.open(path1 / boost::filesystem::path("evals.dat"),std::ios_base::out);
        out << Eigenvalues << std::endl;
        out.close();
        out.open(path1 / boost::filesystem::path("evals_shift.dat"),std::ios_base::out);
        out << __num_format<RealVectorType>(Eigenvalues - RealMatrixType::Identity(Eigenvalues.size(),Eigenvalues.size()).diagonal()*getMinimumEigenvalue()) << std::endl;
        out.close();
        };
    if (Status >= Prepared) {
        out.open(path1 / boost::filesystem::path("evecs.dat"),std::ios_base::out);
        out << H << std::endl;
        out.close();
        };
    out.open(path1 / boost::filesystem::path("info.dat"),std::ios_base::out);
    out << "Quantum numbers: " << QN << std::endl;
    out << "Block number:    " << Block << std::endl;
    out.close();
    return true;
}
#endif

} // end of namespace Pomerol

