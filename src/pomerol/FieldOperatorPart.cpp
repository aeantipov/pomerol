#include "pomerol/FieldOperatorPart.h"

using std::stringstream;

namespace Pomerol{

FieldOperatorPart::FieldOperatorPart(
        const IndexClassification &IndexInfo, const StatesClassification &S, const HamiltonianPart &HFrom,  const HamiltonianPart &HTo, ParticleIndex PIndex) :
        ComputableObject(), IndexInfo(IndexInfo), S(S), HFrom(HFrom), HTo(HTo), PIndex(PIndex),
        MatrixElementTolerance(1e-8)
{}

void FieldOperatorPart::compute()
{
    if ( Status >= Computed ) return;
    BlockNumber to = HTo.getBlockNumber();
    BlockNumber from = HFrom.getBlockNumber();

    const std::vector<FockState>& toStates = S.getFockStates(to);
    const std::vector<FockState>& fromStates = S.getFockStates(from);

    MatrixType RightMat(fromStates.size(), fromStates.size());
    MatrixType LeftMat(toStates.size(), fromStates.size());
    RightMat.setZero();
    LeftMat.setZero();

    /* Rotation is done in the following way:
     * C_{nm} = \sum_{lk} U^{+}_{nl} C_{lk} U_{km} = \sum_{lk} U^{*}_{ln}O_{lk}U_{km},
     * where the actual sum starts from k state. Big letters denote global states, smaller - InnerQuantumStates.
     * We use the fact each column of O_{lk} has only one nonzero elements.
     * */
    for (std::vector<FockState>::const_iterator CurrentState = fromStates.begin();
                                                CurrentState < fromStates.end(); CurrentState++) {
	    FockState K=*CurrentState;
        std::map<FockState, ComplexType> result1 = O->actRight(K);
        if (result1.size()) {
            FockState L=result1.begin()->first;
            int sign = int(std::real(result1.begin()->second));
	        if ( L!=ERROR_FOCK_STATE && std::abs(sign)>std::numeric_limits<RealType>::epsilon() ) {
		        InnerQuantumState l=S.getInnerState(L), k=S.getInnerState(K);

                for (InnerQuantumState n=0; n<toStates.size(); n++) {
                    LeftMat(n,k) = std::conj(HTo.getMatrixElement(l,n));
                }

                for (InnerQuantumState m=0; m<fromStates.size(); m++) {
                    RightMat(k,m) = RealType(sign) * HFrom.getMatrixElement(k,m);
                }
	        }
        }
    }

// Workaround for Eigen issue 1224
// https://gitlab.com/libeigen/eigen/-/issues/1224
//
// Affected versions are some betas of 3.3 but not the 3.3 release
#if EIGEN_VERSION_AT_LEAST(3,2,90) && EIGEN_MAJOR_VERSION<3
    elementsRowMajor = MatrixType(LeftMat * RightMat).sparseView(MatrixElementTolerance);
#else
    elementsRowMajor = (LeftMat * RightMat).sparseView(MatrixElementTolerance);
#endif

    elementsColMajor = elementsRowMajor;
    Status = Computed;
}

const ColMajorMatrixType& FieldOperatorPart::getColMajorValue(void) const
{
    return elementsColMajor;
}

const RowMajorMatrixType& FieldOperatorPart::getRowMajorValue(void) const
{
    return elementsRowMajor;
}

void FieldOperatorPart::print_to_screen() const  //print to screen C and CX
{
    BlockNumber to   = HTo.getBlockNumber();
    BlockNumber from = HFrom.getBlockNumber();
    INFO(S.getQuantumNumbers(from) << "->" << S.getQuantumNumbers(to));
    for (size_t P=0; P<elementsColMajor.outerSize(); ++P)
	for (ColMajorMatrixType::InnerIterator it(elementsColMajor,P); it; ++it) {
	    FockState N = S.getFockState(to, it.row());
	    FockState M = S.getFockState(from, it.col());
	    INFO(N <<" " << M << " : " << it.value());
        }
}

BlockNumber FieldOperatorPart::getLeftIndex(void) const
{
    return HTo.getBlockNumber();
}

BlockNumber FieldOperatorPart::getRightIndex(void) const
{
    return HFrom.getBlockNumber();
}

// Specialized methods

AnnihilationOperatorPart::AnnihilationOperatorPart(const IndexClassification &IndexInfo, const StatesClassification &S,
                                                  const HamiltonianPart &HFrom, const HamiltonianPart &HTo, ParticleIndex PIndex) :
    FieldOperatorPart(IndexInfo,S,HFrom,HTo,PIndex)
{
    O = new Pomerol::OperatorPresets::C(PIndex);
}

CreationOperatorPart::CreationOperatorPart(const IndexClassification &IndexInfo, const StatesClassification &S,
                                                  const HamiltonianPart &HFrom, const HamiltonianPart &HTo, ParticleIndex PIndex) :
    FieldOperatorPart(IndexInfo,S,HFrom,HTo,PIndex)
{
    O = new Pomerol::OperatorPresets::Cdag(PIndex);
}

const CreationOperatorPart& AnnihilationOperatorPart::transpose() const
{
    CreationOperatorPart *CX = new CreationOperatorPart(IndexInfo, S, HTo, HFrom, PIndex); // swapped h_to and h_from
    CX->elementsRowMajor = elementsRowMajor.transpose();
    CX->elementsColMajor = elementsColMajor.transpose();
    return *CX;
}

const AnnihilationOperatorPart& CreationOperatorPart::transpose() const
{
    AnnihilationOperatorPart *C = new AnnihilationOperatorPart(IndexInfo, S, HTo, HFrom, PIndex); // swapped h_to and h_from
    C->elementsRowMajor = elementsRowMajor.transpose();
    C->elementsColMajor = elementsColMajor.transpose();
    return *C;
}

QuadraticOperatorPart::QuadraticOperatorPart(const IndexClassification &IndexInfo, const StatesClassification &S,
                                             const HamiltonianPart &HFrom, const HamiltonianPart &HTo,
                                             ParticleIndex PIndex1, ParticleIndex PIndex2) :
        FieldOperatorPart(IndexInfo,S,HFrom,HTo,9999), Index1(PIndex1), Index2(PIndex2)
        // Index=9999 dummy
{
    O = new Pomerol::OperatorPresets::N_offdiag(PIndex1, PIndex2);
}

} // end of namespace Pomerol
