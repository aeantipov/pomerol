#include "pomerol/Susceptibility.h"

namespace Pomerol {

Susceptibility::Susceptibility(const StatesClassification& S, const Hamiltonian& H,
                               const MonomialOperator& A, const MonomialOperator& B,
                               const DensityMatrix& DM) :
    Thermal(DM.beta), ComputableObject(), S(S), H(H), A(A), B(B), DM(DM)
{
}

Susceptibility::Susceptibility(const Susceptibility& Chi) :
    Thermal(Chi.beta), ComputableObject(Chi), S(Chi.S), H(Chi.H), A(Chi.A), B(Chi.B), DM(Chi.DM),
    Vanishing(Chi.Vanishing), SubtractDisconnected(Chi.SubtractDisconnected), ave_A(Chi.ave_A), ave_B(Chi.ave_B)
{
    for(auto const& p : Chi.parts)
        parts.emplace_back(p);
}

void Susceptibility::prepare()
{
    if(getStatus() >= Prepared) return;

    // Find out non-trivial blocks of A and B.
    MonomialOperator::BlocksBimap const& ANontrivialBlocks = A.getBlockMapping();
    MonomialOperator::BlocksBimap const& BNontrivialBlocks = B.getBlockMapping();

    MonomialOperator::BlocksBimap::left_const_iterator Aiter = ANontrivialBlocks.left.begin();
    MonomialOperator::BlocksBimap::right_const_iterator Biter = BNontrivialBlocks.right.begin();

    while(Aiter != ANontrivialBlocks.left.end() && Biter != BNontrivialBlocks.right.end()){
        // <Aleft|A|Aright><Bleft|B|Bright>
        BlockNumber Aleft = Aiter->first;
        BlockNumber Aright = Aiter->second;
        BlockNumber Bleft = Biter->second;
        BlockNumber Bright = Biter->first;

        // Select a relevant 'world stripe' (sequence of blocks).
        if(Aleft == Bright && Aright == Bleft){
            // check if retained blocks are included. If not, do not push.
            if ( DM.isRetained(Aleft) || DM.isRetained(Aright) )
                parts.emplace_back(
                              (MonomialOperatorPart&)A.getPartFromLeftIndex(Aleft),
                              (MonomialOperatorPart&)B.getPartFromRightIndex(Bright),
                              H.getPart(Aright), H.getPart(Aleft),
                              DM.getPart(Aright), DM.getPart(Aleft));
        }

        unsigned long AleftInt = Aleft;
        unsigned long BrightInt = Bright;

        if(AleftInt <= BrightInt) Aiter++;
        if(AleftInt >= BrightInt) Biter++;
    }

    if(!parts.empty())
        Vanishing = false;

    setStatus(Prepared);
}

void Susceptibility::compute()
{
    if(getStatus() >= Computed) return;
    if(getStatus() < Prepared) prepare();

    if(getStatus() < Computed){
        for(auto & p : parts)
            p.compute();
    }
    setStatus(Computed);
}

void Susceptibility::subtractDisconnected()
{
    EnsembleAverage EA_A(S, H, A, DM);
    EnsembleAverage EA_B(S, H, B, DM);
    subtractDisconnected(EA_A, EA_B);
}

void Susceptibility::subtractDisconnected(ComplexType ave_A, ComplexType ave_B)
{
    SubtractDisconnected = true;
    this->ave_A = ave_A;
    this->ave_B = ave_B;
}

void Susceptibility::subtractDisconnected(EnsembleAverage &EA_A, EnsembleAverage &EA_B)
{
    EA_A.compute();
    EA_B.compute();
    subtractDisconnected(EA_A.getResult(), EA_B.getResult());
}

} // namespace Pomerol
