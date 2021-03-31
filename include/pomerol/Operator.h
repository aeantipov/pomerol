/** \file Operator.h
**  \brief Declaration of the Operator, Operator::Term classes
**
**  \author    Igor Krivenko (Igor.S.Krivenko@gmail.com)
*/


#ifndef POMEROL_OPERATOR_H_uie67d
#define POMEROL_OPERATOR_H_uie67d

#include <ostream>
#include <istream>
#include <tuple>
#include <vector>
#include <map>
#include <limits>
#include <cmath>
#include <iterator>
#include <boost/operators.hpp>

#include "Misc.h"
#include "Index.h"
#include "IndexClassification.h"
#include "Lattice.h"

namespace Pomerol {

class Operator;

namespace OperatorPresets {
    Operator c(ParticleIndex);
    Operator c_dag(ParticleIndex);
    Operator n(ParticleIndex);
    Operator n_offdiag(ParticleIndex, ParticleIndex);
};

class Operator :
    boost::addable<Operator,
    boost::subtractable<Operator,
    boost::multipliable<Operator,
    boost::addable2<Operator, ComplexType,
    boost::subtractable2<Operator, ComplexType,
    boost::multipliable2<Operator, ComplexType
    > > > > > >
{

public:

    Operator(){};
    Operator(Operator const& in):monomials(in.monomials){};
    //Operator(Operator &&);
    Operator& operator=(Operator const & in){monomials = in.monomials; return *this;};
/*
#ifndef TRIQS_WORKAROUND_INTEL_COMPILER_BUGS
    Operator& operator=(Operator && o);
#else
    Operator& operator=(Operator && o) noexcept
        { std::swap(monomials,o.monomials); return *this; }
#endif
*/
    // Type of a fundamental operator
    enum op_type  {creation, annihilation };
    static const int create_annihilate = 0;  // to use std::get<create_annihilate>(...)

    // A composite index labels one fundamental operator. It is LessThanComparable.
    typedef std::tuple<op_type, ParticleIndex> composite_index_t;

    friend std::ostream& operator<<(std::ostream &os, composite_index_t i)
    {
        if(std::get<create_annihilate>(i) == creation) os << "^+";
        os << "(";
        os << std::get<1>(i);
        os << ")";
        return os;
    }

    // Monomial: an ordered set of creation/annihilation operators
    // NB: std::vector supports lexicographical less-than comparison using std::lexicographical_compare
    typedef std::vector<composite_index_t> monomial_t;

    friend std::ostream& operator<<(std::ostream &os, monomial_t const& m)
    {
        for(auto const& c : m)
            os << "C" << c;
        return os;
    }

    friend
    bool operator<(monomial_t const& m1, monomial_t const& m2)
    {
        return m1.size() != m2.size() ? m1.size() < m2.size() :
               std::lexicographical_compare(m1.begin(),m1.end(),m2.begin(),m2.end());
    }

    friend
    bool operator==(monomial_t const& m1, monomial_t const& m2);

    // Map of all monomials with coefficients
    typedef std::map<monomial_t,ComplexType> monomials_map_t;
    // Print Operator itself
    friend std::ostream& operator<<(std::ostream& os, Operator const& op)
    {
        if(op.monomials.size() != 0){
            bool print_plus = false;
            for(auto const& m : op.monomials) {
                os << (print_plus ? " + " : "" ) << m.second;
                if(m.first.size()) os << "*";
                os << m.first;
                print_plus = true;
            }
        } else
            os << "0";
        return os;
    }

    // Iterators (only const!)
    typedef monomials_map_t::const_iterator const_iterator;
    const_iterator begin() const { return monomials.begin(); }
    const_iterator end() const { return monomials.end(); }
    //const_iterator cbegin() const { return monomials.cbegin(); }
    //const_iterator cend() const { return monomials.cend(); }

    // Algebraic operations involving ComplexType constants
    Operator operator-() const
    {
        Operator tmp(*this);
        for(auto & m : tmp.monomials)
            m.second = -m.second;
        return tmp;
    }

    Operator& operator+=(const ComplexType alpha)
    {
        bool is_new_monomial;
        monomials_map_t::iterator it;
        std::tie(it,is_new_monomial) = monomials.insert(std::make_pair(monomial_t(0),alpha));
        if(!is_new_monomial){
            it->second += alpha;
            erase_zero_monomial(monomials,it);
        }
        return *this;
    }

    Operator& operator-=(const ComplexType alpha)
    {
        bool is_new_monomial;
        monomials_map_t::iterator it;
        std::tie(it,is_new_monomial) = monomials.insert(std::make_pair(monomial_t(0),-alpha));
        if(!is_new_monomial){
            it->second -= alpha;
            erase_zero_monomial(monomials,it);
        }
        return *this;
    }

    friend
    Operator operator-(const ComplexType alpha, Operator const& op)
    {
        return -op + alpha;
    }

    Operator& operator*= (const ComplexType alpha)
    {
        if(std::abs(alpha) < 100*std::numeric_limits<RealType>::epsilon()){
            monomials.clear();
        } else {
            for(auto & m : monomials) m.second *= alpha;
        }
        return *this;
    }

    // Algebraic operations
    Operator& operator+=(Operator const& op)
    {
        bool is_new_monomial;
        monomials_map_t::iterator it;
        for(auto const& m : op.monomials) {
            std::tie(it,is_new_monomial) = monomials.insert(m);
            if(!is_new_monomial){
                it->second += m.second;
                erase_zero_monomial(monomials,it);
            }
        }
        return *this;
    }

    Operator& operator-=(Operator const& op)
    {
        bool is_new_monomial;
        monomials_map_t::iterator it;
        for(auto const& m : op.monomials) {
            std::tie(it,is_new_monomial) = monomials.insert(std::make_pair(m.first,-m.second));
            if(!is_new_monomial){
                it->second -= m.second;
                erase_zero_monomial(monomials,it);
            }
        }
        return *this;
    }

    Operator& operator*=(Operator const& op)
    {
        monomials_map_t tmp_map; // product will be stored here
        for(auto const& m : monomials)
            for(auto const& op_m : op.monomials) {
                // prepare an unnormalized product
                monomial_t product_m;
                product_m.reserve(m.first.size() + op_m.first.size());
                std::copy(m.first.begin(), m.first.end(), std::back_inserter(product_m));
                std::copy(op_m.first.begin(), op_m.first.end(), std::back_inserter(product_m));

                normalize_and_insert(product_m, m.second*op_m.second, tmp_map);
            }

        std::swap(monomials, tmp_map);
        return *this;
    }

    bool isEmpty() const {return (monomials.size()==0);};

      /** Returns a matrix element of the operator.
     * \param[in] bra A state to the left of the operator.
     * \param[in] ket A state to the right of the operator.
     * \param[out] Resulting matrix element.
     */
    virtual ComplexType getMatrixElement(const FockState &bra, const FockState &ket) const;

    /** Returns the matrix element of an operator between two states represented by a linear combination of FockState's. */
    virtual ComplexType getMatrixElement( const VectorType & bra, const VectorType &ket, const std::vector<FockState> &states) const;

    /** Returns a result of acting of an operator on a state to the right of the operator.
     * \param[in] ket A state to act on.
     * \param[out] A map of states and corresponding matrix elements, which are the result of an action.
     */
    static std::tuple<FockState,ComplexType> actRight(const monomial_t &in, const FockState &ket);
    virtual std::map<FockState, ComplexType> actRight(const FockState &ket) const;

    /** Returns an operator that is a commutator of the current operator and another one
     * \param[in] rhs An operator to calculate a commutator with.
     * \param[out] Resulting operator.
     */
    Operator getCommutator(const Operator &rhs) const;

    /** Returns an operator that is an anticommutator of the current operator and another one
     * \param[in] rhs An operator to calculate an anticommutator with.
     * \param[out] Resulting operator.
     */
    Operator getAntiCommutator(const Operator &rhs) const;

    /** Checks if current operator commutes with a given one.
     * \param[in] rhs An operator to calculate a commutator with.
    */
    bool commutes(const Operator &rhs) const;


    friend
    bool operator==(const Operator &lhs, const Operator &rhs);


    // Free factory functions
    friend Operator Pomerol::OperatorPresets::c(ParticleIndex);
    friend Operator Pomerol::OperatorPresets::c_dag(ParticleIndex);
    friend Operator Pomerol::OperatorPresets::n(ParticleIndex);
    friend Operator Pomerol::OperatorPresets::n_offdiag(ParticleIndex, ParticleIndex);

protected:

    // Use a template parameter instead of std::complex<double>
    monomials_map_t monomials;

    // Normalize a monomial and insert into a map
    static void normalize_and_insert(monomial_t & m, ComplexType coeff, monomials_map_t & target)
    {
        // The normalization is done by employing a simple bubble sort algorithms.
        // Apart from sorting elements this function keeps track of the sign and
        // recursively calls itself if a permutation of two operators produces a new
        // monomial
        if(m.size() >= 2){
            bool is_swapped;
            do {
                is_swapped = false;
                for (std::size_t n = 1; n < m.size(); ++n){
                    composite_index_t & prev_index = m[n-1];
                    composite_index_t & cur_index = m[n];
                    if(prev_index == cur_index) return;   // The monomial is effectively zero
                    if(prev_index > cur_index){
                        // Are we swapping C and C^+ with the same indices?
                        // A bit ugly ...
                        composite_index_t cur_index_flipped_type(cur_index);
                        std::get<create_annihilate>(cur_index_flipped_type) =
                            op_type(!bool(std::get<create_annihilate>(cur_index_flipped_type)));
                        if(prev_index == cur_index_flipped_type){
                            monomial_t new_m;
                            new_m.reserve(m.size() - 2);
                            std::copy(m.begin(), m.begin() + n-1, std::back_inserter(new_m));
                            std::copy(m.begin() + n+1, m.end(), std::back_inserter(new_m));

                            normalize_and_insert(new_m, coeff, target);
                        }
                        coeff = -coeff;
                        std::swap(prev_index, cur_index);
                        is_swapped = true;
                    }
                }
            } while(is_swapped);
        }

        // Insert the result
        bool is_new_monomial;
        monomials_map_t::iterator it;
        std::tie(it,is_new_monomial) = target.insert(std::make_pair(m, coeff));
        if(!is_new_monomial){
            it->second += coeff;
            erase_zero_monomial(target,it);
        }
    }

    // Erase a monomial with a close-to-zero coefficient.
    static void erase_zero_monomial(monomials_map_t & m,
                                    monomials_map_t::iterator & it)
    {
        if(std::abs(it->second) < 100*std::numeric_limits<RealType>::epsilon())
            m.erase(it);
    }

    public:
    void printAllTerms() const { std::cout << (*this) << std::endl; };
    /** Exception - wrong operation with labels. */
    class exWrongLabel : public std::exception { virtual const char* what() const throw(); };
    /** Exception - Matrix element of term vanishes. */
    class exMelemVanishes : public std::exception { virtual const char* what() const throw(); };

};

// Free functions to make creation/annihilation operators
namespace OperatorPresets {
inline Operator c(ParticleIndex index) {
    typedef Operator c_t;

    c_t tmp;
    c_t::monomial_t m; m.push_back(std::make_tuple(c_t::annihilation, index));
    tmp.monomials.insert(std::make_pair(m,1.0));
    return tmp;
}

inline Operator c_dag(ParticleIndex index) {
    typedef Operator c_dag_t;

    c_dag_t tmp;
    c_dag_t::monomial_t m;
    m.push_back(std::make_tuple(c_dag_t::creation, index));
    tmp.monomials.insert(std::make_pair(m,1.0));
    return tmp;
}

inline Operator n(ParticleIndex index) {
    typedef Operator n_t;

    n_t tmp;
    n_t::monomial_t m;
    m.push_back(std::make_tuple(n_t::creation, index));
    m.push_back(std::make_tuple(n_t::annihilation, index));
    tmp.monomials.insert(std::make_pair(m,1.0));

    return tmp;
}

inline Operator n_offdiag(ParticleIndex index1, ParticleIndex index2) {
    typedef Operator n_t;

    n_t tmp;
    n_t::monomial_t m;
    m.push_back(std::make_tuple(n_t::creation, index1));
    m.push_back(std::make_tuple(n_t::annihilation, index2));
    tmp.monomials.insert(std::make_pair(m,1.0));

    return tmp;
}

} // end of namespace Pomerol::OperatorPresets
} // end of namespace Pomerol

#endif
