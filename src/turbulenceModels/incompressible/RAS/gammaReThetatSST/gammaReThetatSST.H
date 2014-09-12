/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::incompressible::RASModels::gammaReThetatSST

Description
    Implementation of the gamma-ReTheta turbulence model for predicting
    laminar-turbulent transition on arbitrary geometries and grid topo-
    logies.

    The implementation corresponds to the version number CFX-v-1.1 of the
    turbulence model.

    WARNING:
    This implementation is complete and ready to use, but it is not fully
    validated, especially for cases where nonzero pressure gradients occur!
    This implementation is FOR TESTING ONLY and should not be used for
    simulations where prediction of laminar-turbulent transition is critical.

    Please contact me if you find any bugs or errors, I will be grateful to
    further improve this turbulence model: felix.langfeldt@gmail.com

    For information about how to use the turbulence model, what boundary
    conditions to apply etc. please refer to the following literature:

    R.B. LANGTRY and F.R. MENTER (2009): "Correlation-Based Transition
    Modeling for Unstructured Parallelized Computational Fluid Dynamic
    Codes", AIIA Journal 47 (12), pp. 2894-2906

    K. SULUKSNA, P. DECHAUMPHAI, E. JUNTASARO (2009): "Correlations for
    modeling transitional boundary layers under influences of freestream
    turbulence and pressure gradient", Intl. Journal of Heat and Fluid Flow
    30, pp. 66-75

    P. MALAN, K. SULUKSNA, E. JUNTASARO (2009): "Calibrating the gamma-
    ReTheta Transition Model for Commercial CFD", 47th AIAA Aerospace
    Sciences Meeting

    N.N. SORENSEN (2009): "CFD Modelling of Laminar-turbulent Transition
    for Airfoils and Rotors using the gamma-ReTheta Model", Wind Energy 12,
    pp. 715-733


    VERSION HISTORY:
    2011-09-06 : - ensured compatibility with OpenFOAM 2.0.x
    2011-03-10 : - added ReThetac() and Flength() correlations of
                   SORENSEN (2009)
    2011-02-28 : - added missing max() in the return statement of the
                   F1 blending function
		 - reverted ReThetac() to use original correlation of
		   LANGTRY and MENTER
    2011-02-24 : - first public release


    The default model coefficients correspond to the following:
    @verbatim
        gammaReThetatSSTCoeffs
        {
		sigmaThetat     2;
		ca1             2;
		ce1             1;
	        ca2             0.06;
		ce2             50;
		cThetat         0.03;
		sigmaf          1;
		s1              2;
		alphaK1         0.85034;
		alphaK2         1;
		alphaOmega1     0.5;
		alphaOmega2     0.85616;
		gamma1          0.5532;
		gamma2          0.4403;
		beta1           0.075;
		beta2           0.0828;
		betaStar        0.09;
		a1              0.31;
		c1              10;
        }
    @endverbatim

SourceFiles
    gammaReThetatSST.C

\*---------------------------------------------------------------------------*/

#ifndef gammaReThetatSST_H
#define gammaReThetatSST_H

#include "RASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class gammaReThetatSST Declaration
\*---------------------------------------------------------------------------*/

class gammaReThetatSST
:
    public RASModel
{
    // Private data
    
	    // Convergence tolerance for calculating ReThetat
	    // (not ReThetatTilda!!)
	    static const scalar tol_;
	    // Maximum iterations
	    static const int maxIter_;
                        

        // Model coefficients
	
	// gamma-ReThetat-model coefficients
	
	    dimensionedScalar ca1_;
	    dimensionedScalar ce1_;

	    dimensionedScalar ca2_;
	    dimensionedScalar ce2_;

	    dimensionedScalar cThetat_;

	    dimensionedScalar sigmaf_;
	    dimensionedScalar sigmaThetat_;

	    dimensionedScalar s1_;

	// standard k-omega-SST coefficients
	
            dimensionedScalar alphaK1_;
            dimensionedScalar alphaK2_;

            dimensionedScalar alphaOmega1_;
            dimensionedScalar alphaOmega2_;

            dimensionedScalar gamma1_;
            dimensionedScalar gamma2_;

            dimensionedScalar beta1_;
            dimensionedScalar beta2_;

            dimensionedScalar betaStar_;

            dimensionedScalar a1_;
            dimensionedScalar c1_;


        //- Wall distance field
        //  Note: different to wall distance in parent RASModel
        wallDist y_;

        // Fields

            volScalarField gamma_;
            volScalarField ReThetatTilda_;
            volScalarField k_;
            volScalarField omega_;
            volScalarField nut_;


    // Private member functions

        // Empirical correlations
	
	volScalarField Flength() const;
	volScalarField ReThetac() const;

	// Model functions of gamma-ReThetat-model
	
	tmp<volScalarField> Fonset() const;
	tmp<volScalarField> Fonset1() const;
	tmp<volScalarField> Fturb() const;
	tmp<volScalarField> Freattach() const;
	tmp<volScalarField> Fwake() const;
	tmp<volScalarField> FThetat() const;

	// Empirical Functions for calculating the local transition onset
	// momentum thickness Reynolds number
	void ReThetat(volScalarField& ReThetatField) const;
	scalar ReThetatEq(scalar Tu, scalar lambda) const;

	// Intermittency modification for separation-induced transition
	tmp<volScalarField> gammaSep() const;


        tmp<volScalarField> F1(const volScalarField& CDkOmega) const;
        tmp<volScalarField> F2() const;

        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

        tmp<volScalarField> alphaK
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, alphaK1_, alphaK2_);
        }

        tmp<volScalarField> alphaOmega
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, alphaOmega1_, alphaOmega2_);
        }

        tmp<volScalarField> beta
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, beta1_, beta2_);
        }

        tmp<volScalarField> gamma
        (
            const volScalarField& F1
        ) const
        {
            return blend(F1, gamma1_, gamma2_);
        }


public:

    //- Runtime type information
    TypeName("gammaReThetatSST");


    // Constructors

        //- Construct from components
        gammaReThetatSST
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
	    const word& turbulenceModelName = turbulenceModel::typeName,
	    const word& modelName = typeName
        );


    //- Destructor
    virtual ~gammaReThetatSST()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for gamma
        tmp<volScalarField> DgammaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DgammaEff", (nut_/sigmaf_) + nu())
            );
        }

        //- Return the effective diffusivity for ReThetatTilda
        tmp<volScalarField> DReThetatTildaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DReThetatTildaEff", sigmaThetat_*(nut_ + nu()))
            );
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphaK(F1)*nut_ + nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff", alphaOmega(F1)*nut_ + nu())
            );
        }

        //- Return the intermittency
        virtual tmp<volScalarField> gamma() const
        {
            return gamma_;
        }

        //- Return the local transition momentum thickness Reynolds number
        virtual tmp<volScalarField> ReThetatTilda() const
        {
            return ReThetatTilda_;
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    betaStar_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Return the viscosity ratio
        tmp<volScalarField> Rt() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("Rt", k_/(nu()*omega_)) // LANGTRY's definition
            );
        }

	//- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
