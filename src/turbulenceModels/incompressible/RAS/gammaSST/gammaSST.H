/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 Felix Langfeldt
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::incompressible::RASModels::gammaSST

Description
    Implementation of the gamma-SST turbulence model for predicting
    laminar-turbulent transition on arbitrary geometries and grid topo-
    logies.

    The implementation corresponds to the paper  F.R. MENTER et al. (2015):
    "A One-equation Local Correlation-Based Transition Model", Flow,
    Turbulence and Combustion 95(4), pp. 583-619.

    WARNING:
    This implementation is complete and ready to use, but it is not fully
    validated, especially for cases where nonzero pressure gradients occur!
    This implementation is FOR TESTING ONLY and should not be used for
    simulations where prediction of laminar-turbulent transition is critical.

    Please contact me if you find any bugs or errors, I will be grateful to
    further improve this turbulence model: felix.langfeldt@gmail.com


    VERSION HISTORY:
    2016-02-04 : - first public release


    The default model coefficients correspond to the following:
    @verbatim
	gammaSSTCoeffs
	{
	    Flength         100;
	    ca2             0.06;
	    ce2             50;
	    sigmaf          1;
	    CTU1            100;
	    CTU2            1000;
	    CTU3            1;
	    CPG1            14.68;
	    CPG2            -7.34;
	    CPG3            0;
	    CPG1lim         1.5;
	    CPG2lim         3;
	    Ck              1;
	    CSEP            1;
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
	    kInf            0;
	    omegaInf        0;
	}
    @endverbatim

SourceFiles
    gammaSST.C

\*---------------------------------------------------------------------------*/

#ifndef gammaSST_H
#define gammaSST_H

#include "RASModel.H"
#include "wallDist.H"
#include "reflectionVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class gammaSST Declaration
\*---------------------------------------------------------------------------*/

class gammaSST
:
    public RASModel
{
    // Private data

        // Model coefficients
    
        // gamma-ReThetat-model coefficients

        dimensionedScalar Flength_;

        dimensionedScalar ca2_;
        dimensionedScalar ce2_;

        dimensionedScalar sigmaf_;
        
	dimensionedScalar CTU1_;
	dimensionedScalar CTU2_;
	dimensionedScalar CTU3_;
        
	dimensionedScalar CPG1_;
	dimensionedScalar CPG2_;
	dimensionedScalar CPG3_;
	dimensionedScalar CPG1lim_;
	dimensionedScalar CPG2lim_;

	dimensionedScalar Ck_;
	dimensionedScalar CSEP_;

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

        // coefficients for k-omega-SST with controlled decay

        dimensionedScalar kInf_;
        dimensionedScalar omegaInf_;


        //- Wall distance field
        //  Note: different to wall distance in parent RASModel
        wallDist y_;
	reflectionVectors yr_;

        // Fields

        volScalarField gamma_;
        volScalarField k_;
        volScalarField omega_;
        volScalarField nut_;


    // Private member functions

        // Empirical correlations
    
        volScalarField ReThetac() const;

        // Model functions of gamma-ReThetat-model
    
        tmp<volScalarField> Fonset() const;
        tmp<volScalarField> Fonset1() const;
        tmp<volScalarField> Fturb() const;

        tmp<volScalarField> PkLim() const;


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
    TypeName("gammaSST");


    // Constructors

        //- Construct from components
        gammaSST
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~gammaSST()
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

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

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
