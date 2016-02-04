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

\*---------------------------------------------------------------------------*/

#include "gammaSST.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gammaSST, 0);
addToRunTimeSelectionTable(RASModel, gammaSST, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

volScalarField gammaSST::ReThetac() const
{
    volScalarField ReThetac 
        (
            IOobject
            (
                "ReThetac",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
           gamma_ 
    );
    
    volScalarField FPG 
        (
            IOobject
            (
                "FPG",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
           gamma_ 
    );
    
    volScalarField TuL = min(scalar(100)*sqrt(k_/scalar(1.5))/(omega_*y_),scalar(100));
    volScalarField dVdy = fvc::grad(yr_.n() & U_) & yr_.n();
    volScalarField lamThL = min(
        max(scalar(-7.57e-3)*dVdy*sqr(y_)/nu()+scalar(0.0128),scalar(-1)),
        scalar(1)
    );

    forAll(FPG, cellI)
    {
        if(lamThL[cellI]>=scalar(0))
	{
	    FPG[cellI] = min(
	        scalar(1.0)+CPG1_.value()*lamThL[cellI],
	        CPG1lim_.value()
	    );
	}
	else
	{
	    FPG[cellI] = min(
	        scalar(1.0)+CPG2_.value()*lamThL[cellI]+
		CPG3_.value()*min(lamThL[cellI]+scalar(0.0681),0),
	        CPG2lim_.value()
	    );
	}
    }

    FPG = max(FPG, scalar(0));

    ReThetac = CTU1_+CTU2_*exp(-CTU3_*TuL*FPG);
    
    return ReThetac;
}

tmp<volScalarField> gammaSST::Fonset() const
{
    return tmp<volScalarField>
    (
        new volScalarField
            (
                IOobject
                (
                "Fonset",
                    runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
        max
        (
            min(Fonset1(),scalar(2))-max(scalar(1)-pow3(Rt()/scalar(3.5)),scalar(0)),
            scalar(0)
        )
        )
    );
}

tmp<volScalarField> gammaSST::Fonset1() const
{
    return sqr(y_)*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))/(scalar(2.2)*nu()*ReThetac());
}

tmp<volScalarField> gammaSST::PkLim() const
{
    volScalarField S = sqrt(scalar(2))*mag(skew(fvc::grad(U_)));
    volScalarField O = sqrt(scalar(2))*mag(skew(fvc::grad(U_)));
    volScalarField FonLim = min(
        max(sqr(y_)*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))/(scalar(2.2)*nu()*scalar(1100))-scalar(1.0),
	    scalar(0)
	), scalar(3)
    );
    return scalar(5)*Ck_*max(gamma_-scalar(0.2),scalar(0))*(1-gamma_)*FonLim*
	   nu()*max(scalar(3)*CSEP_-nut_/nu(),scalar(0))*S*O;
}

tmp<volScalarField> gammaSST::Fturb() const
{
    return exp(-pow4(Rt()/scalar(2)));
}

tmp<volScalarField> gammaSST::F1(const volScalarField& CDkOmega) const
{
    volScalarField CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    volScalarField arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*nu()/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    // Modified blending function!
    return 
    max(
        tanh(pow4(arg1)),
    exp(-sqr(pow4(y_*sqrt(k_)/(scalar(120)*nu()))))
    );
}

tmp<volScalarField> gammaSST::F2() const
{
    volScalarField arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*nu()/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gammaSST::gammaSST
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, lamTransportModel, turbulenceModelName),

    Flength_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Flength",
            coeffDict_,
            100.0
        )
    ),
    ca2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ca2",
            coeffDict_,
            0.06
        )
    ),
    ce2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce2",
            coeffDict_,
            50.0
        )
    ),
    sigmaf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaf",
            coeffDict_,
            1.0
        )
    ),
    CTU1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTU1",
            coeffDict_,
            100.0
        )
    ),
    CTU2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTU2",
            coeffDict_,
            1000.0
        )
    ),
    CTU3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTU3",
            coeffDict_,
            1.0
        )
    ),
    CPG1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG1",
            coeffDict_,
            14.68
        )
    ),
    CPG2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG2",
            coeffDict_,
            -7.34
        )
    ),
    CPG3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG3",
            coeffDict_,
            0.0
        )
    ),
    CPG1lim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG1lim",
            coeffDict_,
            1.5
        )
    ),
    CPG2lim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CPG2lim",
            coeffDict_,
            3.0
        )
    ),
    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            coeffDict_,
            1.0
        )
    ),
    CSEP_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CSEP",
            coeffDict_,
            1.0
        )
    ),
    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            coeffDict_,
            0.85034
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            coeffDict_,
            0.85616
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            coeffDict_,
            0.5532
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            coeffDict_,
            0.4403
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
        )
    ),
    kInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kInf",
            coeffDict_,
            0.0,
            sqr(dimLength/dimTime)
        )
    ),
    omegaInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "omegaInf",
            coeffDict_,
            0.0,
            dimless/dimTime
        )
    ),

    y_(mesh_),

    yr_(mesh_),

    gamma_
    (
        IOobject
        (
            "gamma",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh_
    ),
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh_
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh_
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    )
{
    nut_ = a1_*k_/max(a1_*omega_, F2()*sqrt(scalar(2))*mag(symm(fvc::grad(U_))));
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> gammaSST::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> gammaSST::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> gammaSST::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> gammaSST::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool gammaSST::read()
{
    if (RASModel::read())
    {
        Flength_.readIfPresent(coeffDict());
        ca2_.readIfPresent(coeffDict());
        ce2_.readIfPresent(coeffDict());
        sigmaf_.readIfPresent(coeffDict());
        CTU1_.readIfPresent(coeffDict());
        CTU2_.readIfPresent(coeffDict());
        CTU3_.readIfPresent(coeffDict());
        Ck_.readIfPresent(coeffDict());
        CSEP_.readIfPresent(coeffDict());
        alphaK1_.readIfPresent(coeffDict());
        alphaK2_.readIfPresent(coeffDict());
        alphaOmega1_.readIfPresent(coeffDict());
        alphaOmega2_.readIfPresent(coeffDict());
        gamma1_.readIfPresent(coeffDict());
        gamma2_.readIfPresent(coeffDict());
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        betaStar_.readIfPresent(coeffDict());
        a1_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());
        kInf_.readIfPresent(coeffDict());
        omegaInf_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void gammaSST::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
        yr_.correct();
    }

    volScalarField S = mag(symm(fvc::grad(U_)));
    volScalarField O = mag(skew(fvc::grad(U_)));
    volScalarField G(GName(), nut_*2*S*O);

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    volScalarField CDkOmega =
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_;

    volScalarField F1 = this->F1(CDkOmega);

    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::laplacian(DomegaEff(F1), omega_)
     ==
        gamma(F1)*2*S*O
      - fvm::Sp(beta(F1)*omega_, omega_)
      - fvm::SuSp
        (
            (F1 - scalar(1))*CDkOmega/omega_,
            omega_
        )
      + beta(F1)*sqr(omegaInf_)
    );

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
        G*gamma_
      + PkLim()
      - fvm::Sp(max(gamma_,scalar(0.1))*betaStar_*omega_, k_)
      + betaStar_*omegaInf_*kInf_
    );
    
    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    nut_ = a1_*k_/max(a1_*omega_, F2()*sqrt(scalar(2))*S);
    nut_.correctBoundaryConditions();


    // Intermittency equation
    
    tmp<fvScalarMatrix> gammaEqn
    (
        fvm::ddt(gamma_)
      + fvm::div(phi_, gamma_)
      - fvm::laplacian(DgammaEff(), gamma_)
     ==
        Flength_*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))*Fonset()*gamma_
      - fvm::Sp
        (
        Flength_*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))*Fonset()*gamma_,
        gamma_
        )
      + ca2_*sqrt(scalar(2))*mag(skew(fvc::grad(U_)))*Fturb()*gamma_
      - fvm::Sp
        (
            ce2_*ca2_*sqrt(scalar(2))*mag(skew(fvc::grad(U_)))*Fturb()*gamma_,
            gamma_
        )
    ); // old equation

    gammaEqn().relax();
    solve(gammaEqn);

    bound(gamma_,scalar(0));

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
