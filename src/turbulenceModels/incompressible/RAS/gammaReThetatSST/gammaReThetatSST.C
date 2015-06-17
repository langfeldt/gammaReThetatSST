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

#include "gammaReThetatSST.H"
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

defineTypeNameAndDebug(gammaReThetatSST, 0);
addToRunTimeSelectionTable(RASModel, gammaReThetatSST, dictionary);

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

// Tolerance and maximum iteration number for calculation of ReThetat
const scalar gammaReThetatSST::tol_ = 1.0e-4;
const int gammaReThetatSST::maxIter_ = 100;

// Selectable correlation names
const word gammaReThetatSST::CORRN_MENTER2009 = "LangtryMenter2009";
const word gammaReThetatSST::CORRN_SULUKSNA2009 = "SuluksnaEtAl2009";
const word gammaReThetatSST::CORRN_MALAN2009 = "MalanEtAl2009";
const word gammaReThetatSST::CORRN_SORENSEN2009 = "Sorensen2009";
const word gammaReThetatSST::CORRN_TOMAC2013 = "TomacEtAl2013";

// Selectable correlation IDs
const int gammaReThetatSST::CORR_MENTER2009 = 0;
const int gammaReThetatSST::CORR_SULUKSNA2009 = 1;
const int gammaReThetatSST::CORR_MALAN2009 = 2;
const int gammaReThetatSST::CORR_SORENSEN2009 = 3;
const int gammaReThetatSST::CORR_TOMAC2013 = 4;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

volScalarField gammaReThetatSST::Flength() const
{
    volScalarField Flength 
        (
            IOobject
            (
                "Flength",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
            ReThetatTilda_ 
    );

    switch (corrID_) {
        case CORR_SULUKSNA2009:
            // CORRELATION: SULUKSNA et al. 2009
            Flength = min(
                scalar(0.1)*exp(scalar(-0.022)*ReThetatTilda_+scalar(12))+scalar(0.45),
                scalar(300)
            );
            break;
        case CORR_MALAN2009:
            // CORRELATION: MALAN et al. 2009
            Flength = min(
                exp(scalar(-0.01173)*ReThetatTilda_+scalar(7.168))+scalar(0.5),
                scalar(300)
            );
            break;
        case CORR_SORENSEN2009: 
            // CORRELATION: SORENSEN 2009
            Flength = min(
                scalar(150)*exp(scalar(-1)*pow(ReThetatTilda_/scalar(120),1.2))+scalar(0.1),
                scalar(30)
            );
            break;
        case CORR_TOMAC2013: 
            // CORRELATION: TOMAC et al. 2013
            Flength = scalar(0.162)+scalar(93.3)*exp(scalar(-1)*sqr(ReThetatTilda_)/scalar(49153))+
		    (scalar(50)/(scalar(260)*sqrt(scalar(6.283))))*exp(scalar(-0.5)*
	            sqr((ReThetatTilda_-scalar(520))/scalar(260)));
            break;
        default:
            // CORRELATION: LANGTRY and MENTER 2009 
            forAll(Flength, cellI)
            {
                if(ReThetatTilda_[cellI] < scalar(400))
                    Flength[cellI] = scalar(398.189e-1)-scalar(119.270e-4)*ReThetatTilda_[cellI]-scalar(132.567e-6)*sqr(ReThetatTilda_[cellI]);
                else if(ReThetatTilda_[cellI] < scalar(596))
                    Flength[cellI] = scalar(263.404)-scalar(123.939e-2)*ReThetatTilda_[cellI]+scalar(194.548e-5)*sqr(ReThetatTilda_[cellI])-scalar(101.695e-8)*pow3(ReThetatTilda_[cellI]);
                else if(ReThetatTilda_[cellI] < scalar(1200))
                    Flength[cellI] = scalar(0.5)-(ReThetatTilda_[cellI]-scalar(596.0))*scalar(3e-4);
                else
                    Flength[cellI] = scalar(0.3188);
            }

            Flength = (Flength*(scalar(1.0)-exp(-sqr(sqr(y_)*omega_/(0.4*500.0*nu()))))+40.0*exp(-sqr(sqr(y_)*omega_/(0.4*500.0*nu()))));
    }

    return Flength;
}

volScalarField gammaReThetatSST::ReThetac() const
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
            ReThetatTilda_ 
    );

    switch (corrID_) {
        case CORR_SULUKSNA2009:
            // CORRELATION: SULUKSNA et al. 2009
            ReThetac = min(
                max(
                scalar(1.47)*ReThetatTilda_-sqr(scalar(0.025)*ReThetatTilda_)-scalar(120),
                scalar(125)
                ),
                ReThetatTilda_
            );
            break;
        case CORR_MALAN2009:
            // CORRELATION: MALAN et al. 2009
            ReThetac = min(
                scalar(0.615)*ReThetatTilda_+scalar(61.5),
                ReThetatTilda_
            );
            break;
        case CORR_SORENSEN2009: 
            // CORRELATION: SORENSEN 2009
            ReThetac = tanh(pow4((ReThetatTilda_-scalar(100))/scalar(400)))*
                   (ReThetatTilda_+scalar(12000))/scalar(25)+
                   (scalar(1)-tanh(pow4((ReThetatTilda_-scalar(100))/scalar(400))))*
                   (scalar(7)*ReThetatTilda_+scalar(100))/scalar(10);
            break;
        case CORR_TOMAC2013: 
            // CORRELATION: TOMAC et al. 2013
            ReThetac = min(
                scalar(0.993)*ReThetatTilda_,
                scalar(0.322)*ReThetatTilda_+(scalar(105900)/(scalar(150)*sqrt(scalar(6.283))))*
		    (exp(scalar(-0.5)*sqr((ReThetatTilda_-scalar(560))/scalar(150)))+
		     exp(scalar(-0.5)*sqr((ReThetatTilda_-scalar(168))/scalar(150))))
            );
            break;
        default:
            // CORRELATION: LANGTRY and MENTER 2009 
            forAll(ReThetac, cellI)
            {
                if(ReThetatTilda_[cellI] > scalar(1870))
                    ReThetac[cellI] = ReThetatTilda_[cellI]-(scalar(593.11)+(ReThetatTilda_[cellI]-scalar(1870.0))*scalar(0.482));
                else
                    ReThetac[cellI] = ReThetatTilda_[cellI]-(scalar(396.035e-2)-scalar(120.656e-4)*ReThetatTilda_[cellI]+scalar(868.230e-6)*sqr(ReThetatTilda_[cellI])-scalar(696.506e-9)*pow3(ReThetatTilda_[cellI])+scalar(174.105e-12)*pow4(ReThetatTilda_[cellI]));
            }
    }
    
    return ReThetac;
    
}

tmp<volScalarField> gammaReThetatSST::Fonset() const
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
            min(max(Fonset1(),pow4(Fonset1())),scalar(2))-max(scalar(1)-pow3(Rt()/scalar(2.5)),scalar(0)),
            scalar(0)
        )
        )
    );
}

tmp<volScalarField> gammaReThetatSST::Fonset1() const
{
    return sqr(y_)*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))/(scalar(2.193)*nu()*ReThetac());
}

tmp<volScalarField> gammaReThetatSST::Fturb() const
{
    return exp(-pow4(Rt()/scalar(4)));
}

tmp<volScalarField> gammaReThetatSST::Freattach() const
{
    return exp(-pow4(Rt()/scalar(20)));
}

tmp<volScalarField> gammaReThetatSST::Fwake() const
{
    return exp(-sqr(sqr(y_)*omega_/(scalar(1.0e5)*nu())));
}

tmp<volScalarField> gammaReThetatSST::FThetat() const
{
        volScalarField magVort = sqrt(scalar(2))*mag(skew(fvc::grad(U_)));
        magVort = max(magVort,dimensionedScalar("smallOmega",magVort.dimensions(),SMALL));
    return min
    (
        max
        (
            Fwake()*exp(-pow4(magSqr(U_)/(scalar(375.0)*nu()*magVort*ReThetatTilda_))),
        scalar(1.0)-sqr((ce2_*gamma_-scalar(1.0))/(ce2_-scalar(1.0)))
        ),
        scalar(1.0)
    );
}


// FULL IMPLEMENTATION according to LANGTRY and MENTER 2009
// SULUKSNA et. al. 2009 suggest fixing lambda = 0 and thus omitting
// pressure gradient influence on ReThetat (to be tested!)
void gammaReThetatSST::ReThetat(volScalarField& ReThetatField) const
{
    scalar Tu, lambda, ReThetatOld, ReThetatNew, ReThetatTol, dUds, K;
    volScalarField U2gradU = (sqr(U_)&&(fvc::grad(U_)));

    forAll(ReThetatField, cellI)
    {
        int iter = 0;
        Tu = max(
            scalar(100)*sqrt(k_[cellI]/scalar(1.5))/max(mag(U_[cellI]),SMALL),
            scalar(0.027)
        );

        dUds = U2gradU[cellI]/(sqr(max(mag(U_[cellI]),SMALL)));

        // Starting value
        ReThetatNew = max(ReThetatEq(Tu, scalar(0), scalar(0)),scalar(20.0));
        ReThetatTol = ReThetatNew*tol_;

        if (dUds_) {
            do
            {
                ReThetatOld = ReThetatNew;
                lambda = max(
                        min(
                        sqr(ReThetatOld)*nu()()[cellI]*dUds/(sqr(max(mag(U_[cellI]),SMALL))),
                    scalar(0.1)
                    ),
                    scalar(-0.1)
                );
                K = max(
                        min(
                        nu()()[cellI]*dUds/(sqr(max(mag(U_[cellI]),SMALL))),
                    scalar(3e-6)
                    ),
                    scalar(-3e-6)
                );
                ReThetatNew = max(ReThetatEq(Tu, lambda, K),scalar(20.0));

                if (iter++ > maxIter_)
                {
                    FatalErrorIn
                    (
                         "gammaReThetatSST::ReThetat(volScalarField& ReThetatField) const"
                    )   << "Maximum number of iterations exceeded"
                        << abort(FatalError);
                }
            } while(mag(ReThetatNew-ReThetatOld) > ReThetatTol);
        }

        ReThetatField[cellI] = ReThetatNew;
    }

}

scalar gammaReThetatSST::ReThetatEq(scalar Tu, scalar lambda, scalar K) const
{
    scalar FTu;
    scalar FlamK;
    switch (corrID_) {
        case CORR_SULUKSNA2009:
        case CORR_SORENSEN2009:
            // "OLD" CORRELATION FROM MENTER ET AL. (2004)
            FTu = scalar(803.73)*pow((Tu+scalar(0.6067)),scalar(-1.027));
            if(lambda > scalar(0)) {
                scalar FK = scalar(0.0962e6)*K+scalar(0.148e12)*sqr(K)+scalar(0.0141e18)*pow3(K);
                FlamK = scalar(1.0)+FK*(scalar(1.0)-exp(-Tu/scalar(1.5)))+scalar(0.556)*(scalar(1.0)-exp(-scalar(23.9)*lambda))*exp(-Tu/scalar(1.5)); 
            }
            else {
                scalar Flam = scalar(10.32)*lambda+scalar(89.47)*sqr(lambda)+265.51*pow3(lambda);
                FlamK = scalar(1.0)+Flam*exp(-Tu/scalar(3.0));
            }
            break;
        default:
            // "NEW" CORRELATION FROM LANGTRY/MENTER (2009)
            if(Tu > scalar(1.3)) {
                FTu = scalar(331.5)*pow((Tu-scalar(0.5658)),scalar(-0.671));
            }
            else {
                FTu = scalar(1173.51)-scalar(589.428)*Tu+scalar(0.2196)/sqr(Tu);
            }
            if(lambda > scalar(0)) {
                FlamK = (scalar(1.0)+scalar(0.275)*(scalar(1.0)-exp(scalar(-35.0)*lambda))*exp(scalar(-2.0)*Tu));
            }
            else {
                FlamK = (scalar(1.0)+(scalar(12.986)*lambda+scalar(123.66)*sqr(lambda)+scalar(405.689)*pow3(lambda))*exp(-pow((Tu/scalar(1.5)),scalar(1.5))));
            }
    }
    return FTu*FlamK;
}

tmp<volScalarField> gammaReThetatSST::gammaSep() const
{
    return FThetat()*min
    (
        s1_*Freattach()*max
        (
            sqr(y_)*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))/(scalar(3.235)*nu()*ReThetac())-scalar(1.0),
        scalar(0.0)
        ),
        scalar(2.0)
    );
}

tmp<volScalarField> gammaReThetatSST::F1(const volScalarField& CDkOmega) const
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

tmp<volScalarField> gammaReThetatSST::F2() const
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

gammaReThetatSST::gammaReThetatSST
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, lamTransportModel, turbulenceModelName),

    ca1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ca1",
            coeffDict_,
            2.0
        )
    ),
    ce1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce1",
            coeffDict_,
            1.0
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
    cThetat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cThetat",
            coeffDict_,
            0.03
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
    sigmaThetat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaThetat",
            coeffDict_,
            2.0
        )
    ),
    s1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "s1",
            coeffDict_,
            2.0
        )
    ),
    dUds_
    (
        Switch::lookupOrAddToDict
        (
            "dUds",
            coeffDict_,
            false
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
    ReThetatTilda_
    (
        IOobject
        (
            "ReThetatTilda",
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
    ),
    convk_
    (
        IOobject
        (
            "convk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", k_.dimensions()/dimTime, 0)
    ),
    diffk_
    (
        IOobject
        (
            "diffk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", k_.dimensions()/dimTime, 0)
    ),
    Pk_
    (
        IOobject
        (
            "Pk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", k_.dimensions()/dimTime, 0)
    ),
    Dk_
    (
        IOobject
        (
            "Dk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", k_.dimensions()/dimTime, 0)
    ),
    convomega_
    (
        IOobject
        (
            "convomega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", omega_.dimensions()/dimTime, 0)
    ),
    diffomega_
    (
        IOobject
        (
            "diffomega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", omega_.dimensions()/dimTime, 0)
    ),
    Pomega_
    (
        IOobject
        (
            "Pomega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", omega_.dimensions()/dimTime, 0)
    ),
    Domega_
    (
        IOobject
        (
            "Domega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", omega_.dimensions()/dimTime, 0)
    ),
    convgamma_
    (
        IOobject
        (
            "convgamma",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", gamma_.dimensions()/dimTime, 0)
    ),
    diffgamma_
    (
        IOobject
        (
            "diffgamma",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", gamma_.dimensions()/dimTime, 0)
    ),
    Pgamma_
    (
        IOobject
        (
            "Pgamma",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", gamma_.dimensions()/dimTime, 0)
    ),
    Dgamma_
    (
        IOobject
        (
            "Dgamma",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", gamma_.dimensions()/dimTime, 0)
    ),
    convReThetatTilda_
    (
        IOobject
        (
            "convReThetatTilda",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", ReThetatTilda_.dimensions()/dimTime, 0)
    ),
    diffReThetatTilda_
    (
        IOobject
        (
            "diffReThetatTilda",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", ReThetatTilda_.dimensions()/dimTime, 0)
    ),
    PReThetatTilda_
    (
        IOobject
        (
            "PReThetatTilda",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", ReThetatTilda_.dimensions()/dimTime, 0)
    )
{
    // get correlations name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word corrName
    (
        IOdictionary
        (
            IOobject
            (
                "RASProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookupOrDefault("gammaReThetatSSTCorrelations",CORRN_MENTER2009)
    );

    // set selected correlation ID
    word corrInfo;
    if (corrName == CORRN_SULUKSNA2009) {
        corrID_ = CORR_SULUKSNA2009;
        corrInfo = "Suluksna et al. (2009)";
    }
    else if (corrName == CORRN_MALAN2009) {
        corrID_ = CORR_MALAN2009;
        corrInfo = "Malan et al. (2009)";
    }
    else if (corrName == CORRN_SORENSEN2009) {
        corrID_ = CORR_SORENSEN2009;
        corrInfo = "Sorensen (2009)";
    }
    else if (corrName == CORRN_TOMAC2013) {
        corrID_ = CORR_TOMAC2013;
        corrInfo = "Tomac et al. (2013)";
    }
    else {
        corrID_ = CORR_MENTER2009;
        corrInfo = "Langtry and Menter (2009)";
    }

    Info << "Using gammaReThetat-correlations by " << corrInfo << endl;

    nut_ = a1_*k_/max(a1_*omega_, F2()*sqrt(scalar(2))*mag(symm(fvc::grad(U_))));
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar gammaReThetatSST::ReThetatTildaInlet(scalar Tu) const
{
    return ReThetatEq(Tu,scalar(0),scalar(0));
}

tmp<volSymmTensorField> gammaReThetatSST::R() const
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


tmp<volSymmTensorField> gammaReThetatSST::devReff() const
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


tmp<fvVectorMatrix> gammaReThetatSST::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> gammaReThetatSST::divDevRhoReff
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


bool gammaReThetatSST::read()
{
    if (RASModel::read())
    {
        ca1_.readIfPresent(coeffDict());
        ce1_.readIfPresent(coeffDict());
        ca2_.readIfPresent(coeffDict());
        ce2_.readIfPresent(coeffDict());
        cThetat_.readIfPresent(coeffDict());
        sigmaf_.readIfPresent(coeffDict());
        sigmaThetat_.readIfPresent(coeffDict());
        s1_.readIfPresent(coeffDict());
        dUds_.readIfPresent("dUds",coeffDict());
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


void gammaReThetatSST::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
    }

    volScalarField S2 = magSqr(symm(fvc::grad(U_)));
    volScalarField G(GName(), nut_*2*S2);

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
        gamma(F1)
       *min(2*S2, (c1_/a1_)*betaStar_*omega_*max(a1_*omega_, F2()*sqrt(scalar(2)*S2)))
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

    convomega_ = fvc::div(phi_, omega_);
    diffomega_ = fvc::laplacian(DomegaEff(F1), omega_);
    Pomega_ = gamma(F1)*min(2*S2, (c1_/a1_)*betaStar_*omega_*max(a1_*omega_, F2()*sqrt(scalar(2)*S2)))
	      - fvc::SuSp
                (
                    (F1 - scalar(1))*CDkOmega/omega_,
                    omega_
                );
    Domega_ = fvc::Sp(beta(F1)*omega_, omega_);


    volScalarField gammaEff = max
    (
        gamma_,
        gammaSep()
    );

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
        min(G, c1_*betaStar_*k_*omega_)*gammaEff
      - fvm::Sp(min(max(gammaEff,scalar(0.1)),scalar(1))*betaStar_*omega_, k_)
      + betaStar_*omegaInf_*kInf_
    );
    
    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);

    convk_ = fvc::div(phi_, k_);
    diffk_ = fvc::laplacian(DkEff(F1), k_);
    Pk_ = min(G, c1_*betaStar_*k_*omega_)*gammaEff;
    Dk_ = fvc::Sp(min(max(gammaEff,scalar(0.1)),scalar(1))*betaStar_*omega_, k_);

    // Re-calculate viscosity
    nut_ = a1_*k_/max(a1_*omega_, F2()*sqrt(scalar(2)*S2));
    nut_.correctBoundaryConditions();


    // local transition onset momentum thickness Reynolds number
    volScalarField ReThetatField
    (
        IOobject
        (
            "ReThetatField",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        ReThetatTilda_
    );
    ReThetat(ReThetatField);
    
    // OUTPUT FUNCTIONS
    if(runTime_.outputTime())
    {
//      this->FThetat(sqr(magU)/max(O,dimensionedScalar("smallOmega",O.dimensions(),SMALL)))().write();
//      ReThetatField.write();
//      Tu.write();
//      dUds.write();
//      Flength().write();
//      Fonset()().write();
//      Fonset1()().write();
//      ReThetac()().write();
//      volScalarField gProd = Flength()*ca1_*mag(symm(fvc::grad(U_)))*sqrt(Fonset()*gamma_)*(scalar(1)-ce1_*gamma_);
//      volScalarField gDest = ca2_*mag(skew(fvc::grad(U_)))*Fturb()*gamma_*(ce2_*gamma_-scalar(1));
//      volScalarField gSrce = gProd-gDest;
//      gProd.write();
//      gDest.write();
//      gSrce.write();
//      mag(symm(fvc::grad(U_)))().write();
//      CDkOmega.write();
//      Fturb()().write();
//      F1.write();     
    }

    // Transition onset momentum thickness Reynolds number equation
    tmp<fvScalarMatrix> ReThetatTildaEqn
    (
        fvm::ddt(ReThetatTilda_)
      + fvm::div(phi_, ReThetatTilda_)
      - fvm::laplacian(DReThetatTildaEff(), ReThetatTilda_)
     ==
        cThetat_*magSqr(U_)*(scalar(1.0)-FThetat())*ReThetatField/(scalar(500.0)*nu())
      - fvm::Sp(cThetat_*magSqr(U_)*(scalar(1.0)-FThetat())/(scalar(500.0)*nu()), ReThetatTilda_)
    );

    ReThetatTildaEqn().relax();
    solve(ReThetatTildaEqn);

    bound(ReThetatTilda_,scalar(20));
    
    convReThetatTilda_ = fvc::div(phi_, ReThetatTilda_);
    diffReThetatTilda_ = fvc::laplacian(DReThetatTildaEff(), ReThetatTilda_);
    PReThetatTilda_ = cThetat_*magSqr(U_)*(scalar(1.0)-FThetat())*ReThetatField/(scalar(500.0)*nu())
                      - fvc::Sp(cThetat_*magSqr(U_)*(scalar(1.0)-FThetat())/(scalar(500.0)*nu()), ReThetatTilda_);
  

    // Intermittency equation
    
    tmp<fvScalarMatrix> gammaEqn
    (
        fvm::ddt(gamma_)
      + fvm::div(phi_, gamma_)
      - fvm::laplacian(DgammaEff(), gamma_)
     ==
        Flength()*ca1_*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))*sqrt(Fonset()*gamma_)
      - fvm::Sp
        (
        Flength()*ca1_*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))*sqrt(Fonset()*gamma_)*ce1_,
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

    convgamma_ = fvc::div(phi_, gamma_);
    diffgamma_ = fvc::laplacian(DgammaEff(), gamma_);
    Pgamma_ = Flength()*ca1_*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))*sqrt(Fonset()*gamma_)
              - fvc::Sp
                (
                Flength()*ca1_*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))*sqrt(Fonset()*gamma_)*ce1_,
                gamma_
                );
    Dgamma_ = - ca2_*sqrt(scalar(2))*mag(skew(fvc::grad(U_)))*Fturb()*gamma_
              + fvc::Sp
              (
                  ce2_*ca2_*sqrt(scalar(2))*mag(skew(fvc::grad(U_)))*Fturb()*gamma_,
                  gamma_
              );


}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
