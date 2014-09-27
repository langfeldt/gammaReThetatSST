/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | Unsupported Contributions for OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Felix Langfeldt
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
    Foam::incompressible::ReThetatTildaInletFvPatchScalarField

Group
    grpIcoRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a inlet boundary condition for the
    ReThetatTilda variable, which is used in the gammaReThetatSST transitional
    turbulence model. The inlet condition is based upon the turbulence
    intensity, defined as a fraction of the mean velocity:

        \f[
            k_p = 1.5 (I |U|)^2
        \f]

    where

    \vartable
        k_p     | turbulence kinetic energy at the patch
        I       | turbulence intensity
        U       | velocity field
    \endvartable

    The inlet value of ReThetatTilda is a function of the turbulence intensity:

        \f[
            ReThetatTilda_p = ReThetat(I)
        \f]

    where

    \vartable
        ReThetatTilda_p     | ReThetatTilda at the patch
        I                   | turbulence intensity
        ReThetat()          | empirical correlation
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        U            | velocity field name     | no          | U
        phi          | flux field name         | no          | phi
        k            | turbulence kinetic energy field name | no | k
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type        ReThetatTildaInlet;
        value       uniform 1;      // placeholder
    }
    \endverbatim

Note
    In the event of reverse flow, a zero-gradient condition is applied

SeeAlso
    Foam::inletOutletFvPatchField

SourceFiles
    ReThetatTildaInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef ReThetatTildaInletFvPatchScalarField_H
#define ReThetatTildaInletFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
   Class ReThetatTildaInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class ReThetatTildaInletFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        // Correlation names
        static const word CORRN_MENTER2009;
        static const word CORRN_SULUKSNA2009;
        static const word CORRN_MALAN2009;
        static const word CORRN_SORENSEN2009;

        // Correlation IDs
        static const int CORR_MENTER2009;
        static const int CORR_SULUKSNA2009;
        static const int CORR_MALAN2009;
        static const int CORR_SORENSEN2009;
        
        // Selected Correlation ID        
        int corrID_;

        //- Name of the velocity field
        word UName_;

        //- Name of the turbulent kinetic energy field
        word kName_;


public:

    //- Runtime type information
    TypeName("ReThetatTildaInlet");


    // Constructors

        //- Construct from patch and internal field
        ReThetatTildaInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ReThetatTildaInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  ReThetatTildaInletFvPatchScalarField
        //  onto a new patch
        ReThetatTildaInletFvPatchScalarField
        (
            const ReThetatTildaInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        ReThetatTildaInletFvPatchScalarField
        (
            const ReThetatTildaInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new ReThetatTildaInletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        ReThetatTildaInletFvPatchScalarField
        (
            const ReThetatTildaInletFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new ReThetatTildaInletFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
