#codeStream{


        code
        #{
            const IOdictionary& d = static_cast<const IOdictionary&>(dict);
            const fvMesh& mesh = refCast<const fvMesh>(d.db());

            const volScalarField& beta = mesh.lookupObject<volScalarField>("beta"); 

            scalar a(1e-5);
            
                       
            //Initialise capacity field
            scalarField omega(mesh.nCells(), scalar(0));
            
            //Create Field K
            volTensorField K
            (
                IOobject
                (
                    "K",
                    mesh.time().constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );
            

            scalarField k((1.0/3.0)*tr(K.internalField()));
            
            scalarField betaM(beta.internalField());
            
            forAll(omega,cellI)
            {
                omega[cellI] =             
                    a*(pow(betaM[cellI],3))
                    /
                    (
                        k[cellI]*pow(1.0-betaM[cellI],2)
                    );
            }            


            os << omega;

        #};

        //- Optional:
        codeInclude
        #{
            #include "fvCFD.H"
            //#include <stdlib.h>
        #};

        //- Optional:
        codeOptions
        #{
            -I$(LIB_SRC)/finiteVolume/lnInclude \
            -I$(LIB_SRC)/meshTools/lnInclude
        #};


};
