#codeStream{


        code
        #{
            const IOdictionary& d = static_cast<const IOdictionary&>(dict);
            const fvMesh& mesh = refCast<const fvMesh>(d.db());

            //Initialise capacity field
            scalarField betaI(mesh.nCells(), scalar(0));
            
            const volScalarField& beta = mesh.lookupObject<volScalarField>("beta"); 

            
            forAll(beta,cellI)
            {
                betaI[cellI] = 
                    (1.0 - beta.internalField()[cellI])
                    /
                    beta.internalField()[cellI];
            }            

            os << betaI;

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
