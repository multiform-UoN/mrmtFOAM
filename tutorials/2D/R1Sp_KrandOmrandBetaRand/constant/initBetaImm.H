#codeStream{


        code
        #{
            const IOdictionary& d = static_cast<const IOdictionary&>(dict);
            const fvMesh& mesh = refCast<const fvMesh>(d.db());

            
            const volScalarField& beta = mesh.lookupObject<volScalarField>("beta"); 

     
            scalarField  betaI
            ( 
                    (1.0 - beta.internalField())
                    /
                    beta.internalField()
            );            

            writeEntry(os,betaI);

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
