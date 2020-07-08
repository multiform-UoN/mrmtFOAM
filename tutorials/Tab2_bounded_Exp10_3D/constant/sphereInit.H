
dimensions      [1 -3 0 0 0 0 0];

 internalField uniform 1.0 ;

boundaryField
{

    inlet
    {
        type            zeroGradient;
    }

    outlet
    {
        type            zeroGradient;
    }

    wall
    {
        type            zeroGradient;
    }

}
