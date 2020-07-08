//- Constrain deltaT based on the velocity-based Courant number 
if (adjustTimeStep)
{
    scalar maxDeltaTFact = maxCo/(CoNum + small);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    scalar deltaTOld(runTime.deltaTValue());
    
    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "Ajusting time step to Courant Number\n"    
        << "Time step corrected from deltaT = " << deltaTOld 
        << "  to deltaT = " << runTime.deltaTValue() << endl;
}

