#include "unsteadyNS.H"
#include "ITHACAstream.H"
#include "ITHACAutilities.H"


int main(int argc, char *argv[])
{


#include "setRootCase.H"
    Foam::Time runTime(Foam::Time::controlDictName, args);
    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
            )
        );
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
            ),
        mesh
        );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
            ),
        mesh
        );

    volVectorField U_rec("uRec", U);
    volScalarField P_rec("pRec", p);

    PtrList<volVectorField> UHF;
    PtrList<volScalarField> PHF;
    PtrList<volVectorField> USUP;
    PtrList<volScalarField> PSUP;
    PtrList<volVectorField> UPPE;
    PtrList<volScalarField> PPPE;
    // PtrList<volVectorField> UNOS;
    // PtrList<volScalarField> PNOS;


    ITHACAstream::read_fields(UHF, U, "./ITHACAoutput/Offline_check/");
    ITHACAstream::read_fields(PHF, p, "./ITHACAoutput/Offline_check/");
    ITHACAstream::read_fields(USUP, U_rec, "./ITHACAoutput/Reconstruction/");
    ITHACAstream::read_fields(PSUP, P_rec, "./ITHACAoutput/Reconstruction/");
    //ITHACAstream::read_fields(UPPE, U_rec, "./ITHACAoutput/ReconstructionPPE/",0,201);
    //ITHACAstream::read_fields(PPPE, P_rec, "./ITHACAoutput/ReconstructionPPE/",0,201);


    Eigen::MatrixXd errorUPPE;
    Eigen::MatrixXd errorPPPE;
    Eigen::MatrixXd errorUSUP;
    Eigen::MatrixXd errorPSUP;
    

    //errorUPPE = ITHACAutilities::error_listfields(UHF, UPPE);
    //errorPPPE = ITHACAutilities::error_listfields(PHF, PPPE);
    errorUSUP = ITHACAutilities::error_listfields(UHF, USUP);
    errorPSUP = ITHACAutilities::error_listfields(PHF, PSUP);
    
    ITHACAstream::exportMatrix(errorUSUP, "errorUSUP", "python", "./ITHACAoutput/postProcessing/");
    ITHACAstream::exportMatrix(errorPSUP, "errorPSUP", "python", "./ITHACAoutput/postProcessing/");


    ITHACAstream::exportMatrix(errorUSUP, "errorUSUP", "matlab", "./ITHACAoutput/postProcessing/");
    ITHACAstream::exportMatrix(errorPSUP, "errorPSUP", "matlab", "./ITHACAoutput/postProcessing/");



    //ITHACAstream::exportMatrix(errorUPPE, "errorUPPE", "python", "./ITHACAoutput/postProcessing/");
    //ITHACAstream::exportMatrix(errorPPPE, "errorPPPE", "python", "./ITHACAoutput/postProcessing/");

}
