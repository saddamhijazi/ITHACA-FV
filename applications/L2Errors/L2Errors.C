#include "unsteadyNS.H"
#include "ITHACAstream.H"
#include "ITHACAutilities.H"


int main(int argc, char *argv[])
{

    //std::string word1 = argv[1];
    //std::string folder1 = argv[2];
    //std::string folder2 = argv[3];


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

    volScalarField nut
    (
        IOobject
        (
            "nut",
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
    volScalarField nut_rec("nutRec", nut);
    std::string word1 = "turbulent";

    PtrList<volVectorField> uFOM;
    PtrList<volScalarField> pFOM;
    PtrList<volScalarField> nutFOM;

    PtrList<volVectorField> uROM;
    PtrList<volScalarField> pROM;
    PtrList<volScalarField> nutROM;


    ITHACAstream::read_fields(uFOM, U, "./ITHACAoutput/Offline_check/");
    ITHACAstream::read_fields(pFOM, p, "./ITHACAoutput/Offline_check/");
    if(word1=="turbulent")
    {
        ITHACAstream::read_fields(nutFOM, nut, "./ITHACAoutput/Offline_check/");
    }


    ITHACAstream::read_fields(uROM, U_rec, "./ITHACAoutput/Reconstruction/");
    ITHACAstream::read_fields(pROM, P_rec, "./ITHACAoutput/Reconstruction/");
    if(word1=="turbulent")
    {
        ITHACAstream::read_fields(nutROM, nut_rec, "./ITHACAoutput/Reconstruction/");
    }



    Eigen::MatrixXd errorU;
    Eigen::MatrixXd errorP;
    Eigen::MatrixXd errorNut;

    

    errorU = ITHACAutilities::error_listfields(uFOM, uROM);
    errorP = ITHACAutilities::error_listfields(pFOM, pROM);
    if(word1=="turbulent")
    {
        errorNut = ITHACAutilities::error_listfields(nutFOM, nutROM);
    }
    
    ITHACAstream::exportMatrix(errorU, "errorU", "matlab", "./ITHACAoutput/postProcessing/");
    ITHACAstream::exportMatrix(errorP, "errorP", "matlab", "./ITHACAoutput/postProcessing/");
    if(word1=="turbulent")
    {
        ITHACAstream::exportMatrix(errorNut, "errorNut", "matlab", "./ITHACAoutput/postProcessing/");
    }

}
