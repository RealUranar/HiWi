#include <json.hpp>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <complex>
#include <set>
#include <Eigen/Core>
#include <Eigen/LU>

#include <limits>
#include <algorithm>
#include <fstream>


#include <H5Cpp.h>
#include <hdf5.h>

//#include <archive.h>
//#include <archive_entry.h>

using json = nlohmann::json;
using namespace std;


/*
FehlerQuellen:
    Einlesen von Matrixen
        -Koordinaten
        -AugCharge
*/


/**
 * Test if the LOBSTER input format is available
 */
bool isInFolder(string filename){
    std::ifstream file(filename);
    std::cout << file.good() << std::endl;
    return file.good();
}

Eigen::ArrayXd _readArray(std::vector<double> input){
	std::vector<double> dat = input;
	Eigen::ArrayXd arr = Eigen::ArrayXd::Map(dat.data(), dat.size());
	return arr;
}





/**
 * Get the valence of the used Pseudo Potentials
 * Formatted as:
 *  map["C"] = (2s, 2p)
 */
map<string, set<string>> getValConfigs(json& LobJSON){
    
    map<string, set<string>> valence;
    for (json PP: LobJSON["paw"]){
        vector<string> s =  PP["valence"];
        set<string> orbitals(s.begin(), s.end());
        valence[PP["element"]] =  orbitals;
    }
    /*Output testen
    std::map<std::string, set<string>>::iterator it;
    for (it = valence.begin(); it != valence.end(); it++){
        string elementSymbol = it->first;
        elementSymbol[0] = toupper(elementSymbol[0]);
        std::cout << elementSymbol << ' ';
        set<string>& orbitals = it->second;
        for (set<string>::iterator itSet=orbitals.begin(); itSet != orbitals.end(); itSet++){
            std::cout << "  " << *itSet<< "\n";
        }
    }*/
    return valence;
}

/**
 * Return the number of spinchannels
 */
int getNumberOfSpins(json& LobJSON){
    return LobJSON["wave_function"]["n_Spins"];
}

/**
 * Return the number of bands used in the calculation
 */
int getNumberOfBands(json& LobJSON){
    return LobJSON["wave_function"]["n_Bands"];
}

int getFractionalKpointsAndWeights(json& LobJSON, Eigen::Matrix3Xd *coordinates, Eigen::VectorXd *weights){
    int nKpoints = LobJSON["wave_function"]["spin_channels"][0]["k_points"].size();
    Eigen::Matrix3Xd kcoord = Eigen::Matrix3Xd::Zero(3,nKpoints);
    Eigen::VectorXd kweights = Eigen::VectorXd::Zero(nKpoints);
    
    int i = 0;
    for (json kPoint : LobJSON["wave_function"]["spin_channels"][0]["k_points"]){//Loop over every Kpoint
            kcoord(0,i) = double(kPoint["coordinates"][0]);
            kcoord(1,i) = double(kPoint["coordinates"][1]);
            kcoord(2,i) = double(kPoint["coordinates"][2]);
            kweights[i] = double(kPoint["weight"]);
            i++;
        }
	(*coordinates) = kcoord;
	(*weights) = kweights;
	return nKpoints;
}


Eigen::Matrix3d getLatticeReal(json& LobJSON){
    Eigen::Matrix3d LatReal;
    for (int x = 0; x < 3; x++){
        for (int y = 0; y < 3; y++){
            LatReal(x,y) = LobJSON["cell"]["lattice_vectors_real"][x][y];
        }
    }
    return LatReal;
}

Eigen::Matrix3d getLatticeRec(Eigen::Matrix3d LatReal){
    Eigen::Matrix3d LatRec = LatReal.inverse().transpose();
    return LatRec;
}

double getCellVolume(json& LobJSON){
    return LobJSON["cell"]["cell_volume"];
}


vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> getFractionalIonPositions(json& LobJSON) {
    vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> ionPos;
    for (json atom: LobJSON["cell"]["atomic_structure"]){
        Eigen::Vector3d coords; coords << atom["coordinates"][0], atom["coordinates"][1], atom["coordinates"][2];
        ionPos.push_back(coords);
    }
    return ionPos;
}

vector<string> getElementNames(json& LobJSON){//ACHTUNG HIER SIND DIE ELEMENTE GROßBUCHSTABEN, VIELLEICHT MÜSSEN ES KLEINE SEIN!!!!
    vector<string> elements;
    for (json atom: LobJSON["cell"]["atomic_structure"]){
        elements.push_back(atom["element"]);
    }
    return elements;
}

double getCutOffEnergy(json& LobJSON){
	return LobJSON["wave_function"]["cutoff_Energy"];
}

double getFermiLevel(json& LobJSON) {
	return LobJSON["wave_function"]["fermi_Energy"];
}


vector<int> getNumberOfPlanewavesAtKpoint(json& LobJSON) {
    vector<int> nPW;
    for (json kPoint : LobJSON["wave_function"]["spin_channels"][0]["k_points"]){
        nPW.push_back(kPoint["n_plane_waves"]);
    }
	return nPW;
}

Eigen::Vector3i getSizeOfGvectorGrid(json& LobJSON){
    json gVec = LobJSON["wave_function"]["sizeOfGvectorGrid"];
    Eigen::Vector3i gSize; gSize << gVec[0], gVec[1], gVec[2];
    return gSize;
}


//Zeile 1007 GVectorgrid Erstellen


double getEnergyMinimumForAnalysis(){
	return -20.00;
}
double getEnergyMaximumForAnalysis(){
	return 20.00;
}
int getEnergyStepsForAnalysis(){
	return 100;
}

void readUPF(json& LobJSON){
    const double wfConversion = pow(0.529177249,-3.0/2.0);

    for (json UPF : LobJSON["paw"]){
        std::string elementSymbol = UPF["element"];
        string PSPtitle = UPF["title"];
        std::cout << PSPtitle << endl;

        std::vector<double> dat = UPF["grid"];
        Eigen::ArrayXd rgrid = Eigen::ArrayXd::Map(dat.data(), dat.size());
        //Eigen::VectorxD rgrid_A = rgrid * AUTOA;

        int nProj = UPF["augmentationCharges"].size();
        Eigen::MatrixXd readMatrix(nProj,nProj);
        for (int x = 0; x < nProj; x++){
            for (int y = 0; y < nProj; y++){
                readMatrix(x,y) = UPF["augmentationCharges"][x][y];
            }
        }
        for (int idummy = 0; idummy < nProj; idummy++){
            Eigen::ArrayXd projektors = _readArray(UPF["projectors"][idummy]);
            projektors *=wfConversion;
        }

    }
    return;
}

void readOccsAndEigen(json& LobJSON){
    int nBands = LobJSON["wave_function"]["n_Bands"];
    for (json spin :LobJSON["wave_function"]["spin_channels"]){
        int ikpt = 0;
        for (json kpoint: spin["k_points"]){
            for (int band = 0; band < nBands ; band++){
                double occval = kpoint["occupations"][band];
                double eval = kpoint["eingen_values"][band];
            }
            ikpt++;
        }
    }
    return;
}

/*std::stringstream _readArchive(char const* archiveName, string filename){
    struct archive *a;
    struct archive_entry *entry;
    int r;
    a = archive_read_new();
    archive_read_support_filter_all(a);
    archive_read_support_format_all(a);
    r = archive_read_open_filename(a, archiveName, 10240); // Note 1
    if (r != ARCHIVE_OK)
        exit(1);
    while (archive_read_next_header(a, &entry) == ARCHIVE_OK) {
        if (string(archive_entry_pathname(entry)) == filename){
            size_t total = archive_entry_size(entry);
            char buff[total];

            ssize_t size = archive_read_data(a, buff, total);
            //ssize_t size = archive_read_data(a, buff, total);
            std::stringstream iss(buff);
            r = archive_read_free(a);
            return iss;
        }
        printf("%s\n",archive_entry_pathname(entry));
        archive_read_data_skip(a);  // Note 2
        }
    r = archive_read_free(a);
      // Note 3
    std::cout << "\n" << "Could not find File!!!" <<"\n";
    exit(1);
}*/



void readPWCoeffs(){
    H5::H5File Kpoint("../kPoint1.hdf5",H5F_ACC_RDONLY);
    H5::DataSet dataSet = Kpoint.openDataSet("PWCoeffs");

	// Read data into buffer matrix of type double before casting into complex type
    H5::DataSpace dataspace = dataSet.getSpace();
    hsize_t dims_out[2];
    int nDims = dataspace.getSimpleExtentDims(dims_out, NULL);

    typedef struct complex_type{
        double r;
        double i;
    } complex_type;

    H5::CompType complex_compound( sizeof(complex_type));
    complex_compound.insertMember("r", HOFFSET(complex_type, r), H5::PredType::NATIVE_DOUBLE);
    complex_compound.insertMember("i", HOFFSET(complex_type, i), H5::PredType::NATIVE_DOUBLE);

	Eigen::MatrixXcd buff = Eigen::MatrixXcd(dims_out[1],dims_out[0]);

    //H5::FloatType datatype(H5::PredType::IEEE_F64LE);
    dataSet.read(buff.data(), complex_compound);

    std::cout << "\n" <<  buff.col(0)[0] << "\n";
	//dataSet.read(buff.data(), datatype, dataspace);

	// Don't forget to close everything...
	dataspace.close();
	dataSet.close();
	Kpoint.close();

    //std::cout << buff;
	// Cast to complex type and fill into coefficients vector
	// TODO: This would be nicer if we could loop over bands in here, but the existing framework does not allow this
	// Refactor?
    return;
}


void readGvectorGrid(){
    H5::H5File Kpoint("../kPoint1.hdf5",H5F_ACC_RDONLY);
    H5::DataSet dataSet = Kpoint.openDataSet("Miller");

	// Read data into buffer matrix of type double before casting into complex type
    H5::DataSpace dataspace = dataSet.getSpace();
    hsize_t dims_out[2];
    int nDims = dataspace.getSimpleExtentDims(dims_out, NULL);

    Eigen::MatrixXd buff(dims_out[1], dims_out[0]);

    dataSet.read(buff.data(), H5::PredType::NATIVE_DOUBLE);

    //Eigen::MatrixXi buffi = buff;

    for (int pw=0; pw < dims_out[0]; pw++){
        Eigen::Vector3i Test = buff.col(pw).cast<int>();
        Eigen::Vector3d gvectorFractional = buff.col(pw);
        Eigen::Vector3d gPlusKFractional = (gvectorFractional+ Eigen::Vector3d(0,0,0));
    }

	// Don't forget to close everything...
	dataspace.close();
	dataSet.close();
	Kpoint.close();

    //std::cout << buff;
	// Cast to complex type and fill into coefficients vector
	// TODO: This would be nicer if we could loop over bands in here, but the existing framework does not allow this
	// Refactor?
    return;
}


int main(){
    json LobJSON;
    if (!isInFolder("TESTING/LobsterInput.json")){
        std::cout << "\n File not found! \n";
        return 0;
    }
    //std::stringstream file = _readArchive("../LobsterInput.tar.gz", "LobsterInput.json");
    std::ifstream file("TESTING/LobsterInput.json");
    file >> LobJSON;
    // Eigen::Matrix3Xd coordinates;
    // Eigen::VectorXd weights;

    // getValConfigs(LobJSON);
    // getNumberOfBands(LobJSON);
    // getFractionalKpointsAndWeights(LobJSON, &coordinates, &weights);
    // Eigen::Matrix3d LatReal = getLatticeReal(LobJSON);
    // Eigen::Matrix3d LatRec = getLatticeRec(LatReal);

    // getCellVolume(LobJSON);
    // getFractionalIonPositions(LobJSON);
    // getElementNames(LobJSON);
    // getCutOffEnergy(LobJSON);
    // double x = getFermiLevel(LobJSON);
    
    // getNumberOfPlanewavesAtKpoint(LobJSON);
    
    readUPF(LobJSON);
    // getSizeOfGvectorGrid(LobJSON);
    // //std::cout << "\n" << x << "\n";
    // readOccsAndEigen(LobJSON);
    // readPWCoeffs();
    // readGvectorGrid();
    return 0;
}