struct LatLonDegrees64 {
	double lat, lon;
	LatLonDegrees64(double lat, double lon) : lat(lat), lon(lon) {}
};

typedef std::vector<LatLonDegrees64> LatLonDegrees64Vector;



// Data Loading
LatLonDegrees64Vector read_latlon_csv(std::string file_path) {
    LatLonDegrees64Vector cornerVector;
    std::ifstream file(file_path);
    std::string row; 
    while (std::getline(file, row)) {
        double lon = stof(row.substr(0, row.find(",")));        
        double lat = stof(row.substr(row.find(",")+1));
        cornerVector.push_back(LatLonDegrees64(lat, lon));       
    }
    file.close();    
    return cornerVector;
}



LatLonDegrees64Vector nodes_latlon = readCSV(file_name);    
ECEF64Vector nodes_ecef = LatLon2ECEFVector(nodes_latlon); 
