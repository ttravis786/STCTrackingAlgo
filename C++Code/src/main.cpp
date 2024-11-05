#include <iostream>
#include <string>
#include <fstream>
#include <cstdint>
#include <array>
#include <cmath>
#include <math.h>    
#include <chrono>

void possiblePointPair(std::array<std::array<float, 3>,2>& pairData, std::array<std::array<float, 2>,2> &coordPair, 
					std::array<std::array<std::array<float,2>, 2>,7> &pp, int &i){
	if (pairData[1][2] == pairData[0][2] ){
		for (int i=0; i<2; ++i){
			for (int j=0; j<0; ++j){
				coordPair[i][j] = 0;
			}
		} 
	}
	else {
		float tRatio = pairData[1][2]/pairData[0][2];
		for (int i=0; i<2; i++) {
			//coordPair[0][i] = pairDaia[0][i]
			//	+ ((pairData[1][i] - pairData[0][i])
			//		/(1+tRatio));
			// this breaks when t_RATIO IS 1 (EDGECASE)
			//coordPair[1][i] = pairData[0][i]
			//	+ ((pairData[1][i] - pairData[0][i])
			//			/(1-tRatio));
			coordPair[0][i] = pairData[1][i]/(1 + tRatio) + pairData[0][i]/(1 + 1/tRatio);
			// t1 - t2
			coordPair[1][i] = ((pairData[1][i] * pairData[0][2]) -  (pairData[0][i] * pairData[1][2]))/(pairData[0][2] - pairData[1][2]);

		}
	}
	
	//std::cout << "CoordPAIR" <<  coordPair[1][0] << ", " << coordPair[1][1] << "    \n";
	pp[i] = coordPair;
}

// edge case for points is same x, same t for this use x =v*t (v unkown)
void potentialPoints(std::array<float, 8>& x_arr, std::array<float, 8>& y_arr, std::array<float, 8>& t_arr, std::array<std::array<std::array<float,2>, 2>,7> &pp,
	std::array<std::array<float,3>, 2> &pairCoords,	std::array<std::array<float, 2>,2> &coordPair) {
	// Assume points are ordered in time.
	for (int i =0; i < 7; i++) {	
		for (int j = 0; j < 2; j++) {
			pairCoords[j][0] = x_arr[i + j];
			pairCoords[j][1] = y_arr[i + j];
			pairCoords[j][2] = t_arr[i + j];
		}
		possiblePointPair(pairCoords, coordPair, pp, i);
	}
}

float calculateGradient(std::array<std::array<std::array<float,2>, 2>,7> &pp ,std::array<int8_t, 7> &pointsSelected, int &pointsAdded){
	float  gradient;
	int prev = 0;
	int skip = pointsAdded -2;
	//std::cout << "\n\n" << "skip " << skip;
	for (int i=0; i<7; ++i){
		if (pointsSelected[i] != -1){
			if (skip!= prev){
				prev += 1;
			}
			else{
				for (int j=i+1; j<7; ++j){
					if (pointsSelected[j] != -1){
						//std::cout <<"i: " << i << "j: " << j << "\n";
						int psi = pointsSelected[i];
						int psj = pointsSelected[j];
						gradient = (pp[i][psi][1] - pp[j][psj][1]) / (pp[i][psi][0] - pp[j][psj][0]);	
						//std::cout << gradient << "\n";
						//std::cout <<"x1,y1: " << pp[i][psi][0] << "," << pp[i][psi][1] << " x2,y2: " << pp[j][psj][0] << "," << pp[j][psj][1]  << "\n";
						break;
					}
				}
				break;
			}
		}
	}
	return gradient;
}


int findNextPoint(std::array<std::array<std::array<float,2>, 2>,7> &pp ,std::array<int8_t, 7> &pointsSelected,
	      int &pointIndex,int &pointsAdded, float &originalGradient)
	   {
	// no points yet
	int np=0;
	// allow up to 5 degree error
	float pError= 10*( 3.1415 / 180);
	if (pointsAdded == 0){
		pointsAdded = 1;
		pointIndex = 1;
		pointsSelected[0] = 1;
		np = findNextPoint(pp, pointsSelected, pointIndex, pointsAdded, originalGradient);
		if (np == 1){
			std::cout <<"a\n";
			return 1;
		}
		// if failure try the other one
		pointsAdded = 1;
		pointIndex = 1;
		pointsSelected[0] = 0;
		np = findNextPoint(pp, pointsSelected, pointIndex, pointsAdded, originalGradient);
		if (np == 1){
			std::cout <<"b\n";
			return 1;
		}
		return 0;
	}
	// one other point
	else if (pointsAdded == 1){
		pointsAdded = 2;
		pointIndex = 2;
		pointsSelected[1] = 1;
		originalGradient = calculateGradient(pp, pointsSelected, pointsAdded);
		std::cout << "\n" << "Grad C" << originalGradient;
		np = findNextPoint(pp, pointsSelected, pointIndex, pointsAdded, originalGradient);
		if (np== 1){
			std::cout <<"c\n";
			return 1;
		}
		// if failure try the other one 0 indicates skip none (use first two points)
		pointsSelected[1] = 0;
		pointsAdded = 2;
		pointIndex = 2;
		originalGradient = calculateGradient(pp, pointsSelected, pointsAdded);
		std::cout << "\n" << "Grad D" << originalGradient;
		np = findNextPoint(pp, pointsSelected, pointIndex, pointsAdded, originalGradient);
		if (np== 1){
			std::cout <<"d\n";
			return 1;
		}
		return 0;
	}
	// return points if there are 4
	else if (pointsAdded == 6){
		//std::cout << "Nearly Works";
		return 1;
	}	
	int oPA = pointsAdded;
	int oPI = pointIndex;
	pointsAdded += 1;
	pointIndex += 1;
	pointsSelected[pointIndex-1] = 1;
	float gradient = calculateGradient(pp, pointsSelected, pointsAdded);
	std::cout << "\n" << gradient << ", " << originalGradient;
	// less than 5% error
	if (fabs(atan(gradient) -atan(originalGradient))<pError){
		np = findNextPoint(pp, pointsSelected, pointIndex, pointsAdded, originalGradient);
		if (np == 1){
			std::cout <<"e\n";
			return 1;
		}
	}
	pointsAdded = oPA + 1;
	pointIndex = oPI + 1;
	pointsSelected[pointIndex-1] = 0; 
	gradient = calculateGradient(pp, pointsSelected, pointsAdded);
	std::cout << " GradB: " << gradient << ", " << originalGradient;
	// less than 5% error
	if (fabs(atan(gradient) -atan(originalGradient))<pError){
		np = findNextPoint(pp, pointsSelected, pointIndex, pointsAdded, originalGradient);
		if (np == 1){
			std::cout <<"f\n";
			return 1;
		}
	pointsAdded = oPA;
	pointIndex = oPI;
	// point check
	std::cout << "Defo get to the end ";
	}
	return np;
}
// pick the inner and outer most point for most distance
void findFit(std::array<std::array<float,2>,6> &ap, float &m, float  &c){
	m = (ap[3][1] - ap[0][1]) / (ap[3][0] -  ap[0][0]);
	c = ap[3][1] - m * ap[3][0];
	//std::cout << m << "\n";
}

void findDistance(float &m, float &c, std::array<float, 8> & y_arr, std::array<float, 8> &x_arr, 
		  std::array<float, 8> &t_arr, std::array<float, 8> &v_arr){
	int normaliser = sqrt((m * m + 1));
	for (int i=0; i <8; ++i){
		//std::cout << "T avlue" << t_arr[i] << "\n";
		v_arr[i] = std::fabs(10000*(((y_arr[i] - m * x_arr[i] - c) / normaliser)/t_arr[i]));
	}
}

void calculate_v(float &v_sum, float &v_sum_square, float &mean_v, float &std_v, std::array<float, 8> &t_arr, std::array<float, 8> &v_arr){
	v_sum = 0;
	v_sum_square = 0;
	float total_weight = 0;
	float min_v = 1000;
	float max_v = 0;
	for (int i=0; i < 8; ++i){
		if (std::isinf(v_arr[i]) or t_arr[i] < 15){
			continue;
		}
		if (v_arr[i] > max_v){
			max_v = v_arr[i];	
		}
		if (v_arr[i] <  min_v){
			max_v = v_arr[i];	
		}
		v_sum += t_arr[i] * v_arr[i];
		total_weight += t_arr[i];
	}
	mean_v = v_sum/total_weight;
	std_v = (max_v-min_v)/2;
	//std_v = sqrt((v_sum_square/size) - (mean_v * mean_v));
}

void get_fit(std::array<std::array<std::array<float,2>, 2>,7> &pp, std::array<std::array<float,2>, 7> &aPoints, float &m, float &c, int &np, float &chi_square, int &numValidPoints, std::array<bool, 7> &validPoints, bool &pairChoice, float &trial_m, float &trial_c, int &combos, float &sumX, float &sumY, float &sumXY, float &sumX2, float &trial_chi_square, int &k){
	np = 0;
	chi_square = 100;
	numValidPoints = 7;
	combos = 1;
	for (int i=0; i<7; ++i){
		if (pp[i][0][0] == 0){
			numValidPoints -= 1;
			validPoints[i] = 0;
		}
		else{
			combos *= 2;
			validPoints[i] = 1;
		}
	}
	// points in track
	int maxPoints = 4;
	if (numValidPoints > maxPoints){
		int pointsToRemove = numValidPoints - maxPoints;
		int j = 0;
		numValidPoints = maxPoints;
		int index;
		while (pointsToRemove > 0) {
			if ((j%2)  == 0){
				index = 3 + std::floor((j+1)/2);
			}
			else{
				index = 3 - std::floor((j+1)/2);
			}
			if (validPoints[index] == 1){
				validPoints[index] = 0;
				pointsToRemove -= 1;
				combos /=2;
			}
			j +=1;
		}
	}
	//std::cout << "\n NEW TRACK ";
	for (int i=0; i<combos; ++i) {
		k = 0;
		sumX = 0;
		sumY = 0;
		sumXY = 0;
		sumX2 = 0;
	//	std::cout << "\n Combo ";
        	for (int j = 0; j <7; ++j) {
			if (!validPoints[j]){
				continue;
			}
           		pairChoice  = (i >> k) & 1;
	//		std::cout <<  pairChoice;
			k += 1;
 		        sumX += pp[j][pairChoice][0];
 		        sumY += pp[j][pairChoice][1];
        		sumXY += pp[j][pairChoice][1] * pp[j][pairChoice][0];
 		        sumX2 += pp[j][pairChoice][0] * pp[j][pairChoice][0];
	        }
		trial_m = (numValidPoints *  sumXY- sumX * sumY) / ( numValidPoints * sumX2 - sumX * sumX);
   		trial_c = (sumY - trial_m * sumX) / numValidPoints;
		k = 0;
		trial_chi_square = 0;
        	for (int j = 0; j <7; ++j) {
			if (!validPoints[j]){
				continue;
			}
           		pairChoice  = (i >> k) & 1;
			k += 1;
	//		std::cout << pp[j][pairChoice][1];
 		        trial_chi_square += (pp[j][pairChoice][1] - (trial_m * pp[j][pairChoice][0]) - trial_c)*(pp[j][pairChoice][1] - (trial_m * pp[j][pairChoice][0]) - trial_c);
		}
	//	std::cout << trial_m;
	///	std::cout << trial_c;	
	//	std::cout << "res" << trial_chi_square;
		if (trial_chi_square < chi_square){
			chi_square = trial_chi_square;
			m = trial_m;
			c = trial_c;
			np = 1;
		}		
		
	}
}


void actualPoints(std::array<std::array<std::array<float,2>, 2>,7> &pp, std::array<std::array<float,2>,6> &aPoints, int &np){
	std::array<int8_t, 7> pointsSelected = {-1,-1,-1,-1,-1,-1,-1};
	int pointIndex = 0;
	int pointsAdded = 0;
	float originalGradient = 0;
	np = findNextPoint(pp, pointsSelected, pointIndex, pointsAdded, originalGradient);
	pointsAdded = 0;
	if (np == 1){
		std::cout << "pointsFound"  << np << "\n";
		for (int i=0; i<7; ++i){
			if(pointsSelected[i] == 1){
				aPoints[pointsAdded] = pp[i][1];
				pointsAdded +=1;
				//std::cout << "Found Point1  " <<  pp[i][1][0] <<  ", " << pp[i][1][1] << "\n";
				//std::cout << pointsAdded << "\n";
			}
			if(pointsSelected[i] == 0){
				aPoints[pointsAdded] = pp[i][0];
				pointsAdded +=1;
				//std::cout << "Found Point0  " <<  pp[i][0][0] <<  ", " << pp[i][0][1] << "\n";
				//std::cout << pointsAdded << "\n";
			}
			}
		for (int i=0; i<6; ++i){
			//continue;
			//std::cout << "Found Point  " <<  aPoints[i][0] <<  ", " << aPoints[i][1] << "\n";
		}
	}
}

void  writeCSV(std::array<float, 8>& x_arr, std::array<float, 8>& y_arr, std::array<float, 8>& t_arr,
		int& tracknum, std::array<std::array<std::array<float,2>, 2>,7> &pp, std::array<std::array<float,2>,7> &ap){
	std::ofstream file("track_"+  std::to_string(tracknum)  + "rawData.csv");
	file << "x,y,t" << '\n';
	for (int i=0; i < 8; ++i){
		file << x_arr[i] << ',' << y_arr[i] << ',' << t_arr[i] << '\n';
	}
	file.close();
	std::ofstream filep("track_"+  std::to_string(tracknum)  + "pairData.csv");
	filep << "x,y" << '\n';
	for (int i=0; i < 7; ++i){
		for (int j=0; j<2; ++j){
			filep << pp[i][j][0] << ',' << pp[i][j][1]  << '\n';
		}
	}
	filep.close();
	std::ofstream fileap("track_"+  std::to_string(tracknum)  + "processedData.csv");
	fileap << "x,y" << '\n';
	for (int i=0; i < 7; ++i){
		fileap << ap[i][0] << ',' << ap[i][1]  << '\n';
	}
	fileap.close();

}
	

int main(int argc, char **argv){
	//TApplication app("app", &argc, argv);
	std::string fileLocation;
	std::cout << "Enter track file name: ";
	std::getline(std::cin, fileLocation);
	std::cout << "Reading in file: " << fileLocation << std::endl;
	std::ifstream file(fileLocation, std::ios::binary);
	if (!file.is_open()) {
		std::cerr << "Unable to open file!" <<std::endl; 
		return 0;
	}
	int tracknum = 0;
	int max_track_num = 1000000;
	const int n = 8;
	std::array<float, n> y_arr;
	std::array<float, n> x_arr;
	std::array<float, n> t_arr;
	std::array<float, n> v_arr;
	std::array<std::array<std::array<float,2>, 2>,7> pp;
	std::array<std::array<float,2>,7> ap;
	std::array<std::array<float,3>, 2> pairCoords;
	std::array<std::array<float,2>, 2> coordPair;
	float v_sum;
	std::ofstream fileParams("firParams.csv");
	fileParams << "v,v_std,beta,beta_std,m,c,np\n";
	float v_sum_square;
	float mean_v;
	float std_v; 
	float m;
	float c;
	float beta;
	int np;
	auto start = std::chrono::high_resolution_clock::now();	

	float chi_square = 1000000;
	int numValidPoints = 7;
	int combos = 1;
	std::array<bool, 7> validPoints;
	bool pairChoice;
	float trial_m;
	float trial_c;
	float sumX;
	float sumY;
	float sumXY;
	float sumX2;
	float trial_chi_square;
	int k;
	while ((!file.eof()) & (tracknum <  max_track_num)){
		//std::cout << "newTrack" << "\n\n\n";
		ap = {{{0,0},{0,0},{0,0},{0,0}, {0,0}, {0,0}, {0,0}}};
		//std::cout << "Track Number" << tracknum;
		for (int i = 0; i < n; ++i) {
			uint16_t value = 0;
			file.read(reinterpret_cast<char*>(&value), sizeof(uint16_t));
			x_arr[i] = ((value) & 0b111); 
			y_arr[i] = (value >> 3) & 0b111; 
			t_arr[i] = (value >> 6) & 0b1111111111;
			if (!(i%2==0)){
				y_arr[i] += 0.5;			
			}
			//std::cout << "Hit" << (i + 1) << " - X: " << x_arr[i]
			//	<< ", Y: " << y_arr[i] << ", Time: " << t_arr[i] << std::endl;
		}
		//pp[0][0][0] = 0; 
		//pp[0][1][0] = 0;
		potentialPoints(x_arr, y_arr, t_arr, pp, pairCoords, coordPair);
		
		//std::cout << pp[0][0][0] << ", "  << pp[0][0][1] << "\n";
		//std::cout << pp[0][1][0] << ", "  << pp[0][1][1] << "\n";
		//actualPoints(pp, ap, np);
		//std::cout << ap[0][0] << ap[0][0] << "\n";
		//findFit(ap, m, c);
		//std::cout << m << "\n";
		get_fit(pp, ap, m, c, np, chi_square, numValidPoints, validPoints, pairChoice, trial_m, trial_c, combos, sumX, sumY, sumXY, sumX2, trial_chi_square, k);
		findDistance(m, c, y_arr, x_arr, t_arr, v_arr);
		calculate_v(v_sum, v_sum_square, mean_v, std_v,t_arr, v_arr);
		//		plotData(x_arr, y_arr, t_arr, tracknum, pp);
		//tracknum += 1;
		beta = atan(m);
		//writeCSV(x_arr, y_arr, t_arr, tracknum, pp, ap);
		fileParams << 2*mean_v << "," << 2*std_v << "," << beta << "," << beta << ","  << m << "," << c << "," << np << '\n';
	}
	auto end = std::chrono::high_resolution_clock::now();
	// Calculate the duration in seconds (with decimal places)
   	std::chrono::duration<float> duration = end - start;
	fileParams.close();
   	// Output time in seconds
   	std::cout << "Time taken: " << duration.count() << " seconds" << std::endl;
	file.close();
	std::cout << "Done";
	//app.Run();
	return 0;
}
