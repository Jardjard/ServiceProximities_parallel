/*
 *	Program:			ServiceProximities_parallel.exe
 *	Module:				ServiceProximities_parallel_main.cpp
 *	Author:				J. Lerner, K. Haldane
 *	Date:				November 20, 2014
 *	Description:		Parallel version of the Service Proximities program.
 *						Calculates the closest distance to a given service (services.dat)
 *						for each residence (residences.dat) and prints statistics on the results.
 */

#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>
#include <iterator>
#include <ctime>

#include <mpi.h>
 
using namespace std;

const string IN_FILE_RESIDENCES = "residences.dat";
const string IN_FILE_SERVICES = "services.dat";
const string IN_FILE_SERVICE_CODES = "service-codes.csv";
const int TAG_DATA = 1, TAG_QUIT = 2; 

// Class used to store an result set and some of its logic
class ResultSet
{
public:
	int count, onekm, twokm, fivekm, morethanfivekm;
	double onekmPercent, twokmPercent, fivekmPercent, morethanfivekmPercent;

	void calcPercentages()
	{
		count = onekm + twokm + fivekm + morethanfivekm;
		onekmPercent = (double)onekm / count * 100.0;
		twokmPercent = (double)twokm / count * 100.0;
		fivekmPercent = (double)fivekm / count * 100.0;
		morethanfivekmPercent = (double)morethanfivekm / count * 100.0;
	}
};

struct Service {
  double lat;
  double lon;
};

MPI_Datatype createRecType()
{
	// Setup the five arguments for a call to MPI_Type_struct()
	int count = 2;	// 2 blocks
	int blocklens[] = { 5, 4 };	// 4 ints, 4 doubles
	MPI_Aint offsets[2];	
	offsets[0] = offsetof(ResultSet, count); // offset in bytes to block 1
	offsets[1] = offsetof(ResultSet, onekmPercent); // offset in bytes to block 2
	MPI_Datatype old_types[] = { MPI_INT, MPI_DOUBLE };	// input data types
	MPI_Datatype t; // output data type
	
	// Call the datatype contructor function
	MPI_Type_struct(count, blocklens, offsets, old_types, &t);

	// Allocate memory for the new type
	MPI_Type_commit(&t);

	return t;
}

void results(string serviceCode, float duration, int count, ResultSet* sepCounts, int numProcs, string locationNamme) {
	int onekm = 0, twokm = 0, fivekm = 0, morethanfivekm = 0;


	cout.imbue(std::locale(""));
	cout << "\n-------------------------------------------------------------\n";
	cout << " Proximities of Residential Addresses to Services in Toronto\n";
	cout << "-------------------------------------------------------------\n\n";
	cout << " Service:" << std::setw(51) << std::right << locationNamme << endl;
	cout << " Service Code:"<< std::setw(46) << std::right << atoi(serviceCode.c_str()) << endl;
	cout << " Number of Service Locations:"<< std::setw(31) << std::right << count << endl;
	cout << " Number of Processes:"<< std::setw(39) << std::right << numProcs << endl;
	cout << " Elapsed Time in Seconds:"<< std::setw(35) << std::right << duration << endl << endl << endl;

	for (int i = 0; i < numProcs; ++i) {
		onekm += sepCounts[i].onekm;
		twokm += sepCounts[i].twokm;
		fivekm += sepCounts[i].fivekm;
		morethanfivekm += sepCounts[i].morethanfivekm;
		
		cout << " Process #" << (i+1) << " Results for " << sepCounts[i].count << " Addresses...\n\n";
		cout << " Nearest Service (km)     # of Addresses      % of Addresses\n";
		cout << "----------------------   ----------------    ----------------\n";
		cout << "                0 - 1" << std::setw(19) << std::right << sepCounts[i].onekm << " " << fixed << setprecision(2) << std::setw(19) << std::right << sepCounts[i].onekmPercent << endl;
		cout << "                1 - 2" << std::setw(19) << std::right << sepCounts[i].twokm << " " << fixed << setprecision(2) << std::setw(19) << std::right << sepCounts[i].twokmPercent << endl;
		cout << "                2 - 5" << std::setw(19) << std::right << sepCounts[i].fivekm << " " << fixed << setprecision(2) << std::setw(19) << std::right << sepCounts[i].fivekmPercent << endl;
		cout << "                  > 5" << std::setw(19) << std::right << sepCounts[i].morethanfivekm << " " << fixed << setprecision(2) << std::setw(19) << std::right << sepCounts[i].morethanfivekmPercent << endl << endl;

	}
	
	int locationcount = onekm + twokm + fivekm + morethanfivekm;

	cout << " Aggregate Results for all " << locationcount << " Addresses...\n\n";
	cout << " Nearest Service (km)     # of Addresses      % of Addresses\n";
	cout << "----------------------   ----------------    ----------------\n";
	cout << "                0 - 1" << std::setw(19) << std::right << onekm << " " << fixed << setprecision(2) << std::setw(19) << std::right << (double)onekm/locationcount * 100 << endl;
	cout << "                1 - 2" << std::setw(19) << std::right << twokm << " " << fixed << setprecision(2) << std::setw(19) << std::right << (double)twokm/locationcount * 100 << endl;
	cout << "                2 - 5" << std::setw(19) << std::right << fivekm << " " << fixed << setprecision(2) << std::setw(19) << std::right << (double)fivekm/locationcount * 100 << endl;
	cout << "                  > 5" << std::setw(19) << std::right << morethanfivekm << " " << fixed << setprecision(2) << std::setw(19) << std::right << (double)morethanfivekm/locationcount * 100 << endl;
	cout << "\n-------------------------------------------------------------\n";
	cout << "-------------------------------------------------------------\n\n";
}

int main(int argc, char* argv[])
{
	if(MPI_Init(&argc, &argv) == MPI_SUCCESS)
	{
		// Get the # of processes and rank of this process
		int numProcs, rank;
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		ResultSet rs;
		
		// Create a derived type for passing the result set
		MPI_Datatype recType = createRecType();

		if (argc <=1) {
			if (rank == 0) { 
				cout << "USAGE ERROR: Include a 6 digit service code on the command line.\n";
			}
			// Free memory allocated for the derived type recType
			MPI_Type_free(&recType);
			MPI_Finalize();
			return EXIT_FAILURE;
		}
		else {
			// Prepare to read and process the input data
			ifstream in_res(IN_FILE_RESIDENCES.c_str());
			ifstream in_serv(IN_FILE_SERVICES.c_str());
			ifstream in_serv_code(IN_FILE_SERVICE_CODES.c_str());

			if (!in_res.is_open() || in_res.fail())
			{
				cerr << "Error opening input file " << IN_FILE_RESIDENCES << endl;
				// Free memory allocated for the derived type recType
				MPI_Type_free(&recType);
				MPI_Finalize();
				return EXIT_FAILURE;
			}
			if (!in_serv.is_open() || in_serv.fail())
			{
				cerr << "Error opening input file " << IN_FILE_SERVICES << endl;
				// Free memory allocated for the derived type recType
				MPI_Type_free(&recType);
				MPI_Finalize();
				return EXIT_FAILURE;
			}
	
			if (!in_serv_code.is_open() || in_serv_code.fail())
			{
				cerr << "Error opening input file " << IN_FILE_SERVICES << endl;
				// Free memory allocated for the derived type recType
				MPI_Type_free(&recType);
				MPI_Finalize();
				return EXIT_FAILURE;
			}

			// get the requested service code from the command line argument
			string serviceCode = argv[1];
			vector<Service> vService;

			string code_in;
			double service_lat, service_long;
			int count = 0;
			int onekm = 0, twokm = 0, fivekm = 0, morethanfivekm = 0;

			clock_t start;
			float duration;
			if (rank == 0) {
				start = clock();
			}

			// Loop through all the service and pick the ones that match the code
			while (in_serv >> code_in >> service_lat >> service_long ) 
			{
				if (code_in.compare(serviceCode) == 0) {
					Service s = {};
					s.lat = service_lat;
					s.lon = service_long;
					vService.push_back( s );
					count++;
				}
			}

			// Determine the number of bytes in a record
			string rec;
			getline(in_res, rec);
			int sinBytes = (int)in_res.tellg();

			// Move the file pointer to MY first record
			in_res.seekg(sinBytes * rank, ios::beg);

			// check each residence location ...
			double residence_lat, residence_long;
			while (in_res >> residence_lat >> residence_long) {

				double shortest_distance = -1;
				// ... against each service location
				for( vector<Service>::iterator it = vService.begin(); it != vService.end(); ++it) {
					double dist, x_dif, y_dif;
					x_dif = it->lat - residence_lat;
					y_dif = it->lon - residence_long;
					dist = sqrt(x_dif * x_dif + y_dif * y_dif);
					// update shortest_distance if a shorter distance is found
					if (shortest_distance == -1 || shortest_distance > dist ) {
						shortest_distance = dist;
					}
				}
				if (shortest_distance <= 1000) ++onekm;
				else if (shortest_distance <= 2000) ++twokm;
				else if (shortest_distance <= 5000) ++fivekm;
				else ++morethanfivekm;

				// seek to the next offset for this process
				in_res.seekg(sinBytes * (numProcs-1), ios::cur);
			}
			// collect results for this process
			rs.onekm = onekm;
			rs.twokm = twokm;
			rs.fivekm = fivekm;
			rs.morethanfivekm = morethanfivekm;
			rs.calcPercentages();
			
			// gather the results from all processes to Process 0
			ResultSet* sepCounts = new ResultSet[numProcs];
			MPI_Gather(&rs, 1, recType, sepCounts, 1, recType, 0, MPI_COMM_WORLD);

			// Process 0 handles the output
			if (rank == 0) {
				duration = (clock() - (float)start) / CLOCKS_PER_SEC;
				string location;
				string locationName = "Unknown Location";
				// get the location name
				while (getline(in_serv_code, location, ',') && getline(in_serv_code, locationName)) {
					if(location == serviceCode) {
						break;
					}
				}
				// print the results
				results(serviceCode, duration, count, sepCounts, numProcs, locationName);
			}
		}
		// Free memory allocated for the derived type recType
		MPI_Type_free(&recType);
		MPI_Finalize();
	}
	return EXIT_SUCCESS;
}





























//
//inline int countDigits(int num) {
//	if (num < 10) {
//		return 1;
//	}
//	else {
//		return countDigits(num / 10) + 1;
//	}
//}
//
//inline void printNumber(int num, int length) {
//	for (int i = countDigits(num); i < length; ++i)
//		cout << " ";
//	cout << num;
//}










// MPI

/*
	Program:			TPV_mpi1.exe
	Module:				TPV_mpi1_main.cpp
	Author:				T. Haworth
	Date:				October 29, 2014
	Description:		This is a parallel version of the program that uses a
						contiguous and homogeneous derived data type based on 
						MPI_DOUBLE. The data type is created using the 
						MPI_Type_contiguous() function since it is designed to
						simplify the creation of this kind of derived type (i.e.
						a new type based on the contiguous replication of a
						homogeneous type).

						Present Value Formula:

						A = R * ( (1 - (1+i)^(-n) ) / i ) * (1+i)^(-k)

						Where:
						A = the present value of the annuity
						R = the amount ($) of each periodic payment of the annuity
						i = the interest rate per payment period
						n = the number of payments (periods) of the annuity
						k = the number of deferment periods before the annuity's
						first payment
						
						Note that you don't really have to understand much at all about
						annuities to grasp the parallel programming concept introduced by
						this set of examples which is using derived data types to pass
						composite data structures between processes.

	Runtime Tip:		Make sure the working directory contains the annuities.dat file.
						if you launch the program using the "wmpiexec" tool make sure to
						click the "more options" checkbox at the bottom left, then browse
						for the correct "working directory" before you click "Execute".

	Status:				Complete.
*/




//
//#include <fstream>
//#include <iostream>
//#include <cmath>
//#include <string>
//#include <iomanip>
//#include <mpi.h>
//using namespace std;
//
//const string IN_FILE = "annuities.dat";
//const int NUM_FLDS = 4;
//const int TAG_DATA = 1, TAG_QUIT = 2;
//
//inline double calcPV(double r[NUM_FLDS])
//{
//	// A = R x(1 –(1 + i) ^ - n) / i x(1 + i) ^ - k
//	return r[0] * (1 - pow(r[1] + 1, -r[2])) / r[1] * pow(1 + r[1], -r[3]);
//}
//
//void terminateSlaves(int numProcs)
//{
//	double dummy[NUM_FLDS];
//	for (int i = 1; i < numProcs; ++i)
//		MPI_Send(&dummy, 4, MPI_DOUBLE, i, TAG_QUIT, MPI_COMM_WORLD);
//}
//
//int main(int argc, char* argv[])
//{
//	if (MPI_Init(&argc, &argv) == MPI_SUCCESS)
//	{
//		// Obtain the rank and # of processes
//		int numProcs, rank;
//		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
//		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//		double totalPV = 0;
//		double rec[NUM_FLDS];
//		int recNum = 0;
//
//		// Create a derived type for passing annuity records
//		MPI_Datatype recType;
//		MPI_Type_contiguous(NUM_FLDS, MPI_DOUBLE, &recType);
//		MPI_Type_commit(&recType);
//
//		if (rank == 0)
//		{
//			// Master process reads and processes the input data
//
//			ifstream in(IN_FILE);
//
//			if (!in.is_open() || in.fail())
//			{
//				cerr << "Error opening input file " << IN_FILE << endl;
//
//				// Exit the program cleanly
//				terminateSlaves(numProcs);
//				MPI_Finalize();
//				return EXIT_FAILURE;
//			}
//
//
//			int sendId;
//			while (in >> rec[0 ]>> rec[1] >> rec[2] >> rec[3])
//			{
//				sendId = recNum % numProcs;
//
//				if (sendId == 0)
//					// Master process calculates this annuity's PV
//					totalPV += calcPV(rec);
//				else
//					// Send the rec to a slave process
//					//MPI_Send(rec, NUM_FLDS, MPI_DOUBLE, sendId, TAG_DATA, MPI_COMM_WORLD);
//					MPI_Send(rec, 1, recType, sendId, TAG_DATA, MPI_COMM_WORLD);
//
//				recNum++;
//			}
//
//			// Tell the slaves to quit
//			terminateSlaves(numProcs);
//		}
//		else
//		{
//			// Slave receives annuity record and calculates the PV
//			MPI_Status status;
//			do
//			{
//				//MPI_Recv(rec, NUM_FLDS, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//				MPI_Recv(rec, 1, recType, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//
//				if (status.MPI_TAG == TAG_DATA)
//					totalPV += calcPV(rec);
//
//			} while ( status.MPI_TAG != TAG_QUIT );
//		}
//
//		// Aggregate the PVs for all processes 
//		double grandTotalPV;
//		MPI_Reduce(&totalPV, &grandTotalPV, 1, MPI_DOUBLE, MPI_SUM, 0, 
//			MPI_COMM_WORLD);
//
//		if (rank == 0)
//		{
//			// Master process reports the total PV for all the annuities
//			cout << "Number of annuities processed: " << recNum << endl;
//			cout << "Number of processes used: " << numProcs << endl;
//			cout << fixed << setprecision(2) << "The total present value is: $"
//				<< grandTotalPV << endl;
//		}
//
//		// Free memory allocated for the derived type recType
//		MPI_Type_free(&recType);
//
//		MPI_Finalize();
//	}
//
//	
//	return EXIT_SUCCESS;
//}