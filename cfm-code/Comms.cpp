/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# comms.cpp
#
# Description: 	Classes for communicating data (e.g. parameters, partial
#				gradients..etc) during parameter update - see
param.cpp.
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Comms.h"
#include "mpi.h"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/set.hpp>
#include <sstream>
#include <string>

Comms::Comms() {
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_nump);
}

void Comms::printWithWorkerId(const char *msg) { std::cout << mpi_rank << ": " << msg << std::endl; }

void Comms::collectGradsInMasterOrigMpi(std::vector<float> &grads) {
	std::vector<float> global_grads(grads.size(), 0.0f);
	MPI_Reduce(&grads[0], &global_grads[0], grads.size(), MPI_FLOAT, MPI_SUM, MASTER, MPI_COMM_WORLD);
	grads = global_grads;
	/*for(auto & grad : grads)
	    grad /= float(mpi_nump);*/
}

int Comms::collectSumInMaster(int partial) {
	int tmp = partial;
	int total;
	MPI_Reduce(&tmp, &total, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
	return total;
}

float Comms::getTimeUsages(float time_used, MPI_Op op) {
	float rev = 0.0;
	MPI_Reduce(&time_used, &rev, 1, MPI_FLOAT, op, MASTER, MPI_COMM_WORLD);
	return rev;
}

void Comms::gatherTimeUsages(float time_used, std::vector<float> &time_used_vector) {
	time_used_vector.resize(mpi_nump);
	std::vector<float> local_time_used(mpi_nump, 0.0f);
	local_time_used[mpi_rank] = time_used;
	MPI_Gather(&time_used_vector[0], time_used_vector.size(), MPI_FLOAT, &time_used_vector[0], 1, MPI_FLOAT, MASTER,
	           MPI_COMM_WORLD);
}

void WorkerComms::setMasterUsedIdxs() {

	MPI_Barrier(MPI_COMM_WORLD); // All threads wait
	num_used = used_idxs.size();

	// Serialize the data for this worker's used idxs
	std::ostringstream str;
	boost::archive::text_oarchive ar(str);
	ar & used_idxs;
	std::string serialized_used_idxs = str.str();

	// Let the master know how many characters will be sent for this worker
	unsigned int data_size = serialized_used_idxs.size();
	MPI_Send(&data_size, 1, MPI_CHAR, MASTER, 0, MPI_COMM_WORLD);

	// Send the master the (serialized) used idxs for this worker
	if (data_size > 0) {
		// std::cout << mpi_rank << "start_sending" << "..." << std::endl;
		MPI_Send(&(serialized_used_idxs[0]), data_size, MPI_CHAR, MASTER, 0, MPI_COMM_WORLD);
		// std::cout << mpi_rank << "end_sending" << "..." << std::endl;
	}
}

void MasterComms::printToMasterOnly(const char *msg) { std::cout << msg << std::endl; }

void MasterComms::setMasterUsedIdxs() {

	MPI_Status status;
	worker_used_idxs.resize(mpi_nump);
	worker_num_used.resize(mpi_nump);

	// Add the master's own used idxs first
	std::set<unsigned int>::iterator it = used_idxs.begin();
	for (; it != used_idxs.end(); ++it) master_used_idxs.insert(*it);

	MPI_Barrier(MPI_COMM_WORLD); // All threads wait

	// Fetch each of the worker used idxs in turn
	for (int i = 1; i < mpi_nump; i++) {

		// Find out the length of the incoming used idxs
		unsigned int data_size = 0;
		MPI_Recv(&data_size, 1, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
		if (data_size == 0) {
			worker_num_used[i] = 0;
			worker_used_idxs[i].clear();
			continue;
		}

		// Receive the (serialized) used idxs for this worker
		std::string serialized_used_idxs;
		serialized_used_idxs.resize(data_size);
		// std::cout << i << "start_recv" << "..." << std::endl;
		MPI_Recv(&(serialized_used_idxs[0]), data_size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
		// std::cout << i << "end_recv" << "..." << std::endl;

		// Unserialize the data and write to worker_used_idxs[i]
		std::istringstream ifs(serialized_used_idxs);
		boost::archive::text_iarchive ar(ifs);
		ar &worker_used_idxs[i];
		worker_num_used[i] = worker_used_idxs[i].size();

		// Populate the master_used_idxs
		std::set<unsigned int>::iterator it = worker_used_idxs[i].begin();
		for (; it != worker_used_idxs[i].end(); ++it) master_used_idxs.insert(*it);
	}
}

float Comms::collectQInMaster(float Q) {

	float Qsum;
	MPI_Barrier(MPI_COMM_WORLD); // All threads wait
	MPI_Reduce(&Q, &Qsum, 1, MPI_FLOAT, MPI_SUM, MASTER, MPI_COMM_WORLD);
	return Qsum; // Note: Only the master has the real Qsum.
}

void Comms::broadcastParamsWeightsOrigMpi(Param *param) {
	std::vector<float> *weights = param->getWeightsPtr();
	MPI_Bcast(&((*weights)[0]), weights->size(), MPI_FLOAT, MASTER, MPI_COMM_WORLD);
}

void Comms::broadcastDropouts(Param *param) {

	auto dropouts = param->getDropoutsPtr();
	if (nullptr != dropouts) MPI_Bcast(&((*dropouts)[0]), dropouts->size(), MPI_CXX_BOOL, MASTER, MPI_COMM_WORLD);
}

void WorkerComms::collectGradsInMaster(std::vector<float> &grads) {

	MPI_Barrier(MPI_COMM_WORLD); // All threads wait
	if (num_used == 0) return;

	std::vector<float> used_grads(num_used);
	std::set<unsigned int>::iterator it = used_idxs.begin();
	for (int i = 0; it != used_idxs.end(); ++it, i++) used_grads[i] = grads[*it];
	MPI_Send(&(used_grads[0]), num_used, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD);
}

void MasterComms::collectGradsInMaster(std::vector<float> &grads) {

	MPI_Status status;
	MPI_Barrier(MPI_COMM_WORLD); // All threads wait

	// Receive and accumulate Gradients
	std::vector<float> used_grads;
	for (int i = 1; i < mpi_nump; i++) {

		if (worker_num_used[i] > 0) {
			used_grads.resize(worker_num_used[i]);
			MPI_Recv(&(used_grads[0]), worker_num_used[i], MPI_FLOAT, i, 0, MPI_COMM_WORLD, &status);

			std::set<unsigned int>::iterator it = worker_used_idxs[i].begin();
			for (int j = 0; it != worker_used_idxs[i].end(); ++it, j++) grads[*it] += used_grads[j];
		}
	}
	// Normarlized By  number of processes
	for (auto &grad : grads) grad /= float(mpi_nump);
}

void WorkerComms::broadcastParamsWeights(Param *param) {

	// Receive updated params from master
	MPI_Status status;
	MPI_Barrier(MPI_COMM_WORLD); // All threads wait
	if (num_used == 0) return;
	std::vector<float> used_params(num_used);
	MPI_Recv(&(used_params[0]), num_used, MPI_FLOAT, MASTER, 0, MPI_COMM_WORLD, &status);

	// Update the params
	std::set<unsigned int>::iterator it = used_idxs.begin();
	for (int i = 0; it != used_idxs.end(); ++it, i++) param->setWeightAtIdx(used_params[i], *it);
}

// careful, this is only a good idea if we only have small amount of process
// MPIBroadcast use tree pattern which could be  more effective
void MasterComms::broadcastParamsWeights(Param *param) {

	MPI_Barrier(MPI_COMM_WORLD); // All threads wait

	// Send each processor the changes of interest to them
	std::vector<float> used_params;
	for (int i = 0; i < mpi_nump; i++) {

		if (worker_num_used[i] == 0) continue;
		used_params.resize(worker_num_used[i]);

		std::set<unsigned int>::iterator it = worker_used_idxs[i].begin();
		for (int j = 0; it != worker_used_idxs[i].end(); ++it, j++) used_params[j] = param->getWeightAtIdx(*it);

		MPI_Send(&(used_params[0]), worker_num_used[i], MPI_FLOAT, i, 0, MPI_COMM_WORLD);
	}
}

int Comms::broadcastIntValue(int value) {

	MPI_Barrier(MPI_COMM_WORLD); // All threads wait
	MPI_Bcast(&value, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	return value;
}

float Comms::broadcastFloatValue(float Q) {

	MPI_Barrier(MPI_COMM_WORLD); // All threads wait
	MPI_Bcast(&Q, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	return Q;
}