/*#########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# message.cpp
#
# Description: 	Utility class for storing and passing around sparse
#               messages.
#
#
# Copyright (c) 2013, Felicity Allen
# All rights reserved.

# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/

#include "Message.h"
#include "Util.h"
#include <cmath>

Message::Message(unsigned int size) { reset(size); }

double Message::getIdx(unsigned int i) {

	if (!normalized) normalize();
	logdbl dbl = message[i];
	return double(dbl);
}

void Message::addWeightedMessage(Message &msg, double weight) {
	Message::const_iterator it = msg.begin();
	for (; it != msg.end(); it++) addToIdx(it.index(), (double)*it + log(weight));
}

void Message::addToIdx(unsigned int i, double value) {
	logdbl dbl = message[i];
	message[i] = logdbl(logAdd((double)dbl, value));
	empty      = 0;
	normalized = 0;
}

void Message::multIdx(unsigned int i, double value) {
	logdbl dbl = message[i];
	if ((double)dbl > -A_BIG_DBL)
		message[i] = logdbl((double)dbl + value);
	else
		message[i] = logdbl(value);
	empty      = 0;
	normalized = 0;
}

void Message::normalize() {
	log_sum     = -A_BIG_DBL;
	iterator it = message.begin();
	for (; it != message.end(); ++it) log_sum = logAdd(log_sum, (double)*it);
	it = message.begin();
	for (; it != message.end(); ++it) *it = logdbl((double)*it - log_sum);
	log_sum    = 0;
	normalized = 1;
}

void Message::reset(unsigned int size) {
	empty      = 1;
	normalized = 1;
	message.resize(size);
	message.clear();
	log_sum = -A_BIG_DBL;
}

Message &Message::operator=(const Message &rhs) {

	log_sum    = rhs.log_sum;
	empty      = rhs.empty;
	message    = rhs.message;
	normalized = rhs.normalized;
	return *this;
}

void Message::print() {
	if (!normalized) normalize();
	iterator it = message.begin();
	for (; it != message.end(); ++it)
		if (exp((double)*it) > 0.001) std::cout << it.index() << ": " << exp((double)*it) << ", ";
	std::cout << std::endl;
}
