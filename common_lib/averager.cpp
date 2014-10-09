#include "averager.h"
/*
class CAverager
{
protected:
	int numToRemember;
	int first, numRemembered;
	double *buffer;
	double currentAverage;

public:
	CAverager(int _numToRemember);
	virtual ~CAverager();

	void add(double x);
	double getAverage();
};
*/

namespace CommonLib {

CAverager::CAverager(int _numToRemember)
{
	numToRemember = _numToRemember;
	first = -1;
	numRemembered = 0;
	buffer = new double[numToRemember];
	currentAverage = 0;
}

CAverager::~CAverager()
{
	if (buffer)
		delete []buffer;
}

void CAverager::reset()
{
	first = -1;
	numRemembered = 0;
	currentAverage = 0;
}

void CAverager::add(double x)
{
	if (numRemembered == 0) {
		buffer[0] = x;
		first = 0;
		numRemembered = 1;
		currentAverage = x;
	} else if (numRemembered < numToRemember) {
		buffer[numRemembered] = x;
		currentAverage = (currentAverage*numRemembered + x)/(numRemembered+1);
		numRemembered++;
	} else {
		currentAverage -= buffer[first]/numToRemember;
		buffer[first] = x;
		currentAverage += buffer[first]/numToRemember;
		first = (first+1) % numToRemember;
	}
}

double CAverager::getAverage()
{
	return currentAverage;
}

}; //end namespace CommonLib

