#ifndef __AVERAGER_H__
#define __AVERAGER_H__

/* This file implements:

	CAverager: class for averaging last N numbers

	usage:

		CAverager averager(3); //number of numbers to remember

		averager.add(1);
		averager.add(2);
		averager.add(6);
		
		cout << "Average of first 3 is " << averager.getAverage(); //=3

		averager.add(1);
		averager.add(1);
		averager.add(1);

		cout << "Average of last 3 is " << averager.getAverage(); //=1

	Marjan Sterk, 2002
*/

namespace CommonLib {

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

	void reset();
	void add(double x);
	double getAverage();
};

//--------------------------------------------------------------------------------------------
template <class T>
class CFirstNList
//remembers pointers of first N of those T's that are inserted
//T's should have operator < 
{
protected:
	int numToRemember, numRemembered;
	T **buffer;

public:
	CFirstNList(int _numToRemember);
	virtual ~CFirstNList();

	void reset();
	void add(T *item);
	int getLength()
	{
		return numRemembered;
	}
	T *operator[](int i)
	{
		return buffer[i];
	}
};

template <class T>
CFirstNList<T>::CFirstNList(int _numToRemember)
{
	numToRemember = _numToRemember;
	numRemembered = 0;
	buffer = new T*[numToRemember];
}

template <class T>
CFirstNList<T>::~CFirstNList()
{
	for (int i = 0; i < numRemembered; i++)
		delete buffer[i];
	if (buffer)
		delete []buffer;
}

template <class T>
void CFirstNList<T>::reset()
{
	for (int i = 0; i < numRemembered; i++)
		delete buffer[i];
	numRemembered = 0;
}

template <class T>
void CFirstNList<T>::add(T *item)
{
	if (numRemembered == 0) {
		buffer[0] = item;
		numRemembered = 1;
	} else {
		int placeToPut;
		for (placeToPut = 0; placeToPut < numRemembered; placeToPut++)
			if (*item < *buffer[placeToPut])
				break;
		if (numRemembered == numToRemember) {
			if (placeToPut == numRemembered) {
				delete item;
				return;
			}
			delete buffer[numRemembered-1];
		}
		else
			numRemembered++;
		for (int i = numRemembered-1; i > placeToPut; i--)
			buffer[i] = buffer[i-1];
		buffer[placeToPut] = item;
	}
}

};

#endif 


