#ifndef __INIFILER_H__
#define __INIFILER_H__

/* This file implements:

	CIniSetting<type> template for setting: can read & write itself
	CIniFiler: ini file manager

	usage:

		CIniFiler iniFiler;
		int windowSize;
		char windowName[100];

		ini.addSetting(new CIniSetting<int>("windowSize", windowSize, 100)); //100 is default value
		ini.addSetting(new CIniSetting<char*>("windowName", windowName, "noname"));
		ini.readIniFile(filename);
		if (ini.getNumUnfound() > 0)
			cerr << "Missing settings.\n";
		if (ini.getNumMultipleDefs() > 0)
			cerr << "Multiply defined settings.\n";
		
		...

		windowSize++;
		ini.writeIniFile(filename);

	Marjan Sterk, 2000-2002
*/

#include <string.h>
#include <fstream>

namespace CommonLib {

class IIniSetting //abstract setting class
{
protected:
	char *name;
	IIniSetting *next;
	bool set;

	virtual void readSetting(std::istream &stream) = 0;
	virtual void writeSetting(std::ostream &stream) = 0;
	void insert(IIniSetting *setting);
	IIniSetting *find(char *_name);

public:
	IIniSetting(char *_name);
	virtual ~IIniSetting();

	bool isSet() { return set; }

	friend class CIniFiler;
};


template<class T> //descendant class template for simple types
class CIniSetting: public IIniSetting
{
protected:
	T &ref;
	T defValue;
	virtual void readSetting(std::istream &stream)
	{
		stream >> ref;
		set = true;
	}
	virtual void writeSetting(std::ostream &stream)
	{ 
		stream << name << " = " << ref << std::endl;
		set = true;
	}
public:
	CIniSetting(char *_name, T &_ref, T _defValue)
		: IIniSetting(_name), ref(_ref)
	{ ref = defValue = _defValue; }
	friend class CIniFiler;
};

template <>
class CIniSetting<char*> : public IIniSetting //specialization for strings
{
protected:
	char *p, *defValue;
	int maxLen;
	virtual void readSetting(std::istream &stream)
	{ 
		stream.get(p, maxLen);
		set = true;
	}
	virtual void writeSetting(std::ostream &stream)
	{ 
		stream << name << "=" << p << std::endl;
		set = true;
	}
public:
	CIniSetting(char *_name, char *_p, int _maxLen, char *_defValue)
		: IIniSetting(_name), p(_p), maxLen(_maxLen)
	{ 
		defValue = new char[strlen(_defValue)+1];
		strcpy(defValue, _defValue);
		strcpy(p, defValue);
	}
	virtual ~CIniSetting()
	{
		if (defValue) delete defValue;
	}
	friend class CIniFiler;
};
/*
class CIniSetting<class T*> : public IIniSetting //specialization for arrays

	- all share the same default
	- p must have enough space
	- iniSetting len is a param that chooses array length
	- len must be a part of iniFiler eventually used with this
	- len must come in ini file before this; otherwise, default len will be used

{
protected:
	T *p, defValue;
	CIniSetting<int> &len;
	virtual void readSetting(istream &stream)
	{
		for (int i = 0; i < len.ref; i++)
		{
			stream >> p[i];
	}
	virtual void writeSetting(istream &stream)
	{
		stream << p[0] << endl;
		for (int i = 0; i < len.ref; i++)
		{
	}
public:
	CIniSetting(char *_name, T *_p, CIniSetting<int> &_len, T _defValue)
		: IInitSetting(_name), p(_p), len(_len)
	{
		defValue = T;
	}
	friend class CIniFiler;
}
*/
class CIniFiler //class for managing ini files
{
	int numMultipleDefs;
	IIniSetting *first;

public:
	CIniFiler();
	virtual ~CIniFiler();

	void addSetting(IIniSetting *setting);
	void readIniFile(char *filename);
	void writeIniFile(char *filename);
	int getNumUnfound();
	int getNumFound();
	int getNumMultipleDefs() { return numMultipleDefs; }
};

}; //end namespace CommonLib

#endif 


