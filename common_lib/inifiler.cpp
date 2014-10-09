#include "inifiler.h"
#include <ctype.h>
#include <stdio.h>

namespace CommonLib {

IIniSetting::IIniSetting(char *_name)
{
	next = NULL;
	set = false;
	name = new char[strlen(_name)+1];
	strcpy(name, _name);
}

IIniSetting::~IIniSetting()
{
	delete []name;
	if (next)
		delete next;
}

void IIniSetting::insert(IIniSetting *setting)
{
	setting->next = next;
	next = setting;
}

IIniSetting *IIniSetting::find(char *_name)
{
	if (strcmp(name, _name) == 0)
		return this;
	if (next)
		return next->find(_name);
	return NULL;
}

CIniFiler::CIniFiler()
{
	first = NULL;
}

CIniFiler::~CIniFiler()
{
	if (first)
		delete first;
}

void CIniFiler::addSetting(IIniSetting *setting)
{
	if (first)
		first->insert(setting);
	else
		first = setting;
}

void streamEatwhite(std::istream &s)
{
	char c;
	do
		s >> c;
	while (isspace(c));
	s.putback(c);
}

void CIniFiler::readIniFile(char *filename)
{
	numMultipleDefs = 0;
	char settingName[100];
	IIniSetting *setting;
	std::ifstream stream(filename, std::ios::in);
	if (!stream)
		return;
	while (!stream.eof())
	{
		streamEatwhite(stream);
		if (stream.peek() == '%')
		{
			stream.ignore(1000, '\n');
			continue;
		}
		stream.getline(settingName, 100, '=');
		while (isspace(settingName[strlen(settingName)-1]))
			settingName[strlen(settingName)-1] = '\0';
		setting = first->find(settingName);
		if (setting)
		{
			if (setting->isSet())
				numMultipleDefs ++;
			setting->readSetting(stream);
		}
		stream.ignore(100, '\n');
	}
}

void CIniFiler::writeIniFile(char *filename)
{
	char *tempFilename;
	char c, settingName[100];
	IIniSetting *setting;

	tempFilename = new char[strlen(filename)+5];
	strcpy(tempFilename, filename);
	strcat(tempFilename, "_t");

	rename(filename, tempFilename);
	std::ifstream oldStream(tempFilename, std::ios::in);
	std::ofstream newStream(filename, std::ios::out);

	if (!oldStream || !newStream)
		return;

	while (!oldStream.eof())
	{
		c = oldStream.get();
		if (c == -1) //no idea why get() returns -1 without the loop terminating
			continue;
		while (isspace(c))
		{
			newStream << (char)c;
			c = oldStream.get();
		}
		oldStream.putback(c);
		if (oldStream.peek() == '%')
		{
			do {
				c = oldStream.get();
				newStream << (char)c;
			} while (c != '\n');
			continue;
		}
		oldStream.getline(settingName, 100, '=');
		while (isspace(settingName[strlen(settingName)-1]))
			settingName[strlen(settingName)-1] = '\0';
		setting = first->find(settingName);
		if (setting)
		{
			setting->writeSetting(newStream);
			oldStream.ignore(100, '\n');
		} else {
			newStream << settingName;
			do {
				c = oldStream.get();
				newStream << (char)c;
			} while (c != '\n' && c != -1);
		}
	}

	oldStream.close();
	newStream.close();
	remove(tempFilename);
}

int CIniFiler::getNumUnfound()
{
	int numFound = 0, numUnfound = 0;
	for (IIniSetting *setting = first; setting != NULL; setting = setting->next)
		(setting->set ? numFound : numUnfound) ++;
	return numUnfound;
}

int CIniFiler::getNumFound()
{
	int numFound = 0, numUnfound = 0;
	for (IIniSetting *setting = first; setting != NULL; setting = setting->next)
		(setting->set ? numFound : numUnfound) ++;
	return numFound;
}

}; //end namespace CommonLib

